#include <iostream>


#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "hs71.hpp"
#include "Uno.hpp"
#include "optimization/ModelFactory.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"
#include "tools/Timer.hpp"

size_t memory_allocation_amount = 0;
extern void find_preset(const std::string& preset_name, Options& options);

Statistics create_statistics(const Model& model, const Options& options) {
   Statistics statistics(options);
   statistics.add_column("iters", Statistics::int_width, options.get_int("statistics_major_column_order"));
   statistics.add_column("step norm", Statistics::double_width, options.get_int("statistics_step_norm_column_order"));
   statistics.add_column("objective", Statistics::double_width, options.get_int("statistics_objective_column_order"));
   if (model.is_constrained()) {
      statistics.add_column("primal infeas.", Statistics::double_width, options.get_int("statistics_primal_infeasibility_column_order"));
   }
   statistics.add_column("complementarity", Statistics::double_width, options.get_int("statistics_complementarity_column_order"));
   statistics.add_column("stationarity", Statistics::double_width, options.get_int("statistics_stationarity_column_order"));
   return statistics;
}

void run_uno(const Options& options) {
   // Model
   std::unique_ptr<Model> hs71 = std::make_unique<HS71>();

   // initialize initial primal and dual points
   Iterate initial_iterate(hs71->number_variables, hs71->number_constraints);
   hs71->get_initial_primal_point(initial_iterate.primals);
   hs71->get_initial_dual_point(initial_iterate.multipliers.constraints);
   hs71->project_primals_onto_bounds(initial_iterate.primals);

   // reformulate (scale, add slacks, relax the bounds, ...) if necessary
   std::unique_ptr<Model> model = ModelFactory::reformulate(std::move(hs71), initial_iterate, options);

   // create the statistics
   Statistics statistics = create_statistics(*model, options);

   // create the constraint relaxation strategy
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(statistics, *model, options);

   // create the globalization mechanism
   auto globalization_mechanism = GlobalizationMechanismFactory::create(statistics, *constraint_relaxation_strategy, options);

   // instantiate the combination of ingredients and solve the problem
   Uno uno = Uno(*globalization_mechanism, options);
   try {
      Result result = uno.solve(statistics, *model, initial_iterate, options);

      // print the optimization summary
      std::string combination = options.get_string("globalization_mechanism") + " " + options.get_string("constraint_relaxation_strategy") + " " +
                                options.get_string("globalization_strategy") + " " + options.get_string("subproblem");
      std::cout << "\nUno (" << combination << ")\n";
      std::cout << Timer::get_current_date();
      std::cout << "────────────────────────────────────────\n";
      const bool print_solution = options.get_bool("print_solution");
      result.print(print_solution);
      std::cout << "memory_allocation_amount = " << memory_allocation_amount << '\n';
   }
   catch (const std::exception& e) {
      std::cout << "Uno terminated with an error\n";
   }
}

Level Logger::level = Level::WARNING;

void print_uno_version() {
   std::cout << "Welcome in Uno 1.0\n";
   std::cout << "To solve an AMPL model, type ./uno_ampl path_to_file/file.nl\n";
   std::cout << "To choose a constraint relaxation strategy, use the argument -constraint_relaxation_strategy "
                "[feasibility_restoration|l1_relaxation]\n";
   std::cout << "To choose a subproblem method, use the argument -subproblem [QP|LP|primal_dual_interior_point]\n";
   std::cout << "To choose a globalization mechanism, use the argument -globalization_mechanism [LS|TR]\n";
   std::cout << "To choose a globalization strategy, use the argument -globalization_strategy "
                "[l1_merit|leyffer_filter_method|waechter_filter_method]\n";
   std::cout << "To choose a preset, use the argument -preset [filtersqp|ipopt|byrd]\n";
   std::cout << "The options can be combined in the same command line. Autocompletion is possible (see README).\n";
}



int main(int argc, char* argv[]) {

    // if (1 < argc) {
      // get the default options
      Options options = get_default_options("uno.options");
      // override them with the command line options
      
      find_preset("ipopt", options);
      options["linear_solver"] = "MA27";

      Logger::set_logger(options.get_string("logger"));

      //if (std::string(argv[1]) == "-v") {
      //   print_uno_version();
      //}
      //else if (std::string(argv[1]) == "--strategies") {
      //   Uno::print_available_strategies();
      //}
      //else {
         options.print();
         // run Uno on the .nl file (last command line argument)
    //  std::string model_name = std::string(argv[argc - 1]);
         run_uno(options);
      //}
   //}
   //else {
   //   print_uno_version();
   //}
   return EXIT_SUCCESS;
}