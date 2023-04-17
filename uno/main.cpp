// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "preprocessing/Preprocessing.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "interfaces/AMPL/AMPLModel.hpp"
#include "Uno.hpp"
#include "optimization/ModelFactory.hpp"
#include "optimization/ScaledModel.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"
#include "tools/Timer.hpp"

size_t memory_allocation_amount = 0;

void* operator new(size_t size) {
   memory_allocation_amount += size;
   return malloc(size);
}

void run_uno_ampl(const std::string& model_name, const Options& options) {
   // AMPL model
   std::unique_ptr<Model> ampl_model = std::make_unique<AMPLModel>(model_name);

   // initialize initial primal and dual points
   Iterate first_iterate(ampl_model->number_variables, ampl_model->number_constraints);
   ampl_model->get_initial_primal_point(first_iterate.primals);
   ampl_model->get_initial_dual_point(first_iterate.multipliers.constraints);
   ampl_model->project_primals_onto_bounds(first_iterate.primals);

   // reformulate (scale, add slacks, relax the bounds, ...) if necessary
   std::unique_ptr<Model> model = ModelFactory::reformulate(std::move(ampl_model), first_iterate, options);

   // enforce linear constraints at initial point
   if (options.get_bool("enforce_linear_constraints")) {
      Preprocessing::enforce_linear_constraints(options, *model, first_iterate.primals, first_iterate.multipliers);
   }

   // create the constraint relaxation strategy
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(*model, options);

   // create the globalization mechanism
   auto mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);

   // instantiate the combination of ingredients and solve the problem
   Uno uno = Uno(*mechanism, options);
   Result result = uno.solve(*model, first_iterate, options);

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

Level Logger::logger_level = INFO;

int main(int argc, char* argv[]) {
   if (1 < argc) {
      // get the default options
      Options options = get_default_options("uno.options");
      // override them with the command line options
      get_command_line_options(argc, argv, options);
      set_logger(options.get_string("logger"));

      if (std::string(argv[1]) == "-v") {
         std::cout << "Welcome in Uno 1.0\n";
         std::cout << "To solve an AMPL model, type ./uno_ampl path_to_file/file.nl\n";
         std::cout << "To choose a globalization mechanism, use the argument -globalization_mechanism [LS|TR]\n";
         std::cout << "To choose a constraint relaxation strategy, use the argument -constraint_relaxation [feasibility-restoration|l1-relaxation]\n";
         std::cout << "To choose a globalization strategy, use the argument -globalization_strategy [l1-merit|leyffer-filter-strategy|waechter-filter-strategy]\n";
         std::cout << "To choose a subproblem method, use the argument -subproblem [QP|LP|primal_dual_interior_point]\n";
         std::cout << "To choose a preset, use the argument -preset [filtersqp|ipopt|byrd]\n";
         std::cout << "The options can be combined in the same command line. Autocompletion is possible.\n";
      }
      else if (std::string(argv[1]) == "--strategies") {
         Uno::print_available_strategies();
      }
      else {
         options.print();
         // run Uno on the .nl file (last command line argument)
         std::string model_name = std::string(argv[argc - 1]);
         run_uno_ampl(model_name, options);
      }
   }
   return EXIT_SUCCESS;
}
