// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "Uno.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
      globalization_mechanism(globalization_mechanism),
      tolerance(options.get_double("tolerance")),
      max_iterations(options.get_unsigned_int("max_iterations")),
      small_step_factor(options.get_double("small_step_factor")) {
}

Result Uno::solve(const Model& model, Iterate& current_iterate, const Options& options) {
   Timer timer{};
   timer.start();
   size_t major_iterations = 0;

   std::cout << "\nProblem " << model.name << '\n';
   std::cout << model.number_variables << " variables, " << model.number_constraints << " constraints\n";
   std::cout << "Problem type: " << Model::type_to_string[model.problem_type] << "\n\n";

   Statistics statistics = Uno::create_statistics(model, options);

   // use the current point to initialize the strategies and generate the initial iterate
   this->globalization_mechanism.initialize(statistics, current_iterate);

   TerminationStatus termination_status = NOT_OPTIMAL;
   try {
      // check for termination
      while (!this->termination_criterion(termination_status, major_iterations)) {
         statistics.new_line();
         major_iterations++;
         DEBUG << "### Outer iteration " << major_iterations << '\n';

         // compute an acceptable iterate by solving a subproblem at the current point
         auto [new_iterate, step_norm] = this->globalization_mechanism.compute_acceptable_iterate(statistics, current_iterate);

         // compute the status of the new iterate
         termination_status = this->check_termination(model, new_iterate, step_norm);
         Uno::add_statistics(statistics, model, new_iterate, major_iterations);
         if (Logger::logger_level == INFO) statistics.print_current_line();

         current_iterate = std::move(new_iterate);
      }
   }
   catch (std::exception& exception) {
      ERROR << exception.what();
   }
   if (Logger::logger_level == INFO) statistics.print_footer();
   timer.stop();

   const size_t number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
   const size_t hessian_evaluation_count = this->globalization_mechanism.get_hessian_evaluation_count();
   Result result = {termination_status, std::move(current_iterate), model.number_variables, model.number_constraints, major_iterations,
         timer.get_duration(), Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_jacobian, hessian_evaluation_count,
          number_subproblems_solved};
   return result;
}

Statistics Uno::create_statistics(const Model& model, const Options& options) {
   Statistics statistics(options);
   statistics.add_column("major", Statistics::int_width, options.get_int("statistics_major_column_order"));
   statistics.add_column("minor", Statistics::int_width, options.get_int("statistics_minor_column_order"));
   statistics.add_column("step norm", Statistics::double_width, options.get_int("statistics_step_norm_column_order"));
   statistics.add_column("objective", Statistics::double_width, options.get_int("statistics_objective_column_order"));
   if (model.is_constrained()) {
      statistics.add_column("primal infeas.", Statistics::double_width, options.get_int("statistics_primal_infeasibility_column_order"));
   }
   statistics.add_column("dual infeas.", Statistics::double_width, options.get_int("statistics_dual_infeasibility_column_order"));
   statistics.add_column("complementarity", Statistics::double_width, options.get_int("statistics_complementarity_column_order"));
   statistics.add_column("stationarity", Statistics::double_width, options.get_int("statistics_stationarity_column_order"));
   return statistics;
}

void Uno::add_statistics(Statistics& statistics, const Model& model, const Iterate& iterate, size_t major_iterations) {
   statistics.add_statistic(std::string("major"), major_iterations);
   statistics.add_statistic("objective", iterate.model_evaluations.objective);
   if (model.is_constrained()) {
      statistics.add_statistic("primal infeas.", iterate.primal_constraint_violation);
   }
   statistics.add_statistic("dual infeas.", 0.); // TODO
   statistics.add_statistic("complementarity", iterate.complementarity_error);
   statistics.add_statistic("stationarity", iterate.stationarity_error);
}

bool Uno::termination_criterion(TerminationStatus current_status, size_t iteration) const {
   return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

bool not_all_zero_multipliers(const Model& model, const Multipliers& multipliers, double tolerance) {
   // constraint multipliers
   for (double multiplier_j: multipliers.constraints) {
      if (tolerance < std::abs(multiplier_j)) {
         return true;
      }
   }
   // bound multipliers
   for (size_t i = 0; i < model.number_variables; i++) {
      if (tolerance < std::abs(multipliers.lower_bounds[i] + multipliers.upper_bounds[i])) {
         return true;
      }
   }
   return false;
}

TerminationStatus Uno::check_termination(const Model& model, Iterate& current_iterate, double step_norm) const {
   // evaluate termination conditions based on optimality conditions
   const bool stationarity = (current_iterate.stationarity_error <= this->tolerance * std::sqrt(model.number_variables));
   const bool complementarity = (current_iterate.complementarity_error <= this->tolerance * static_cast<double>(model.number_variables + model.number_constraints));
   const bool primal_feasibility = (current_iterate.primal_constraint_violation <= this->tolerance * static_cast<double>(model.number_variables));
   // TODO dual feasibility

   if (stationarity && complementarity) {
      if (primal_feasibility) {
         if (0. < current_iterate.multipliers.objective) {
            // feasible regular stationary point
            return FEASIBLE_KKT_POINT;
         }
         else if (current_iterate.multipliers.objective == 0. && not_all_zero_multipliers(model, current_iterate.multipliers, this->tolerance)) {
            // feasible but CQ failure
            return FJ_POINT;
         }
      }
      else if (current_iterate.multipliers.objective == 0. && not_all_zero_multipliers(model, current_iterate.multipliers, this->tolerance)) {
         // no primal feasibility, minimum of constraint violation
        return INFEASIBLE_KKT_POINT;
      }
   }
   // stationarity & complementarity not achieved, but we can terminate with a small step
   if (step_norm <= this->tolerance / this->small_step_factor) {
      if (primal_feasibility) {
         return FEASIBLE_SMALL_STEP;
      }
      else {
         return INFEASIBLE_SMALL_STEP;
      }
   }
   return NOT_OPTIMAL;
}


