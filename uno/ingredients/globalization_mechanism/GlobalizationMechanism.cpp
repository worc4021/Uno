// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options) :
      constraint_relaxation_strategy(constraint_relaxation_strategy),
      direction(this->constraint_relaxation_strategy.maximum_number_variables(), this->constraint_relaxation_strategy.maximum_number_constraints()),
      tight_tolerance(options.get_double("tolerance")),
      loose_tolerance(options.get_double("loose_tolerance")),
      loose_tolerance_consecutive_iteration_threshold(options.get_unsigned_int("loose_tolerance_consecutive_iteration_threshold")),
      progress_norm(norm_from_string(options.get_string("progress_norm"))),
      unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")) {
}

void GlobalizationMechanism::assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
      double primal_step_length, double dual_step_length, double bound_dual_step_length) {
   trial_iterate.set_number_variables(current_iterate.primals.size());
   // take primal step
   add_vectors(current_iterate.primals, direction.primals, primal_step_length, trial_iterate.primals);
   // project the trial iterate onto the bounds to avoid numerical errors
   model.project_onto_variable_bounds(trial_iterate.primals);
   // take dual step: line-search carried out only on constraint multipliers. Bound multipliers updated with full bound dual step length
   add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, dual_step_length, trial_iterate.multipliers.constraints);
   add_vectors(current_iterate.multipliers.lower_bounds, direction.multipliers.lower_bounds, bound_dual_step_length, trial_iterate.multipliers.lower_bounds);
   add_vectors(current_iterate.multipliers.upper_bounds, direction.multipliers.upper_bounds, bound_dual_step_length, trial_iterate.multipliers.upper_bounds);
   trial_iterate.progress.reset();
   trial_iterate.is_objective_computed = false;
   trial_iterate.is_objective_gradient_computed = false;
   trial_iterate.are_constraints_computed = false;
   trial_iterate.is_constraint_jacobian_computed = false;
   trial_iterate.status = TerminationStatus::NOT_OPTIMAL;
}

TerminationStatus GlobalizationMechanism::check_termination(const Model& model, Iterate& current_iterate) {
   // test convergence wrt the tight tolerance
   const TerminationStatus status_tight_tolerance = this->check_convergence_with_given_tolerance(model, current_iterate, this->tight_tolerance);
   if (status_tight_tolerance != TerminationStatus::NOT_OPTIMAL || this->loose_tolerance <= this->tight_tolerance) {
      return status_tight_tolerance;
   }

   // if not converged, check convergence wrt loose tolerance (provided it is strictly looser than the tight tolerance)
   const TerminationStatus status_loose_tolerance = this->check_convergence_with_given_tolerance(model, current_iterate, this->loose_tolerance);
   // if converged, keep track of the number of consecutive iterations
   if (status_loose_tolerance != TerminationStatus::NOT_OPTIMAL) {
      this->loose_tolerance_consecutive_iterations++;
   }
   else {
      this->loose_tolerance_consecutive_iterations = 0;
      return TerminationStatus::NOT_OPTIMAL;
   }
   // check if loose tolerance achieved for enough consecutive iterations
   if (this->loose_tolerance_consecutive_iteration_threshold <= this->loose_tolerance_consecutive_iterations) {
      return status_loose_tolerance;
   }
   else {
      return TerminationStatus::NOT_OPTIMAL;
   }
}

TerminationStatus GlobalizationMechanism::check_convergence_with_given_tolerance(const Model& model, Iterate& current_iterate, double tolerance) const {
   // evaluate termination conditions based on optimality conditions
   const bool optimality_stationarity =
         (current_iterate.residuals.optimality_stationarity / current_iterate.residuals.stationarity_scaling <= tolerance);
   const bool feasibility_stationarity =
         (current_iterate.residuals.feasibility_stationarity / current_iterate.residuals.stationarity_scaling <= tolerance);
   const bool optimality_complementarity =
         (current_iterate.residuals.optimality_complementarity / current_iterate.residuals.complementarity_scaling <= tolerance);
   const bool feasibility_complementarity =
         (current_iterate.residuals.feasibility_complementarity / current_iterate.residuals.complementarity_scaling <= tolerance);
   const bool primal_feasibility = (current_iterate.residuals.infeasibility <= tolerance);
   const bool no_trivial_duals = current_iterate.multipliers.not_all_zero(model.number_variables, tolerance);

   DEBUG << "Termination criteria for tolerance = " << tolerance << ":\n";
   DEBUG << "Stationarity (optimality): " << std::boolalpha << optimality_stationarity << '\n';
   DEBUG << "Stationarity (feasibility): " << std::boolalpha << feasibility_stationarity << '\n';
   DEBUG << "Complementarity (optimality): " << std::boolalpha << optimality_complementarity << '\n';
   DEBUG << "Complementarity (feasibility): " << std::boolalpha << feasibility_complementarity << '\n';
   DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';
   DEBUG << "Not all zero multipliers: " << std::boolalpha << no_trivial_duals << "\n\n";

   if (current_iterate.is_objective_computed && current_iterate.evaluations.objective < this->unbounded_objective_threshold) {
      return TerminationStatus::UNBOUNDED;
   }
   else if (optimality_complementarity && primal_feasibility) {
      if (0. < current_iterate.objective_multiplier && optimality_stationarity) {
         // feasible regular stationary point
         return TerminationStatus::FEASIBLE_KKT_POINT;
      }
      else if (feasibility_stationarity && no_trivial_duals) {
         // feasible but violation of CQ
         return TerminationStatus::FEASIBLE_FJ_POINT;
      }
   }
   else if (feasibility_complementarity && feasibility_stationarity) {
      // no primal feasibility, stationary point of constraint violation
      return TerminationStatus::INFEASIBLE_STATIONARY_POINT;
   }
   return TerminationStatus::NOT_OPTIMAL;
}

size_t GlobalizationMechanism::get_hessian_evaluation_count() const {
   return this->constraint_relaxation_strategy.get_hessian_evaluation_count();
}

size_t GlobalizationMechanism::get_number_subproblems_solved() const {
   return this->constraint_relaxation_strategy.get_number_subproblems_solved();
}
