#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations) :
      relaxation_strategy(constraint_relaxation_strategy),
      max_iterations(max_iterations) {
}

Iterate GlobalizationMechanism::assemble_trial_iterate(Iterate& current_iterate, Direction& direction, double step_length) {
   if (direction.norm == 0.) {
      add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, step_length, current_iterate.multipliers.constraints);
      copy_from(current_iterate.multipliers.lower_bounds, direction.multipliers.lower_bounds);
      copy_from(current_iterate.multipliers.upper_bounds, direction.multipliers.upper_bounds);
      current_iterate.multipliers.objective = direction.objective_multiplier;
      current_iterate.progress = {0., 0.};
      return current_iterate;
   }
   else {
      Iterate trial_iterate(direction.x.size(), direction.multipliers.constraints.size());
      add_vectors(current_iterate.x, direction.x, step_length, trial_iterate.x);
      add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, step_length, trial_iterate.multipliers.constraints);
      copy_from(trial_iterate.multipliers.lower_bounds, direction.multipliers.lower_bounds);
      copy_from(trial_iterate.multipliers.upper_bounds, direction.multipliers.upper_bounds);
      trial_iterate.multipliers.objective = direction.objective_multiplier;
      return trial_iterate;
   }
}

int GlobalizationMechanism::get_hessian_evaluation_count() const {
   return this->relaxation_strategy.get_hessian_evaluation_count();
}

int GlobalizationMechanism::get_number_subproblems_solved() const {
   return this->relaxation_strategy.get_number_subproblems_solved();
}

void GlobalizationMechanism::print_warning(const char* message) {
   WARNING << RED << message << RESET << "\n";
}