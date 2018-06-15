#include "Problem.hpp"

Problem::Problem(std::string name): name(name) {
}

//! compute || c(feasible_constraints) ||_1 and || c(infeasible_constraints) ||_1 for index sets
// feasible_constraints, infeasible_constraints (overall length m)
double Problem::feasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints) {
	double feasible_residual = 0.;
	
	/* compute residuals for infeasible constraints */
	// TODO useful?
	for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
		int j = constraint_partition.infeasible_set[k];
		if (constraint_partition.status[j] == INFEASIBLE_LOWER) {
			feasible_residual += std::max(0., constraints[j] - this->constraint_ub[j]);
		}
		else {
			feasible_residual += std::max(0., this->constraint_lb[j] - constraints[j]);
		}
	}

	/* compute residuals for feasible constraints */
	for (unsigned int k = 0; k < constraint_partition.feasible_set.size(); k++) {
		int j = constraint_partition.feasible_set[k];
		feasible_residual += std::max(0., this->constraint_lb[j] - constraints[j]);
		feasible_residual += std::max(0., constraints[j] - this->constraint_ub[j]);
	}
	return feasible_residual;
}

double Problem::infeasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints) {
	double infeasible_residual = 0.;
	
	/* compute residuals for infeasible constraints */
	for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
		int j = constraint_partition.infeasible_set[k];
		if (constraint_partition.status[j] == INFEASIBLE_LOWER) {
			infeasible_residual += std::max(0., this->constraint_lb[j] - constraints[j]);
		}
		else {
			infeasible_residual += std::max(0., constraints[j] - this->constraint_ub[j]);
		}
	}
	return infeasible_residual;
}

/* compute ||c||_1 */
double Problem::l1_inf_norm(std::vector<double> constraints) {
	double norm = 0.;

	for (int j = 0; j < this->number_constraints; j++) {
		double residual = std::max(constraints[j] - this->constraint_ub[j], this->constraint_lb[j] - constraints[j]);
		norm += std::max(0., residual);
	}
	return norm;
}