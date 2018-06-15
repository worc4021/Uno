#ifndef PROBLEM_H
#define PROBLEM_H

#include <string>
#include <vector>
#include <map>
#include "Constraint.hpp"
#include "Matrix.hpp"

enum FunctionType {
	LINEAR = 0, /*!< Linear function */
	QUADRATIC, /*!< Quadratic function */
	NONLINEAR /*!< Nonlinear function */
};

/*! \class Problem
* \brief Optimization problem
*
*  Description of an optimization problem
*/
class Problem {
	public:
		Problem(std::string name);
	
		int number_variables; /*!< Number of variables */
		int number_constraints; /*!< Number of constraints */
		
		std::string name;
		
		/* objective */
		double obj_sign; /*!< Sign of the objective function (1: minimization, -1: maximization) */
		std::string objective_name;
		FunctionType objective_type; /*!< Type of the objective (LINEAR, QUADRATIC, NONLINEAR) */
		std::map<int,double> objective_variables;
		virtual double objective(std::vector<double> x) = 0;
		virtual std::vector<double> objective_dense_gradient(std::vector<double> x) = 0;
		virtual std::map<int,double> objective_sparse_gradient(std::vector<double> x) = 0;
		
		/* variables */
		std::vector<std::string> variable_name;
		std::vector<bool> variable_discrete;
		std::vector<double> variable_lb;
		std::vector<double> variable_ub;
		
		/* constraints */
		std::vector<std::string> constraint_name;
		std::vector<std::map<int,double> > constraint_variables;
		std::vector<double> constraint_lb;
		std::vector<double> constraint_ub;
		std::vector<FunctionType> constraints_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
		virtual double evaluate_constraint(int j, std::vector<double> x) = 0;
		virtual std::vector<double> evaluate_constraints(std::vector<double> x) = 0;
		virtual std::vector<double> constraint_dense_gradient(int j, std::vector<double> x) = 0;
		virtual std::map<int,double> constraint_sparse_gradient(int j, std::vector<double> x) = 0;
		
		std::vector<int> jacobian_sparsity;
		virtual std::vector<std::vector<double> > constraints_jacobian_dense(std::vector<double> x) = 0;
		virtual void create_jacobian_sparsity() = 0;
		
		/* Hessian */
		int hessian_maximum_number_nonzero;  /*!< Number of nonzero elements in the Hessian */
		std::vector<int> hessian_column_start; /*!< Column description of sparse Hessian */
		std::vector<int> hessian_row_number; /*!< Row description of sparse Hessian */
		virtual Matrix lagrangian_hessian(std::vector<double> x, double objective_multiplier, std::vector<double> multipliers) = 0;
		
		virtual std::vector<double> primal_initial_solution() = 0;
		virtual std::vector<double> dual_initial_solution() = 0;
		
		double feasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints);
		double infeasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints);
		
		double l1_inf_norm(std::vector<double> constraints);
		
		int number_eval_objective;
		int number_eval_constraints;
		int number_eval_hessian;
};

#endif // PROBLEM_H