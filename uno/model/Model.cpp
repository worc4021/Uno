// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <iostream>
#include <utility>
#include "Model.hpp"
#include "symbolic/VectorExpression.hpp"
#include "linear_algebra/Vector.hpp"

// abstract Problem class
Model::Model(std::string name, size_t number_variables, size_t number_constraints, double objective_sign) :
      name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints), objective_sign(objective_sign) {
}

void Model::project_onto_variable_bounds(std::vector<double>& x) const {
   for (size_t variable_index: Range(this->number_variables)) {
      x[variable_index] = std::max(std::min(x[variable_index], this->variable_upper_bound(variable_index)), this->variable_lower_bound(variable_index));
   }
}

bool Model::is_constrained() const {
   return (0 < this->number_constraints);
}

// individual constraint violation
double Model::constraint_violation(double constraint_value, size_t constraint_index) const {
   const double lower_bound_violation = std::max(0., this->constraint_lower_bound(constraint_index) - constraint_value);
   const double upper_bound_violation = std::max(0., constraint_value - this->constraint_upper_bound(constraint_index));
   return std::max(lower_bound_violation, upper_bound_violation);
}

// compute ||c||
double Model::constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const {
   const VectorExpression constraint_violation(Range(constraints.size()), [&](size_t constraint_index) {
      return this->constraint_violation(constraints[constraint_index], constraint_index);
   });
   return norm(residual_norm, constraint_violation);
}

double Model::linearized_constraint_violation(const std::vector<double>& primal_direction, const std::vector<double>& constraints,
      const RectangularMatrix<double>& constraint_jacobian, double step_length, Norm residual_norm) const {
   // determine the linearized constraint violation term: ||c(x_k) + α ∇c(x_k)^T d||
   const VectorExpression linearized_constraints(Range(this->number_constraints), [&](size_t constraint_index) {
      const double linearized_constraint_j = constraints[constraint_index] + step_length * dot(primal_direction, constraint_jacobian[constraint_index]);
      return this->constraint_violation(linearized_constraint_j, constraint_index);
   });
   return norm(residual_norm, linearized_constraints);
}