// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <iostream>
#include <cassert>
#include <utility>
#include <algorithm>
#include "Model.hpp"
#include "linear_algebra/VectorExpression.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

// abstract Problem class
Model::Model(std::string name, size_t number_variables, size_t number_constraints) :
      name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints), slacks(number_constraints) {
   this->equality_constraints.reserve(number_constraints);
   this->inequality_constraints.reserve(number_constraints);
   this->lower_bounded_variables.reserve(number_variables);
   this->upper_bounded_variables.reserve(number_variables);
   this->single_lower_bounded_variables.reserve(number_variables);
   this->single_upper_bounded_variables.reserve(number_variables);
}

void Model::determine_bounds_types(std::vector<Interval>& bounds, std::vector<BoundType>& status) {
   assert(bounds.size() == status.size());
   // build the "status" vector as a mapping (map/transform operation) of the "bounds" vector
   std::transform(begin(bounds), end(bounds), begin(status), [](const Interval& bounds_i) {
      if (bounds_i.lb == bounds_i.ub) {
         return EQUAL_BOUNDS;
      }
      else if (is_finite(bounds_i.lb) && is_finite(bounds_i.ub)) {
         return BOUNDED_BOTH_SIDES;
      }
      else if (is_finite(bounds_i.lb)) {
         return BOUNDED_LOWER;
      }
      else if (is_finite(bounds_i.ub)) {
         return BOUNDED_UPPER;
      }
      else {
         return UNBOUNDED;
      }
   });
}

void Model::project_primals_onto_bounds(std::vector<double>& x) const {
   for (size_t variable_index: Range(this->number_variables)) {
      x[variable_index] = std::max(std::min(x[variable_index], this->get_variable_upper_bound(variable_index)), this->get_variable_lower_bound(variable_index));
   }
}

bool Model::is_constrained() const {
   return (0 < this->number_constraints);
}

double Model::compute_constraint_violation(double constraint_value, size_t constraint_index) const {
   const double lower_bound_violation = std::max(0., this->get_constraint_lower_bound(constraint_index) - constraint_value);
   const double upper_bound_violation = std::max(0., constraint_value - this->get_constraint_upper_bound(constraint_index));
   return std::max(lower_bound_violation, upper_bound_violation);
}

// compute ||c||
double Model::compute_constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const {
   VectorExpression<double> constraint_violation(constraints.size(), [&](size_t constraint_index) {
      return this->compute_constraint_violation(constraints[constraint_index], constraint_index);
   });
   return norm(residual_norm, constraint_violation);
}

double Model::compute_linearized_constraint_violation(const std::vector<double>& primal_direction, const std::vector<double>& constraints,
      const RectangularMatrix<double>& constraint_jacobian, double step_length, Norm residual_norm) const {
   // determine the linearized constraint violation term: ||c(x_k) + α ∇c(x_k)^T d||
   VectorExpression<double> linearized_constraints(this->number_constraints, [&](size_t constraint_index) {
      const double linearized_constraint_j = constraints[constraint_index] + step_length * dot(primal_direction, constraint_jacobian[constraint_index]);
      return this->compute_constraint_violation(linearized_constraint_j, constraint_index);
   });
   return norm(residual_norm, linearized_constraints);
}