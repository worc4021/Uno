// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModel.hpp"

namespace uno {
   // zero Hessian
   class ZeroHessian : public HessianModel {
   public:
      ZeroHessian() = default;

      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void evaluate_hessian(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
            const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const OptimizationProblem& problem, const Vector<double>& vector,
         const Vector<double>& constraint_multipliers, Vector<double>& result) override;
   };
} // namespace