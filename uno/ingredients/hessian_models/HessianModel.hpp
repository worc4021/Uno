// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODEL_H
#define UNO_HESSIANMODEL_H

#include <cstddef>

namespace uno {
   // forward declarations
   class OptimizationProblem;
   class Options;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;

   class HessianModel {
   public:
      HessianModel() = default;
      virtual ~HessianModel();

      size_t evaluation_count{0};

      virtual void initialize_statistics(Statistics& statistics, const Options& options) const = 0;
      virtual void evaluate_hessian(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
            const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) = 0;
      virtual void compute_hessian_vector_product(const OptimizationProblem& problem, const Vector<double>& vector,
         const Vector<double>& constraint_multipliers, Vector<double>& result) = 0;

   };
} // namespace

#endif // UNO_HESSIANMODEL_H