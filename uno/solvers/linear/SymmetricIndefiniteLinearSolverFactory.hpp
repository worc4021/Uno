// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERFACTORY_H
#define UNO_LINEARSOLVERFACTORY_H

#include <memory>
#include "SymmetricIndefiniteLinearSolver.hpp"

#ifdef HAS_MA57
#include "MA57Solver.hpp"
#endif

namespace uno {
   class SymmetricIndefiniteLinearSolverFactory {
   public:
      static std::unique_ptr<SymmetricIndefiniteLinearSolver<double>> create([[maybe_unused]] const std::string& linear_solver_name, [[maybe_unused]] size_t dimension, [[maybe_unused]] size_t number_nonzeros) {
   #ifdef HAS_MA57
         if (linear_solver_name == "MA57") {
            return std::make_unique<MA57Solver>(dimension, number_nonzeros);
         }
   #endif
         throw std::invalid_argument("Linear solver name is unknown");
      }

      // return the list of available QP solvers
      static std::vector<std::string> available_solvers() {
         std::vector<std::string> solvers{};
         #ifdef HAS_MA57
         solvers.emplace_back("MA57");
         #endif
         return solvers;
      }
   };
} // namespace

#endif // UNO_LINEARSOLVERFACTORY_H
