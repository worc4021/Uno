// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <cassert>
#include "ma27.hpp"
#include "MA27Solver.hpp"
#include "linear_algebra/Vector.hpp"

MA27Solver::MA27Solver(size_t max_dimension, size_t max_number_nonzeros) : SymmetricIndefiniteLinearSolver<double>(max_dimension),
   iwork(5 * max_dimension),
   lwork(static_cast<int>(1.2 * static_cast<double>(max_dimension))),
   work(static_cast<size_t>(this->lwork)), residuals(max_dimension) {
   this->row_indices.reserve(max_number_nonzeros);
   this->column_indices.reserve(max_number_nonzeros);
   // set the default values of the controlling parameters
   FC_ma27id(this->cntl.data(), this->icntl.data());
   // suppress warning messages
   this->icntl[4] = 0;
   // iterative refinement enabled
   this->icntl[8] = 1;
}

void MA27Solver::factorize(const SymmetricMatrix<double>& matrix) {
   // general factorization method: symbolic factorization and numerical factorization
   this->do_symbolic_factorization(matrix);
   this->do_numerical_factorization(matrix);
}

void MA27Solver::do_symbolic_factorization(const SymmetricMatrix<double>& matrix) {
   assert(matrix.dimension <= this->max_dimension && "MA27Solver: the dimension of the matrix is larger than the preallocated size");
   assert(matrix.number_nonzeros <= this->row_indices.capacity() &&
      "MA27Solver: the number of nonzeros of the matrix is larger than the preallocated size");

   // build the internal matrix representation
   this->save_matrix_to_local_format(matrix);

   const int n = static_cast<int>(matrix.dimension);
   const int nnz = static_cast<int>(matrix.number_nonzeros);

   // sparsity pattern
   const int lkeep = 5 * n + nnz + std::max(n, nnz) + 42;
   std::vector<int> keep(static_cast<size_t>(lkeep));

   // symbolic factorization
   FC_ma27ad(/* const */ &n,
         /* const */ &nnz,
         /* const */ this->row_indices.data(),
         /* const */ this->column_indices.data(),
         /* const */ &lkeep,
         /* const */ keep.data(),
         /* out */ this->iwork.data(),
         /* const */ this->icntl.data(),
         /* out */ this->info.data(),
         /* out */ this->rinfo.data());

   assert(0 <= info[0] && "MA27: the symbolic factorization failed");
   if (0 < info[0]) {
      WARNING << "MA27 has issued a warning: info(1) = " << info[0] << '\n';
   }
   int lfact = 2 * this->info[8];
   std::vector<double> fact(static_cast<size_t>(lfact));
   int lifact = 2 * this->info[9];
   std::vector<int> ifact(static_cast<size_t>(lifact));

   // store the symbolic factorization
   this->factorization = {n, nnz, std::move(fact), lfact, std::move(ifact), lifact, lkeep, std::move(keep)};
}

void MA27Solver::do_numerical_factorization(const SymmetricMatrix<double>& matrix) {
   assert(matrix.dimension <= this->max_dimension && "MA27Solver: the dimension of the matrix is larger than the preallocated size");
   assert(this->factorization.nnz == static_cast<int>(matrix.number_nonzeros) && "MA27Solver: the numbers of nonzeros do not match");

   const int n = static_cast<int>(matrix.dimension);
   // numerical factorization
   FC_ma27bd(&n,
         &this->factorization.nnz,
         /* const */ matrix.data_raw_pointer(),
         /* out */ this->factorization.fact.data(),
         /* const */ &this->factorization.lfact,
         /* out*/ this->factorization.ifact.data(),
         /* const */ &this->factorization.lifact,
         /* const */ &this->factorization.lkeep,
         /* const */ this->factorization.keep.data(), this->iwork.data(), this->icntl.data(), this->cntl.data(),
         /* out */ this->info.data(),
         /* out */ this->rinfo.data());
}

void MA27Solver::solve_indefinite_system(const SymmetricMatrix<double>& matrix, const std::vector<double>& rhs, std::vector<double>& result) {
   // solve
   const int n = static_cast<int>(matrix.dimension);
   const int lrhs = n; // integer, length of rhs

   // solve the linear system
   if (this->use_iterative_refinement) {
      FC_ma27dd(&this->job, &n, &this->factorization.nnz, matrix.data_raw_pointer(), this->row_indices.data(), this->column_indices.data(),
            this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(), &this->factorization.lifact,
            rhs.data(), result.data(), this->residuals.data(), this->work.data(), this->iwork.data(), this->icntl.data(),
            this->cntl.data(), this->info.data(), this->rinfo.data());
   }
   else {
      // copy rhs into result (overwritten by MA27)
      copy_from(result, rhs);

      FC_ma27cd(&this->job, &n, this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(),
            &this->factorization.lifact, &this->nrhs, result.data(), &lrhs, this->work.data(), &this->lwork, this->iwork.data(),
            this->icntl.data(), this->info.data());
   }
}

std::tuple<size_t, size_t, size_t> MA27Solver::get_inertia() const {
   // rank = number_positive_eigenvalues + number_negative_eigenvalues
   // n = rank + number_zero_eigenvalues
   const size_t rank = this->rank();
   const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
   const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
   const size_t number_zero_eigenvalues = static_cast<size_t>(this->factorization.n) - rank;
   return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
}

size_t MA27Solver::number_negative_eigenvalues() const {
   return static_cast<size_t>(this->info[23]);
}

/*
bool MA27Solver::matrix_is_positive_definite() const {
   // positive definite = non-singular and no negative eigenvalues
   return not this->matrix_is_singular() && this->number_negative_eigenvalues() == 0;
}
*/

bool MA27Solver::matrix_is_singular() const {
   return (this->info[0] == 4);
}

size_t MA27Solver::rank() const {
   return static_cast<size_t>(this->info[24]);
}

void MA27Solver::save_matrix_to_local_format(const SymmetricMatrix<double>& matrix) {
   // build the internal matrix representation
   this->row_indices.clear();
   this->column_indices.clear();
   matrix.for_each([&](size_t row_index, size_t column_index, double /*entry*/) {
      this->row_indices.push_back(static_cast<int>(row_index + this->fortran_shift));
      this->column_indices.push_back(static_cast<int>(column_index + this->fortran_shift));
   });
}
