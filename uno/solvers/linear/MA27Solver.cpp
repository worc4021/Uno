// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <cassert>
#include <algorithm>
#include "ma27.hpp"
#include "MA27Solver.hpp"
#include "linear_algebra/Vector.hpp"

MA27Solver::MA27Solver(size_t max_dimension, size_t max_number_nonzeros)
    : SymmetricIndefiniteLinearSolver<double>(max_dimension)
    , nz_max(static_cast<int>(max_number_nonzeros))
    , n(static_cast<int>(max_dimension))
    , nnz(static_cast<int>(max_number_nonzeros))
{
   irn.resize(max_number_nonzeros);
   icn.resize(max_number_nonzeros);
   iw.resize(3 * (max_dimension + max_number_nonzeros));
   ikeep.resize(3 * (max_dimension + max_number_nonzeros));
   iw1.resize(max_dimension);

   iflag = 0;
   // set the default values of the controlling parameters
   FC_ma27id(icntl.data(), cntl.data());
   // suppress warning messages
   icntl[1 - 1] = 0;
   icntl[2 - 1] = 0;
   icntl[3 - 1] = 0;
}

void MA27Solver::factorize(const SymmetricMatrix<double> &matrix)
{
   // // general factorization method: symbolic factorization and numerical factorization
   do_symbolic_factorization(matrix);
   do_numerical_factorization(matrix);
}

void MA27Solver::do_symbolic_factorization(const SymmetricMatrix<double> &matrix)
{
   assert(matrix.dimension <= max_dimension && "MA27Solver: the dimension of the matrix is larger than the preallocated size");
   assert(matrix.number_nonzeros <= irn.capacity() &&
          "MA27Solver: the number of nonzeros of the matrix is larger than the preallocated size");

   // build the internal matrix representation
   save_matrix_to_local_format(matrix);

   n = static_cast<int>(matrix.dimension);
   nnz = static_cast<int>(matrix.number_nonzeros);

   // symbolic factorization
   int liw = static_cast<int>(iw.size());
   FC_ma27ad(&n, &nnz, irn.data(), icn.data(), iw.data(), &liw, ikeep.data(), iw1.data(), &nsteps,
             &iflag, icntl.data(), cntl.data(), info.data(), &ops);

   factor.resize(3 * static_cast<std::size_t>(info[4]) / 2);

   std::copy(matrix.data_raw_pointer(), matrix.data_raw_pointer() + matrix.number_nonzeros, factor.begin());

   assert(0 <= info[0] && "MA27: the symbolic factorization failed");
   if (0 < info[0])
   {
      WARNING << "MA27 has issued a warning: info(1) = " << info[0] << '\n';
   }
}

void MA27Solver::do_numerical_factorization([[maybe_unused]]const SymmetricMatrix<double> &matrix)
{
   assert(matrix.dimension <= this->max_dimension && "MA27Solver: the dimension of the matrix is larger than the preallocated size");
   assert(nnz == static_cast<int>(matrix.number_nonzeros) && "MA27Solver: the numbers of nonzeros do not match");

   // numerical factorization
   int la = static_cast<int>(factor.size());
   int liw = static_cast<int>(iw.size());
   FC_ma27bd(&n, &nnz, irn.data(), icn.data(), factor.data(), &la, iw.data(), &liw,
             ikeep.data(), &nsteps, &maxfrt, iw1.data(), icntl.data(), cntl.data(), info.data());
}

void MA27Solver::solve_indefinite_system([[maybe_unused]]const SymmetricMatrix<double> &matrix, const std::vector<double> &rhs, std::vector<double> &result)
{
   // solve
   std::vector<double> w(maxfrt); // double workspace
   int la = static_cast<int>(factor.size());
   int liw = static_cast<int>(iw.size());

   std::copy(rhs.cbegin(), rhs.cend(), result.begin());

   FC_ma27cd(&n, factor.data(), &la, iw.data(), &liw, w.data(), &maxfrt, result.data(), iw1.data(),
             &nsteps, icntl.data(), info.data());
}

std::tuple<size_t, size_t, size_t> MA27Solver::get_inertia() const
{
   // rank = number_positive_eigenvalues + number_negative_eigenvalues
   // n = rank + number_zero_eigenvalues
   const size_t rankA = rank();
   const size_t num_negative_eigenvalues = number_negative_eigenvalues();
   const size_t num_positive_eigenvalues = rankA - num_negative_eigenvalues;
   const size_t num_zero_eigenvalues = static_cast<size_t>(n) - rankA;
   return std::make_tuple(num_positive_eigenvalues, num_negative_eigenvalues, num_zero_eigenvalues);
}

size_t MA27Solver::number_negative_eigenvalues() const
{
   return static_cast<size_t>(info[15 - fortran_shift]);
}

// bool MA27Solver::matrix_is_positive_definite() const {
//    // positive definite = non-singular and no negative eigenvalues
//    return (!matrix_is_singular() && (number_negative_eigenvalues() == 0));
// }

bool MA27Solver::matrix_is_singular() const
{
   return (info[1 - fortran_shift] == -5);
}

size_t MA27Solver::rank() const
{
   return info[1 - fortran_shift] == 3 ? static_cast<size_t>(info[2 - fortran_shift]) : static_cast<size_t>(n);
}

void MA27Solver::save_matrix_to_local_format(const SymmetricMatrix<double> &matrix)
{
   // build the internal matrix representation
   irn.clear();
   icn.clear();

   matrix.for_each([&](size_t row_index, size_t column_index, double /*entry*/)
                   {
      irn.push_back(static_cast<int>(row_index + fortran_shift));
      icn.push_back(static_cast<int>(column_index + fortran_shift)); });
}
