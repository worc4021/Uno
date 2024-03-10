// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA27SOLVER_H
#define UNO_MA27SOLVER_H

#include <array>
#include <vector>
#include "SymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"

/*! \class MA27Solver
 * \brief Interface for MA27
 * see https://github.com/YimingYAN/linSolve
 *
 *  Interface to the symmetric indefinite linear solver MA27
 */
class MA27Solver : public SymmetricIndefiniteLinearSolver<double>
{
public:
   explicit MA27Solver(size_t max_dimension, size_t max_number_nonzeros);
   ~MA27Solver() override {
      SymmetricIndefiniteLinearSolver<double>::~SymmetricIndefiniteLinearSolver();
   };

   void factorize(const SymmetricMatrix<double> &matrix) override;
   void do_symbolic_factorization(const SymmetricMatrix<double> &matrix) override;
   void do_numerical_factorization(const SymmetricMatrix<double> &matrix) override;
   void solve_indefinite_system(const SymmetricMatrix<double> &matrix, const std::vector<double> &rhs, std::vector<double> &result) override;

   [[nodiscard]] std::tuple<size_t, size_t, size_t> get_inertia() const override;
   [[nodiscard]] size_t number_negative_eigenvalues() const override;
   // [[nodiscard]] bool matrix_is_positive_definite() const override;
   [[nodiscard]] bool matrix_is_singular() const override;
   [[nodiscard]] size_t rank() const override;

private:
   int nz_max{};         // maximal number of nonzeros entries
   int n{};			   // dimension of current factorisation (maximal value here is <= max_dimension)
   int nnz{};			   // number of nonzeros of current factorisation
   std::array<int,30> icntl{};  // integer array of length 30; integer control values
   std::array<double,5> cntl{}; // double array of length 5; double control values

   std::vector<int> irn{}; // row index of input
   std::vector<int> icn{}; // col index of input

   std::vector<int> iw{};    // integer workspace of length liw
   std::vector<int> ikeep{}; // integer array of 3*n; pivot sequence
   std::vector<int> iw1{};   // integer workspace array of length n
   int nsteps{};             // integer, to be set by ma27
   int iflag{};              // integer; 0 if pivot order chosen automatically; 1 if the pivot order set by ikeep
   std::array<int,20> info{};           // integer array of length 20
   double ops{};             // double, operations count

   std::vector<double> factor{}; // data array of length la;
   int maxfrt{};                 // integer, to be set by ma27

   const size_t fortran_shift{1};

   // bool use_iterative_refinement{false}; // Not sure how to do this with ma27
   void save_matrix_to_local_format(const SymmetricMatrix<double> &row_index);
};

#endif // UNO_MA27SOLVER_H
