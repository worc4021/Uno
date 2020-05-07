#include <exception>
#include "HessianEvaluation.hpp"
#include "MA57Solver.hpp"
#include "Utils.hpp"

HessianEvaluation::HessianEvaluation(int dimension) : dimension(dimension), convexify(false) {
}

HessianEvaluation::~HessianEvaluation() {
}

CSCMatrix HessianEvaluation::modify_inertia(CSCMatrix& hessian) {
    MA57Solver solver;
    double beta = 1e-4;
    
    // Nocedal and Wright, p51
    double smallest_diagonal_entry = hessian.smallest_diagonal_entry();
    DEBUG << "The minimal diagonal entry of the Hessian is " << hessian.smallest_diagonal_entry() << "\n";
    
    double inertia = 0.;
    double previous_inertia = 0.;
    if (smallest_diagonal_entry <= 0.) {
        inertia = beta - smallest_diagonal_entry;
    }
    
    if (0. < inertia) {
        hessian = hessian.add_identity_multiple(inertia - previous_inertia);
    }
    COOMatrix coo_hessian = hessian.to_COO();
    DEBUG << "Testing factorization with inertia term " << inertia << "\n";
    MA57Factorization factorization = solver.factorize(coo_hessian);
    
    bool good_inertia = false;
    while (!good_inertia) {
        DEBUG << factorization.number_negative_eigenvalues() << " negative eigenvalues\n";
        if (!factorization.matrix_is_singular() && factorization.number_negative_eigenvalues() == 0) {
            good_inertia = true;
            DEBUG << "Factorization was a success with inertia " << inertia << "\n";
        }
        else {
            if (inertia == 0.) {
                inertia = beta;
                previous_inertia = 0.;
            }
            else {
                previous_inertia = inertia;
                inertia *= 2.;
            }
            hessian = hessian.add_identity_multiple(inertia - previous_inertia);
            coo_hessian = hessian.to_COO();
            DEBUG << "Testing factorization with inertia term " << inertia << "\n";
            factorization = solver.factorize(coo_hessian);
        }
    }
    return hessian;
}

/* Exact Hessian */

ExactHessianEvaluation::ExactHessianEvaluation(int dimension): HessianEvaluation(dimension) {
}

void ExactHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    /* compute Hessian */
    iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    
    if (this->convexify) {
        //std::cout.precision(17);
        DEBUG << "hessian before convexification: " << iterate.hessian;
        /* modify the inertia to make the problem strictly convex */
        iterate.hessian = this->modify_inertia(iterate.hessian);
    }
    return;
}

/* BFGS Hessian */

BFGSHessianEvaluation::BFGSHessianEvaluation(int dimension): HessianEvaluation(dimension), previous_hessian(dimension, 1), previous_x(dimension) {
}

void BFGSHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    // the BFGS Hessian is already positive definite, do not convexify
    iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    return;
}

/* Factory */

std::shared_ptr<HessianEvaluation> HessianEvaluationFactory::create(std::string hessian_evaluation_method, int dimension) {
    if (hessian_evaluation_method == "exact") {
        return std::make_shared<ExactHessianEvaluation>(dimension);
    }
//    else if (hessian_evaluation_method == "BFGS") {
//        return std::make_shared<BFGSHessianEvaluation>(dimension);
//    }
    else {
        throw std::invalid_argument("Hessian evaluation method " + hessian_evaluation_method + " does not exist");
    }
}