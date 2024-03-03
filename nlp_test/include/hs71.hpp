#ifndef NLPTEST_HS71_HPP
#define NLPTEST_HS71_HPP

#include <vector>
#include <limits>
#include "optimization/Model.hpp"
#include "linear_algebra/RectangularMatrix.hpp"

class HS71
    : public Model 
{
private:
    std::vector<double> x_l {1.0, 1.0, 1.0, 1.0};
    std::vector<double> x_u {5.0, 5.0, 5.0, 5.0};
    std::vector<double> g_l {25.0, 40.0};
    std::vector<double> g_u {std::numeric_limits<double>::infinity(), 40.0};
    std::vector<double> x0 {3., 3., 3., 3.};
public:
    HS71() 
        : Model("HS71", 4, 2)
    {
        equality_constraints.resize(1);
        equality_constraints[0] = 1;
        
        inequality_constraints.resize(1);
        inequality_constraints[0] = 0;
        
        lower_bounded_variables.resize(x0.size());
        
        for (std::size_t i = 0; i < lower_bounded_variables.size(); i++)
            lower_bounded_variables[i] = i;
        
        upper_bounded_variables.resize(x0.size());
        for (std::size_t i = 0; i < upper_bounded_variables.size(); i++)
            upper_bounded_variables[i] = i;

        single_lower_bounded_variables.resize(0);
        single_upper_bounded_variables.resize(0);

        number_objective_gradient_nonzeros = get_number_objective_gradient_nonzeros();
        number_jacobian_nonzeros = get_number_jacobian_nonzeros();
        number_hessian_nonzeros = get_number_hessian_nonzeros();
    }

   double get_variable_lower_bound(size_t variable_index) const override {
        return x_l[variable_index];
   }
    double get_variable_upper_bound(size_t variable_index) const override {
          return x_u[variable_index];
    }
    double get_constraint_lower_bound(size_t constraint_index) const override {
          return g_l[constraint_index];
    }
    double get_constraint_upper_bound(size_t constraint_index) const override {
          return g_u[constraint_index];
    }
    BoundType get_variable_bound_type(size_t variable_index) const override {
        if (x_l[variable_index] == x_u[variable_index]) {
            return BoundType::EQUAL_BOUNDS;
        }
        else if (std::isfinite(x_l[variable_index]) && std::isfinite(x_u[variable_index])) {
            return BoundType::BOUNDED_BOTH_SIDES;
        }
        else if (std::isfinite(x_l[variable_index])) {
            return BoundType::BOUNDED_LOWER;
        }
        else if (std::isfinite(x_u[variable_index])) {
            return BoundType::BOUNDED_UPPER;
        }
        else {
            return BoundType::UNBOUNDED;
        }
    }
    FunctionType get_constraint_type(size_t constraint_index) const override {
        return FunctionType::NONLINEAR;
    }

    BoundType get_constraint_bound_type(size_t constraint_index) const override {
        if (g_l[constraint_index] == g_u[constraint_index]) {
            return BoundType::EQUAL_BOUNDS;
        }
        else if (std::isfinite(g_l[constraint_index]) && std::isfinite(g_u[constraint_index])) {
            return BoundType::BOUNDED_BOTH_SIDES;
        }
        else if (std::isfinite(g_l[constraint_index])) {
            return BoundType::BOUNDED_LOWER;
        }
        else if (std::isfinite(g_u[constraint_index])) {
            return BoundType::BOUNDED_UPPER;
        }
        else {
            return BoundType::UNBOUNDED;
        }
    }
    size_t get_number_objective_gradient_nonzeros() const override {
        return x0.size();
    }

    size_t get_number_jacobian_nonzeros() const override {
        return 2*x0.size();
    }

    size_t get_number_hessian_nonzeros() const override {
        return x0.size()*(x0.size()+1)/2;
    }
double evaluate_objective(const std::vector<double>& x) const override {
    return x[0]*x[3]*(x[0]+x[1]+x[2])+x[2];
}
void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override {
    gradient.insert(0, x[0]*x[3] + x[3]*(x[0]+x[1]+x[2]));
    gradient.insert(1, x[0]*x[3]);
    gradient.insert(2, x[0]*x[3] + 1);
    gradient.insert(3, x[0]*(x[0]+x[1]+x[2]));
}
void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override {
    constraints[0] = x[0]*x[1]*x[2]*x[3];
    constraints[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
}
void evaluate_constraint_gradient(const std::vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override {
    if (constraint_index == 0) {
        gradient.insert(0, x[1]*x[2]*x[3]);
        gradient.insert(1, x[0]*x[2]*x[3]);
        gradient.insert(2, x[0]*x[1]*x[3]);
        gradient.insert(3, x[0]*x[1]*x[2]);
    }
    else if (constraint_index == 1) {
        gradient.insert(0, 2*x[0]);
        gradient.insert(1, 2*x[1]);
        gradient.insert(2, 2*x[2]);
        gradient.insert(3, 2*x[3]);
    }
}
void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override {
    evaluate_constraint_gradient(x, 0, constraint_jacobian[0]);
    evaluate_constraint_gradient(x, 1, constraint_jacobian[1]);
}
void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, 
    const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const override {
    double H[4][4] = {{2*x[3]           ,0    ,0    ,0},
                      {x[3]             ,0    ,0    ,0},
                      {x[3]             ,0    ,0    ,0},
                      {2*x[0]+x[1]+x[2] ,x[0] ,x[0] ,0}};

    double cH[4][4] = {{0,0,0,0},
                       {x[2]*x[3],0,0,0},
                       {x[1]*x[3],x[0]*x[3],0,0},
                       {x[1]*x[2],x[0]*x[2],x[0]*x[1],0}};
    const double cH2[4][4] = {{2,0,0,0},
                              {0,2,0,0},
                              {0,0,2,0},
                              {0,0,0,2}};
    for (size_t iRow = 0; iRow < 4; iRow++) {
        for (size_t jCol = 0; jCol <= iRow; jCol++) {
            hessian.insert(objective_multiplier*H[iRow][jCol] + multipliers[0]*cH[iRow][jCol] + multipliers[1]*cH2[iRow][jCol], iRow, jCol);
        }
    }
    }
    void get_initial_primal_point(std::vector<double>& x) const override {
        std::copy(x0.cbegin(), x0.cend(), x.begin());
    }
    void get_initial_dual_point(std::vector<double>& multipliers) const override {
        std::fill(multipliers.begin(), multipliers.end(), 0.0);
    }
    void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override {
        // do nothing
    }
    const std::vector<size_t>& get_linear_constraints() const override {
        return equality_constraints;
    }
}; // class HS71


#endif // NLPTEST_HS71_HPP