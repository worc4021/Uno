#ifndef PENALTYSTRATEGY_H
#define PENALTYSTRATEGY_H

#include <vector>
#include "GlobalizationStrategy.hpp"
#include "Constraint.hpp"

/*! \class PenaltyStrategy
 * \brief Step acceptance strategy based on a penalty method
 *
 *  Strategy that accepts or declines a trial step
 */
class PenaltyStrategy : public GlobalizationStrategy {
    public:
        /*!
         *  Constructor that takes an optimization problem and a set of constants
         */
        PenaltyStrategy(Subproblem& subproblem, double tolerance);

        LocalSolution compute_step(Problem& problem, Iterate& current_iterate, double radius) override;

        /*!
         *  Check the validity of a step
         *  Implements the purely virtual method of the superclass
         */
        bool check_step(Problem& problem, Iterate& current_iterate, LocalSolution& solution, double step_length) override;
        Iterate initialize(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers, bool use_trust_region) override;
        void compute_measures(Problem& problem, Iterate& iterate) override;

        double penalty_parameter; /*!< Penalty */

    private:
        PenaltyDimensions penalty_dimensions;
        double tau;
        double eta;
        double epsilon1;
        double epsilon2;

        double compute_linear_model(Problem& problem, LocalSolution& solution);
        std::vector<double> compute_bound_multipliers(Problem& problem, LocalSolution& solution);
        std::vector<double> compute_constraint_multipliers(Problem& problem, LocalSolution& solution);
        double compute_error(Problem& problem, Iterate& current_iterate, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers, double penalty_parameter);
        OptimalityStatus compute_status(Problem& problem, Iterate& current_iterate, double step_norm);
};

#endif // PENALTYSTRATEGY_H