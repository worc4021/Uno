// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WAECHTERFILTERSTRATEGY_H
#define UNO_WAECHTERFILTERSTRATEGY_H

#include "FilterStrategy.hpp"

class WaechterFilterStrategy : public FilterStrategy {
public:
   explicit WaechterFilterStrategy(const Options& options);

   void initialize(const Iterate& first_iterate) override;
   [[nodiscard]] bool is_iterate_acceptable(const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures,
         const ProgressMeasures& predicted_reduction, double objective_multiplier) override;

protected:
   double initial_infeasibility{INF<double>};
};

#endif // UNO_WAECHTERFILTERSTRATEGY_H