#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "l1Relaxation.hpp"

std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(const Model& model, const Options& options) {
   const std::string constraint_relaxation_type = options.at("constraint-relaxation");
   if (constraint_relaxation_type == "feasibility-restoration") {
      return std::make_unique<FeasibilityRestoration>(model, options);
   }
   else if (constraint_relaxation_type == "l1-relaxation") {
      return std::make_unique<l1Relaxation>(model, options);
   }
   throw std::invalid_argument("ConstraintRelaxationStrategy " + constraint_relaxation_type + " is not supported");
}