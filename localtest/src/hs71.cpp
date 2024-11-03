#include <iostream>
#include <string>
#include <stdexcept>
#include "ingredients/globalization_mechanisms/GlobalizationMechanism.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "hs71.hpp"
#include "Uno.hpp"
#include "model/ModelFactory.hpp"
#include "options/Options.hpp"
#include "options/DefaultOptions.hpp"
#include "tools/Logger.hpp"



int main() {

   uno::Options options = uno::DefaultOptions::load();
   uno::Options solvers_options = uno::DefaultOptions::determine_solvers_and_preset();
   // uno::Options::set_preset(solvers_options, "ipopt");
   options.overwrite_with(solvers_options);

   std::unique_ptr<uno::Model> hs_model = std::make_unique<local::HS71>();
   std::unique_ptr<uno::Model> model = uno::ModelFactory::reformulate(std::move(hs_model), options);
   
   uno::Iterate initial_iterate(model->number_variables, model->number_constraints);
   model->initial_primal_point(initial_iterate.primals);
   model->project_onto_variable_bounds(initial_iterate.primals);
   model->initial_dual_point(initial_iterate.multipliers.constraints);
   initial_iterate.feasibility_multipliers.reset();

   auto constraint_relaxation_strategy = uno::ConstraintRelaxationStrategyFactory::create(*model, options);
   auto globalization_mechanism = uno::GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
   uno::Uno uno = uno::Uno(*globalization_mechanism, options);

   // solve the instance
   uno.solve(*model, initial_iterate, options);
   
   return EXIT_SUCCESS;
}
