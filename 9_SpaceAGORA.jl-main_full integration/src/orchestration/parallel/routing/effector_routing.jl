module EffectorRouting

using ..DynamicEffectors: NBodyGravityModel, SolarRadiationPressureModel,
    InverseSquaredGravityModel, InverseSquaredJ2GravityModel,
    GravitationalHarmonicsModel, AerodynamicCoefficientfM
import ..DynamicEffectors

# One eligibility contract for both runtime dispatch and cost accounting.
# These predicates identify existing pre-passes; they do not choose a route.
@inline batchable_effector(::Any)::Bool = false
@inline batchable_effector(::NBodyGravityModel)::Bool = true
@inline batchable_effector(::SolarRadiationPressureModel)::Bool = true
@inline batchable_effector(e::InverseSquaredGravityModel)::Bool = !e.gravity_gradient
@inline batchable_effector(e::InverseSquaredJ2GravityModel)::Bool = !e.gravity_gradient

# Harmonics has its own pre-pass using per-satellite frame/workspace state.
@inline harmonics_prepass_effector(::Any)::Bool = false
@inline harmonics_prepass_effector(::GravitationalHarmonicsModel)::Bool = true

# Per-link atmosphere stays on the queue, as in the existing engine path.
@inline aero_prepass_effector(::Any)::Bool = false
@inline aero_prepass_effector(e::AerodynamicCoefficientfM)::Bool =
    !DynamicEffectors.AerodynamicEffectors._per_link_enabled(e)

end # module EffectorRouting
