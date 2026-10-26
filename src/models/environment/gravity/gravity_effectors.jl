# Compatibility facade only. Physical field kernels are in this directory;
# spacecraft wrench adapters are owned by DynamicEffectors.GravityEffectors.
module GravityEffectors

using ..DynamicEffectors.GravityEffectors: ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
using ..DynamicEffectors.GravityEffectors: GravityGradientTorqueModel
using ..DynamicEffectors.GravityEffectors: aerobraking_gravity_force_ii, j2_secular_rates

export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
export GravityGradientTorqueModel
export aerobraking_gravity_force_ii, j2_secular_rates

    using ..DynamicEffectors.GravityEffectors: _inverse_squared_gravity_accel, _inverse_squared_j2_gravity_accel
    using ..DynamicEffectors.PerturbationEffectors: NBodyGravityModel, nbody_acceleration_ii_at_epoch
end # module GravityEffectors
