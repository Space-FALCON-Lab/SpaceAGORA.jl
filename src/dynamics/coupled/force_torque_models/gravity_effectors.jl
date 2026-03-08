module GravityEffectors
    using ...ConfigTypes: ODEParams
    using ...AbstractTypes: AbstractForceTorqueModel
    import ..DynamicEffectors: calcForceTorque

    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    export aerobraking_gravity_force_ii, j2_secular_rates

    include(joinpath(@__DIR__, "..", "..", "..", "environment", "gravity", "gravity_models.jl"))
end
