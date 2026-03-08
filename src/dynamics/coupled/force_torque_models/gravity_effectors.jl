module GravityEffectors
    using ...ConfigTypes: ODEParams
    using ...AbstractTypes: AbstractForceTorqueModel
    using ...EffectorSampling: StateSample, EnvironmentSample, EffectorEnvironmentRequirements
    using ...SimulationModel: rot
    import ..DynamicEffectors: calcForceTorque, wrench, environment_requirements

    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    export aerobraking_gravity_force_ii, j2_secular_rates

    include(joinpath(@__DIR__, "..", "..", "..", "environment", "gravity", "gravity_models.jl"))
end
