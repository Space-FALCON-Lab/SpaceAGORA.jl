module AerodynamicEffectors
    using ...Structure
    using ...ConfigTypes: ODEParams, AeroScratchWorkspace
    using ...AbstractTypes: AbstractForceTorqueModel
    using ...EffectorSampling: StateSample, EnvironmentSample, EffectorEnvironmentRequirements
    using ...ParallelPolicy
    using ...Kinematics
    import ...SimulationModel
    using ...SimulationModel: rot
    using LinearAlgebra
    using StaticArrays
    import ..DynamicEffectors: calcForceTorque, wrench, environment_requirements

    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight

    include(joinpath(@__DIR__, "..", "..", "..", "core", "interfaces", "reference_system.jl"))
    include(joinpath(@__DIR__, "..", "aerodynamic_wrench_models.jl"))
end
