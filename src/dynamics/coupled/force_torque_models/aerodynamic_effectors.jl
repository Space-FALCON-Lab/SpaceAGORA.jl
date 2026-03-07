module AerodynamicEffectors
    using ...Structure
    using ...ConfigTypes: ODEParams, AeroScratchWorkspace
    using ...AbstractTypes: AbstractForceTorqueModel
    using ...ParallelPolicy
    using ...Kinematics
    using ...SimulationModel: rot
    using LinearAlgebra
    using StaticArrays
    import ..DynamicEffectors: calcForceTorque

    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight

    include(joinpath(@__DIR__, "..", "..", "..", "core", "interfaces", "reference_system.jl"))
    include(joinpath(@__DIR__, "..", "aerodynamic_wrench_models.jl"))
end
