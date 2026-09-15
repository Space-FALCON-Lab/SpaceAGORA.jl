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
    import ..DynamicEffectors: calcForceTorque, wrench, wrench_caching!, environment_requirements, solver_partition

    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
    export MeshAeroPanels, MeshAeroSurrogate, AerodynamicCoefficientMeshSurrogate, mesh_aero_panels, panel_aero_coefficients, panel_aero_coefficients_split, panel_shadow_mask, panel_projected_area, fit_mesh_aero_surrogate, mesh_aero_coefficients, write_mesh_aero_surrogate, read_mesh_aero_surrogate, fibonacci_directions

    include(joinpath(@__DIR__, "..", "..", "..", "core", "interfaces", "reference_system.jl"))
    include(joinpath(@__DIR__, "..", "aerodynamic_wrench_models.jl"))
    include(joinpath(@__DIR__, "..", "aerodynamic_mesh_surrogate.jl"))
end
