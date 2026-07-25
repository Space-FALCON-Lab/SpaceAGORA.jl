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
    using SpecialFunctions: erfc, erfcx
    using StaticArrays
    import ..DynamicEffectors: calcForceTorque, wrench, environment_requirements, solver_partition

    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
    export AerodynamicSurfaceModel
    export modified_newtonian_cp_max
    export regular_newtonian_sphere_cone_cn_ca, modified_newtonian_sphere_cone_cn_ca
    export sphere_cone_cn_ca, cn_ca_to_cl_cd, sphere_cone_cl_cd
    export NewtonianSurfaceQuadrature, NewtonianAerodynamicGeometry, NewtonianAerodynamicCoefficients
    export newtonian_surface_quadrature, combine_newtonian_surfaces
    export sphere_cone_newtonian_geometry, newtonian_aerodynamic_coefficients
    export newtonian_stability_derivatives
    export AerodynamicSurfaceQuadrature, AerodynamicGeometry, SurfaceAerodynamicCoefficients
    export aerodynamic_surface_quadrature, aerodynamic_plate_surface, combine_aerodynamic_surfaces
    export sphere_cone_aerodynamic_geometry
    export free_molecular_surface_coefficients, free_molecular_aerodynamic_coefficients
    export gas_mean_free_path, gas_knudsen_number, transitional_free_molecular_weight
    export AerodynamicRegimeResult, aerodynamic_regime_coefficients

    include(joinpath(@__DIR__, "..", "..", "..", "core", "interfaces", "reference_system.jl"))
    include(joinpath(@__DIR__, "newtonian_aerodynamics.jl"))
    include(joinpath(@__DIR__, "rarefied_surface_aerodynamics.jl"))
    include(joinpath(@__DIR__, "..", "aerodynamic_wrench_models.jl"))
end
