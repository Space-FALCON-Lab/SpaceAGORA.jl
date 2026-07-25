"""
    Wrapper module for all dynamic effector models (all forces/torques)
"""
module DynamicEffectors
    function calcForceTorque end
    using ..EffectorSampling: wrench, environment_requirements, solver_partition

    include(joinpath(@__DIR__, "force_torque_models", "gravity_effectors.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "aerodynamic_effectors.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "perturbation_effectors.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "thruster_models.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "guidance_models.jl"))

    # Gravity helpers still rely on perturbation calculations for legacy aerobraking paths.
    @eval GravityEffectors using ..PerturbationEffectors

    using .GravityEffectors: ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    using .GravityEffectors: aerobraking_gravity_force_ii
    using .AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
    using .AerodynamicEffectors: AerodynamicSurfaceModel
    using .AerodynamicEffectors: modified_newtonian_cp_max
    using .AerodynamicEffectors: regular_newtonian_sphere_cone_cn_ca, modified_newtonian_sphere_cone_cn_ca
    using .AerodynamicEffectors: sphere_cone_cn_ca, cn_ca_to_cl_cd, sphere_cone_cl_cd
    using .AerodynamicEffectors: NewtonianSurfaceQuadrature, NewtonianAerodynamicGeometry, NewtonianAerodynamicCoefficients
    using .AerodynamicEffectors: newtonian_surface_quadrature, combine_newtonian_surfaces
    using .AerodynamicEffectors: sphere_cone_newtonian_geometry, newtonian_aerodynamic_coefficients
    using .AerodynamicEffectors: newtonian_stability_derivatives
    using .AerodynamicEffectors: AerodynamicSurfaceQuadrature, AerodynamicGeometry, SurfaceAerodynamicCoefficients
    using .AerodynamicEffectors: aerodynamic_surface_quadrature, aerodynamic_plate_surface, combine_aerodynamic_surfaces
    using .AerodynamicEffectors: sphere_cone_aerodynamic_geometry
    using .AerodynamicEffectors: free_molecular_surface_coefficients, free_molecular_aerodynamic_coefficients
    using .AerodynamicEffectors: gas_mean_free_path, gas_knudsen_number, transitional_free_molecular_weight
    using .AerodynamicEffectors: AerodynamicRegimeResult, aerodynamic_regime_coefficients
    using .AerodynamicEffectors: _parse_bool_env, _multibody_outer_parallel_hint, collect_and_reset_link_wrenches!
    using .AerodynamicEffectors: _multibody_parallel_mode, _multibody_thread_threshold, _multibody_max_threads
    using .AerodynamicEffectors: _threadid_capacity, _multibody_use_threads, _multibody_thread_decision
    using .AerodynamicEffectors: _make_aero_scratch_workspace, _ensure_aero_workspace_capacity!, _aero_workspace_for_sat!
    using .PerturbationEffectors: NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
    using .PerturbationEffectors: srp, srp_cannonball_accel, _spice_query_name
    using .PerturbationEffectors: planetary_albedo_accel, planetary_ir_accel
    using .PerturbationEffectors: _make_nbody_scratch_workspace, _ensure_nbody_workspace_capacity!, _nbody_workspace_for_sat!
    using .PerturbationEffectors: _make_harmonics_scratch_workspace, _harmonics_workspace_for_sat!
    using .PerturbationEffectors: _nbody_body_position_from_cache_j2000_m, _srp_sun_position_from_cache_j2000_m
    using .PerturbationEffectors: eclipse_area_calc
    using .ThrusterModels: BaseThrusterModel
    using .GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel

    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    export NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
    export aerobraking_gravity_force_ii, srp, srp_cannonball_accel, planetary_albedo_accel, planetary_ir_accel
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
    export calcForceTorque
    export wrench, environment_requirements, solver_partition
    export BaseThrusterModel
    export AerobrakingCampaignPropulsiveManeuverGuidanceModel
end
