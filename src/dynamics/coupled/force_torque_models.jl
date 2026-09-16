"""
    Wrapper module for all dynamic effector models (all forces/torques)
"""
module DynamicEffectors
    function calcForceTorque end
    using ..EffectorSampling: wrench, wrench_caching!, environment_requirements, solver_partition

    include(joinpath(@__DIR__, "force_torque_models", "gravity_effectors.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "aerodynamic_effectors.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "perturbation_effectors.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "thruster_models.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "guidance_models.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "robot_arm_reaction_effector.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "plume_surface_interaction.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "regolith_erosion.jl"))
    include(joinpath(@__DIR__, "force_torque_models", "ejecta_transport.jl"))

    # Gravity helpers still rely on perturbation calculations for legacy aerobraking paths.
    @eval GravityEffectors using ..PerturbationEffectors

    using .GravityEffectors: ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    using .GravityEffectors: aerobraking_gravity_force_ii
    using .AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
    using .AerodynamicEffectors: MeshAeroPanels, MeshAeroSurrogate, AerodynamicCoefficientMeshSurrogate, mesh_aero_panels, panel_aero_coefficients, panel_aero_coefficients_split, panel_shadow_mask, panel_projected_area, fit_mesh_aero_surrogate, mesh_aero_coefficients, write_mesh_aero_surrogate, read_mesh_aero_surrogate, fibonacci_directions
    using .AerodynamicEffectors: _parse_bool_env, _multibody_outer_parallel_hint, collect_and_reset_link_wrenches!
    using .AerodynamicEffectors: _multibody_parallel_mode, _multibody_thread_threshold, _multibody_max_threads, refresh_multibody_parallel_mode!
    using .AerodynamicEffectors: _threadid_capacity, _multibody_use_threads, _multibody_thread_decision
    using .AerodynamicEffectors: _make_aero_scratch_workspace, _ensure_aero_workspace_capacity!, _aero_workspace_for_sat!
    using .PerturbationEffectors: NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
    using .PerturbationEffectors: MagneticTorqueRodModel, get_magnetic_field_dipole, calculate_magnetic_torque
    using .PerturbationEffectors: EddyCurrentDampingModel, eddy_damping_torque
    using .PerturbationEffectors: LVLHCascadeAttitudeControlModel
    using .PerturbationEffectors: srp, srp_cannonball_accel, _spice_query_name
    using .PerturbationEffectors: planetary_albedo_accel, planetary_ir_accel
    using .PerturbationEffectors: _make_nbody_scratch_workspace, _ensure_nbody_workspace_capacity!, _nbody_workspace_for_sat!
    using .PerturbationEffectors: _make_harmonics_scratch_workspace, _harmonics_workspace_for_sat!
    using .PerturbationEffectors: _nbody_body_position_from_cache_j2000_m, _srp_sun_position_from_cache_j2000_m
    using .PerturbationEffectors: eclipse_area_calc
    using .ThrusterModels: BaseThrusterModel
    using .GuidanceModels: AerobrakingCampaignPropulsiveManeuverGuidanceModel
    using .RobotArmReactionEffectors: RobotArmReactionEffector
    using .PlumeSurfaceInteraction: PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState
    using .PlumeSurfaceInteraction: plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities, plume_engine_axis
    using .RegolithErosion: RegolithProperties, lunar_mare_regolith, ErosionEnvironment, erosion_environment
    using .RegolithErosion: ErosionRegimeKind, NoErosion, ViscousErosion, DiffusionDrivenFlowRegime, BearingCapacityFailureRegime
    using .RegolithErosion: AbstractErosionRegime, ViscousErosionRoberts, ViscousErosionEnergyFlux
    using .RegolithErosion: DiffusionDrivenFlow, BearingCapacityFailure
    using .RegolithErosion: erosion_rate, erosion_onset, regime_kind, default_erosion_regimes, regolith_erosion_rate
    using .RegolithErosion: shields_threshold_shear_pa, energy_flux_threshold_shear_pa, soil_bearing_capacity_pa
    using .RegolithErosion: mean_thermal_speed_mps, mean_lift_height_m, gas_dynamic_viscosity_pa_s
    using .RegolithErosion: pressure_diffusion_depth_m, soil_tensile_strength_pa
    using .EjectaTransport: EjectaSoil, EjectaTransportConfig, EjectaReferenceGasField
    using .EjectaTransport: ejecta_gas_state, ejecta_drag_coefficient, ejecta_flow_regime
    using .EjectaTransport: ejecta_launch_speed, ejecta_trajectory, ejecta_distribution, ejecta_escape_speed

    export ConstantGravityModel, InverseSquaredGravityModel, InverseSquaredJ2GravityModel
    export NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
    export aerobraking_gravity_force_ii, srp, srp_cannonball_accel, planetary_albedo_accel, planetary_ir_accel
    export MagneticTorqueRodModel, get_magnetic_field_dipole, calculate_magnetic_torque
    export EddyCurrentDampingModel, eddy_damping_torque
    export LVLHCascadeAttitudeControlModel
    export AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
    export MeshAeroPanels, MeshAeroSurrogate, AerodynamicCoefficientMeshSurrogate, mesh_aero_panels, panel_aero_coefficients, panel_aero_coefficients_split, panel_shadow_mask, panel_projected_area, fit_mesh_aero_surrogate, mesh_aero_coefficients, write_mesh_aero_surrogate, read_mesh_aero_surrogate, fibonacci_directions
    export calcForceTorque
    export wrench, wrench_caching!, environment_requirements, solver_partition
    export BaseThrusterModel
    export AerobrakingCampaignPropulsiveManeuverGuidanceModel
    export RobotArmReactionEffector
    export PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState
    export plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities
    export RegolithProperties, lunar_mare_regolith, ErosionEnvironment, erosion_environment
    export ErosionRegimeKind, NoErosion, ViscousErosion, DiffusionDrivenFlowRegime, BearingCapacityFailureRegime
    export AbstractErosionRegime, ViscousErosionRoberts, ViscousErosionEnergyFlux
    export DiffusionDrivenFlow, BearingCapacityFailure
    export erosion_rate, erosion_onset, regime_kind, default_erosion_regimes, regolith_erosion_rate
    export shields_threshold_shear_pa, energy_flux_threshold_shear_pa, soil_bearing_capacity_pa
    export mean_thermal_speed_mps, mean_lift_height_m, gas_dynamic_viscosity_pa_s
    export pressure_diffusion_depth_m, soil_tensile_strength_pa
    export EjectaSoil, EjectaTransportConfig, EjectaReferenceGasField
    export ejecta_gas_state, ejecta_drag_coefficient, ejecta_flow_regime
    export ejecta_launch_speed, ejecta_trajectory, ejecta_distribution, ejecta_escape_speed
end
