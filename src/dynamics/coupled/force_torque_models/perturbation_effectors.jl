module PerturbationEffectors
    using ...ConfigTypes: ODEParams, NBodyScratchWorkspace, HarmonicsScratchWorkspace
    using ...AbstractTypes: AbstractPlanet, AbstractForceTorqueModel
    using ...EffectorSampling: StateSample, EnvironmentSample, ThirdBodyEphemerisSample, EffectorEnvironmentRequirements
    using ...Planets: Earth, Mars, Venus, Moon, Titan
    using ...EphemeridesModels: spice_position_j2000_m
    using ...ParallelPolicy
    using ....RuntimeServices: SPICE_LOCK
    using ...SimulationModel: SRPSunEphemerisCache, NBodyEphemerisCache, SpiceRhsMemo
    using ...SimulationModel: ephemerides_cache_key, ephemerides_requires_spice, planet_frame_lpi
    using ..AerodynamicEffectors: _multibody_thread_decision
    # Geodetic conversion for the IGRF field option; AerodynamicEffectors owns
    # a reference_system.jl include, so reuse its binding rather than adding
    # another include copy (the module-replacement bug class).
    using ..AerodynamicEffectors: rtolatlong
    using StaticArrays
    import ..DynamicEffectors: calcForceTorque, wrench, environment_requirements
    import ...EffectorSampling: gravity_backbone_structure, gravity_backbone_acceleration_ii
    import ...EffectorSampling: gravity_backbone_kick_structure, gravity_backbone_kick_acceleration_ii

    export NBodyGravityModel, GravitationalHarmonicsModel, SolarRadiationPressureModel
    export srp, srp_cannonball_accel, planetary_albedo_accel, planetary_ir_accel
    export MagneticTorqueRodModel, get_magnetic_field_dipole, get_magnetic_field, calculate_magnetic_torque
    export EddyCurrentDampingModel, eddy_damping_torque
    export LVLHCascadeAttitudeControlModel

    # Read once at module load; avoids a syscall + string allocation on every harmonics call.
    const _DEBUG_COMPARE_J2 = Ref{Bool}(
        lowercase(strip(get(ENV, "SPACEAGORA_DEBUG_COMPARE_J2", "0"))) in ("1", "true", "yes", "on")
    )

    include(joinpath(@__DIR__, "..", "perturbations.jl"))
end
