
module ConfigTypes
# import .SpacecraftModel
# include("core/simulation_model.jl")
using ..SpacecraftModels: SpacecraftModel
using ..SimConfig: SimulationConfiguration
using ..EnvironmentModels: GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate
using ..LegacyModelCodes:
    LegacyGravityModelCode,
    LegacyDensityModelCode,
    LegacyAerodynamicModelCode,
    LegacyThermalModelCode,
    LegacyThrustControlCode,
    LegacyGravityConstant,
    LegacyDensityConstant,
    LegacyAeroCdClConstant,
    LegacyThermalConvectiveRadiative,
    LegacyThrustNone,
    _compat_enum_parse
using ..CommandTypes: PropulsiveManeuverCommand, PropulsiveBurnPlan
using ..EffectorSampling: StateSample
using StaticArrays
using AstroTime
using OrdinaryDiffEq
using Reexport

export Initial_condition, Aerodynamics, Engines, Model, Cnf, Solution, ODEParams, IntermediateSolution, Mission, InitialParameters
export SaveCache, SaveData
export SRPSunEphemerisCache, NBodyEphemerisCache, PlanetFrameEphemerisCache, SpiceRuntimeCounters, SpiceRhsMemo
export GramTrackCache, VacuumPredictedGRAMCache, AeroScratchWorkspace, NBodyScratchWorkspace, HarmonicsScratchWorkspace
export PolicyDecisionEnvConfig, GramTrackCacheConfig, CallbackEnvConfig, RhsPlanEnvConfig
export RhsEffectorDecision, RhsExecutionPlan
    @kwdef struct Mission
        e::Int64 = 0
        d::Int64 = 0
        l::Int64 = 0
        a::Int64 = 0
        planet::Int64 = 0
    end

    @kwdef struct InitialParameters
        M::Mission = Mission()
        gm::LegacyGravityModelCode = LegacyGravityConstant
        dm::LegacyDensityModelCode = LegacyDensityConstant
        wm::Int64 = 0
        am::LegacyAerodynamicModelCode = LegacyAeroCdClConstant
        tm::LegacyThermalModelCode = LegacyThermalConvectiveRadiative
        cm::Int64 = 0
        tc::LegacyThrustControlCode = LegacyThrustNone
        mc::Int64 = 0
    end

    # Backward-compatible constructor for legacy integer-coded callers.
    function InitialParameters(
        M::Mission,
        gm::Union{LegacyGravityModelCode, Integer},
        dm::Union{LegacyDensityModelCode, Integer},
        wm::Integer,
        am::Union{LegacyAerodynamicModelCode, Integer},
        tm::Union{LegacyThermalModelCode, Integer},
        cm::Integer,
        tc::Union{LegacyThrustControlCode, Integer},
        mc::Integer
    )
        return InitialParameters(
            M=M,
            gm=_compat_enum_parse(LegacyGravityModelCode, gm),
            dm=_compat_enum_parse(LegacyDensityModelCode, dm),
            wm=Int64(wm),
            am=_compat_enum_parse(LegacyAerodynamicModelCode, am),
            tm=_compat_enum_parse(LegacyThermalModelCode, tm),
            cm=Int64(cm),
            tc=_compat_enum_parse(LegacyThrustControlCode, tc),
            mc=Int64(mc)
        )
    end

    # @kwdef mutable struct Planet
    #     Rp_e::Float64 = 0.0
    #     Rp_p::Float64 = 0.0
    #     Rp_m::Float64 = 0.0
    #     mass::Float64  = 0.0
    #     p::Float64 = 0.0
    #     k::Float64 = 0.0
    #     ω::Vector{Float64} = [0.0, 0.0, 0.0]
    #     g_ref::Float64 = 0.0
    #     ρ_ref::Float64 = 0.0
    #     h_ref::Float64 = 0.0
    #     H::Float64 = 0.0
    #     R::Float64 = 0.0
    #     γ::Float64 = 0.0
    #     T::Float64 = 0.0
    #     J2::Float64 = 0.0
    #     μ::Float64 = 0.0
    #     μ_fluid::Float64 = 0.0
    #     Lz::Float64 = 0.0
    #     α::Float64 = 0.0
    #     δ::Float64 = 0.0
    #     inertial_frame::Symbol = :J2000
    #     L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
    #     Clm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     Slm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     Clm_topo::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     Slm_topo::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     name::String = ""
    #     A_grav::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     A_topo::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     VR01::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     VR11::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     N1::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     N2::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     A::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    #     Re::Vector{Float64} = [0.0, 0.0, 0.0]
    #     Im::Vector{Float64} = [0.0, 0.0, 0.0]
    #     polyfit_coeffs::Vector{Float64} = [0.0, 0.0, 0.0]
    #     topography_function::Function = (args, Clm, Slm, latitude, longitude) -> 0.0
    # end

    @kwdef mutable struct Aerodynamics
        δ::Float64 = 0.0
        α::Float64 = 0.0
        thermal_accomodation_factor::Float64 = 0.0
        reflection_coefficient::Float64 = 0.0
        thermal_contact::Float64 = 0.0
        heat_rate_limit::Float64 = 0.0
        heat_load_limit::Float64 = 0.0
    end

    # Thermal

    @kwdef mutable struct Engines
        ϕ::Float64 = 0.0
        g_e::Float64 = 0.0
        T::Float64 = 0.0
        Isp::Float64 = 0.0
    end

    @kwdef struct Initial_condition
        a::Float64 = 0.0
        e::Float64 = 0.0
        i::Float64 = 0.0
        Ω::Float64 = 0.0
        ω::Float64 = 0.0
        vi::Float64 = 0.0
        m::Float64 = 0.0
        year::Int64 = 0
        month::Int64 = 0
        day::Int64 = 0
        hour::Int64 = 0
        minute::Int64 = 0
        second::Float64 = 0.0
        time_rot::Float64 = 0.0
        el_time::Float64 = 0.0 # Elapsed time in seconds
        DateTimeIC::Epoch = from_utc(2000, 1, 1, 12, 0, 0)
        DateTimeJ2000::Epoch = from_utc(2000, 1, 1, 12, 0, 0)
    end

    @kwdef struct Model 
        body::SpacecraftModel = SpacecraftModel()
        aerodynamics::Aerodynamics = Aerodynamics()
        engines::Engines = Engines()
        initial_condition::Initial_condition = Initial_condition()
    end

    # Struct to store simulation results at each time step
    @kwdef struct IntermediateSolution
        time::Float64 = 0.0
        year::Int64 = 2000
        month::Int64 = 1
        day::Int64 = 1
        hour::Int64 = 12
        minute::Int64 = 0
        second::Float64 = 0.0
        number_of_passage::Int64 = 0
        pos_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0) # Inertial position vector
        pos_ii_mag::Float64 = 0.0
        vel_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0) # Inertial velocity vector
        vel_ii_mag::Float64 = 0.0
        pos_pp::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0) # Planet-relative position vector
        pos_pp_mag::Float64 = 0.0
        vel_pp::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0) # Planet-relative velocity vector
        vel_pp_mag::Float64 = 0.0
        oe::SVector{6,Float64} = SVector{6,Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # Orbital elements
        lat::Float64 = 0.0
        lon::Float64 = 0.0
        alt::Float64 = 0.0
        γ_ii::Float64 = 0.0 # Inertial flight path angle
        γ_pp::Float64 = 0.0 # Planet-relative flight path angle
        h_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0) # Inertial angular momentum vector
        h_ii_mag::Float64 = 0.0
        h_pp::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0) # Planet-relative angular momentum vector
        h_pp_mag::Float64 = 0.0
        uD::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        uE::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        uN::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        vN::Float64 = 0.0
        vE::Float64 = 0.0
        azi_pp::Float64 = 0.0
        ρ::Float64 = 0.0
        T::Float64 = 0.0
        p::Float64 = 0.0
        wind::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        cL::Float64 = 0.0
        cD::Float64 = 0.0
        S::Float64 = 0.0 # Molecular speed ratio
        mass::Float64 = 0.0
        T_r::Float64 = 0.0 # Not sure what this is, seem to only be set to 0 in complete passage
        dynamic_pressure::Float64 = 0.0
        gravity_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        drag_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        drag_pp::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        lift_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        lift_pp::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        force_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        τ_body::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        energy::Float64 = 0.0
        MC_index::Int64 = 0
        drag_state::Int64 = 0
        quaternion::SVector{4,Float64} = SVector{4,Float64}(0.0, 0.0, 0.0, 0.0)
        ω::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        α_control::Float64 = 0.0 # Angle of attack control, should probably move to another struct for control results
        inertia_tensor::SVector{9, Float64} = SVector{9, Float64}(zeros(9)) # inertia tensor components
        τ_rw::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0) # total torque applied by all reaction wheels, move to another struct for control results
        α::Vector{Float64} = zeros(3) # Angle of attack
        β::Vector{Float64} = zeros(3) # Sideslip angle
        heat_rate::Vector{Float64} = [0.0, 0.0, 0.0]
        heat_load::Vector{Float64} = [0.0, 0.0, 0.0]
        rw_h::Vector{Float64} = [0.0, 0.0, 0.0] # angular momentum magnitudes of each reaction wheel, move to another struct for control results
        rw_τ::Vector{Float64} = [0.0, 0.0, 0.0] # torque magnitude applied by each reaction wheel, move to another struct for control results
        thruster_forces::Vector{Float64} = [0.0, 0.0, 0.0, 0.0] # forces applied by each thruster, move to another struct for control results
    end

    @kwdef mutable struct Cnf
        impact::Bool = false
        altitude_periapsis::Vector{Float64} = []
        latitude_periapsis::Vector{Float64} = []
        longitude_periapsis::Vector{Float64} = []
        max_heatrate::Vector{Float64} = []
        solution_intermediate::Vector{IntermediateSolution} = []
        atmospheric_data::Dict{String,Float64} = Dict()
        previous_atmospheric_data::Dict{String,Float64} = Dict()
        drag_state::Bool = false
        ascending_phase::Bool = false
        evaluate_switch_heat_load::Bool = false
        security_mode::Bool = false
        time_IEI::Float64 = 0.0
        time_OEI::Float64 = 0.0
        time_switch_1::Float64 = 0.0
        time_switch_2::Float64 = 0.0
        state_inner_boundary_atmosphere::Vector{Float64} = []
        count_aerobraking::Int64 = 0
        count_dori::Int64 = 0
        count_phase::Int64 = 0
        count_numberofpassage::Int64 = 0
        count_overcome_hr::Int64 = 0
        counter_random::Int64 = 0
        save_index_heat::Int64 = 0
        index_warning_alt::Int64 = 0
        index_warning_flow::Int64 = 0
        index_Mars_Gram_call::Int64 = 0
        index_MonteCarlo::Int64 = 1
        index_propellant_mass::Int64 = 1
        T_w::Float64 = 4.0
        Δv_man::Float64 = 0.0
        closed_form_solution_off::Int64 = 1
        α::Float64 = pi/2
        α_past::Float64 = pi/2
        raise_man_orbit::Vector{Float64} = []
        lower_man_orbit::Vector{Float64} = []
        et::Float64 = 0.0

        # Results to delete
        periapsis_list::Vector{Float64} = []
        Δv_list::Vector{Float64} = []
        orbit_number_list::Vector{Int64} = []
        heat_load_past::Vector{Float64} = []
        heat_load_ppast::Vector{Float64} = []
        state_flesh1::Vector{Vector{Float64}} = [[]]
        α_list::Vector{Float64} = []
        initial_position_closed_form::SVector{7,Float64} = SVector{7,Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        continue_simulation::Bool = true
        timer_revaluation::Float64 = 0.0
        MarsGram_recall::Int64 = 0
        heat_rate_prev::Vector{Float64} = []
        sensible_loads::Bool = false
        counter_integrator::Int64 = 0
        prev_step_integrator::Float64 = 0.0
        initial_time_saved::Float64 = 0.0
        prev_timestep::Float64 = 0.0
        ω_wheel_derivatives::Vector{Vector{Float64}} = []

        # Extra variables missing in python version
        counter::Int64 = 0
        heat_rate_limit::Float64 = 0.0
        time_OP::Float64 = 0.0
        time_IP::Float64 = 0.0
        Gram_justrecalled::Int64 = 0
        Gram_directory::String = ""
        heat_rate_list::Vector{Float64} = []
        stop_simulation::Bool = false
        results_save::Int64 = 0
        count_eventfirststep::Int64 = 0
        eventfirststep_periapsis::Int64 = 0
        count_eventsecondstep::Int64 = 0
        count_reached_EI::Int64 = 0
        count_reached_AE::Int64 = 0
        count_out_drag_passage::Int64 = 0
        count_in_drag_passage::Int64 = 0
        count_in_drag_passage_nt::Int64 = 0
        count_apoapsispoint::Int64 = 0
        count_periapsispoint::Int64 = 0
        count_impact::Int64 = 0
        count_apoapsisgreaterperiapsis::Int64 = 0
        count_stop_firing::Int64 = 0
        count_guidance::Int64 = 0
        count_heat_rate_check::Int64 = 0
        count_heat_load_check_exit::Int64 = 0
        count_final_entry_altitude_reached::Int64 = 0
        time_termination::Bool = false
        t_out_drag_passage::Float64 = 0.0
        t_time_switch_func::Vector{Float64} = []
        t_time_switch_targ::Vector{Float64} = []
        ts_targ_1::Float64 = 0.0
        ts_targ_2::Float64 = 0.0
        
        prob::ODEProblem = ODEProblem((u, p, t) -> u, [0.0], (0.0, 1.0))
        prob_set::Bool = false
        P::Matrix{Float64} = zeros(3,3)

        DU::Float64 = 0.0
        TU::Float64 = 0.0
        MU::Float64 = 0.0   

        targeting::Int64 = 0
        Vf::Float64 = 0.0
        hf::Float64 = 0.0
        γf::Float64 = 0.0

        lambda_switch_list::Vector{Float64} = []
        time_switch_list::Vector{Float64} = []
        time_list::Vector{Float64} = []
        lamv_list::Vector{Float64} = []

        t_switch_targeting::Float64 = 0.0

        drag_pp::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        lift_pp::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        drag_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        lift_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        gravity_cent_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        gravity_nbody_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)
        gravity_harmonics_ii::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)

        CL_current::Float64 = 0.0
        CD_current::Float64 = 0.0

        vel_pp_rw::SVector{3,Float64} = SVector{3,Float64}(0.0, 0.0, 0.0)

        β_body::Vector{Float64} = []
        α_body::Vector{Float64} = []

        T_p::Float64 = 0.0
        S::Float64 = 0.0
        q::Float64 = 0.0
        rot_body_to_inertial::Matrix{Float64} = zeros(3,3)
    end

    # cnf = Cnf()

    @kwdef mutable struct Controller
        guidance_t_eval::Vector{Float64} = []
        count_controller::Int64 = 1
        count_prev_controller::Int64 = 0
        stored_state::Int64 = 1
        prev_time::Float64 = 0.0
        t::Float64 = 0.0
    end

    # controller = Controller()

    @kwdef mutable struct Orientation
        time::Vector{Float64} = []
        year::Vector{Float64} = []
        month::Vector{Float64} = []
        day::Vector{Float64} = []
        hour::Vector{Float64} = []
        minute::Vector{Float64} = []
        second::Vector{Float64} = []
        number_of_passage::Vector{Int64} = []
        pos_ii::Vector{Vector{Float64}} = [[], [], []]
        vel_ii::Vector{Vector{Float64}} = [[], [], []]
        pos_ii_mag::Vector{Float64} = []
        vel_ii_mag::Vector{Float64} = []
        quaternion::Vector{Vector{Float64}} = [[], [], [], []]
        ω::Vector{Vector{Float64}} = [[], [], []]
        pos_pp::Vector{Vector{Float64}} = [[], [], []]
        pos_pp_mag::Vector{Float64} = []
        vel_pp::Vector{Vector{Float64}} = [[], [], []]
        vel_pp_mag::Vector{Float64} = []
        oe::Vector{Vector{Float64}} = [[], [], [], [], [], []]
        lat::Vector{Float64} = []
        lon::Vector{Float64} = []
        alt::Vector{Float64} = []
        γ_ii::Vector{Float64} = []
        γ_pp::Vector{Float64} = []
        h_ii::Vector{Vector{Float64}} = [[], [], []]
        h_pp::Vector{Vector{Float64}} = [[], [], []]
        h_ii_mag::Vector{Float64} = []
        h_pp_mag::Vector{Float64} = []
        uD::Vector{Vector{Float64}} = [[], [], []]
        uE::Vector{Vector{Float64}} = [[], [], []]
        uN::Vector{Vector{Float64}} = [[], [], []]
        vN::Vector{Float64} = []
        vE::Vector{Float64} = []
        azi_pp::Vector{Float64} = []
    end

    @kwdef mutable struct Physical_properties
        ρ::Vector{Float64} = []
        T::Vector{Float64} = []
        p::Vector{Float64} = []
        wind::Vector{Vector{Float64}} = [[], [], []] 
        cL::Vector{Float64} = []
        cD::Vector{Float64} = []
        α::Vector{Vector{Float64}} = []
        β::Vector{Vector{Float64}} = []
        S::Vector{Float64} = []
        inertia_tensor::Vector{Vector{Float64}} = [[], [], [], [], [], [], [], [], []] # inertia tensor components
        α_control::Vector{Float64} = []
        rw_h::Vector{Vector{Float64}} = [] # angular momentum magnitudes of each reaction wheel
        rw_τ::Vector{Vector{Float64}} = [] # torque magnitude applied by each reaction wheel
        thruster_forces::Vector{Vector{Float64}} = [] # forces applied by each thruster
        τ_rw::Vector{Vector{Float64}} = [[], [], []] # total torque applied by all reaction wheels
    end

    @kwdef mutable struct Performance
        mass::Vector{Float64} = []
        heat_rate::Vector{Vector{Float64}} = []
        heat_load::Vector{Vector{Float64}} = []
        T_r::Vector{Float64} = []
        q::Vector{Float64} = []
    end

    @kwdef mutable struct Forces
        gravity_ii::Vector{Vector{Float64}} = [[], [], []]
        drag_pp::Vector{Vector{Float64}} = [[], [], []]
        drag_ii::Vector{Vector{Float64}} = [[], [], []]
        lift_pp::Vector{Vector{Float64}} = [[], [], []]
        lift_ii::Vector{Vector{Float64}} = [[], [], []]
        force_ii::Vector{Vector{Float64}} = [[], [], []]
        τ_body::Vector{Vector{Float64}} = [[], [], []] # total torque applied by all forces
        energy::Vector{Float64} = []
    end

    @kwdef mutable struct Simulation
        MC_seed::Vector{Int64} = []
        drag_passage::Vector{Int64} = []
        solution_states::Int64 = 0
    end

    @kwdef mutable struct Closed_form
        t_cf::Vector{Float64} = []
        h_cf::Vector{Float64} = []
        γ_cf::Vector{Float64} = []
        v_cf::Vector{Float64} = []
    end

    struct SRPSunEphemerisCache
        ets::Vector{Float64}
        positions_j2000_m::Vector{SVector{3, Float64}}
    end

    struct NBodyEphemerisCache
        primary_body_name::String
        body_query_names::Vector{String}
        body_index_by_name::Dict{String, Int}
        ets::Vector{Float64}
        positions_j2000_m::Matrix{SVector{3, Float64}}
    end

    struct PlanetFrameEphemerisCache
        ets::Vector{Float64}
        quaternions::Vector{SVector{4, Float64}}
    end

    @kwdef struct SpiceRuntimeCounters
        nbody_spkpos_runtime_calls::Base.Threads.Atomic{Int64} = Base.Threads.Atomic{Int64}(0)
        nbody_spkpos_cache_build_calls::Base.Threads.Atomic{Int64} = Base.Threads.Atomic{Int64}(0)
        srp_spkpos_runtime_calls::Base.Threads.Atomic{Int64} = Base.Threads.Atomic{Int64}(0)
        srp_spkpos_cache_build_calls::Base.Threads.Atomic{Int64} = Base.Threads.Atomic{Int64}(0)
        planet_pxform_runtime_calls::Base.Threads.Atomic{Int64} = Base.Threads.Atomic{Int64}(0)
        planet_pxform_cache_build_calls::Base.Threads.Atomic{Int64} = Base.Threads.Atomic{Int64}(0)
    end

    @kwdef mutable struct SpiceRhsMemo
        lock::ReentrantLock = ReentrantLock()
        et::Float64 = NaN
        primary_body_name::String = ""
        body_positions_j2000_m::Dict{String, SVector{3, Float64}} = Dict{String, SVector{3, Float64}}()
    end

    @kwdef mutable struct Solution
        orientation::Orientation = Orientation()
        physical_properties::Physical_properties = Physical_properties()
        performance::Performance = Performance()
        forces::Forces = Forces()
        simulation::Simulation = Simulation()
        closed_form::Closed_form = Closed_form()
    end

    mutable struct GramTrackCache
        valid::Bool
        t0::Float64
        t1::Float64
        index_hint::Int
        times::Vector{Float64}
        alts::Vector{Float64}
        lats::Vector{Float64}
        lons::Vector{Float64}
        rhos::Vector{Float64}
        Ts::Vector{Float64}
        winds::Vector{SVector{3, Float64}}
    end

    # Vacuum-predicted GRAM density cache.
    # Propagates a drag-free (two-body + J2) trajectory, samples the density model
    # at N evenly-spaced knots, and fits natural cubic splines to log(ρ) and T vs
    # time.  Subsequent density queries interpolate from the splines when the actual
    # trajectory stays within a configurable position deviation of the vacuum prediction.
    mutable struct VacuumPredictedGRAMCache
        valid::Bool
        t0::Float64                         # time of first spline knot
        t1::Float64                         # time of last spline knot
        h::Float64                          # uniform knot spacing in seconds
        log_rhos::Vector{Float64}           # log(ρ) at knots (spline a-coefficients)
        Ms_rho::Vector{Float64}             # natural cubic spline second derivatives for log(ρ)
        Ts::Vector{Float64}                 # temperature at knots
        Ms_T::Vector{Float64}              # natural cubic spline second derivatives for T
        winds::Vector{SVector{3, Float64}}  # wind vectors at knots (linear interpolation)
        vac_alts::Vector{Float64}           # vacuum-predicted altitude at each knot
        vac_positions::Vector{SVector{3, Float64}} # vacuum-predicted inertial position at each knot
    end

    struct AeroScratchWorkspace
        thread_force::Vector{MVector{3, Float64}}
        thread_cl::Vector{Float64}
        thread_cd::Vector{Float64}
        thread_area::Vector{Float64}
    end

    struct NBodyScratchWorkspace
        pos_primary_k_all::Vector{SVector{3, Float64}}
        thread_force::Vector{MVector{3, Float64}}
    end

    struct HarmonicsScratchWorkspace
        A::Matrix{Float64}
        R::Vector{Float64}
        I::Vector{Float64}
    end

    # ── Run-scoped env-config snapshots ─────────────────────────────────────────
    # Typed snapshots of SPACEAGORA_* environment configuration, resolved once at
    # run_simulation setup (inside any active SimulationEngineConfig override
    # scope) and stored on SharedBuffers so RHS/callback hot paths read plain
    # struct fields instead of re-parsing ENV per evaluation.  ENV access is
    # process-global in Julia, so per-call parsing serializes concurrent
    # in-process ensemble members.  Env-override semantics: values are captured
    # at run start; changing ENV mid-run has no effect on an active run.
    # The structs are pure data; the builders live with the env parsers
    # (ParallelPolicy, SimulationCallbacks, SimulationEngine).

    # Feeds ParallelPolicy.thread_policy_decision and budget-dependent routing.
    #
    # The adaptive-controller knobs below (persistent_hints through
    # adaptive_trim_quanta) were previously read straight from ENV on every
    # adaptive decision *and* every policy observation. Each read is a Dict
    # lookup plus strip/lowercase/parse, all of which allocate, and the
    # observation path fires once per RHS call on the flat-effector route --
    # so an L50 constellation solve paid several thousand of them per second
    # for values that cannot change inside a run. They are snapshotted here
    # for the same reason inner_thread_budget already was.
    struct PolicyDecisionEnvConfig
        inner_thread_budget::Int
        outer_parallel_active::Bool
        auto_min_budget_default::Int
        auto_min_budget_density::Int
        auto_min_budget_density_lockfree::Int
        auto_min_budget_thermal::Int
        adaptive_enabled::Bool
        persistent_hints::Bool
        adaptive_measured_reward::Bool
        adaptive_bootstrap_threads::Bool
        adaptive_control_tail_guard::Bool
        adaptive_rho::Float64
        adaptive_delta::Float64
        adaptive_window::Int
        adaptive_trim_quanta::Int
    end

    # GRAM along-track density cache tolerances (built by
    # SimulationCallbacks._gram_track_cache_config from SPACEAGORA_GRAM_* vars).
    struct GramTrackCacheConfig
        mode::Symbol
        entry_horizon_s::Float64
        entry_alt_tol_m::Float64
        entry_ang_tol_rad::Float64
        entry_points::Int
        orbit_horizon_s::Float64
        orbit_alt_tol_m::Float64
        orbit_ang_tol_rad::Float64
        orbit_points::Int
        transition_band_m::Float64
    end

    # Per-invocation knobs consulted by density/control/thermal callbacks and
    # RHS-side atmosphere sampling.
    struct CallbackEnvConfig
        gram_track_cache::GramTrackCacheConfig
        gram_runtime_stats_enabled::Bool
        gram_track_cache_ignore_time_window::Bool
        gram_track_cache_target_use_j2::Bool
        # When true, sample_buffered_atmosphere reuses the density/temperature/wind
        # sampled once per accepted step (by the density DiscreteCallback) for every
        # RHS stage evaluation within that step, instead of requiring an exact time
        # match. Real (non-cached) GRAM density includes small-scale perturbation
        # noise that is not smooth in time; resampling it at every RK stage makes an
        # adaptive solver's step-size controller see that noise as stiffness and
        # collapse dt (measured: 2 sats / 1s mission -> 2.4M steps, 604s wall time).
        # Freezing density for the whole step removes that noise from the local
        # truncation-error estimate; altitude (the dominant driver of the smooth mean
        # density) changes negligibly over one integration step for a LEO orbit, so
        # this is a small, standard approximation -- the same idea the vacuum-predicted
        # cache already uses at a coarser (look-ahead-window) grain.
        density_freeze_per_step::Bool
        vacuum_gram_cache_enabled::Bool
        vacuum_gram_cache_npoints::Int
        vacuum_gram_cache_horizon_s::Float64
        vacuum_gram_cache_deviation_m::Float64
        density_parallel_mode::Symbol
        density_thread_threshold::Int
        density_allow_with_outer::Bool
        density_assume_threadsafe::Bool
        density_batch_mode::Symbol
        density_batch_threshold::Int
        gram_isolated_pool_mode::Symbol
        gram_isolated_pool_threshold::Int
        gram_isolated_pool_max_workers::Int
        control_parallel_mode::Symbol
        control_thread_threshold::Int
        control_allow_with_outer::Bool
        control_assume_threadsafe::Bool
        thermal_parallel_mode::Symbol
        thermal_thread_threshold::Int
        thermal_allow_with_outer::Bool
    end

    # Knobs consulted by the per-RHS-call execution-plan routing chain in
    # SimulationEngine (setup.jl / dynamics_rhs.jl).
    struct RhsPlanEnvConfig
        execution_mode::Symbol
        profile_forces_serial::Bool
        batch_parallel_mode::Symbol
        batch_thread_threshold::Int
        effector_parallel_mode::Symbol
        effector_thread_threshold::Int
        effector_max_threads::Int
        effector_allow_with_outer::Bool
        effector_heavy_only::Bool
        effector_cost_ns_per_item_default::Float64
        effector_cost_min_samples::Int
        effector_cost_ema_alpha::Float64
        effector_work_ns_per_worker_threshold::Float64
        effector_outer_work_scale::Float64
        flat_min_sats::Int
        flat_min_effectors::Int
        flat_work_ns_threshold::Float64
        flat_work_per_worker_ns_threshold::Float64
        flat_cost_heterogeneity_threshold::Float64
        flat_min_thread_budget::Int
        harmonics_batch_enabled::Bool
        harmonics_batch_min_sats_per_worker::Int
        harmonics_batch_spin_barrier::Bool
        harmonics_batch_allow_with_outer::Bool
        rhs_effector_cost_min_samples::Int
        flat_packet_target_min_ns::Float64
        flat_packet_scheduler_mode::Symbol
        flat_packet_min_items::Int
        flat_packet_work_ns_threshold::Float64
        flat_packet_heterogeneity_threshold::Float64
        flat_packet_overhead_disable_ratio::Float64
        flat_packet_overhead_min_samples::Int
    end

    # Concrete shape of the per-RHS-call execution plan returned by
    # SimulationEngine._rhs_execution_plan and stored in rhs_plan_override by the
    # calibration sweep.  isbits, so plans are stack-allocated in the RHS.
    const RhsEffectorDecision = @NamedTuple{use_threads::Bool, allotment::Int64, mode::Symbol, policy_applied::Bool}
    const RhsExecutionPlan = @NamedTuple{
        mode::Symbol,
        allotment::Int64,
        scheduler::Symbol,
        dominant_axis::Symbol,
        policy_applied::Bool,
        effector_decision::RhsEffectorDecision,
    }

    const _PerSatDensityModel = Union{GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate}
    const _HarmonicsWorkspaceMap = Dict{UInt, HarmonicsScratchWorkspace}

    @inline function _typed_nothing_vector(::Type{T}, n::Int) where {T}
        out = Vector{Union{Nothing, T}}(undef, n)
        fill!(out, nothing)
        return out
    end

    # A struct to hold the data shared between the callback and the integrator.
    # n_sats is a plain runtime field (not a type parameter): none of these
    # buffers are StaticArrays sized by it, so making it part of the type
    # bought no runtime performance and instead forced a fresh JIT
    # specialization of the whole RHS/effector call graph per distinct
    # satellite count -- ruinous for constellation-size sweeps (see
    # benchmarks/studies/gram_mars_fix_and_constellation_scaling).
    @kwdef struct SharedBuffers
        n_sats::Int
        densities::Vector{Float64} = zeros(Float64, n_sats)
        temperatures::Vector{Float64} = ones(Float64, n_sats)
        winds::Vector{SVector{3,Float64}} = [SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:n_sats]
        density_sample_t::Vector{Float64} = fill(NaN, n_sats)
        density_batch_altitudes::Vector{Float64} = zeros(Float64, n_sats)
        density_batch_latitudes::Vector{Float64} = zeros(Float64, n_sats)
        density_batch_longitudes::Vector{Float64} = zeros(Float64, n_sats)
        heat_rates::Vector{Vector{Float64}} = [Float64[] for _ in 1:n_sats]
        density_models::Vector{_PerSatDensityModel} = _PerSatDensityModel[]
        gram_density_cache::Vector{Union{Nothing, GramTrackCache}} = _typed_nothing_vector(GramTrackCache, n_sats)
        vacuum_gram_caches::Vector{Union{Nothing, VacuumPredictedGRAMCache}} = _typed_nothing_vector(VacuumPredictedGRAMCache, n_sats)
        gram_isolated_pool_models::Vector{GRAMAtmosphereModel} = GRAMAtmosphereModel[]
        gram_isolated_pool_locks::Vector{ReentrantLock} = ReentrantLock[]
        harmonics_workspaces::Vector{Union{Nothing, _HarmonicsWorkspaceMap}} = _typed_nothing_vector(_HarmonicsWorkspaceMap, n_sats)
        nbody_workspaces::Vector{Union{Nothing, NBodyScratchWorkspace}} = _typed_nothing_vector(NBodyScratchWorkspace, n_sats)
        aero_workspaces::Vector{Union{Nothing, AeroScratchWorkspace}} = _typed_nothing_vector(AeroScratchWorkspace, n_sats)
        nbody_ephemeris_cache::Base.RefValue{Union{Nothing, NBodyEphemerisCache}} = Ref{Union{Nothing, NBodyEphemerisCache}}(nothing)
        srp_sun_ephemeris_cache::Base.RefValue{Union{Nothing, SRPSunEphemerisCache}} = Ref{Union{Nothing, SRPSunEphemerisCache}}(nothing)
        planet_frame_ephemeris_cache::Base.RefValue{Union{Nothing, PlanetFrameEphemerisCache}} = Ref{Union{Nothing, PlanetFrameEphemerisCache}}(nothing)
        harmonics_lpi_lock::ReentrantLock = ReentrantLock()
        harmonics_lpi_key::Base.RefValue{Any} = Ref{Any}(nothing)
        harmonics_lpi::Base.RefValue{SMatrix{3,3,Float64,9}} = Ref(SMatrix{3,3,Float64,9}((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)))
        maneuver_commands::Vector{PropulsiveManeuverCommand} = [PropulsiveManeuverCommand() for _ in 1:n_sats]
        maneuver_burn_plans::Vector{PropulsiveBurnPlan} = [PropulsiveBurnPlan() for _ in 1:n_sats]
        spice_runtime_counters::SpiceRuntimeCounters = SpiceRuntimeCounters()
        spice_rhs_memo_enabled::Base.RefValue{Bool} = Ref(true)
        spice_rhs_memo::SpiceRhsMemo = SpiceRhsMemo()
        current_time::Base.RefValue{Float64} = Ref(0.0)
        et_start::Base.RefValue{Float64} = Ref(0.0)
        solve_segment_end_time::Base.RefValue{Float64} = Ref(NaN)
        debug_control::Base.RefValue{Bool} = Ref(false)
        debug_initial_derivative::Base.RefValue{Bool} = Ref(false)
        effector_cost_ns_per_item::Base.RefValue{Float64} = Ref(NaN)
        effector_cost_samples::Base.RefValue{Int64} = Ref(Int64(0))
        rhs_effector_cost_ns::Base.RefValue{Vector{Float64}} = Ref(Float64[])
        rhs_effector_cost_samples::Base.RefValue{Vector{Int64}} = Ref(Int64[])
        rhs_flat_effector_partials::Base.RefValue{Array{Float64, 3}} = Ref(Array{Float64, 3}(undef, 0, 0, 0))
        rhs_flat_effector_totals::Base.RefValue{Matrix{Float64}} = Ref(Matrix{Float64}(undef, 0, 0))
        rhs_flat_state_samples::Base.RefValue{Vector{Union{Nothing, StateSample}}} = Ref(Vector{Union{Nothing, StateSample}}())
        rhs_flat_state_pos_ii::Base.RefValue{Vector{SVector{3, Float64}}} = Ref(SVector{3, Float64}[])
        rhs_flat_state_vel_ii::Base.RefValue{Vector{SVector{3, Float64}}} = Ref(SVector{3, Float64}[])
        rhs_flat_state_mass_kg::Base.RefValue{Vector{Float64}} = Ref(Float64[])
        rhs_flat_state_q_ib::Base.RefValue{Vector{SVector{4, Float64}}} = Ref(SVector{4, Float64}[])
        rhs_flat_state_omega_body::Base.RefValue{Vector{SVector{3, Float64}}} = Ref(SVector{3, Float64}[])
        rhs_flat_planet_lpi::Base.RefValue{SMatrix{3, 3, Float64, 9}} = Ref(SMatrix{3,3,Float64,9}((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)))
        rhs_flat_planet_pos_pp::Base.RefValue{Vector{SVector{3, Float64}}} = Ref(SVector{3, Float64}[])
        rhs_flat_planet_vel_pp::Base.RefValue{Vector{SVector{3, Float64}}} = Ref(SVector{3, Float64}[])
        rhs_flat_planet_alt_m::Base.RefValue{Vector{Float64}} = Ref(Float64[])
        rhs_flat_planet_lat_rad::Base.RefValue{Vector{Float64}} = Ref(Float64[])
        rhs_flat_planet_lon_rad::Base.RefValue{Vector{Float64}} = Ref(Float64[])
        rhs_flat_solar_pos_ii::Base.RefValue{SVector{3, Float64}} = Ref(SVector{3, Float64}(0.0, 0.0, 0.0))
        rhs_flat_solar_t::Base.RefValue{Float64} = Ref(NaN)
        rhs_flat_work_items::Base.RefValue{Vector{Int}} = Ref(Int[])
        rhs_flat_packet_starts::Base.RefValue{Vector{Int}} = Ref(Int[])
        rhs_flat_packet_ends::Base.RefValue{Vector{Int}} = Ref(Int[])
        rhs_flat_packet_costs::Base.RefValue{Vector{Float64}} = Ref(Float64[])
        rhs_flat_packet_elapsed_ns::Base.RefValue{Vector{Int64}} = Ref(Int64[])
        rhs_flat_packet_overhead_ema::Base.RefValue{Float64} = Ref(NaN)
        rhs_flat_packet_overhead_samples::Base.RefValue{Int64} = Ref(Int64(0))
        rhs_flat_packet_disabled::Base.RefValue{Bool} = Ref(false)
        rhs_planet_frame_prefilled::Base.RefValue{Bool} = Ref(false)
        rhs_atmosphere_prefilled::Base.RefValue{Bool} = Ref(false)
        rhs_solar_prefilled::Base.RefValue{Bool} = Ref(false)
        rhs_harmonics_batch_pool::Base.RefValue{Any} = Ref{Any}(nothing)
        # Per-satellite atmosphere presence flag, maintained by get_drag_state_callback.
        # The timestamp is NaN until the callback has staged a value for a known
        # integrator time, so RHS code can distinguish current state from defaults.
        in_atmosphere::Vector{Bool} = fill(false, n_sats)
        in_atmosphere_sample_t::Vector{Float64} = fill(NaN, n_sats)
        # Pre-solve calibration override: when non-nothing, _rhs_execution_plan returns
        # this plan directly, bypassing all heuristic routing logic.  Set by the
        # auto-calibration sweep in rhs_calibration.jl and cleared after the solve.
        # Concretely typed (RhsExecutionPlan): a Ref{Any} here made the plan infer
        # as Any in the RHS, boxing every plan access and re-boxing ODEParams
        # (~1.6 KB inline) on each dynamics call — ~2/3 of solve-path allocations.
        rhs_plan_override::Base.RefValue{Union{Nothing, RhsExecutionPlan}} = Ref{Union{Nothing, RhsExecutionPlan}}(nothing)
        # Per-accepted-step cache of _rhs_execution_plan's routing decision, active
        # only when SPACEAGORA_RHS_PLAN_STEP_CACHE=on (_rhs_plan_step_cache_enabled,
        # setup.jl). N_sats/thread policy don't change mid-step, so re-deriving the
        # plan (active-satellite count, effector cost decision, routing heuristics)
        # on every solver STAGE call is pure waste when a solve takes several stages
        # per accepted step. Invalidated (set to `nothing`) once per accepted step by
        # update_planet_frame_callback's affect! (planet_frame.jl), which -- like this
        # cache -- already runs exactly once per accepted step via a DiscreteCallback
        # with an always-true condition (OrdinaryDiffEq only invokes DiscreteCallback
        # affect! on accepted steps, never per stage, never on rejected steps).
        # Distinct from rhs_plan_override above: that one is a permanent, solve-wide
        # pin from calibration; this one is refreshed every accepted step.
        rhs_plan_step_cache::Base.RefValue{Union{Nothing, RhsExecutionPlan}} = Ref{Union{Nothing, RhsExecutionPlan}}(nothing)
        # Run-scoped env-config snapshots (see struct docs above).  `nothing`
        # until _initialize_runtime_env_config! runs at run_simulation setup;
        # hot-path accessors fall back to live ENV parsing when unset so
        # hand-constructed ODEParams (unit tests, withenv probes) keep the
        # read-ENV-at-use behavior.
        # Lazily-computed "any effector carries a robot-arm plan" flag.
        # _apply_coupled_robot_arm_rhs! runs for every satellite on every RHS
        # evaluation, and scanning the control + dynamic effector tuples per
        # satellite is pure waste in the (overwhelmingly common) no-robot-arm
        # case. Effector tuples are fixed for a run, so no invalidation needed.
        robot_arm_present::Base.RefValue{Union{Nothing, Bool}} = Ref{Union{Nothing, Bool}}(nothing)
        policy_env_config::Base.RefValue{Union{Nothing, PolicyDecisionEnvConfig}} = Ref{Union{Nothing, PolicyDecisionEnvConfig}}(nothing)
        rhs_env_config::Base.RefValue{Union{Nothing, RhsPlanEnvConfig}} = Ref{Union{Nothing, RhsPlanEnvConfig}}(nothing)
        callback_env_config::Base.RefValue{Union{Nothing, CallbackEnvConfig}} = Ref{Union{Nothing, CallbackEnvConfig}}(nothing)
    end

    # SaveData is an output/persistence boundary and intentionally remains heterogeneous.
    const SaveData = Dict{Symbol, Any}
    @kwdef struct SaveCache
        # Define fields to store cached data for saving results during the simulation for expensive computations, e.g., drag, lift, density, etc.
        # Vector stores the values for each satellite if simulating multiple satellites in parallel, at a given time step, then passes to the SaveData struct at the end of the time step for saving results
        ρ_cache::Vector{Float64} = []
        heat_rate_cache::Vector{Float64} = []
        drag_cache::Vector{SVector{3,Float64}} = []
        lift_cache::Vector{SVector{3,Float64}} = []
        cross_cache::Vector{SVector{3,Float64}} = []
        # Add more fields as needed to store the relevant data for saving results
    end

    
    # n_sats is a plain runtime field, not a type parameter -- see the
    # SharedBuffers docs above for why. `A` remains a genuine type parameter:
    # it varies by density-model/config type, a legitimate multiple-dispatch
    # axis (unlike n_sats, nothing dispatches on a *specific* n_sats value).
    #
    # Deliberately NOT @kwdef: with only one type parameter left, a @kwdef
    # auto-generated keyword constructor (`ODEParams{A}(; ...) where A`, zero
    # positional args, A free) has the exact same callable signature as a
    # hand-written outer constructor that also leaves A free to be inferred
    # from `args` -- Julia treats those as the same method and silently
    # overwrites one, which is a hard error under module precompilation. So
    # this struct declares plain fields (no inline defaults) and gets its
    # positional inner constructor from Julia's normal default; the keyword
    # constructor below (with defaults/validation) is the only public API and
    # has a genuinely different signature (0 positional args vs. 6), so the
    # two coexist without conflict.
    struct ODEParams{A <: SimulationConfiguration}
        # m::Model = Model()                      # Model struct
        # cnf::Cnf = Cnf()            # Configuration parameters
        # solution::Solution = Solution() # Solution struct
        # index_phase_aerobraking::Float64 = 0.0  # Phase index for aerobraking control
        # ip::InitialParameters                     # Input parameters struct
        # aerobraking_phase::Int64 = 0          # Aerobraking phase
        # t_prev::Float64 = 0.0                 # Previous time for Gram calls
        # date_initial::Any          # Initial date time
        # time_0::Float64 = 0.0                  # Initial time
        # initial_state::Initial_condition         # Initial state struct
        # gram_atmosphere::Any = nothing   # GRAM atmosphere data
        # gram::Any = nothing              # GRAM object
        # numberofpassage::Int64 = 0       # Current passage number
        # orientation_sim::Bool = false    # Flag for orientation simulation
        n_sats::Int
        args::A # Arguments dictionary
        shared_buffers::SharedBuffers # Shared buffers for callback and integrator
        is_active::Vector{Bool} # Vector to track which satellites are still active in the simulation
        orbit_counter::Vector{Int64} # Counter for the number of orbits completed
        save_cache::SaveCache # Cache for saving results
    end

    function ODEParams(;
        n_sats::Int,
        args=SimulationConfiguration(),
        shared_buffers=SharedBuffers(n_sats=n_sats),
        is_active::Vector{Bool}=[true for _ in 1:n_sats],
        orbit_counter::Vector{Int64}=ones(Int64, n_sats),
        save_cache::SaveCache=SaveCache(),
    )
        args isa SimulationConfiguration || throw(ArgumentError("ODEParams args must be a SimulationConfiguration."))
        shared_buffers isa SharedBuffers || throw(ArgumentError("ODEParams shared_buffers must be a SharedBuffers."))
        return ODEParams{typeof(args)}(
            n_sats,
            args,
            shared_buffers,
            is_active,
            orbit_counter,
            save_cache,
        )
    end
    # solution = Solution()

    # function reset_config()
    #     global model = Model()
    #     global cnf = Cnf()
    #     global controller = Controller()
    #     global solution = Solution()
    # end
end
