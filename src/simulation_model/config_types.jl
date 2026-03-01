
module ConfigTypes
# import .SpacecraftModel
# include("simulation_model/SimulationModel.jl")
using ..PhysicalModel: SpacecraftModel
using ..SimConfig: SimulationConfiguration
using StaticArrays
using AstroTime
using OrdinaryDiffEq
using Reexport

export Initial_condition, Aerodynamics, Engines, Model, Cnf, Solution, ODEParams, IntermediateSolution, Mission, InitialParameters
export SaveCache, SaveData
export SRPSunEphemerisCache, NBodyEphemerisCache, PlanetFrameEphemerisCache, SpiceRuntimeCounters, SpiceRhsMemo
export LegacyGravityModelCode, LegacyDensityModelCode, LegacyAerodynamicModelCode, LegacyThermalModelCode, LegacyThrustControlCode
export LegacyGravityConstant, LegacyGravityInverseSquared, LegacyGravityInverseSquaredJ2, LegacyGravityGRAM
export LegacyDensityConstant, LegacyDensityExponential, LegacyDensityNoDensity, LegacyDensityGRAM, LegacyDensityNRLMSISE
export LegacyAeroCdClConstant, LegacyAeroDiffusive, LegacyAeroNoBallisticAxial
export LegacyThermalConvectiveRadiative, LegacyThermalMaxwellian
export LegacyThrustNone, LegacyThrustAerobrakingManeuver, LegacyThrustDragPassageFiring

    @enum LegacyGravityModelCode::Int8 begin
        LegacyGravityConstant = 0
        LegacyGravityInverseSquared = 1
        LegacyGravityInverseSquaredJ2 = 2
        LegacyGravityGRAM = 3
    end

    @enum LegacyDensityModelCode::Int8 begin
        LegacyDensityConstant = 0
        LegacyDensityExponential = 1
        LegacyDensityNoDensity = 2
        LegacyDensityGRAM = 3
        LegacyDensityNRLMSISE = 4
    end

    @enum LegacyAerodynamicModelCode::Int8 begin
        LegacyAeroCdClConstant = 0
        LegacyAeroDiffusive = 1
        LegacyAeroNoBallisticAxial = 2
    end

    @enum LegacyThermalModelCode::Int8 begin
        LegacyThermalConvectiveRadiative = 1
        LegacyThermalMaxwellian = 2
    end

    @enum LegacyThrustControlCode::Int8 begin
        LegacyThrustNone = 0
        LegacyThrustAerobrakingManeuver = 1
        LegacyThrustDragPassageFiring = 2
    end

    @inline function _legacy_enum_parse(::Type{T}, x::T) where {T <: Enum}
        return x
    end
    @inline function _legacy_enum_parse(::Type{LegacyGravityModelCode}, x::Integer)
        x == 0 && return LegacyGravityConstant
        x == 1 && return LegacyGravityInverseSquared
        x == 2 && return LegacyGravityInverseSquaredJ2
        x == 3 && return LegacyGravityGRAM
        throw(ArgumentError("Invalid gravity model code $x."))
    end
    @inline function _legacy_enum_parse(::Type{LegacyDensityModelCode}, x::Integer)
        x == 0 && return LegacyDensityConstant
        x == 1 && return LegacyDensityExponential
        x == 2 && return LegacyDensityNoDensity
        x == 3 && return LegacyDensityGRAM
        x == 4 && return LegacyDensityNRLMSISE
        throw(ArgumentError("Invalid density model code $x."))
    end
    @inline function _legacy_enum_parse(::Type{LegacyAerodynamicModelCode}, x::Integer)
        x == 0 && return LegacyAeroCdClConstant
        x == 1 && return LegacyAeroDiffusive
        x == 2 && return LegacyAeroNoBallisticAxial
        throw(ArgumentError("Invalid aerodynamic model code $x."))
    end
    @inline function _legacy_enum_parse(::Type{LegacyThermalModelCode}, x::Integer)
        x == 1 && return LegacyThermalConvectiveRadiative
        x == 2 && return LegacyThermalMaxwellian
        throw(ArgumentError("Invalid thermal model code $x."))
    end
    @inline function _legacy_enum_parse(::Type{LegacyThrustControlCode}, x::Integer)
        x == 0 && return LegacyThrustNone
        x == 1 && return LegacyThrustAerobrakingManeuver
        x == 2 && return LegacyThrustDragPassageFiring
        throw(ArgumentError("Invalid thrust control code $x."))
    end

    # Keep legacy numeric comparisons (`ip.gm == 1`) working while fields are typed enums.
    @inline Base.:(==)(lhs::LegacyGravityModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
    @inline Base.:(==)(lhs::Integer, rhs::LegacyGravityModelCode) = rhs == lhs
    @inline Base.:(==)(lhs::LegacyDensityModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
    @inline Base.:(==)(lhs::Integer, rhs::LegacyDensityModelCode) = rhs == lhs
    @inline Base.:(==)(lhs::LegacyAerodynamicModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
    @inline Base.:(==)(lhs::Integer, rhs::LegacyAerodynamicModelCode) = rhs == lhs
    @inline Base.:(==)(lhs::LegacyThermalModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
    @inline Base.:(==)(lhs::Integer, rhs::LegacyThermalModelCode) = rhs == lhs
    @inline Base.:(==)(lhs::LegacyThrustControlCode, rhs::Integer) = Int(lhs) == Int(rhs)
    @inline Base.:(==)(lhs::Integer, rhs::LegacyThrustControlCode) = rhs == lhs

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
            gm=_legacy_enum_parse(LegacyGravityModelCode, gm),
            dm=_legacy_enum_parse(LegacyDensityModelCode, dm),
            wm=Int64(wm),
            am=_legacy_enum_parse(LegacyAerodynamicModelCode, am),
            tm=_legacy_enum_parse(LegacyThermalModelCode, tm),
            cm=Int64(cm),
            tc=_legacy_enum_parse(LegacyThrustControlCode, tc),
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
    #     J2000_to_pci::SMatrix{3, 3, Float64} = SMatrix{3, 3, Float64}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
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
        positions_j2000::Vector{SVector{3, Float64}}
    end

    struct NBodyEphemerisCache
        primary_body_name::String
        body_query_names::Vector{String}
        body_index_by_name::Dict{String, Int}
        ets::Vector{Float64}
        positions_j2000::Matrix{SVector{3, Float64}}
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
        body_positions_j2000::Dict{String, SVector{3, Float64}} = Dict{String, SVector{3, Float64}}()
    end

    @kwdef mutable struct Solution
        orientation::Orientation = Orientation()
        physical_properties::Physical_properties = Physical_properties()
        performance::Performance = Performance()
        forces::Forces = Forces()
        simulation::Simulation = Simulation()
        closed_form::Closed_form = Closed_form()
    end

    # A struct to hold the data shared between the callback and the integrator
    @kwdef struct SharedBuffers{N_sats}
        densities::Vector{Float64} = zeros(Float64, N_sats)
        temperatures::Vector{Float64} = ones(Float64, N_sats)
        winds::Vector{SVector{3,Float64}} = [SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:N_sats]
        heat_rates::Vector{Vector{Float64}} = [Float64[] for _ in 1:N_sats]
        density_models::Vector{Any} = Any[]
        gram_density_cache::Vector{Any} = Any[]
        harmonics_workspaces::Vector{Any} = [nothing for _ in 1:N_sats]
        nbody_workspaces::Vector{Any} = [nothing for _ in 1:N_sats]
        aero_workspaces::Vector{Any} = [nothing for _ in 1:N_sats]
        nbody_ephemeris_cache::Base.RefValue{Union{Nothing, NBodyEphemerisCache}} = Ref{Union{Nothing, NBodyEphemerisCache}}(nothing)
        srp_sun_ephemeris_cache::Base.RefValue{Union{Nothing, SRPSunEphemerisCache}} = Ref{Union{Nothing, SRPSunEphemerisCache}}(nothing)
        planet_frame_ephemeris_cache::Base.RefValue{Union{Nothing, PlanetFrameEphemerisCache}} = Ref{Union{Nothing, PlanetFrameEphemerisCache}}(nothing)
        spice_runtime_counters::SpiceRuntimeCounters = SpiceRuntimeCounters()
        spice_rhs_memo_enabled::Base.RefValue{Bool} = Ref(true)
        spice_rhs_memo::SpiceRhsMemo = SpiceRhsMemo()
        current_time::Base.RefValue{Float64} = Ref(0.0)
        et_start::Base.RefValue{Float64} = Ref(0.0)
        debug_control::Base.RefValue{Bool} = Ref(false)
        debug_initial_derivative::Base.RefValue{Bool} = Ref(false)
    end

    const SaveData = Dict{Symbol, Any}
    @kwdef struct SaveCache
        # Define fields to store cached data for saving results during the simulation for expensive computations, e.g., drag, lift, density, etc.
        # Vector stores the values for each satellite if simulating multiple satellites in parallel, at a given time step, then passes to the SaveData struct at the end of the time step for saving results
        ρ_cache::Vector{Float64} = []
        heat_rate_cache::Vector{Float64} = []
        drag_cache::Vector{SVector{3,Float64}} = []
        lift_cache::Vector{SVector{3,Float64}} = []
        # Add more fields as needed to store the relevant data for saving results
    end

    
    @kwdef struct ODEParams{N_sats}
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
        args::SimulationConfiguration = SimulationConfiguration() # Arguments dictionary
        shared_buffers::SharedBuffers = SharedBuffers{N_sats}() # Shared buffers for callback and integrator
        is_active::Vector{Bool} = [true for _ in 1:N_sats] # Vector to track which satellites are still active in the simulation
        orbit_counter::Vector{Int64} = ones(Int64, N_sats) # Counter for the number of orbits completed
        save_cache::SaveCache = SaveCache() # Cache for saving results
    end
    # solution = Solution()

    # function reset_config()
    #     global model = Model()
    #     global cnf = Cnf()
    #     global controller = Controller()
    #     global solution = Solution()
    # end
end
