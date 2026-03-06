using Test
using CSV
using DataFrames
using Dates
using LinearAlgebra
using Logging
using Serialization
using StaticArrays
using ComponentArrays
using DiffEqBase
using DiffEqCallbacks
using OrdinaryDiffEq
using Quaternions
using SPICE
using TOML
using JET
using Aqua

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel
include(joinpath(REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))

# SimulationEngine uses SimulationModel and provides canonical runtime entrypoints.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end
if !isdefined(@__MODULE__, :_solver_policy_mode)
    const build_initial_conditions = SimulationEngine.build_initial_conditions
    const _build_solver_tolerances = SimulationEngine._build_solver_tolerances
    const _solve_with_solver_policy = SimulationEngine._solve_with_solver_policy
    const _split_subproblem = SimulationEngine._split_subproblem
    const _resolve_component_tolerance = SimulationEngine._resolve_component_tolerance
    const spacecraft_dynamics! = SimulationEngine.spacecraft_dynamics!
    const _debug_print_nan_parameter_paths! = SimulationEngine._debug_print_nan_parameter_paths!
    const _load_checkpoint = SimulationEngine._load_checkpoint
    const _normalize_warning_emitted = SimulationEngine._normalize_warning_emitted
    const _solver_policy_mode = SimulationEngine._solver_policy_mode
    const _gram_per_sat_instances_enabled = SimulationEngine._gram_per_sat_instances_enabled
    const _solver_maxiters = SimulationEngine._solver_maxiters
    const _effector_parallel_mode = SimulationEngine._effector_parallel_mode
    const _effector_thread_threshold = SimulationEngine._effector_thread_threshold
    const _effector_max_threads = SimulationEngine._effector_max_threads
    const _effector_long_mission_threshold_s = SimulationEngine._effector_long_mission_threshold_s
    const _effector_cost_ns_per_item_default = SimulationEngine._effector_cost_ns_per_item_default
    const _effector_cost_ema_alpha = SimulationEngine._effector_cost_ema_alpha
    const _effector_work_ns_per_worker_threshold = SimulationEngine._effector_work_ns_per_worker_threshold
    const _dynamic_effector_thread_decision = SimulationEngine._dynamic_effector_thread_decision
    const _retcode_is_stiff_symptom = SimulationEngine._retcode_is_stiff_symptom
    const _split_imex_solver_spec = SimulationEngine._split_imex_solver_spec
    const _auto_stiff_switched = SimulationEngine._auto_stiff_switched
    const _solve_with_explicit_solver = SimulationEngine._solve_with_explicit_solver
    const _solve_with_multirate_solver = SimulationEngine._solve_with_multirate_solver
    const _multirate_fast_substeps = SimulationEngine._multirate_fast_substeps
    const _multirate_slow_dt_s = SimulationEngine._multirate_slow_dt_s
    const _multirate_slow_solver_spec = SimulationEngine._multirate_slow_solver_spec
    const _multirate_fast_solver_spec = SimulationEngine._multirate_fast_solver_spec
    const _append_series_columns! = SimulationEngine._append_series_columns!
    const _atomic_write_file = SimulationEngine._atomic_write_file
    const _checkpoint_directory = SimulationEngine._checkpoint_directory
    const _checkpoint_paths = SimulationEngine._checkpoint_paths
    const _clear_ephemeris_reuse_cache! = SimulationEngine._clear_ephemeris_reuse_cache!
    const _effector_long_orbit_threshold = SimulationEngine._effector_long_orbit_threshold
    const _ephemeris_reuse_enabled = SimulationEngine._ephemeris_reuse_enabled
    const _ephemeris_reuse_max_entries = SimulationEngine._ephemeris_reuse_max_entries
    const _ephemeris_reuse_store! = SimulationEngine._ephemeris_reuse_store!
    const _find_sample_value = SimulationEngine._find_sample_value
    const _has_active_srp_effector = SimulationEngine._has_active_srp_effector
    const _initialize_aero_workspace_buffers! = SimulationEngine._initialize_aero_workspace_buffers!
    const _initialize_harmonics_workspace_buffers! = SimulationEngine._initialize_harmonics_workspace_buffers!
    const _initialize_heat_rate_buffers! = SimulationEngine._initialize_heat_rate_buffers!
    const _initialize_nbody_ephemeris_cache! = SimulationEngine._initialize_nbody_ephemeris_cache!
    const _initialize_nbody_ephemeris_cache_buffer! = SimulationEngine._initialize_nbody_ephemeris_cache_buffer!
    const _initialize_nbody_workspace_buffers! = SimulationEngine._initialize_nbody_workspace_buffers!
    const _initialize_planet_frame_cache_buffer! = SimulationEngine._initialize_planet_frame_cache_buffer!
    const _initialize_planet_frame_ephemeris_cache! = SimulationEngine._initialize_planet_frame_ephemeris_cache!
    const _initialize_srp_sun_cache_buffer! = SimulationEngine._initialize_srp_sun_cache_buffer!
    const _initialize_srp_sun_ephemeris_cache! = SimulationEngine._initialize_srp_sun_ephemeris_cache!
    const _mission_is_long_for_effector_threads = SimulationEngine._mission_is_long_for_effector_threads
    const _nbody_ephemeris_cache_dt_s = SimulationEngine._nbody_ephemeris_cache_dt_s
    const _nbody_ephemeris_cache_max_samples = SimulationEngine._nbody_ephemeris_cache_max_samples
    const _planet_frame_cache_dt_s = SimulationEngine._planet_frame_cache_dt_s
    const _planet_frame_cache_max_samples = SimulationEngine._planet_frame_cache_max_samples
    const _reset_spice_rhs_memo! = SimulationEngine._reset_spice_rhs_memo!
    const _reset_spice_runtime_counters! = SimulationEngine._reset_spice_runtime_counters!
    const _srp_ephemeris_cache_dt_s = SimulationEngine._srp_ephemeris_cache_dt_s
    const _srp_ephemeris_cache_max_samples = SimulationEngine._srp_ephemeris_cache_max_samples
    const _validate_orientation_inertia! = SimulationEngine._validate_orientation_inertia!
    const _validate_thermal_model_support! = SimulationEngine._validate_thermal_model_support!
end
# Simulation campaign entrypoint.
include(joinpath(REPO_ROOT, "src", "mission", "campaigns", "run.jl"))

struct ConstantDensityModel <: SimulationModel.AbstractDensityModel
    rho::Float64
    temp::Float64
end

struct TimedTangentialThrusterModel <: SimulationModel.AbstractControlEffectorModel
    thrust::Float64
    direction_sign::Float64 # +1.0 prograde, -1.0 retrograde
    start_time::Float64
    stop_time::Float64
end

mutable struct CountingGuidanceModel <: SimulationModel.AbstractTypes.AbstractGuidanceModel
    hits::Vector{Int}
end

mutable struct CountingNavigationModel
    hits::Vector{Int}
end

mutable struct CountingControlModel <: SimulationModel.AbstractControlEffectorModel
    hits::Vector{Int}
end

struct ThrowingForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct NaNForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct NaNParamForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct ConstantTorqueModel <: SimulationModel.AbstractForceTorqueModel
    torque::SVector{3, Float64}
end

struct ConstantForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
end

struct ThrowingOrbitPlanet <: SimulationModel.AbstractPlanet
    Rp_e::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlForceTorque(
    model::TimedTangentialThrusterModel,
    u::AbstractVector{Float64},
    p::ODEParams,
    i::Int64,
    t::Float64
)
    if t < model.start_time || t > model.stop_time
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    v = SVector{3, Float64}(u.vel)
    vm = norm(v)
    if vm == 0.0 || !isfinite(vm)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    force = model.thrust * model.direction_sign * (v / vm)
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlEffect!(
    model::TimedTangentialThrusterModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

function SimulationModel.calcGuidanceEffect!(
    model::CountingGuidanceModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcNavigationEffect!(
    model::CountingNavigationModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcControlEffect!(
    model::CountingControlModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcForceTorque(
    model::ThrowingForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    error("intentional derivative failure")
end

function SimulationModel.calcForceTorque(
    model::NaNForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(NaN, NaN, NaN), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::NaNParamForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    p.shared_buffers.current_time[] = NaN
    return SVector{3, Float64}(NaN, NaN, NaN), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::ConstantTorqueModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), model.torque
end

function SimulationModel.calcForceTorque(
    model::ConstantForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return model.force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.ControlHooks.rvtoorbitalelement(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet::ThrowingOrbitPlanet
)
    throw(ErrorException("forced-orbital-element-conversion-failure"))
end

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

angle_distance(a::Float64, b::Float64) = abs(mod(a - b + pi, 2pi) - pi)

function specific_energy(df::DataFrame, mu::Float64)
    r = sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    v2 = df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2
    return 0.5 .* v2 .- mu ./ r
end

function make_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    panel = Link{0}(root=false, m=30.0, ref_area=6.0, r=MVector{3, Float64}(0.0, 1.2, 0.0))

    if isnothing(orientation_state)
        ic = InitialCondition(
            ra=EARTH.Rp_e + ra_alt_m,
            rp=EARTH.Rp_e + rp_alt_m,
            i=i_deg,
            ω=ω_deg,
            Ω=Ω_deg,
            ν=ν_deg
        )
    else
        q0, w0 = orientation_state
        ra = EARTH.Rp_e + ra_alt_m
        rp = EARTH.Rp_e + rp_alt_m
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        ic = InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = root.m + panel.m
    return SpacecraftModel(
        Joint[],
        [root, panel],
        root,
        true,
        dry_mass,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function make_single_link_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=EARTH.Rp_e + ra_alt_m,
        rp=EARTH.Rp_e + rp_alt_m,
        i=i_deg,
        ω=ω_deg,
        Ω=Ω_deg,
        ν=ν_deg
    )

    return SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function make_base_thruster_model(;
    thrust::Float64,
    direction::Float64=0.0,
    Δv::Float64=0.0,
    start_burn_time::Float64=0.0,
    stop_burn_time::Float64=0.0,
    Isp::Float64=300.0
)
    return BaseThrusterModel(
        thrust=[thrust],
        direction=[direction],
        Δv=[Δv],
        start_burn_time=[start_burn_time],
        stop_burn_time=[stop_burn_time],
        Isp=[Isp]
    )
end

function build_config(;
    spacecraft::SpacecraftModel,
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::SimulationModel.InitialTime=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    planet=EARTH
)
    environment_model = EnvironmentModel(
        planet=planet,
        EI=EI_km,
        density_model=density_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    )

    return SimulationConfiguration(
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

function build_config_multi(;
    spacecraft::Vector{SpacecraftModel},
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::SimulationModel.InitialTime=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    planet=EARTH
)
    environment_model = EnvironmentModel(
        planet=planet,
        EI=EI_km,
        density_model=density_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    )

    return SimulationConfiguration(
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

function interp_linear(times::AbstractVector{<:Real}, values::AbstractVector{<:Real}, t::Float64)
    idx = searchsortedlast(times, t)
    if idx <= 0
        return Float64(values[1])
    elseif idx >= length(times)
        return Float64(values[end])
    end
    t1 = Float64(times[idx])
    t2 = Float64(times[idx + 1])
    y1 = Float64(values[idx])
    y2 = Float64(values[idx + 1])
    α = (t - t1) / (t2 - t1)
    return (1.0 - α) * y1 + α * y2
end

function make_agora_earth_spacecraft()
    main_bus = Link{0}(root=true, m=620.0, ref_area=2.05 * 2.8)
    left_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, -2.05 / 2.0 - 5.7 / 4.0, 0.0))
    right_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, 2.05 / 2.0 + 5.7 / 4.0, 0.0))
    ic = InitialCondition(
        ra=56_378.7978559e3,
        rp=EARTH.Rp_e + 200_590.0,
        i=89.876,
        ω=75.505,
        Ω=104.115,
        ν=175.0
    )

    return SpacecraftModel(
        Joint[],
        [main_bus, left_panel, right_panel],
        main_bus,
        true,
        main_bus.m + left_panel.m + right_panel.m,
        200.0,
        main_bus.inertia,
        0,
        0,
        ic,
        1
    )
end

function run_case(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            run_simulation(args; isolate_state=isolate_state)
            @test isfile(results_csv_path)
            return CSV.read(results_csv_path, DataFrame)
        end
    end
end

function run_case_silent(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                run_simulation(args; isolate_state=isolate_state)
            end
            @test isfile(results_csv_path)
            return CSV.read(results_csv_path, DataFrame)
        end
    end
end

function run_case_capture_stdout(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            output = ""
            mktemp() do path, io
                redirect_stdout(io) do
                    run_simulation(args; isolate_state=isolate_state)
                end
                flush(io)
                seekstart(io)
                output = read(io, String)
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame), output
            else
                @test !isfile(results_csv_path)
                return DataFrame(), output
            end
        end
    end
end

function run_case_via_execute_analysis(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                execute_analysis(args; isolate_state=isolate_state)
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame)
            else
                @test !isfile(results_csv_path)
                return DataFrame()
            end
        end
    end
end

function run_case_via_campaign(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true, state=nothing)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                if state === nothing
                    execute_campaign(args; isolate_state=isolate_state)
                else
                    execute_campaign(args, state; isolate_state=isolate_state)
                end
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame)
            else
                @test !isfile(results_csv_path)
                return DataFrame()
            end
        end
    end
end

function run_case_via_execute_campaign(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true, state=nothing)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                if state === nothing
                    execute_campaign(args; isolate_state=isolate_state)
                else
                    execute_campaign(args, state; isolate_state=isolate_state)
                end
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame)
            else
                @test !isfile(results_csv_path)
                return DataFrame()
            end
        end
    end
end

function run_case_via_execute_orbital_elements_campaign(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                execute_orbital_elements_campaign(args; isolate_state=isolate_state)
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame)
            else
                @test !isfile(results_csv_path)
                return DataFrame()
            end
        end
    end
end

function run_case_via_execute_vgamma_campaign(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                execute_vgamma_campaign(args; isolate_state=isolate_state)
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame)
            else
                @test !isfile(results_csv_path)
                return DataFrame()
            end
        end
    end
end

function seed_solution_for_save_csv!(
    solution::Solution;
    n_bodies::Int,
    n_reaction_wheels::Int,
    n_thrusters::Int,
    base::Float64=1.0,
    closed_form::Bool=false
)
    solution.physical_properties.α = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.β = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_rate = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_load = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.rw_h = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.rw_τ = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.thruster_forces = [Float64[] for _ in 1:n_thrusters]

    function _push_sample!(field)
        if field isa Vector{Float64}
            push!(field, base)
        elseif field isa Vector{Int64}
            push!(field, 1)
        elseif field isa Vector{Vector{Float64}}
            for subfield in field
                push!(subfield, base)
            end
        end
        return nothing
    end

    for group in (solution.orientation, solution.physical_properties, solution.performance, solution.forces)
        for fname in fieldnames(typeof(group))
            _push_sample!(getfield(group, fname))
        end
    end

    if closed_form
        solution.closed_form.t_cf = [base + 10.0]
        solution.closed_form.h_cf = [base + 20.0]
        solution.closed_form.γ_cf = [base + 30.0]
        solution.closed_form.v_cf = [base + 40.0]
    end

    return solution
end


module IncludeOrderSandbox
end
const INCLUDE_ORDER_SANDBOX = IncludeOrderSandbox

module ExportImportSandbox
end
const EXPORT_IMPORT_SANDBOX = ExportImportSandbox

const GUIDANCE_SANDBOX_LOADED = Ref(false)
module GuidanceSandbox
using ..SimulationModel
using ..SimulationModel.AbstractTypes: AbstractGuidanceModel
using ComponentArrays
end
const GUIDANCE_SANDBOX = GuidanceSandbox

@testset "SimulationModel Export Contract" begin
    sandbox = EXPORT_IMPORT_SANDBOX
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
    @test_nowarn Core.eval(sandbox, :(using .SimulationModel))

    required_public_names = [
        :SimulationConfiguration,
        :InitialCondition,
        :InitialTime,
        :MissionConfiguration,
        :EnvironmentModel,
        :IntegrationTolerances,
        :ControlModel,
        :GuidanceModel,
        :NavigationModel,
        :DynamicsModel
    ]

    for sym in required_public_names
        @test Core.eval(sandbox, :(isdefined(@__MODULE__, $(QuoteNode(sym)))))
    end
end

@testset "Include-Order + Name Ambiguity Smoke" begin
    sandbox = INCLUDE_ORDER_SANDBOX

    Base.include_string(sandbox, """
    module ConflictingExports
        export SimulationConfiguration
        struct SimulationConfiguration end
    end
    using .ConflictingExports
    """)

    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
    Core.eval(sandbox, :(const quat_mult = SimulationModel.quat_mult))
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
    @test isdefined(sandbox, :SimulationEngine)
    @test Core.eval(sandbox, :(isdefined(SimulationEngine, :run_simulation)))
end

@testset "Simulation Filename Canonical Contract" begin
    engine_path = joinpath(REPO_ROOT, "src", "simulation", "engine", "execution.jl")
    legacy_exec_dir = joinpath(REPO_ROOT, "src", "simulation", "execution")
    legacy_run_path = joinpath(legacy_exec_dir, "run_simulation.jl")
    legacy_elements_path = joinpath(legacy_exec_dir, "simulation_elements.jl")
    legacy_execution_path = joinpath(legacy_exec_dir, "simulation_execution.jl")
    legacy_aerobraking_path = joinpath(REPO_ROOT, "src", "simulation", "Aerobraking.jl")
    legacy_complete_path = joinpath(REPO_ROOT, "src", "simulation", "Complete_passage.jl")

    @test isfile(engine_path)
    @test !isfile(legacy_run_path)
    @test !isfile(legacy_elements_path)
    @test !isfile(legacy_execution_path)
    @test !isfile(legacy_aerobraking_path)
    @test !isfile(legacy_complete_path)
end











function ensure_guidance_sandbox_loaded!()
    if GUIDANCE_SANDBOX_LOADED[]
        return
    end

    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_models.jl"))))
    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_functions.jl"))))
    GUIDANCE_SANDBOX_LOADED[] = true
end









@testset "API Convenience Constructors" begin
    @testset "Mission/Environment Config Validation" begin
        mc_default = MissionConfiguration()
        @test mc_default.mission_type == MissionTime

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_str = MissionConfiguration(mission_type="Time", mission_time=600.0, number_of_orbits=1, num_steps_to_save=10)
            @test mc_str.mission_type == MissionTime
            @test mc_str.mission_type == "Time"
            @test "Time" == mc_str.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_sym = MissionConfiguration(mission_type=:orbits, mission_time=600.0, number_of_orbits=2, num_steps_to_save=10)
            @test mc_sym.mission_type == MissionOrbits
            @test mc_sym.mission_type == :orbits
            @test :orbits == mc_sym.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_enum = MissionConfiguration(mission_type=MissionOrbits, mission_time=600.0, number_of_orbits=3, num_steps_to_save=10)
            @test mc_enum.mission_type == MissionOrbits
            @test mc_enum.mission_type == "Orbits"
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "1") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs (:warn, r"mission_type=.*deprecated") MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == true
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == false
        end

        @test_throws ArgumentError MissionConfiguration(mission_type="invalid")
        @test_throws ArgumentError MissionConfiguration(mission_time=0.0)
        @test_throws ArgumentError MissionConfiguration(number_of_orbits=0)
        @test_throws ArgumentError MissionConfiguration(num_steps_to_save=0)

        env_ok = EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=8,
            wind=false
        )
        @test env_ok.EI == 120.0

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=-1.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=false,
            topo_degree=8,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=-1,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=-1,
            wind=false
        )
    end

    ic = InitialCondition()
    @test ic isa InitialCondition
    @test ic.a == 0.0
    @test ic.e == 0.0

    link = Link()
    @test link isa Link{0}

    joint = Joint()
    @test joint isa Joint

    sc = SpacecraftModel()
    @test sc isa SpacecraftModel
    @test sc.root.root

    custom_root = Link{0}(root=true, m=100.0)
    sc_custom = SpacecraftModel(root=custom_root, id=42)
    @test sc_custom.id == 42
    @test sc_custom.root === custom_root
    @test any(link -> link === custom_root, sc_custom.links)

    nbody = NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)
    @test nbody isa NBodyGravityModel
    @test nbody.primary_body_name == "Earth"
    @test nbody.body_names == ("Sun",)

    nbody_mars = NBodyGravityModel(["Sun"], "Mars", SPICE_PATH)
    nbody_venus = NBodyGravityModel(["Sun"], "Venus", SPICE_PATH)
    nbody_titan = NBodyGravityModel(["Sun"], "Titan", SPICE_PATH)
    @test lowercase(nbody_mars.planet.name) == "mars"
    @test lowercase(nbody_venus.planet.name) == "venus"
    @test lowercase(nbody_titan.planet.name) == "titan"
    @test_throws ArgumentError NBodyGravityModel(["Sun"], "Pluto", SPICE_PATH)

    nbody_jupiter = NBodyGravityModel(["Jupiter"], "Earth", SPICE_PATH)
    nbody_state = [EARTH.Rp_e + 500e3, 0.0, 0.0, 0.0, 0.0, 0.0, 500.0]
    sc_nbody = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=175.0)
    args_nbody = build_config(
        spacecraft=sc_nbody,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    force_nbody, torque_nbody = calcForceTorque(nbody_jupiter, nbody_state, ODEParams{1}(args=args_nbody), 1)
    @test all(isfinite, force_nbody)
    @test torque_nbody == SVector{3, Float64}(0.0, 0.0, 0.0)

    harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics_l20 = GravitationalHarmonicsModel(20, 20, harmonics_file, EARTH)
    @test size(harmonics_l20.C) == (21, 21)
    @test size(harmonics_l20.S) == (21, 21)
    @test_throws ArgumentError GravitationalHarmonicsModel(10, 11, harmonics_file, EARTH)

    child_link = Link(root=false, q=MVector{4, Float64}(sin(pi / 4), 0.0, 0.0, cos(pi / 4)))
    rot_child = rotate_to_body(child_link)
    @test size(rot_child) == (3, 3)
    @test isapprox(det(Matrix(rot_child)), 1.0; atol=1e-12)
    @test norm(Matrix(rot_child) - Matrix{Float64}(I, 3, 3)) > 0.1

    @testset "Quaternion DCM Conversion Negative-Trace Branch" begin
        dcm_180_x = SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
        q_neg_trace = SimulationModel.dcm_to_quaternion(dcm_180_x)
        @test isapprox(norm(q_neg_trace), 1.0; atol=1e-12, rtol=0.0)

        dcm_roundtrip = SimulationModel.rot(q_neg_trace)
        @test isapprox(Matrix(dcm_roundtrip), Matrix(dcm_180_x); atol=1e-12, rtol=0.0)
    end

    @testset "Effector Rate Validation" begin
        @test GuidanceModel((), Float64[]) isa GuidanceModel
        @test NavigationModel((), Float64[]) isa NavigationModel
        @test ControlModel((), Float64[]) isa ControlModel

        @test_throws ArgumentError GuidanceModel((:g1,), Float64[])
        @test_throws ArgumentError NavigationModel((:n1,), Float64[])
        @test_throws ArgumentError ControlModel((:c1,), Float64[])

        @test_throws ArgumentError GuidanceModel((:g1,), [0.0])
        @test_throws ArgumentError NavigationModel((:n1,), [-1.0])
        @test_throws ArgumentError ControlModel((:c1,), [Inf])

        sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
        sc2 = make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        args_bad_slots = build_config_multi(
            spacecraft=[sc1, sc2],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=60.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(make_base_thruster_model(thrust=1.0, Δv=1.0, start_burn_time=-1.0, stop_burn_time=-1.0),),
            control_rates=[1.0],
            keplerian=true
        )
        @test_throws ArgumentError run_case_silent(args_bad_slots)
    end
end

@testset "Run Simulation Metadata Return" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        keplerian=true
    )

    metadata = mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; return_solution=true, return_solver_metadata=true)
        end
    end

    @test metadata isa NamedTuple
    @test hasproperty(metadata, :solution)
    @test hasproperty(metadata, :solver_mode)
    @test hasproperty(metadata, :solver_trace)
    @test hasproperty(metadata, :parallel_policy)
    @test hasproperty(metadata, :spice_counters)
end

@testset "RHS Completeness: Mass Derivative" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    u = build_initial_conditions(args)
    du = copy(u)
    du.sc[1].mass = 789.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du, u, p, 0.0)
    @test du.sc[1].mass == 0.0

    du_inactive = copy(u)
    du_inactive.sc[1].mass = 123.0
    p_inactive = ODEParams{1}(args=args, is_active=[false])
    spacecraft_dynamics!(du_inactive, u, p_inactive, 0.0)
    @test du_inactive.sc[1].mass == 0.0
end







@testset "Guidance Thruster Campaign Helpers" begin
    ensure_guidance_sandbox_loaded!()
    sandbox = GUIDANCE_SANDBOX

    guidance_model = sandbox.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=[2, 3],
        maneuver_Δv=[-5.0, 20.0]
    )

    thruster = BaseThrusterModel(
        thrust=[500.0],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )
    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p = ODEParams{1}(args=args)
    u = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])

    p.orbit_counter[1] = 3
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 20.0
    @test isapprox(thruster.direction[1], π; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 2
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 5.0
    @test isapprox(thruster.direction[1], 0.0; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 1
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 0.0
end

@testset "Callback Internal Branch Coverage" begin
    mutable struct MockCallbackOpts
        dtmax::Float64
        reltol::Float64
        abstol::Float64
    end

    mutable struct MockCallbackIntegrator{P, U, O}
        p::P
        u::U
        t::Float64
        opts::O
        tdir::Int
        tstop_max::Float64
    end
    DiffEqBase.get_tstops_max(integrator::MockCallbackIntegrator) = integrator.tstop_max
    DiffEqBase.add_tstop!(integrator::MockCallbackIntegrator, tstop) = begin
        integrator.tstop_max = max(integrator.tstop_max, Float64(tstop))
        nothing
    end

    args_base = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false)
    )
    args_orient = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false)
    )

    mission_orbits = MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=true,
        number_of_orbits=1,
        mission_time=120.0,
        orientation_sim=false,
        num_steps_to_save=100
    )
    args_orbits = SimulationConfiguration(
        file_paths=args_base.file_paths,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false),
        mission_configuration=mission_orbits,
        environment_model=args_base.environment_model,
        dynamics_model=args_base.dynamics_model,
        guidance_model=args_base.guidance_model,
        navigation_model=args_base.navigation_model,
        control_model=args_base.control_model,
        initial_time=args_base.initial_time,
        integration_tolerances=args_base.integration_tolerances
    )

    _ = SimulationModel.SimulationCallbacks.get_callbacks(
        1,
        args_orbits.dynamics_model.dynamic_effectors,
        args_orbits
    )

    orbit_cb = SimulationModel.SimulationCallbacks.get_orbit_end_callback(1)
    p_orbit = ODEParams{1}(args=args_orbits)
    u_orbit = build_initial_conditions(args_orbits)
    integrator_orbit = MockCallbackIntegrator(
        p_orbit,
        u_orbit,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    out_orbit = zeros(1)
    orbit_cb.condition(out_orbit, u_orbit, 0.0, integrator_orbit)
    @test isfinite(out_orbit[1])
    orbit_count_before = p_orbit.orbit_counter[1]
    orbit_cb.affect!(integrator_orbit, 1)
    @test p_orbit.orbit_counter[1] == orbit_count_before + 1

    args_impact = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    impact_cb = SimulationModel.SimulationCallbacks.get_impact_callback(2)
    p_impact = ODEParams{2}(args=args_impact)
    u_impact = build_initial_conditions(args_impact)
    integrator_impact = MockCallbackIntegrator(
        p_impact,
        u_impact,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    impact_out = zeros(2)
    impact_cb.condition(impact_out, u_impact, 0.0, integrator_impact)
    @test all(impact_out .> 0.0)
    @test impact_cb.affect! === nothing
    impact_cb.affect_neg!(integrator_impact, 1)
    @test p_impact.is_active[1] == false
    @test p_impact.is_active[2] == true

    args_drag = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    drag_cb = SimulationModel.SimulationCallbacks.get_drag_state_callback(1)
    p_drag = ODEParams{1}(args=args_drag)
    u_drag = build_initial_conditions(args_drag)
    integrator_drag = MockCallbackIntegrator(
        p_drag,
        u_drag,
        12.0,
        MockCallbackOpts(
            args_drag.integration_tolerances.dt_max_atmosphere,
            args_drag.integration_tolerances.reltol_atmosphere,
            args_drag.integration_tolerances.abstol_atmosphere
        ),
        1,
        Inf
    )
    drag_output = ""
    mktemp() do path, io
        redirect_stdout(io) do
            drag_cb.affect!(integrator_drag, 1)
        end
        flush(io)
        seekstart(io)
        drag_output = read(io, String)
    end
    @test occursin("Switching to space integration", drag_output)
    @test integrator_drag.opts.dtmax == args_drag.integration_tolerances.dt_max_orbit
    @test integrator_drag.opts.reltol == args_drag.integration_tolerances.reltol_orbit
    @test integrator_drag.opts.abstol == args_drag.integration_tolerances.abstol_orbit

    quat_proj_cb = SimulationModel.SimulationCallbacks.get_quaternion_projection_callback(1, args_orient)
    p_orient = ODEParams{1}(args=args_orient)
    u_orient = build_initial_conditions(args_orient)
    u_orient.sc[1].q .= [0.0, 0.0, 0.0, 2.0]
    integrator_orient = MockCallbackIntegrator(
        p_orient,
        u_orient,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    @test quat_proj_cb.condition(u_orient, 0.0, integrator_orient) == true
    quat_proj_cb.affect!(integrator_orient)
    @test isapprox(norm(integrator_orient.u.sc[1].q), 1.0; atol=1e-12, rtol=0.0)

    counting_navigation = CountingNavigationModel([0])
    args_navigation = SimulationConfiguration(
        file_paths=args_base.file_paths,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=args_base.mission_configuration,
        environment_model=args_base.environment_model,
        dynamics_model=args_base.dynamics_model,
        guidance_model=args_base.guidance_model,
        navigation_model=NavigationModel(navigation_effectors=(counting_navigation,), navigation_rates=[1.0]),
        control_model=args_base.control_model,
        initial_time=args_base.initial_time,
        integration_tolerances=args_base.integration_tolerances
    )
    navigation_cbs = SimulationModel.SimulationCallbacks.get_navigation_callbacks(1, args_navigation)
    p_navigation = ODEParams{1}(args=args_navigation)
    u_navigation = build_initial_conditions(args_navigation)
    integrator_navigation = MockCallbackIntegrator(
        p_navigation,
        u_navigation,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    navigation_cbs[1].affect!.affect!(integrator_navigation)
    @test counting_navigation.hits == [1]
    withenv("SPACEAGORA_DEV_HOT_RELOAD" => "1") do
        navigation_cbs_hot = SimulationModel.SimulationCallbacks.get_navigation_callbacks(1, args_navigation)
        navigation_cbs_hot[1].affect!.affect!(integrator_navigation)
    end
    @test counting_navigation.hits == [2]

    counting_guidance = CountingGuidanceModel([0])
    args_guidance = SimulationConfiguration(
        file_paths=args_base.file_paths,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=args_base.mission_configuration,
        environment_model=args_base.environment_model,
        dynamics_model=args_base.dynamics_model,
        guidance_model=GuidanceModel(guidance_effectors=(counting_guidance,), guidance_rates=[1.0]),
        navigation_model=args_base.navigation_model,
        control_model=args_base.control_model,
        initial_time=args_base.initial_time,
        integration_tolerances=args_base.integration_tolerances
    )
    guidance_cbs = SimulationModel.SimulationCallbacks.get_guidance_callbacks(1, args_guidance)
    p_guidance = ODEParams{1}(args=args_guidance)
    u_guidance = build_initial_conditions(args_guidance)
    integrator_guidance = MockCallbackIntegrator(
        p_guidance,
        u_guidance,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    guidance_cbs[1].affect!.affect!(integrator_guidance)
    @test counting_guidance.hits == [1]
    withenv("SPACEAGORA_DEV_HOT_RELOAD" => "1") do
        guidance_cbs_hot = SimulationModel.SimulationCallbacks.get_guidance_callbacks(1, args_guidance)
        guidance_cbs_hot[1].affect!.affect!(integrator_guidance)
    end
    @test counting_guidance.hits == [2]

    control_model = CountingControlModel([0, 0])
    args_control = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(control_model,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_control = ODEParams{2}(args=args_control)
    u_control = build_initial_conditions(args_control)
    integrator_control = MockCallbackIntegrator(
        p_control,
        u_control,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv("SPACEAGORA_DEV_HOT_RELOAD" => "1") do
        control_cbs = SimulationModel.SimulationCallbacks.get_control_callbacks(2, args_control)
        control_cbs[1].affect!.affect!(integrator_control)
    end
    @test control_model.hits == [1, 1]

    requires_density = SimulationModel.SimulationCallbacks._requires_density_callback
    requires_orbit_end = SimulationModel.SimulationCallbacks._requires_orbit_end_callback
    requires_drag_state = SimulationModel.SimulationCallbacks._requires_drag_state_callback
    requires_quat_projection = SimulationModel.SimulationCallbacks._requires_quaternion_projection_callback
    density_use_threads = SimulationModel.SimulationCallbacks._density_callback_use_threads
    control_use_threads = SimulationModel.SimulationCallbacks._control_callback_use_threads

    @test requires_density((InverseSquaredGravityModel(),), args_base) == false
    @test requires_density((InverseSquaredGravityModel(), AerodynamicCoefficientfM()), args_base) == true
    @test requires_density((InverseSquaredGravityModel(),), args_drag) == true
    @test requires_orbit_end(args_base) == false
    @test requires_orbit_end(args_orbits) == true
    @test requires_drag_state((InverseSquaredGravityModel(),), args_base) == false
    @test requires_drag_state((InverseSquaredGravityModel(), AerodynamicCoefficientfM()), args_drag) == true
    @test requires_quat_projection(args_base) == false
    @test requires_quat_projection(args_orient) == true

    has_worker_threads = Threads.nthreads() > 1
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        @test density_use_threads(args_control, 8) == false
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        @test density_use_threads(args_control, 8) == has_worker_threads
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "64"
    ) do
        @test density_use_threads(args_control, 2) == has_worker_threads
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test density_use_threads(args_control, 8) == false
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1"
    ) do
        @test density_use_threads(args_control, 8) == has_worker_threads
    end

    control_thruster = BaseThrusterModel(
        thrust=[0.5, 0.6],
        direction=[0.0, π],
        Δv=[0.0, 0.0],
        start_burn_time=[0.0, 0.0],
        stop_burn_time=[10.0, 10.0],
        Isp=[300.0, 300.0]
    )
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        @test control_use_threads(control_thruster, 8, false) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        @test control_use_threads(control_thruster, 8, false) == has_worker_threads
        @test control_use_threads(control_thruster, 8, true) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test control_use_threads(control_thruster, 8, false) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1"
    ) do
        @test control_use_threads(control_thruster, 8, false) == has_worker_threads
    end

    n_parallel_sats = 4
    threaded_spacecraft = [
        make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0 - 5.0 * (k - 1))
        for k in 1:n_parallel_sats
    ]
    threaded_control = CountingControlModel(zeros(Int, n_parallel_sats))
    args_control_parallel = build_config_multi(
        spacecraft=threaded_spacecraft,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(threaded_control,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_control_parallel = ODEParams{n_parallel_sats}(args=args_control_parallel)
    u_control_parallel = build_initial_conditions(args_control_parallel)
    integrator_control_parallel = MockCallbackIntegrator(
        p_control_parallel,
        u_control_parallel,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv(
        "SPACEAGORA_DEV_HOT_RELOAD" => "0",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1"
    ) do
        control_cbs_parallel = SimulationModel.SimulationCallbacks.get_control_callbacks(
            n_parallel_sats,
            args_control_parallel
        )
        control_cbs_parallel[1].affect!.affect!(integrator_control_parallel)
    end
    @test threaded_control.hits == ones(Int, n_parallel_sats)

    args_density_parallel = build_config_multi(
        spacecraft=threaded_spacecraft,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_density_parallel = ODEParams{n_parallel_sats}(args=args_density_parallel)
    u_density_parallel = build_initial_conditions(args_density_parallel)
    integrator_density_parallel = MockCallbackIntegrator(
        p_density_parallel,
        u_density_parallel,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        density_cb_parallel = SimulationModel.SimulationCallbacks.get_density_callback(
            n_parallel_sats,
            args_density_parallel
        )
        density_cb_parallel.affect!(integrator_density_parallel)
    end
    @test all(isfinite, p_density_parallel.shared_buffers.densities)
    @test all(ρ -> ρ >= 0.0, p_density_parallel.shared_buffers.densities)
end

@testset "Callback Env Helper Branch Coverage" begin
    callbacks = SimulationModel.SimulationCallbacks

    @test callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", false) == false
    withenv("SPACEAGORA_CB_TEST_BOOL" => "yes") do
        @test callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", false) == true
    end
    withenv("SPACEAGORA_CB_TEST_BOOL" => "off") do
        @test callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", true) == false
    end
    withenv("SPACEAGORA_CB_TEST_BOOL" => "invalid") do
        @test_throws ArgumentError callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", false)
    end

    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off") do
        @test callbacks._density_callback_parallel_mode() == :off
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on") do
        @test callbacks._density_callback_parallel_mode() == :on
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto") do
        @test callbacks._density_callback_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "invalid") do
        @test_throws ArgumentError callbacks._density_callback_parallel_mode()
    end

    withenv("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "4") do
        @test callbacks._density_callback_thread_threshold() == 4
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "0") do
        @test callbacks._density_callback_thread_threshold() == 1
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError callbacks._density_callback_thread_threshold()
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1") do
        @test callbacks._density_callback_allow_with_outer() == true
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "invalid") do
        @test_throws ArgumentError callbacks._density_callback_allow_with_outer()
    end

    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "off") do
        @test callbacks._control_callback_parallel_mode() == :off
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "on") do
        @test callbacks._control_callback_parallel_mode() == :on
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto") do
        @test callbacks._control_callback_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "invalid") do
        @test_throws ArgumentError callbacks._control_callback_parallel_mode()
    end

    withenv("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "4") do
        @test callbacks._control_callback_thread_threshold() == 4
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "0") do
        @test callbacks._control_callback_thread_threshold() == 1
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError callbacks._control_callback_thread_threshold()
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1") do
        @test callbacks._control_callback_allow_with_outer() == true
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "invalid") do
        @test_throws ArgumentError callbacks._control_callback_allow_with_outer()
    end

    @test callbacks.density_model_threadsafe(NoAtmosphereModel()) == true
    @test callbacks.density_model_threadsafe(ConstantDensityModel(1e-6, 220.0)) == false

    thruster_model = make_base_thruster_model(thrust=0.1, Δv=0.0, start_burn_time=0.0, stop_burn_time=1.0)
    @test callbacks.control_model_threadsafe(thruster_model) == true
    @test callbacks.control_model_threadsafe(CountingControlModel([0])) == false

    callbacks._gram_runtime_stats_reset!()
    @test callbacks._gram_runtime_stats_update!(s -> begin
        s.density_calls += 1
        s.direct_calls += 1
    end) === nothing
    stats_snap = callbacks._gram_runtime_stats_snapshot()
    @test stats_snap.density_calls == 1
    @test stats_snap.direct_calls == 1

    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "off") do
        @test callbacks._gram_track_cache_ignore_time_window() == false
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_TARGET_USE_J2" => "on") do
        @test callbacks._gram_track_cache_target_use_j2() == true
    end

    withenv("SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "on") do
        @test callbacks._thermal_callback_parallel_mode() == :on
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => nothing
    ) do
        @test callbacks._thermal_callback_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => "5") do
        @test callbacks._thermal_callback_thread_threshold() == 5
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "6",
        "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => nothing
    ) do
        @test callbacks._thermal_callback_thread_threshold() == 6
    end
    withenv("SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1") do
        @test callbacks._thermal_callback_allow_with_outer() == true
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1",
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => nothing
    ) do
        @test callbacks._thermal_callback_allow_with_outer() == true
    end
    withenv("SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "invalid") do
        @test_throws ArgumentError callbacks._thermal_callback_allow_with_outer()
    end

    @test callbacks._is_gram_density_model(NoAtmosphereModel()) == false

    args_density_lookup = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_density_lookup = ODEParams{1}(args=args_density_lookup)
    resize!(p_density_lookup.shared_buffers.density_models, 1)
    p_density_lookup.shared_buffers.density_models[1] = ConstantDensityModel(1e-6, 220.0)
    @test callbacks._density_model_for_sat(p_density_lookup, 1) isa ConstantDensityModel
    @test callbacks._density_model_for_sat(p_density_lookup, 2) isa NoAtmosphereModel

    withenv("SPACEAGORA_CB_TEST_FLOAT" => "oops") do
        @test_throws ArgumentError callbacks._parse_float_env("SPACEAGORA_CB_TEST_FLOAT", 1.0)
    end
    withenv("SPACEAGORA_CB_TEST_FLOAT_OPT" => "oops") do
        @test_throws ArgumentError callbacks._parse_float_env_optional("SPACEAGORA_CB_TEST_FLOAT_OPT")
    end
    withenv("SPACEAGORA_CB_TEST_INT_OPT" => "oops") do
        @test_throws ArgumentError callbacks._parse_int_env_optional("SPACEAGORA_CB_TEST_INT_OPT")
    end

    withenv("SPACEAGORA_GRAM_TRACK_CACHE" => "on") do
        @test callbacks._gram_track_cache_mode() == :on
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE" => "auto") do
        @test callbacks._gram_track_cache_mode() == :auto
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE" => "invalid") do
        @test_throws ArgumentError callbacks._gram_track_cache_mode()
    end

    cfg = callbacks.GramTrackCacheConfig(
        :on,
        2.0,
        100.0,
        deg2rad(0.5),
        8,
        10.0,
        1000.0,
        deg2rad(4.0),
        32,
        20_000.0
    )
    @test callbacks._gram_track_cache_enabled(cfg, NoAtmosphereModel()) == false
    entry_profile = callbacks._gram_track_cache_profile(cfg, p_density_lookup, 130e3)
    orbit_profile = callbacks._gram_track_cache_profile(cfg, p_density_lookup, 200e3)
    @test entry_profile[1] == cfg.entry_horizon_s
    @test orbit_profile[1] == cfg.orbit_horizon_s

    @test callbacks._lerp(0.0, 10.0, 0.25) == 2.5
    @test isapprox(callbacks._angdiff_rad(0.0, 2pi - 0.1), 0.1; atol=1e-12, rtol=0.0)
    @test isapprox(callbacks._lerp_angle_rad(0.0, Float64(pi), 0.5), Float64(pi) / 2; atol=1e-12, rtol=0.0)

    cache = callbacks.GramTrackCache()
    cache.valid = true
    cache.t0 = 0.0
    cache.t1 = 10.0
    cache.index_hint = 1
    cache.times = [0.0, 5.0, 10.0]
    cache.alts = [1000.0, 2000.0, 3000.0]
    cache.lats = [0.0, 0.1, 0.2]
    cache.lons = [0.0, 0.1, 0.2]
    cache.rhos = [1.0, 2.0, 3.0]
    cache.Ts = [200.0, 220.0, 240.0]
    cache.winds = [SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(1.0, 0.0, 0.0), SVector{3, Float64}(2.0, 0.0, 0.0)]

    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "0") do
        @test callbacks._gram_track_cache_segment(cache, -1.0) === nothing
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "1") do
        seg_clamped = callbacks._gram_track_cache_segment(cache, -1.0)
        @test seg_clamped == (1, 0.0)
    end

    seg_hit = callbacks._gram_track_cache_ready(cache, 2.5, 1500.0, 0.05, 0.05, 1e-9, 1e-9)
    @test seg_hit == (1, 0.5)
    @test callbacks._gram_track_cache_ready(cache, 2.5, 2000.0, 0.05, 0.05, 1e-9, 1e-9) === nothing
    ρ_eval, T_eval, wind_eval = callbacks._gram_track_cache_eval(cache, 1, 0.25)
    @test isapprox(ρ_eval, 1.25; atol=1e-12, rtol=0.0)
    @test isapprox(T_eval, 205.0; atol=1e-12, rtol=0.0)
    @test isapprox(wind_eval[1], 0.25; atol=1e-12, rtol=0.0)

    cache_flat = callbacks.GramTrackCache()
    cache_flat.valid = true
    cache_flat.t0 = 2.0
    cache_flat.t1 = 2.0
    cache_flat.index_hint = 1
    cache_flat.times = [2.0, 2.0]
    cache_flat.alts = [1000.0, 1000.0]
    cache_flat.lats = [0.0, 0.0]
    cache_flat.lons = [0.0, 0.0]
    cache_flat.rhos = [1.0, 1.0]
    cache_flat.Ts = [200.0, 200.0]
    cache_flat.winds = [SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)]
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "1") do
        @test callbacks._gram_track_cache_segment(cache_flat, 2.0) == (1, 0.0)
    end

    @test callbacks._uses_j2_gravity_effector((InverseSquaredGravityModel(), InverseSquaredJ2GravityModel())) == true

    tol_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_quaternion=1e-9,
            abstol_quaternion=1e-9,
            reltol_angular_rate=4e-7,
            abstol_angular_rate=5e-8
        )
    )
    u_tol = build_initial_conditions(tol_args)
    template_reltol, template_abstol = _build_solver_tolerances(u_tol, tol_args)
    reltol_phase, abstol_phase = callbacks._callback_tolerances_for_phase(template_reltol, template_abstol, tol_args, false)
    @test all(isapprox.(reltol_phase.sc[1].q, tol_args.integration_tolerances.reltol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_phase.sc[1].q, tol_args.integration_tolerances.abstol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_phase.sc[1].ω, tol_args.integration_tolerances.reltol_angular_rate; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_phase.sc[1].ω, tol_args.integration_tolerances.abstol_angular_rate; atol=0.0, rtol=0.0))

    q_id = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    cache_lpi = PlanetFrameEphemerisCache(
        [0.0, 5.0, 10.0],
        [q_id, q_id, q_id]
    )
    @test callbacks._planet_lpi_from_cache(cache_lpi, -1.0) === nothing
    @test callbacks._planet_lpi_from_cache(cache_lpi, NaN) isa SMatrix{3, 3, Float64}
    @test callbacks._planet_lpi_from_cache(cache_lpi, 10.0) isa SMatrix{3, 3, Float64}
    cache_lpi_flip = PlanetFrameEphemerisCache(
        [0.0, 10.0],
        [q_id, -q_id]
    )
    @test callbacks._planet_lpi_from_cache(cache_lpi_flip, 5.0) isa SMatrix{3, 3, Float64}

    args_density_stats = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_density_stats = ODEParams{1}(args=args_density_stats)
    u_density_stats = build_initial_conditions(args_density_stats)
    withenv(
        "SPACEAGORA_GRAM_PROFILE" => "1",
        "SPACEAGORA_GRAM_TRACK_CACHE" => "off"
    ) do
        callbacks._gram_runtime_stats_reset!()
        density_cb_stats = callbacks.get_density_callback(1, args_density_stats)
        density_cb_stats.affect!((p=p_density_stats, u=u_density_stats, t=0.0))
        stats_after_density = callbacks._gram_runtime_stats_snapshot()
        @test stats_after_density.density_calls >= 1
        @test stats_after_density.direct_calls >= 1
    end

    args_thermal_branches = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_thermal_branches = ODEParams{1}(args=args_thermal_branches)
    u_thermal_branches = build_initial_conditions(args_thermal_branches)
    thermal_cb_branches = callbacks.get_thermal_callback(1, args_thermal_branches)

    p_thermal_branches.shared_buffers.heat_rates[1] = [1.0, 2.0, 3.0]
    p_thermal_branches.shared_buffers.densities[1] = 0.0
    p_thermal_branches.shared_buffers.temperatures[1] = 250.0
    thermal_cb_branches.affect!((p=p_thermal_branches, u=u_thermal_branches, t=0.0))
    @test length(p_thermal_branches.shared_buffers.heat_rates[1]) == length(args_thermal_branches.dynamics_model.spacecraft[1].links)

    p_thermal_branches.shared_buffers.densities[1] = 1e-6
    p_thermal_branches.shared_buffers.temperatures[1] = 250.0
    p_thermal_branches.shared_buffers.winds[1] = SVector{3, Float64}(NaN, NaN, NaN)
    thermal_cb_branches.affect!((p=p_thermal_branches, u=u_thermal_branches, t=0.0))

    p_thermal_branches.shared_buffers.winds[1] = SVector{3, Float64}(0.0, 0.0, 0.0)
    p_thermal_branches.shared_buffers.densities[1] = 1e-6
    p_thermal_branches.shared_buffers.temperatures[1] = 250.0
    args_thermal_branches.dynamics_model.spacecraft[1].links[1].α = NaN
    thermal_cb_branches.affect!((p=p_thermal_branches, u=u_thermal_branches, t=0.0))
    @test p_thermal_branches.shared_buffers.heat_rates[1][1] == 0.0

    mission_orbits = MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=true,
        number_of_orbits=1,
        mission_time=120.0,
        orientation_sim=false,
        num_steps_to_save=100
    )
    args_orbit_multi_base = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    args_orbit_multi = SimulationConfiguration(
        file_paths=args_orbit_multi_base.file_paths,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=mission_orbits,
        environment_model=args_orbit_multi_base.environment_model,
        dynamics_model=args_orbit_multi_base.dynamics_model,
        guidance_model=args_orbit_multi_base.guidance_model,
        navigation_model=args_orbit_multi_base.navigation_model,
        control_model=args_orbit_multi_base.control_model,
        initial_time=args_orbit_multi_base.initial_time,
        integration_tolerances=args_orbit_multi_base.integration_tolerances
    )
    p_orbit_multi = ODEParams{2}(args=args_orbit_multi)
    p_orbit_multi.orbit_counter .= [2, 1]
    p_orbit_multi.is_active .= [true, true]
    orbit_cb_multi = callbacks.get_orbit_end_callback(2)
    mutable struct OrbitStopIntegrator{P, U}
        p::P
        u::U
        t::Float64
        terminated::Bool
    end
    DiffEqBase.terminate!(integrator::OrbitStopIntegrator) = begin
        integrator.terminated = true
        nothing
    end
    orbit_integrator = OrbitStopIntegrator(
        p_orbit_multi,
        build_initial_conditions(args_orbit_multi),
        0.0,
        false
    )
    orbit_cb_multi.affect!(orbit_integrator, 1)
    @test orbit_integrator.terminated == false
    p_orbit_multi.orbit_counter .= [2, 2]
    orbit_cb_multi.affect!(orbit_integrator, 1)
    @test orbit_integrator.terminated == true
end

@testset "Coverage Threaded Probe Driver" begin
    if Base.JLOptions().code_coverage == 0
        @test true
    else
        probe_script = joinpath(REPO_ROOT, "test", "coverage_threaded_probes.jl")
        cmd = `$(Base.julia_cmd()) --startup-file=no --depwarn=error --project=$(REPO_ROOT) --code-coverage=user --threads=2 $(probe_script)`
        cmd = addenv(
            cmd,
            "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
            "SPACEAGORA_WARN_NORMALIZE" => "0"
        )

        output = IOBuffer()
        proc = run(pipeline(ignorestatus(cmd), stdout=output, stderr=output))
        text = String(take!(output))
        if !success(proc)
            println(text)
        end

        @test success(proc)
        @test occursin("coverage_threaded_probes_ok", text)
    end
end

@testset "Multibody Parallel Policy Gates" begin
    use_threads = SimulationModel.DynamicEffectors._multibody_use_threads
    has_worker_threads = Threads.nthreads() > 1

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"
    ) do
        @test use_threads(64) == false
    end

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0"
    ) do
        @test use_threads(64) == has_worker_threads
    end

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "on",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"
    ) do
        @test use_threads(64) == has_worker_threads
    end
end

@testset "Parallel Policy Adaptive Controller" begin
    policy = SimulationModel.ParallelPolicy
    policy.reset_policy_telemetry!()

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
        "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "1",
        "SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA" => "1",
        "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
        "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1"
    ) do
        @test policy.use_threads_policy(
            8;
            mode=:auto,
            threshold=1,
            source=:density_callback
        ) == false

        policy.record_policy_observation!(
            :density_callback;
            mode=:auto,
            num_items=1,
            use_threads=false,
            elapsed_ns=10
        )
        snap1 = policy.policy_telemetry_snapshot()
        @test snap1.last_classification == "efficient_satisfied"
        @test snap1.adaptation_updates_total >= 1
        @test snap1.last_desire >= 2

        policy.record_policy_observation!(
            :density_callback;
            mode=:auto,
            num_items=1,
            use_threads=false,
            elapsed_ns=11
        )
        snap2 = policy.policy_telemetry_snapshot()
        @test snap2.last_classification == "efficient_deprived"

        policy.record_policy_observation!(
            :density_callback;
            mode=:auto,
            num_items=0,
            use_threads=false,
            elapsed_ns=12
        )
        snap3 = policy.policy_telemetry_snapshot()
        @test snap3.last_classification == "inefficient"
        @test snap3.last_utilization <= 0.1
        @test snap3.serial_elapsed_ns_total >= 33
        @test snap3.quantum_length == 1
        @test snap3.trim_quanta_budget == 1
        @test snap3.quantums_total >= 3
        @test snap3.quantums_inefficient >= 1
        @test snap3.quantums_efficient_satisfied >= 1
        @test snap3.quantums_efficient_deprived >= 1
        @test snap3.accounted_fraction_proxy >= 0.0
        @test snap3.trimmed_accounted_fraction_proxy >= 0.0
    end

    if Threads.nthreads() > 1
        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "1",
            "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
            "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            use_history = Bool[]
            for _ in 1:8
                decision = policy.thread_policy_decision(
                    6;
                    mode=:auto,
                    threshold=1,
                    source=:other_source
                )
                push!(use_history, decision.use_threads)
                policy.record_policy_observation!(
                    :other_source;
                    mode=:auto,
                    num_items=6,
                    use_threads=decision.use_threads,
                    elapsed_ns=1
                )
            end
            @test any(use_history)
            @test use_history[end] == true
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "1",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == true
            @test decision.allotment >= 2
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "0",
            "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "0",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == false
            @test decision.allotment == 1
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "0",
            "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "1",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == true
            @test decision.allotment == min(Threads.nthreads(), 8)
        end
    end

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1"
    ) do
        @test policy.use_threads_policy(4; mode=:on, threshold=1, source=:other_source) == false
        policy.record_policy_observation!(
            :other_source;
            mode=:on,
            num_items=4,
            use_threads=true,
            elapsed_ns=UInt(5)
        )
        snap = policy.policy_telemetry_snapshot()
        @test snap.threaded_elapsed_ns_total >= 5
        @test snap.other_decisions >= 1
    end

    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "oops") do
        @test_throws ArgumentError policy.effective_inner_thread_budget()
    end
end

@testset "Parallel Policy threaded_reduce" begin
    policy = SimulationModel.ParallelPolicy
    budget = max(1, Threads.nthreads())
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(budget)) do
        reduced = policy.threaded_reduce(
            16,
            budget,
            () -> MVector{2, Int}(0, 0),
            (local_acc, idx) -> begin
                local_acc[1] += idx
                local_acc[2] += 1
                return nothing
            end,
            (dest, src) -> begin
                dest[1] += src[1]
                dest[2] += src[2]
                return nothing
            end
        )
        @test reduced[1] == sum(1:16)
        @test reduced[2] == 16

        reduced_empty = policy.threaded_reduce(
            0,
            budget,
            () -> Ref(0),
            (local_acc, idx) -> begin
                local_acc[] += idx
                return nothing
            end,
            (dest, src) -> begin
                dest[] += src[]
                return nothing
            end
        )
        @test reduced_empty[] == 0
    end
end

@testset "Parallel Policy threaded_foreach_persistent" begin
    policy = SimulationModel.ParallelPolicy
    budget = max(1, Threads.nthreads())

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1"
    ) do
        acc = Base.Threads.Atomic{Int}(0)
        policy.with_policy_context() do
            policy.threaded_foreach_persistent(:density_callback, 16, budget) do idx
                Base.Threads.atomic_add!(acc, idx)
            end
            policy.threaded_foreach_persistent(:density_callback, 16, budget) do idx
                Base.Threads.atomic_add!(acc, idx)
            end
        end
        @test acc[] == 2 * sum(1:16)
    end

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "0"
    ) do
        acc = Base.Threads.Atomic{Int}(0)
        policy.threaded_foreach_persistent(:density_callback, 8, budget) do idx
            Base.Threads.atomic_add!(acc, idx)
        end
        @test acc[] == sum(1:8)
    end

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1"
    ) do
        err = try
            policy.with_policy_context() do
                policy.threaded_foreach_persistent(:thermal_callback, 8, budget) do idx
                    if idx == 3
                        error("threaded_foreach_persistent_probe")
                    end
                end
            end
            nothing
        catch e
            e
        end
        @test err !== nothing
    end
end

@testset "Aerodynamic Helper Branch Coverage" begin
    dynamic_effectors = SimulationModel.DynamicEffectors

    @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false) == false
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "yes") do
        @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false) == true
    end
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "off") do
        @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", true) == false
    end
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "invalid") do
        @test_throws ArgumentError dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false)
    end

    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "off") do
        @test dynamic_effectors._multibody_parallel_mode() == :off
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
        @test dynamic_effectors._multibody_parallel_mode() == :on
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "auto") do
        @test dynamic_effectors._multibody_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "invalid") do
        @test_throws ArgumentError dynamic_effectors._multibody_parallel_mode()
    end

    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "4") do
        @test dynamic_effectors._multibody_thread_threshold() == 4
    end
    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "0") do
        @test dynamic_effectors._multibody_thread_threshold() == 1
    end
    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_thread_threshold()
    end
    withenv("SPACEAGORA_MULTIBODY_MAX_THREADS" => "3") do
        @test dynamic_effectors._multibody_max_threads() == 3
    end
    withenv("SPACEAGORA_MULTIBODY_MAX_THREADS" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_max_threads()
    end

    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
        @test dynamic_effectors._multibody_outer_parallel_hint() == true
    end
    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0") do
        @test dynamic_effectors._multibody_outer_parallel_hint() == false
    end
    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_outer_parallel_hint()
    end

    @test dynamic_effectors._multibody_use_threads(1) == false
    if Threads.nthreads() > 1
        withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
            @test dynamic_effectors._multibody_use_threads(64) == true
        end
        withenv(
            "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
            "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
            "SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER" => "0"
        ) do
            @test dynamic_effectors._multibody_use_threads(64) == false
        end
        withenv(
            "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
            "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
            "SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY" => "1"
        ) do
            @test dynamic_effectors._multibody_use_threads(64; heavy_work=false) == false
        end
    end

    @test dynamic_effectors._threadid_capacity() >= Threads.maxthreadid()

    body_a = Link{0}(root=true)
    body_b = Link{0}(root=false)
    body_a.net_force .= SVector{3, Float64}(1.0, 2.0, 3.0)
    body_a.net_torque .= SVector{3, Float64}(4.0, 5.0, 6.0)
    body_b.net_force .= SVector{3, Float64}(-0.5, 0.0, 0.5)
    body_b.net_torque .= SVector{3, Float64}(1.0, -1.0, 0.0)
    force_sum, torque_sum = dynamic_effectors.collect_and_reset_link_wrenches!([body_a, body_b])

    @test force_sum == SVector{3, Float64}(0.5, 2.0, 3.5)
    @test torque_sum == SVector{3, Float64}(5.0, 4.0, 6.0)
    @test body_a.net_force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_a.net_torque == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_b.net_force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_b.net_torque == SVector{3, Float64}(0.0, 0.0, 0.0)
end




@testset "Typed Planet Constructors + Topography Workspace" begin
    mars = Mars("", SPICE_PATH)
    venus = Venus("", SPICE_PATH)
    titan = Titan("", SPICE_PATH)

    @test mars.name == "Mars"
    @test venus.name == "Venus"
    @test titan.name == "Titan"
    @test mars.μ > 0.0
    @test venus.μ > 0.0
    @test titan.μ > 0.0

    mktempdir() do tmp
        topo_file = joinpath(tmp, "topo_harmonics.csv")
        write(topo_file, "degree,order,C,S\n0,0,1.0,0.0\n1,0,0.1,0.0\n1,1,0.05,0.02\n")
        earth = Earth("", SPICE_PATH)
        SimulationModel.Planets.TopographyHarmonicsWorkspace!(topo_file, earth)
        @test size(earth.topography_workspace.Clm) == (2, 2)
        @test size(earth.topography_workspace.Slm) == (2, 2)
        @test isapprox(earth.topography_workspace.Clm[2, 2], 0.05; atol=0.0, rtol=0.0)
        @test isapprox(earth.topography_workspace.Slm[2, 2], 0.02; atol=0.0, rtol=0.0)
    end
end

@testset "Run Simulation State Isolation" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=45.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    l_pi_before = args.environment_model.planet.L_PI
    @test_nowarn run_simulation(args)
    @test args.environment_model.planet.L_PI == l_pi_before

    args_no_isolation = deepcopy(args)
    @test_nowarn run_simulation(args_no_isolation; isolate_state=false)

    if Threads.nthreads() >= 2
        args_parallel = deepcopy(args)
        t1 = Threads.@spawn run_simulation(args_parallel)
        t2 = Threads.@spawn run_simulation(args_parallel)
        @test_nowarn fetch(t1)
        @test_nowarn fetch(t2)
    end
end


@testset "Deterministic Smoke + No-Drag Energy Invariant" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=1200.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df = run_case(args)
    @test nrow(df) > 100
    @test all(isfinite, df.sc1_pos_1)
    @test all(isfinite, df.sc1_vel_1)

    eps = specific_energy(df, EARTH.μ)
    energy_drift = maximum(abs.(eps .- first(eps))) / abs(first(eps))
    @test energy_drift < 1e-8
end

@testset "Deterministic Replay (No-Drag)" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df1 = run_case_silent(args)
    df2 = run_case_silent(args)
    @test nrow(df1) == nrow(df2)
    @test nrow(df1) > 10

    sample_idxs = round.(Int, range(1, nrow(df1), length=8))
    for idx in sample_idxs
        t = Float64(df1.time[idx])
        p1 = SVector{3, Float64}(Float64(df1.sc1_pos_1[idx]), Float64(df1.sc1_pos_2[idx]), Float64(df1.sc1_pos_3[idx]))
        v1 = SVector{3, Float64}(Float64(df1.sc1_vel_1[idx]), Float64(df1.sc1_vel_2[idx]), Float64(df1.sc1_vel_3[idx]))
        p2 = SVector{3, Float64}(
            interp_linear(df2.time, df2.sc1_pos_1, t),
            interp_linear(df2.time, df2.sc1_pos_2, t),
            interp_linear(df2.time, df2.sc1_pos_3, t)
        )
        v2 = SVector{3, Float64}(
            interp_linear(df2.time, df2.sc1_vel_1, t),
            interp_linear(df2.time, df2.sc1_vel_2, t),
            interp_linear(df2.time, df2.sc1_vel_3, t)
        )
        @test norm(p1 - p2) < 0.1
        @test norm(v1 - v2) < 1e-4
    end
end


@testset "Typed Results Bundle Persistence" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    mktempdir() do tmp
        cd(tmp) do
            run_simulation(args)
            @test isfile(joinpath("output", "simulation_results.csv"))

            bundle_prefix = joinpath("output", "simulation_results")
            feather_path = bundle_prefix * ".feather"
            bundle_csv_path = bundle_prefix * ".csv"
            manifest_path = bundle_prefix * ".manifest.toml"

            @test isfile(feather_path)
            @test isfile(bundle_csv_path)
            @test isfile(manifest_path)

            manifest = TOML.parsefile(manifest_path)
            @test get(manifest, "schema_version", "") == "1"
            @test get(manifest, "steps", 0) > 10
            @test haskey(manifest, "files")

            files = manifest["files"]
            @test haskey(files, "feather")
            @test haskey(files, "csv")
            @test get(files["feather"], "size_bytes", 0) > 0
            @test length(get(files["feather"], "sha256", "")) == 64
        end
    end
end

@testset "IOConfig Owner Coverage" begin
    mktempdir() do tmp
        results_dir = joinpath(tmp, "results")
        mkpath(results_dir)

        args_default_ckpt = (; simulation_settings=(; results_directory=results_dir, checkpoint_directory=""))
        args_custom_ckpt = (; simulation_settings=(; results_directory=results_dir, checkpoint_directory=joinpath(tmp, "custom_ckpt")))

        bundle_prefix = SimulationModel.IOConfig._results_bundle_prefix(args_default_ckpt)
        @test bundle_prefix == joinpath(results_dir, "simulation_results")

        compat_csv = SimulationModel.IOConfig._compat_results_csv_path(args_default_ckpt)
        @test compat_csv == joinpath(results_dir, "simulation_results.csv")

        collision_csv = SimulationModel.IOConfig._collision_results_csv_path(args_default_ckpt)
        @test startswith(collision_csv, joinpath(results_dir, "simulation_results."))
        @test endswith(collision_csv, ".csv")

        @test SimulationModel.IOConfig._checkpoint_directory(args_default_ckpt) == joinpath(results_dir, "checkpoints")
        @test SimulationModel.IOConfig._checkpoint_directory(args_custom_ckpt) == joinpath(tmp, "custom_ckpt")

        paths_default = SimulationModel.IOConfig._checkpoint_paths(args_default_ckpt)
        @test paths_default.data == joinpath(results_dir, "checkpoints", "simulation_checkpoint.bin")
        @test paths_default.manifest == joinpath(results_dir, "checkpoints", "simulation_checkpoint.manifest.toml")

        paths_custom = SimulationModel.IOConfig._checkpoint_paths(args_custom_ckpt)
        @test paths_custom.data == joinpath(tmp, "custom_ckpt", "simulation_checkpoint.bin")
        @test paths_custom.manifest == joinpath(tmp, "custom_ckpt", "simulation_checkpoint.manifest.toml")
    end
end

@testset "Checkpoint Resume Parity" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    baseline_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    baseline_args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=240.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=baseline_settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    df_full = run_case_silent(baseline_args)

    checkpoint_dir = joinpath("output", "checkpoints")
    checkpoint_data_path = joinpath(checkpoint_dir, "simulation_checkpoint.bin")
    checkpoint_manifest_path = joinpath(checkpoint_dir, "simulation_checkpoint.manifest.toml")

    phase1_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false,
        checkpoint_enabled=true,
        checkpoint_interval_s=40.0,
        checkpoint_directory=checkpoint_dir,
        resume_from_checkpoint=false
    )
    phase1_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=phase1_settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    resume_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false,
        checkpoint_enabled=true,
        checkpoint_interval_s=40.0,
        checkpoint_directory=checkpoint_dir,
        resume_from_checkpoint=true
    )
    resume_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=240.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=resume_settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    mktempdir() do tmp
        cd(tmp) do
            run_simulation(phase1_args)
            @test isfile(checkpoint_data_path)
            @test isfile(checkpoint_manifest_path)

            run_simulation(resume_args)
            @test isfile(checkpoint_data_path)
            @test isfile(checkpoint_manifest_path)
            @test isfile(joinpath("output", "simulation_results.csv"))

            df_resume = CSV.read(joinpath("output", "simulation_results.csv"), DataFrame)
            @test nrow(df_resume) >= 2
            @test issorted(df_resume.time)
            @test Float64(df_resume.time[1]) > 0.0
            @test Float64(df_resume.time[end]) > Float64(df_resume.time[1])
            @test abs(Float64(df_resume.time[end]) - Float64(df_full.time[end])) < 1e-8

            p_full = SVector{3, Float64}(Float64(df_full.sc1_pos_1[end]), Float64(df_full.sc1_pos_2[end]), Float64(df_full.sc1_pos_3[end]))
            v_full = SVector{3, Float64}(Float64(df_full.sc1_vel_1[end]), Float64(df_full.sc1_vel_2[end]), Float64(df_full.sc1_vel_3[end]))
            p_resume = SVector{3, Float64}(Float64(df_resume.sc1_pos_1[end]), Float64(df_resume.sc1_pos_2[end]), Float64(df_resume.sc1_pos_3[end]))
            v_resume = SVector{3, Float64}(Float64(df_resume.sc1_vel_1[end]), Float64(df_resume.sc1_vel_2[end]), Float64(df_resume.sc1_vel_3[end]))

            @test norm(p_resume - p_full) < 1.0
            @test norm(v_resume - v_full) < 1e-3
            @test abs(Float64(df_resume.sc1_mass[end]) - Float64(df_full.sc1_mass[end])) < 1e-8
        end
    end
end

@testset "Checkpoint Guards + Missing Resume + Bundle Toggle" begin
    base_spacecraft() = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    base_settings(; kwargs...) = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false;
        kwargs...
    )

    args_bad_interval = build_config(
        spacecraft=base_spacecraft(),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=base_settings(
            checkpoint_enabled=true,
            checkpoint_interval_s=0.0,
            checkpoint_directory=joinpath("output", "checkpoints")
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    @test_throws ArgumentError run_simulation(args_bad_interval)

    args_resume_missing = build_config(
        spacecraft=base_spacecraft(),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=base_settings(
            checkpoint_enabled=true,
            checkpoint_interval_s=20.0,
            checkpoint_directory=joinpath("output", "checkpoints"),
            resume_from_checkpoint=true
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    mktempdir() do tmp
        cd(tmp) do
            @test_logs (:warn, r"resume_from_checkpoint=true but no checkpoint file was found") run_simulation(args_resume_missing)
            @test isfile(joinpath("output", "simulation_results.csv"))

            df = CSV.read(joinpath("output", "simulation_results.csv"), DataFrame)
            @test nrow(df) > 10
            @test abs(Float64(df.time[end]) - 60.0) < 1e-8

            ckpt_manifest_path = joinpath("output", "checkpoints", "simulation_checkpoint.manifest.toml")
            @test isfile(ckpt_manifest_path)
            ckpt_manifest = TOML.parsefile(ckpt_manifest_path)
            @test get(ckpt_manifest, "schema_version", "") == "1"
            @test get(ckpt_manifest, "time_s", 0.0) > 0.0
            @test get(ckpt_manifest, "data_size_bytes", 0) > 0
            @test length(get(ckpt_manifest, "data_sha256", "")) == 64
        end
    end

    args_bundle_disabled = build_config(
        spacecraft=base_spacecraft(),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=base_settings(),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    mktempdir() do tmp
        cd(tmp) do
            withenv("SPACEAGORA_SAVE_BUNDLE" => "0") do
                run_simulation(args_bundle_disabled)
            end
            @test isfile(joinpath("output", "simulation_results.csv"))
            @test !isfile(joinpath("output", "simulation_results.feather"))
            @test !isfile(joinpath("output", "simulation_results.manifest.toml"))
        end
    end
end


@testset "Units/Normalization Consistency Audit" begin
    function strip_comments(src::String)
        no_block = replace(src, r"#=.*?=#"s => "")
        no_line = map(line -> first(split(line, '#'; limit=2)), split(no_block, '\n'; keepempty=true))
        return join(no_line, "\n")
    end

    run_src = strip_comments(read(joinpath(REPO_ROOT, "src", "simulation", "engine", "execution.jl"), String))
    legacy_elements_path = joinpath(REPO_ROOT, "src", "simulation", "execution", "simulation_elements.jl")

    # Canonical engine path should stay SI-native.
    @test !occursin("cnf.DU", run_src)
    @test !occursin("cnf.TU", run_src)
    @test !occursin("cnf.MU", run_src)
    @test !isfile(legacy_elements_path)

    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    settings_norm_true = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=true
    )
    settings_norm_false = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=false
    )

    args_norm_true = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_norm_true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_norm_false = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_norm_false,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    # Build-state sanity: typed initial conditions are SI-scale, not O(1) normalized values.
    u0 = build_initial_conditions(args_norm_true)
    @test norm(SVector{3, Float64}(u0.sc[1].pos)) > 1e6
    @test norm(SVector{3, Float64}(u0.sc[1].vel)) > 1e3
    @test u0.sc[1].mass > 1.0

    @test_throws ArgumentError run_case_silent(args_norm_true)

    df_norm_true = withenv(
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "1",
        "SPACEAGORA_WARN_NORMALIZE" => "0"
    ) do
        run_case_silent(args_norm_true)
    end
    df_norm_false = run_case_silent(args_norm_false)

    @test nrow(df_norm_true) == nrow(df_norm_false)
    @test nrow(df_norm_true) > 10
    sample_idxs = round.(Int, range(1, nrow(df_norm_true), length=8))
    for idx in sample_idxs
        t = Float64(df_norm_true.time[idx])
        p_true = SVector{3, Float64}(Float64(df_norm_true.sc1_pos_1[idx]), Float64(df_norm_true.sc1_pos_2[idx]), Float64(df_norm_true.sc1_pos_3[idx]))
        v_true = SVector{3, Float64}(Float64(df_norm_true.sc1_vel_1[idx]), Float64(df_norm_true.sc1_vel_2[idx]), Float64(df_norm_true.sc1_vel_3[idx]))
        p_false = SVector{3, Float64}(
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_1, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_2, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_3, t)
        )
        v_false = SVector{3, Float64}(
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_1, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_2, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_3, t)
        )
        @test norm(p_true - p_false) < 0.1
        @test norm(v_true - v_false) < 1e-4
    end
end

@testset "Normalize Flag Runtime Policy" begin
    args_warn = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=true),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    withenv("SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "0") do
        @test_throws ArgumentError run_simulation(args_warn)
    end

    withenv(
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "1",
        "SPACEAGORA_WARN_NORMALIZE" => "1"
    ) do
        _normalize_warning_emitted[] = false
        @test_logs (:warn, r"normalize=true is legacy-only") run_simulation(args_warn)
        @test _normalize_warning_emitted[] == true
        @test_logs run_simulation(args_warn)
    end

    withenv(
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "1",
        "SPACEAGORA_WARN_NORMALIZE" => "0"
    ) do
        _normalize_warning_emitted[] = false
        @test_logs run_simulation(args_warn)
        @test _normalize_warning_emitted[] == false
    end
end

@testset "Run Simulation Debug Branches" begin
    debug_thruster = TimedTangentialThrusterModel(5.0, 1.0, 0.0, 120.0)
    args_debug_control = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(debug_thruster,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    withenv("SPACEAGORA_DEBUG_CONTROL" => "1") do
        _, output = run_case_capture_stdout(args_debug_control)
        @test occursin("Applying control effect for spacecraft 1", output)
        @test occursin("Control force:", output)
    end

    args_debug_throw = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(ThrowingForceModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    withenv("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE" => "1") do
        @test_logs (:error, r"The derivative function itself crashed!") begin
            @test_throws ErrorException run_simulation(args_debug_throw)
        end
    end

    args_debug_nan = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(NaNForceModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    nan_output = ""
    withenv("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE" => "1") do
        mktemp() do path, io
            redirect_stdout(io) do
                @test_throws ErrorException run_simulation(args_debug_nan)
            end
            flush(io)
            seekstart(io)
            nan_output = read(io, String)
        end
    end
    @test occursin("--- INITIAL NaN DETECTED ---", nan_output)
    @test occursin("NaN found in Satellite 1 derivative!", nan_output)
    @test occursin("Pos:", nan_output)
    @test occursin("Vel:", nan_output)

    args_debug_nan_param = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(NaNParamForceModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    nan_param_output = ""
    withenv("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE" => "1") do
        mktemp() do path, io
            redirect_stdout(io) do
                @test_throws ErrorException run_simulation(args_debug_nan_param)
            end
            flush(io)
            seekstart(io)
            nan_param_output = read(io, String)
        end
    end
    @test occursin("NaN found in parameter: p.shared_buffers.current_time[]", nan_param_output)

    nan_scan_output = ""
    mktemp() do path, io
        redirect_stdout(io) do
            _debug_print_nan_parameter_paths!(NaN, "p.scalar_test")
            _debug_print_nan_parameter_paths!([1.0, NaN], "p.array_test")
        end
        flush(io)
        seekstart(io)
        nan_scan_output = read(io, String)
    end
    @test occursin("NaN found in parameter: p.scalar_test", nan_scan_output)
    @test occursin("NaN found in parameter: p.array_test[2]", nan_scan_output)
end

@testset "Verbose Gating for Callback/Runtime Logs" begin
    settings_quiet = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    settings_verbose = SimulationSettings(
        results=true,
        verbose=true,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )

    args_orbit_quiet = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_quiet,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    _, out_orbit_quiet = run_case_capture_stdout(args_orbit_quiet)
    @test !occursin("Initial conditions:", out_orbit_quiet)
    @test !occursin("Orbit ", out_orbit_quiet)

    args_orbit_verbose = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_verbose,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    _, out_orbit_verbose = run_case_capture_stdout(args_orbit_verbose)
    @test occursin("Initial conditions:", out_orbit_verbose)

    args_drag_quiet = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=settings_quiet,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    _, out_drag_quiet = run_case_capture_stdout(args_drag_quiet)
    @test !occursin("Switching to atmosphere integration", out_drag_quiet)
    @test !occursin("Impact detected", out_drag_quiet)

    args_drag_verbose = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=settings_verbose,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    _, out_drag_verbose = run_case_capture_stdout(args_drag_verbose)
    @test occursin("Switching to atmosphere integration", out_drag_verbose)
    @test occursin("Impact detected", out_drag_verbose)
end

@testset "Orbital Elements Round-Trip Invariant" begin
    ic = InitialCondition(
        ra=EARTH.Rp_e + 550e3,
        rp=EARTH.Rp_e + 300e3,
        i=37.0,
        ω=55.0,
        Ω=25.0,
        ν=130.0
    )
    r, v = orbitalelemtorv(ic, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)

    @test isapprox(oe[1], ic.a; rtol=1e-10, atol=1e-6)
    @test isapprox(oe[2], ic.e; rtol=1e-10, atol=1e-10)
    @test angle_distance(oe[3], ic.i) < 1e-8
    @test angle_distance(oe[4], ic.Ω) < 1e-8
    @test angle_distance(oe[5], ic.ω) < 1e-8
    @test angle_distance(oe[6], ic.ν) < 1e-8
end

@testset "Reference System Branch Coverage" begin
    ic_arg_wrap = InitialCondition(
        ra=EARTH.Rp_e + 700e3,
        rp=EARTH.Rp_e + 300e3,
        i=55.0,
        ω=300.0,
        Ω=40.0,
        ν=25.0
    )
    r_wrap, v_wrap = orbitalelemtorv(ic_arg_wrap, EARTH)
    oe_wrap = rvtoorbitalelement(SVector{3, Float64}(r_wrap), SVector{3, Float64}(v_wrap), EARTH)
    @test oe_wrap[5] > pi

    ic_circ_incl = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 500e3,
        i=40.0,
        ω=0.0,
        Ω=15.0,
        ν=250.0
    )
    r_circ, v_circ = orbitalelemtorv(ic_circ_incl, EARTH)
    oe_circ = rvtoorbitalelement(SVector{3, Float64}(r_circ), SVector{3, Float64}(v_circ), EARTH)
    @test abs(oe_circ[2]) < 1e-10
    @test oe_circ[6] > pi

    planet_topo = (
        Rp_e=10.0,
        Rp_p=9.0,
        topography_function=(a, Clm, Slm, lat, lon, A) -> 9.5,
        Clm_topo=zeros(1, 1),
        Slm_topo=zeros(1, 1),
        A_topo=0.0
    )
    rp_topo = SVector{3, Float64}(0.1, 0.5, -8.944723618090451)
    lla = rtolatlong(rp_topo, planet_topo)
    @test all(isfinite, lla)

    topo_module_name = gensym(:ReferenceSystemTopoSandbox)
    Core.eval(Main, :(module $topo_module_name
        using LinearAlgebra
        using StaticArrays
        args = :legacy_topography_args
        include(joinpath(Main.REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))
    end))
    topo_sandbox = getfield(Main, topo_module_name)
    lla_topo = topo_sandbox.rtolatlong(rp_topo, planet_topo, true)
    @test isapprox(lla_topo[1], norm(rp_topo) - 9.5; atol=1e-12, rtol=0.0)

    q_branch2 = orbital_elements_to_lvlh_quaternion(
        1.186823891356144,
        0.13305160284932352,
        0.0,
        2.0245819323134224
    )
    q_branch3 = orbital_elements_to_lvlh_quaternion(
        0.0,
        0.1,
        0.0,
        0.10471975511965977
    )
    @test isapprox(norm(q_branch2), 1.0; atol=1e-12, rtol=0.0)
    @test isapprox(norm(q_branch3), 1.0; atol=1e-12, rtol=0.0)
end

@testset "Circular Orbit Robustness" begin
    ic_eq = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 500e3,
        i=0.0,
        ω=0.0,
        Ω=0.0,
        ν=30.0
    )
    r, v = orbitalelemtorv(ic_eq, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)
    @test all(isfinite, oe)
    @test abs(oe[2]) < 1e-10

    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=0.0,
        ω_deg=0.0,
        Ω_deg=0.0,
        ν_deg=30.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=300.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    df = run_case(args)
    @test nrow(df) > 10
end

@testset "Quaternion Norm Invariant" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=420e3,
        i_deg=40.0,
        ω_deg=10.0,
        Ω_deg=20.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            reltol_quaternion=1e-9,
            abstol_quaternion=1e-9,
            dt_max_orbit=5.0
        )
    )

    df = run_case(args)
    qnorm = sqrt.(df.sc1_q_1.^2 .+ df.sc1_q_2.^2 .+ df.sc1_q_3.^2 .+ df.sc1_q_4.^2)
    @test maximum(abs.(qnorm .- 1.0)) < 1e-6
end

@testset "DynamicsRotational Core Contracts" begin
    rotational = SimulationModel.DynamicsRotational

    q = normalize(SVector{4, Float64}(0.15, -0.25, 0.35, 0.85))
    ω = SVector{3, Float64}(0.02, -0.03, 0.01)

    qdot_new = rotational.quaternion_derivative(ω, q)
    qdot_legacy = 0.5 * SVector{4, Float64}(SimulationModel.quat_mult(SVector{4, Float64}(ω..., 0.0), q))
    @test isapprox(qdot_new, qdot_legacy; atol=1e-12, rtol=1e-12)
    @test abs(dot(q, qdot_new)) < 1e-12

    J = SMatrix{3, 3, Float64}(2.0, 0.1, 0.0, 0.1, 3.0, 0.0, 0.0, 0.0, 4.0)
    τ = SVector{3, Float64}(0.12, -0.08, 0.05)
    ωdot_new = rotational.angular_acceleration(ω, J, τ; include_gyroscopic=true)
    ωdot_legacy = J \ (τ - cross(ω, J * ω))
    @test isapprox(ωdot_new, ωdot_legacy; atol=1e-12, rtol=1e-12)

    ωdot_no_gyro = rotational.angular_acceleration(ω, J, τ; include_gyroscopic=false)
    @test isapprox(ωdot_no_gyro, J \ τ; atol=1e-12, rtol=1e-12)

    @test rotational.body_torque(MVector{3, Float64}(τ)) == τ
    @test rotational.body_angular_velocity(MVector{3, Float64}(ω)) == ω
    @test rotational.combine_torques(τ, -τ) == SVector{3, Float64}(0.0, 0.0, 0.0)

    J_spherical = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0)
    ω_any = SVector{3, Float64}(0.3, -0.25, 0.11)
    ωdot_torque_free = rotational.angular_acceleration(
        ω_any,
        J_spherical,
        SVector{3, Float64}(0.0, 0.0, 0.0);
        include_gyroscopic=true
    )
    @test isapprox(ωdot_torque_free, SVector{3, Float64}(0.0, 0.0, 0.0); atol=1e-12, rtol=0.0)
end

@testset "Solver Tolerances Apply Quaternion Overrides" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=420e3,
        i_deg=40.0,
        ω_deg=10.0,
        Ω_deg=20.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    tols = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_quaternion=3e-7,
        abstol_quaternion=4e-8,
        reltol_mass=5e-7,
        abstol_mass=6e-8,
        reltol_heat_load=7e-7,
        abstol_heat_load=8e-8,
        reltol_angular_rate=9e-7,
        abstol_angular_rate=1e-7,
        dt_max_orbit=5.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols
    )

    u0 = build_initial_conditions(args)
    reltol_vec, abstol_vec = _build_solver_tolerances(u0, args)
    @test all(isapprox.(reltol_vec.sc[1].pos, tols.reltol_orbit; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].pos, tols.abstol_orbit; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_vec.sc[1].q, tols.reltol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].q, tols.abstol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_vec.sc[1].ω, tols.reltol_angular_rate; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].ω, tols.abstol_angular_rate; atol=0.0, rtol=0.0))
    @test isapprox(reltol_vec.sc[1].mass, tols.reltol_mass; atol=0.0, rtol=0.0)
    @test isapprox(abstol_vec.sc[1].mass, tols.abstol_mass; atol=0.0, rtol=0.0)
    @test all(isapprox.(reltol_vec.sc[1].heat_loads, tols.reltol_heat_load; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].heat_loads, tols.abstol_heat_load; atol=0.0, rtol=0.0))

    tols_no_orient = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=0.0,
        abstol_mass=0.0,
        reltol_heat_load=0.0,
        abstol_heat_load=0.0,
        reltol_angular_rate=0.0,
        abstol_angular_rate=0.0,
        dt_max_orbit=5.0
    )
    args_no_orient = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_no_orient
    )
    u0_no_orient = build_initial_conditions(args_no_orient)
    reltol_scalar, abstol_scalar = _build_solver_tolerances(u0_no_orient, args_no_orient)
    @test reltol_scalar == tols_no_orient.reltol_orbit
    @test abstol_scalar == tols_no_orient.abstol_orbit

    tols_no_orient_component = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=3e-7,
        abstol_mass=4e-8,
        reltol_heat_load=5e-7,
        abstol_heat_load=6e-8,
        dt_max_orbit=5.0
    )
    args_no_orient_component = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_no_orient_component
    )
    u0_no_orient_component = build_initial_conditions(args_no_orient_component)
    reltol_vec_no_orient, abstol_vec_no_orient = _build_solver_tolerances(u0_no_orient_component, args_no_orient_component)
    @test isapprox(reltol_vec_no_orient.sc[1].mass, tols_no_orient_component.reltol_mass; atol=0.0, rtol=0.0)
    @test isapprox(abstol_vec_no_orient.sc[1].mass, tols_no_orient_component.abstol_mass; atol=0.0, rtol=0.0)
    @test all(isapprox.(reltol_vec_no_orient.sc[1].heat_loads, tols_no_orient_component.reltol_heat_load; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec_no_orient.sc[1].heat_loads, tols_no_orient_component.abstol_heat_load; atol=0.0, rtol=0.0))

    tols_ω_no_orient = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=0.0,
        abstol_mass=0.0,
        reltol_heat_load=0.0,
        abstol_heat_load=0.0,
        reltol_angular_rate=3e-7,
        dt_max_orbit=5.0
    )
    args_ω_no_orient = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_ω_no_orient
    )
    u0_ω_no_orient = build_initial_conditions(args_ω_no_orient)
    reltol_ω_no_orient, abstol_ω_no_orient = _build_solver_tolerances(u0_ω_no_orient, args_ω_no_orient)
    @test reltol_ω_no_orient == tols_ω_no_orient.reltol_orbit
    @test abstol_ω_no_orient == tols_ω_no_orient.abstol_orbit
end

@testset "Solver/Env Helper Parsing Coverage" begin
    withenv("SPACEAGORA_SOLVER_MODE" => nothing) do
        @test _solver_policy_mode() == :tsit5
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "auto") do
        @test _solver_policy_mode() == :auto_stiff
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "rodas") do
        @test _solver_policy_mode() == :rodas5p
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "split_imex") do
        @test _solver_policy_mode() == :split_imex
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "unsupported-mode") do
        @test_throws ArgumentError _solver_policy_mode()
    end

    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "on") do
        @test _gram_per_sat_instances_enabled() == true
    end
    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "off") do
        @test _gram_per_sat_instances_enabled() == false
    end
    withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES" => "maybe") do
        @test_throws ArgumentError _gram_per_sat_instances_enabled()
    end

    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "off") do
        @test _effector_parallel_mode() == :off
    end
    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "on") do
        @test _effector_parallel_mode() == :on
    end
    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "auto") do
        @test _effector_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "invalid") do
        @test_throws ArgumentError _effector_parallel_mode()
    end
    withenv("SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "3") do
        @test _effector_thread_threshold() == 3
    end
    withenv("SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError _effector_thread_threshold()
    end
    withenv("SPACEAGORA_EFFECTOR_MAX_THREADS" => "3") do
        @test _effector_max_threads() == 3
    end
    withenv("SPACEAGORA_EFFECTOR_MAX_THREADS" => "oops") do
        @test_throws ArgumentError _effector_max_threads()
    end
    withenv("SPACEAGORA_EFFECTOR_LONG_MISSION_THRESHOLD_S" => "10.0") do
        @test _effector_long_mission_threshold_s() == 10.0
    end
    withenv("SPACEAGORA_EFFECTOR_LONG_MISSION_THRESHOLD_S" => "0.0") do
        @test_throws ArgumentError _effector_long_mission_threshold_s()
    end
    withenv("SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "12345.0") do
        @test _effector_cost_ns_per_item_default() == 12345.0
    end
    withenv("SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "0.0") do
        @test_throws ArgumentError _effector_cost_ns_per_item_default()
    end
    withenv("SPACEAGORA_EFFECTOR_COST_EMA_ALPHA" => "0.3") do
        @test _effector_cost_ema_alpha() == 0.3
    end
    withenv("SPACEAGORA_EFFECTOR_COST_EMA_ALPHA" => "1.5") do
        @test_throws ArgumentError _effector_cost_ema_alpha()
    end
    withenv("SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "25000.0") do
        @test _effector_work_ns_per_worker_threshold() == 25000.0
    end
    withenv("SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "0.0") do
        @test_throws ArgumentError _effector_work_ns_per_worker_threshold()
    end

    args_eff_single = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), InverseSquaredJ2GravityModel()),
        keplerian=true
    )
    p_eff_single = ODEParams{1}(args=args_eff_single)
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "1",
        "SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "120000.0",
        "SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES" => "999",
        "SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "1000.0",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0"
    ) do
        decision_single = _dynamic_effector_thread_decision(args_eff_single, p_eff_single, args_eff_single.dynamics_model.dynamic_effectors, 1)
        if Threads.nthreads() > 1
            @test decision_single.use_threads == true
            @test decision_single.allotment >= 2
            @test decision_single.allotment <= min(Threads.nthreads(), 4)
        else
            @test decision_single.use_threads == false
        end
    end

    args_eff_multi = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), InverseSquaredJ2GravityModel()),
        keplerian=true
    )
    p_eff_multi = ODEParams{1}(args=args_eff_multi)
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "1",
        "SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "1.0",
        "SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES" => "999",
        "SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "1.0e9",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0"
    ) do
        decision_multi = _dynamic_effector_thread_decision(args_eff_multi, p_eff_multi, args_eff_multi.dynamics_model.dynamic_effectors, 1)
        @test decision_multi.use_threads == false
    end

    args_eff_constellation = args_eff_single
    p_eff_constellation = ODEParams{1}(args=args_eff_constellation)
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "1",
        "SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT" => "120000.0",
        "SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES" => "999",
        "SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD" => "1000.0",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_EFFECTOR_MAX_THREADS" => string(Threads.nthreads())
    ) do
        decision_constellation = _dynamic_effector_thread_decision(
            args_eff_constellation,
            p_eff_constellation,
            args_eff_constellation.dynamics_model.dynamic_effectors,
            4
        )
        share_budget = max(1, fld(max(1, Threads.nthreads()), min(4, max(1, Threads.nthreads()))))
        inner_floor = Threads.nthreads() > 1 ? min(2, Threads.nthreads()) : 1
        expected_cap = min(Threads.nthreads(), max(share_budget, inner_floor))
        if expected_cap > 1
            @test decision_constellation.use_threads == true
            @test decision_constellation.allotment <= expected_cap
        else
            @test decision_constellation.use_threads == false
            @test decision_constellation.allotment == 1
        end
    end

    @test _retcode_is_stiff_symptom(:Unstable)
    @test _retcode_is_stiff_symptom("DtLessThanMin")
    @test _retcode_is_stiff_symptom(:InitialFailure)
    @test !_retcode_is_stiff_symptom(:Success)

    withenv("SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        @test _solver_maxiters() === nothing
    end
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "2500") do
        @test _solver_maxiters() == 2500
    end
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "0") do
        @test_throws ArgumentError _solver_maxiters()
    end
    withenv("SPACEAGORA_SOLVER_MAXITERS" => "not-an-int") do
        @test_throws ArgumentError _solver_maxiters()
    end

    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => nothing) do
        split_solver = _split_imex_solver_spec()
        @test split_solver.label == "KenCarp4"
    end
    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp47") do
        split_solver = _split_imex_solver_spec()
        @test split_solver.label == "KenCarp47"
    end
    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp58") do
        split_solver = _split_imex_solver_spec()
        @test split_solver.label == "KenCarp58"
    end
    withenv("SPACEAGORA_SPLIT_IMEX_SOLVER" => "unsupported-split-solver") do
        @test_throws ArgumentError _split_imex_solver_spec()
    end

    struct DummyNoAlgChoice end
    struct DummyAlgChoice
        alg_choice::Vector{Int}
    end
    @test _auto_stiff_switched(DummyNoAlgChoice()) == false
    @test _auto_stiff_switched(DummyAlgChoice(Int[])) == false
    @test _auto_stiff_switched(DummyAlgChoice([1, 1, 1])) == false
    @test _auto_stiff_switched(DummyAlgChoice([1, 2, 1])) == true

    solver_args = args_eff_single
    prob_simple = ODEProblem(
        (du, u, p, t) -> begin
            du[1] = -u[1]
        end,
        [1.0],
        (0.0, 1.0)
    )

    withenv("SPACEAGORA_SOLVER_MODE" => "tsit5", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "Tsit5"
        @test meta.fallback_used == false
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "rodas5p", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "Rodas5P"
        @test meta.fallback_used == false
    end
    withenv("SPACEAGORA_SOLVER_MODE" => "auto_stiff", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "AutoTsit5(Rodas5P)"
        @test meta.initial_solver == "AutoTsit5"
    end
    split_prob_simple = SplitODEProblem(
        (du, u, p, t) -> begin
            du[1] = -u[1]
        end,
        (du, u, p, t) -> begin
            du[1] = -2u[1]
        end,
        [1.0],
        (0.0, 1.0)
    )
    withenv("SPACEAGORA_SOLVER_MODE" => "split_imex", "SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol, meta = _solve_with_solver_policy(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "KenCarp4(IMEX)"
        @test meta.initial_solver == "KenCarp4"
    end
    withenv(
        "SPACEAGORA_SOLVER_MODE" => "split_imex",
        "SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp47",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing
    ) do
        sol, meta = _solve_with_solver_policy(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test meta.solver == "KenCarp47(IMEX)"
        @test meta.initial_solver == "KenCarp47"
    end

    withenv("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "6") do
        @test _multirate_fast_substeps() == 6
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "bad") do
        @test_throws ArgumentError _multirate_fast_substeps()
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "0") do
        @test_throws ArgumentError _multirate_fast_substeps()
    end

    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => nothing) do
        @test _multirate_slow_dt_s(solver_args) == min(solver_args.integration_tolerances.dt_max_orbit, 2.0)
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.4") do
        @test _multirate_slow_dt_s(solver_args) == 0.4
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => "bad") do
        @test_throws ArgumentError _multirate_slow_dt_s(solver_args)
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.0") do
        @test_throws ArgumentError _multirate_slow_dt_s(solver_args)
    end

    withenv("SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "auto") do
        spec = _multirate_slow_solver_spec()
        @test spec.label == "AutoTsit5(Rodas5P)"
        @test spec.auto_switch_capable == true
    end
    withenv("SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "rodas5p") do
        spec = _multirate_slow_solver_spec()
        @test spec.label == "Rodas5P"
        @test spec.auto_switch_capable == false
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SOLVER" => "kencarp4") do
        spec = _multirate_fast_solver_spec()
        @test spec.label == "KenCarp4"
        @test spec.auto_switch_capable == false
    end
    withenv("SPACEAGORA_MULTIRATE_FAST_SOLVER" => "unsupported") do
        @test_throws ArgumentError _multirate_fast_solver_spec()
    end

    multirate_subprob = _split_subproblem(split_prob_simple, split_prob_simple.f.f1, [1.0], (0.0, 0.5))
    @test multirate_subprob.tspan == (0.0, 0.5)

    @test_throws ArgumentError _solve_with_multirate_solver(prob_simple, solver_args, 1e-8, 1e-8)

    split_prob_zero = SplitODEProblem(
        (du, u, p, t) -> begin
            du[1] = -u[1]
        end,
        (du, u, p, t) -> begin
            du[1] = -2u[1]
        end,
        [1.0],
        (1.0, 1.0)
    )
    withenv("SPACEAGORA_SOLVER_MAXITERS" => nothing) do
        sol_zero, meta_zero = _solve_with_multirate_solver(split_prob_zero, solver_args, 1e-8, 1e-8)
        @test string(sol_zero.retcode) == "Success"
        @test meta_zero.macro_steps == 0
        @test meta_zero.fast_substeps == 0
        @test meta_zero.slow_solver == "Tsit5"
        @test meta_zero.fast_solver == "Tsit5"
    end

    withenv(
        "SPACEAGORA_SOLVER_MAXITERS" => nothing,
        "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "4",
        "SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.2"
    ) do
        sol_mr, meta_mr = _solve_with_multirate_solver(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol_mr.retcode) == "Success"
        @test meta_mr.macro_steps >= 1
        @test meta_mr.fast_substeps == 4
        @test meta_mr.slow_solver == "Tsit5"
        @test meta_mr.fast_solver == "Tsit5"
        @test isapprox(meta_mr.slow_dt_s, 0.2; atol=0.0, rtol=0.0)
        @test isapprox(meta_mr.fast_dt_s, 0.05; atol=0.0, rtol=0.0)
        @test meta_mr.auto_switch_events == 0
    end

    withenv(
        "SPACEAGORA_SOLVER_MODE" => "multirate",
        "SPACEAGORA_SOLVER_MAXITERS" => nothing,
        "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "4",
        "SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.2"
    ) do
        sol, meta = _solve_with_solver_policy(split_prob_simple, solver_args, 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test occursin("Multirate(Strang;", meta.solver)
        @test meta.initial_solver == "Tsit5"
    end

    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "oops") do
        @test_throws ArgumentError _ephemeris_reuse_max_entries()
    end

    reuse_cache = Dict{Any, SRPSunEphemerisCache}()
    reuse_value_a = SRPSunEphemerisCache([0.0], SVector{3, Float64}[SVector{3, Float64}(1.0, 0.0, 0.0)])
    reuse_value_b = SRPSunEphemerisCache([0.0], SVector{3, Float64}[SVector{3, Float64}(2.0, 0.0, 0.0)])
    @test _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_a, 0) === reuse_value_a
    @test !haskey(reuse_cache, :k1)
    _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_a, 2)
    @test _ephemeris_reuse_store!(reuse_cache, :k1, reuse_value_b, 2) === reuse_value_a
    _ephemeris_reuse_store!(reuse_cache, :k2, reuse_value_b, 1)
    @test !haskey(reuse_cache, :k1)
    @test haskey(reuse_cache, :k2)

    withenv("SPACEAGORA_SOLVER_MAXITERS" => "1000") do
        sol = _solve_with_explicit_solver(prob_simple, solver_args, Tsit5(), 1e-8, 1e-8)
        @test string(sol.retcode) == "Success"
        @test !isempty(sol.t)
    end

    sc_nan_inertia = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    sc_nan_inertia.inertia_tensor = SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, NaN, 0.0, 0.0, 0.0, 1.0)
    args_nan_inertia = build_config(
        spacecraft=sc_nan_inertia,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    @test_throws ArgumentError _validate_orientation_inertia!(args_nan_inertia)

    sc_nonsym_inertia = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    sc_nonsym_inertia.inertia_tensor = SMatrix{3, 3, Float64}(2.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0)
    args_nonsym_inertia = build_config(
        spacecraft=sc_nonsym_inertia,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    @test_throws ArgumentError _validate_orientation_inertia!(args_nonsym_inertia)

    struct ThermalModelNoHeatRate <: SimulationModel.AbstractThermalModel end
    @test_throws ArgumentError _validate_thermal_model_support!((environment_model=(thermal_model=ThermalModelNoHeatRate(),),))

    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "12.5") do
        @test _srp_ephemeris_cache_dt_s() == 12.5
    end
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "bad") do
        @test_throws ArgumentError _srp_ephemeris_cache_dt_s()
    end
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "1234") do
        @test _srp_ephemeris_cache_max_samples() == 1234
    end
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "bad") do
        @test_throws ArgumentError _srp_ephemeris_cache_max_samples()
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "15.0") do
        @test _nbody_ephemeris_cache_dt_s() == 15.0
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "bad") do
        @test_throws ArgumentError _nbody_ephemeris_cache_dt_s()
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "1234") do
        @test _nbody_ephemeris_cache_max_samples() == 1234
    end
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "bad") do
        @test_throws ArgumentError _nbody_ephemeris_cache_max_samples()
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_DT_S" => "2.5") do
        @test _planet_frame_cache_dt_s() == 2.5
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_DT_S" => "bad") do
        @test_throws ArgumentError _planet_frame_cache_dt_s()
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES" => "4321") do
        @test _planet_frame_cache_max_samples() == 4321
    end
    withenv("SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES" => "bad") do
        @test_throws ArgumentError _planet_frame_cache_max_samples()
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "1") do
        @test _ephemeris_reuse_enabled() == true
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "0") do
        @test _ephemeris_reuse_enabled() == false
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "bad") do
        @test_throws ArgumentError _ephemeris_reuse_enabled()
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "7") do
        @test _ephemeris_reuse_max_entries() == 7
    end
    withenv("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "-1") do
        @test_throws ArgumentError _ephemeris_reuse_max_entries()
    end
    withenv("SPACEAGORA_EFFECTOR_LONG_ORBIT_THRESHOLD" => "9") do
        @test _effector_long_orbit_threshold() == 9
        args_orbit_mission = (
            mission_configuration=(
                mission_type=SimulationModel.MissionOrbits,
                number_of_orbits=9,
                mission_time=1.0
            ),
        )
        @test _mission_is_long_for_effector_threads(args_orbit_mission)
    end
    @test _has_active_srp_effector((SolarRadiationPressureModel(1.2, 10.0),))

    args_eff_unsupported = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), ConstantForceModel(SVector{3, Float64}(1.0, 0.0, 0.0))),
        keplerian=true
    )
    decision_unsupported = _dynamic_effector_thread_decision(args_eff_unsupported, args_eff_unsupported.dynamics_model.dynamic_effectors, 1)
    @test decision_unsupported.use_threads == false
    @test decision_unsupported.policy_applied == false

    p_workspace_resize = ODEParams{1}(args=args_eff_single)
    resize!(p_workspace_resize.shared_buffers.harmonics_workspaces, 0)
    resize!(p_workspace_resize.shared_buffers.nbody_workspaces, 0)
    resize!(p_workspace_resize.shared_buffers.aero_workspaces, 0)
    _initialize_harmonics_workspace_buffers!(p_workspace_resize)
    _initialize_nbody_workspace_buffers!(p_workspace_resize)
    _initialize_aero_workspace_buffers!(p_workspace_resize)
    @test length(p_workspace_resize.shared_buffers.harmonics_workspaces) == 1
    @test length(p_workspace_resize.shared_buffers.nbody_workspaces) == 1
    @test length(p_workspace_resize.shared_buffers.aero_workspaces) == 1

    args_srp = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), SolarRadiationPressureModel(1.2, 12.0)),
        keplerian=true
    )
    p_srp = ODEParams{1}(args=args_srp)
    _initialize_srp_sun_cache_buffer!(p_srp)
    withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE" => "1") do
        _initialize_srp_sun_ephemeris_cache!(p_srp, 0.0, 0.0)
    end
    @test p_srp.shared_buffers.srp_sun_ephemeris_cache[] === nothing
    withenv(
        "SPACEAGORA_SRP_EPHEMERIS_CACHE" => "1",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "1.0",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "2"
    ) do
        @test_logs (:warn, r"SRP ephemeris cache disabled") _initialize_srp_sun_ephemeris_cache!(p_srp, 0.0, 10.0)
    end
    @test p_srp.shared_buffers.srp_sun_ephemeris_cache[] === nothing

    _clear_ephemeris_reuse_cache!()
    p_srp_reuse_a = ODEParams{1}(args=args_srp)
    p_srp_reuse_b = ODEParams{1}(args=args_srp)
    _initialize_srp_sun_cache_buffer!(p_srp_reuse_a)
    _initialize_srp_sun_cache_buffer!(p_srp_reuse_b)
    _reset_spice_runtime_counters!(p_srp_reuse_a)
    _reset_spice_runtime_counters!(p_srp_reuse_b)
    withenv(
        "SPACEAGORA_SRP_EPHEMERIS_CACHE" => "1",
        "SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "1",
        "SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "16",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S" => "10.0",
        "SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES" => "100"
    ) do
        _initialize_srp_sun_ephemeris_cache!(p_srp_reuse_a, 0.0, 10.0)
        cache_a = p_srp_reuse_a.shared_buffers.srp_sun_ephemeris_cache[]
        @test cache_a isa SRPSunEphemerisCache
        @test p_srp_reuse_a.shared_buffers.spice_runtime_counters.srp_spkpos_cache_build_calls[] == 2
        _initialize_srp_sun_ephemeris_cache!(p_srp_reuse_b, 0.0, 10.0)
        cache_b = p_srp_reuse_b.shared_buffers.srp_sun_ephemeris_cache[]
        @test cache_b === cache_a
        @test p_srp_reuse_b.shared_buffers.spice_runtime_counters.srp_spkpos_cache_build_calls[] == 0
    end
    _clear_ephemeris_reuse_cache!()

    args_nbody = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), NBodyGravityModel(body_names=("moon",), primary_body_name="Earth")),
        keplerian=true
    )
    p_nbody = ODEParams{1}(args=args_nbody)
    _initialize_nbody_ephemeris_cache_buffer!(p_nbody)
    withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "1") do
        _initialize_nbody_ephemeris_cache!(p_nbody, 0.0, 0.0)
    end
    @test p_nbody.shared_buffers.nbody_ephemeris_cache[] === nothing
    withenv(
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "1",
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "1.0",
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "2"
    ) do
        @test_logs (:warn, r"N-body ephemeris cache disabled") _initialize_nbody_ephemeris_cache!(p_nbody, 0.0, 10.0)
    end
    @test p_nbody.shared_buffers.nbody_ephemeris_cache[] === nothing

    p_planet_frame = ODEParams{1}(args=args_srp)
    _initialize_planet_frame_cache_buffer!(p_planet_frame)
    withenv("SPACEAGORA_PLANET_FRAME_CACHE" => "1") do
        _initialize_planet_frame_ephemeris_cache!(p_planet_frame, 0.0, 0.0)
    end
    @test p_planet_frame.shared_buffers.planet_frame_ephemeris_cache[] === nothing
    withenv(
        "SPACEAGORA_PLANET_FRAME_CACHE" => "1",
        "SPACEAGORA_PLANET_FRAME_CACHE_DT_S" => "1.0",
        "SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES" => "2"
    ) do
        @test_logs (:warn, r"Planet frame cache disabled") _initialize_planet_frame_ephemeris_cache!(p_planet_frame, 0.0, 10.0)
    end
    @test p_planet_frame.shared_buffers.planet_frame_ephemeris_cache[] === nothing

    args_nbody_srp = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(body_names=("moon", "sun"), primary_body_name="Earth"),
            SolarRadiationPressureModel(1.2, 12.0)
        ),
        keplerian=true
    )
    p_nbody_srp = ODEParams{1}(args=args_nbody_srp)
    _initialize_nbody_ephemeris_cache_buffer!(p_nbody_srp)
    _initialize_srp_sun_cache_buffer!(p_nbody_srp)
    _reset_spice_runtime_counters!(p_nbody_srp)
    _reset_spice_rhs_memo!(p_nbody_srp)
    p_nbody_srp.shared_buffers.et_start[] = 0.0
    p_nbody_srp.shared_buffers.current_time[] = 123.0
    nbody_model = args_nbody_srp.dynamics_model.dynamic_effectors[2]
    srp_model = args_nbody_srp.dynamics_model.dynamic_effectors[3]
    x_a = Float64[7000e3, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]
    x_b = Float64[6999e3, 10.0, 0.0, 0.0, 0.0, 0.0, 200.0]
    SimulationModel.calcForceTorque(nbody_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 2
    SimulationModel.calcForceTorque(nbody_model, x_b, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 2
    SimulationModel.calcForceTorque(srp_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls[] == 0
    p_nbody_srp.shared_buffers.current_time[] = 124.0
    SimulationModel.calcForceTorque(srp_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls[] == 1
    p_nbody_srp.shared_buffers.spice_rhs_memo_enabled[] = false
    _reset_spice_runtime_counters!(p_nbody_srp)
    _reset_spice_rhs_memo!(p_nbody_srp)
    p_nbody_srp.shared_buffers.current_time[] = 223.0
    SimulationModel.calcForceTorque(nbody_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 2
    SimulationModel.calcForceTorque(nbody_model, x_b, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls[] == 4
    SimulationModel.calcForceTorque(srp_model, x_a, p_nbody_srp, 1)
    @test p_nbody_srp.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls[] == 1

    dyn = SimulationModel.DynamicEffectors
    harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
    @test_throws ArgumentError GravitationalHarmonicsModel(-1, 0, harmonics_file, EARTH)
    @test_throws ArgumentError GravitationalHarmonicsModel(10_000, 0, harmonics_file, EARTH)

    nbody_ws = dyn._make_nbody_scratch_workspace(1)
    dyn._ensure_nbody_workspace_capacity!(nbody_ws, 3, 4)
    @test length(nbody_ws.pos_primary_k_all) == 3
    @test length(nbody_ws.thread_force) == 4
    nbody_ws_oob = dyn._nbody_workspace_for_sat!(p_nbody_srp, 5, 2, 2)
    @test length(nbody_ws_oob.pos_primary_k_all) == 2
    @test length(nbody_ws_oob.thread_force) == 2

    nbody_positions = reshape(
        SVector{3, Float64}[
            SVector{3, Float64}(1.0, 0.0, 0.0),
            SVector{3, Float64}(2.0, 0.0, 0.0),
            SVector{3, Float64}(3.0, 0.0, 0.0),
            SVector{3, Float64}(4.0, 0.0, 0.0)
        ],
        4,
        1
    )
    nbody_cache = NBodyEphemerisCache(
        "earth_barycenter",
        ["moon_barycenter"],
        Dict("moon_barycenter" => 1),
        [0.0, 5.0, 10.0, 15.0],
        nbody_positions
    )
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, -1.0, "moon_barycenter", "earth_barycenter") === nothing
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, NaN, "moon_barycenter", "earth_barycenter") isa SVector{3, Float64}
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, 15.0, "moon_barycenter", "earth_barycenter") == nbody_positions[4, 1]
    @test dyn._nbody_body_position_from_cache_j2000(nbody_cache, 2.5, "moon_barycenter", "earth_barycenter") == SVector{3, Float64}(1.5, 0.0, 0.0)
    nbody_interp_catmull = dyn._nbody_body_position_from_cache_j2000(nbody_cache, 7.5, "moon_barycenter", "earth_barycenter")
    @test nbody_interp_catmull isa SVector{3, Float64}

    srp_cache = SRPSunEphemerisCache(
        [0.0, 5.0, 10.0, 15.0],
        SVector{3, Float64}[
            SVector{3, Float64}(1.0, 0.0, 0.0),
            SVector{3, Float64}(2.0, 0.0, 0.0),
            SVector{3, Float64}(3.0, 0.0, 0.0),
            SVector{3, Float64}(4.0, 0.0, 0.0)
        ]
    )
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, -1.0) === nothing
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, NaN) isa SVector{3, Float64}
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, 15.0) == srp_cache.positions_j2000[end]
    @test dyn._srp_sun_position_from_cache_j2000(srp_cache, 2.5) == SVector{3, Float64}(1.5, 0.0, 0.0)
    srp_interp_catmull = dyn._srp_sun_position_from_cache_j2000(srp_cache, 7.5)
    @test srp_interp_catmull isa SVector{3, Float64}

    force_zero_srp, torque_zero_srp = calcForceTorque(SolarRadiationPressureModel(1.2, 0.0), x_a, p_nbody_srp, 1)
    @test force_zero_srp == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_zero_srp == SVector{3, Float64}(0.0, 0.0, 0.0)

    p_srp.shared_buffers.srp_sun_ephemeris_cache[] = SRPSunEphemerisCache(
        [0.0, 10.0],
        SVector{3, Float64}[SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)]
    )
    p_srp.shared_buffers.et_start[] = 0.0
    p_srp.shared_buffers.current_time[] = 5.0
    x_zero_dist = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0]
    force_zero_dist, torque_zero_dist = calcForceTorque(SolarRadiationPressureModel(1.2, 12.0), x_zero_dist, p_srp, 1)
    @test force_zero_dist == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_zero_dist == SVector{3, Float64}(0.0, 0.0, 0.0)

    @test dyn.eclipse_area_calc(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(1.0, 0.0, 0.0), EARTH.Rp_e) == 1.0
    @test dyn.eclipse_area_calc(
        SVector{3, Float64}(-5.271937279754128e6, -4.185218527153555e6, -1.0434271238606143e6),
        SVector{3, Float64}(-1.1107937751409042e10, 1.6143690900498734e11, 8.75872644289004e10),
        EARTH.Rp_e
    ) == 0.0
    annular_ratio = dyn.eclipse_area_calc(
        SVector{3, Float64}(0.0, 0.0, 1.0e10),
        SVector{3, Float64}(0.0, 0.0, -1.5e11),
        EARTH.Rp_e
    )
    @test 0.0 < annular_ratio < 1.0
    partial_ratio = dyn.eclipse_area_calc(
        SVector{3, Float64}(6.930027129188876e6, -2.6352977471555886e6, -3.363004422597388e6),
        SVector{3, Float64}(-9.438128429326639e10, -6.696979722657439e10, 1.7822072441008075e11),
        EARTH.Rp_e
    )
    @test 0.0 < partial_ratio < 1.0
    none_ratio = dyn.eclipse_area_calc(
        SVector{3, Float64}(1.2249535847697716e7, -5.782145543435082e6, -7.299266925237677e6),
        SVector{3, Float64}(-4.937007687846062e10, -7.233778731136734e10, 1.7733640551002377e11),
        EARTH.Rp_e
    )
    @test isapprox(none_ratio, 1.0; atol=1e-12, rtol=0.0)

    @testset "SRP Regression Contracts" begin
        # On/off delta: SRP-enabled model must produce non-zero force while area=0 is zero force.
        p_srp.shared_buffers.srp_sun_ephemeris_cache[] = SRPSunEphemerisCache(
            [0.0, 10.0],
            SVector{3, Float64}[
                SVector{3, Float64}(149_597_870.7, 0.0, 0.0),
                SVector{3, Float64}(149_597_870.7, 0.0, 0.0)
            ]
        )
        p_srp.shared_buffers.et_start[] = 0.0
        p_srp.shared_buffers.current_time[] = 5.0
        x_srp_regression = Float64[7_000_000.0, 100.0, 50.0, 0.0, 0.0, 0.0, 200.0]

        force_on, _ = calcForceTorque(SolarRadiationPressureModel(1.2, 12.0), x_srp_regression, p_srp, 1)
        force_off, _ = calcForceTorque(SolarRadiationPressureModel(1.2, 0.0), x_srp_regression, p_srp, 1)
        @test norm(force_on) > 0.0
        @test force_off == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test norm(force_on - force_off) > 0.0

        # Eclipse attenuation: partial eclipse force norm must scale by eclipse ratio.
        pos_partial = SVector{3, Float64}(6.930027129188876e6, -2.6352977471555886e6, -3.363004422597388e6)
        sun_partial = SVector{3, Float64}(-9.438128429326639e10, -6.696979722657439e10, 1.7822072441008075e11)
        partial_ratio_regression = dyn.eclipse_area_calc(pos_partial, sun_partial, EARTH.Rp_e)
        force_partial = dyn.srp_cannonball_accel(pos_partial, sun_partial, EARTH.Rp_e, 4.56e-6, 1.2, 12.0, 200.0)
        force_no_eclipse = dyn.srp_cannonball_accel(pos_partial, sun_partial, 0.0, 4.56e-6, 1.2, 12.0, 200.0)
        @test 0.0 < partial_ratio_regression < 1.0
        @test norm(force_no_eclipse) > 0.0
        @test isapprox(norm(force_partial), partial_ratio_regression * norm(force_no_eclipse); atol=0.0, rtol=1e-12)

        # Solver parity: legacy solver helper (`srp`) must match canonical SRP kernel for same geometry.
        planet = args_srp.environment_model.planet
        pos_solver = SVector{3, Float64}(x_a[1:3])
        et_solver = 123.0
        primary_body_name = dyn._spice_query_name(planet.name)
        sun_j2000 = lock(SimulationModel.SPICE_LOCK) do
            SVector{3, Float64}(spkpos("sun", et_solver, "J2000", "none", primary_body_name)[1])
        end
        sun_pci = SVector{3, Float64}(planet.J2000_to_pci * sun_j2000 * 1e3)
        force_kernel = dyn.srp_cannonball_accel(pos_solver, sun_pci, planet.Rp_e, 4.56e-6, 1.2, 12.0, 200.0)
        force_solver = dyn.srp(planet, 4.56e-6, 1.2, 12.0, 200.0, pos_solver, et_solver)
        @test isapprox(force_solver, force_kernel; atol=1e-12, rtol=1e-10)
    end

    default_ckpt_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    args_default_ckpt = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=default_ckpt_settings
    )
    @test _checkpoint_directory(args_default_ckpt) == joinpath("output", "checkpoints")

    mktempdir() do tmp
        bad_ckpt_settings = SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false,
            checkpoint_directory=tmp
        )
        args_bad_ckpt = build_config(
            spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=10.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true,
            simulation_settings=bad_ckpt_settings
        )
        paths_bad_ckpt = _checkpoint_paths(args_bad_ckpt)
        mkpath(dirname(paths_bad_ckpt.data))
        open(paths_bad_ckpt.data, "w") do io
            serialize(io, Dict(:invalid => true))
        end
        @test_throws ArgumentError _load_checkpoint(args_bad_ckpt)
    end

    @test _find_sample_value([nothing, nothing]) === nothing
    results_df = DataFrame()
    _append_series_columns!(results_df, "meta", [(a=1, b=2), (a=3, b=4)])
    _append_series_columns!(results_df, "dictmeta", [Dict(:z => 1, :a => 2), Dict(:z => 3, :a => 4)])
    struct BoxedValue
        value::Int
    end
    _append_series_columns!(results_df, "boxed", [BoxedValue(1), BoxedValue(2)])
    @test results_df.meta_a == [1, 3]
    @test results_df.meta_b == [2, 4]
    @test results_df.dictmeta_a == [2, 4]
    @test results_df.dictmeta_z == [1, 3]
    @test length(results_df.boxed) == 2

    mktempdir() do tmp
        out_path = joinpath(tmp, "artifact.txt")
        @test_throws ErrorException _atomic_write_file(out_path, tmp_path -> begin
            write(tmp_path, "tmp-data")
            throw(ErrorException("forced writer failure"))
        end)
        @test !isfile(out_path)
        @test isempty(readdir(tmp))
    end

    args_heat_copy = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    u_heat_copy = build_initial_conditions(args_heat_copy)
    du_heat_copy = copy(u_heat_copy)
    du_heat_copy .= 0.0
    p_heat_copy = ODEParams{1}(args=args_heat_copy)
    _initialize_heat_rate_buffers!(p_heat_copy)
    p_heat_copy.shared_buffers.heat_rates[1] = Float64[1.25, 9.0]
    spacecraft_dynamics!(du_heat_copy, u_heat_copy, p_heat_copy, 0.0)
    @test du_heat_copy.sc[1].heat_loads[1] == 1.25

    @test_throws ArgumentError _resolve_component_tolerance(-1.0, 1.0, "unit_test_tol")
end

@testset "Rigid-Body Angular Dynamics Uses Inertia Tensor" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.02, -0.03, 0.015)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)
    sc.inertia_tensor = inertia_tensor
    applied_torque = SVector{3, Float64}(0.12, -0.08, 0.05)

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(ConstantTorqueModel(applied_torque),),
        keplerian=true
    )

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    ω = SVector{3, Float64}(u0.sc[1].ω)
    qdot_expected = SimulationModel.DynamicsRotational.quaternion_derivative(ω, SVector{4, Float64}(u0.sc[1].q))
    ωdot_expected = SimulationModel.DynamicsRotational.angular_acceleration(
        ω,
        inertia_tensor,
        applied_torque;
        include_gyroscopic=true
    )
    ωdot_legacy = inertia_tensor \ (applied_torque - cross(ω, inertia_tensor * ω))
    @test isapprox(SVector{4, Float64}(du0.sc[1].q), qdot_expected; atol=1e-12, rtol=1e-10)
    @test isapprox(ωdot_expected, ωdot_legacy; atol=1e-12, rtol=1e-10)
    @test isapprox(SVector{3, Float64}(du0.sc[1].ω), ωdot_expected; atol=1e-12, rtol=1e-10)
end

@testset "Orientation Simulation Rejects Invalid Inertia Tensor" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    sc.inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    @test_throws ArgumentError run_simulation(args)
end

@testset "Heat Loads Are Not Coupled To Force Magnitude" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(ConstantForceModel(SVector{3, Float64}(1.0e5, -2.0e5, 3.0e5)),),
        keplerian=true
    )
    args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    @test norm(SVector{3, Float64}(du0.sc[1].vel)) > 0.0
    @test all(==(0.0), du0.sc[1].heat_loads)
end

@testset "Heat Load Derivative Uses Thermal Model" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=175.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(I(3))

    u0 = build_initial_conditions(args)
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams{1}(args=args)
    _initialize_heat_rate_buffers!(p)
    p.shared_buffers.densities[1] = 1.0e-6
    p.shared_buffers.temperatures[1] = 250.0
    p.shared_buffers.winds[1] = SVector{3, Float64}(0.0, 0.0, 0.0)
    thermal_cb = SimulationModel.SimulationCallbacks.get_thermal_callback(1, args)
    thermal_cb.affect!((p=p, u=u0, t=0.0))
    spacecraft_dynamics!(du0, u0, p, 0.0)

    @test all(isfinite, p.shared_buffers.heat_rates[1])
    @test any(>(0.0), p.shared_buffers.heat_rates[1])
    @test all(isfinite, du0.sc[1].heat_loads)
    @test any(>(0.0), du0.sc[1].heat_loads)
end

@testset "Drag Dissipates Specific Orbital Energy" begin
    sc = make_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=180.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )

    df = run_case(args)
    eps = specific_energy(df, EARTH.μ)
    @test last(eps) < first(eps) - 1e5
end

@testset "AGORA Earth Regression (Golden)" begin
    golden_path = joinpath(@__DIR__, "golden", "agora_earth_regression.csv")
    @test isfile(golden_path)
    golden = CSV.read(golden_path, DataFrame)

    sc = make_agora_earth_spacecraft()
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=6.0 * 3600.0,
        EI_km=300.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        initial_time=SimulationModel.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=15.0
        )
    )

    df = run_case(args)
    @test nrow(df) > 1000
    times = Vector{Float64}(df.time)

    for row in eachrow(golden)
        t = Float64(row.time)
        pos_atol = Float64(row.pos_atol_m)
        vel_atol = Float64(row.vel_atol_mps)

        pos1 = interp_linear(times, df.sc1_pos_1, t)
        pos2 = interp_linear(times, df.sc1_pos_2, t)
        pos3 = interp_linear(times, df.sc1_pos_3, t)
        vel1 = interp_linear(times, df.sc1_vel_1, t)
        vel2 = interp_linear(times, df.sc1_vel_2, t)
        vel3 = interp_linear(times, df.sc1_vel_3, t)

        @test isapprox(pos1, Float64(row.pos1); atol=pos_atol, rtol=0.0)
        @test isapprox(pos2, Float64(row.pos2); atol=pos_atol, rtol=0.0)
        @test isapprox(pos3, Float64(row.pos3); atol=pos_atol, rtol=0.0)
        @test isapprox(vel1, Float64(row.vel1); atol=vel_atol, rtol=0.0)
        @test isapprox(vel2, Float64(row.vel2); atol=vel_atol, rtol=0.0)
        @test isapprox(vel3, Float64(row.vel3); atol=vel_atol, rtol=0.0)
    end
end

@testset "N-Body Gravity Typed API Smoke" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df = run_case(args)
    @test nrow(df) > 10
    @test all(isfinite, df.sc1_pos_1)
    @test all(isfinite, df.sc1_vel_1)
end

@testset "Two-Spacecraft Isolation vs Single-Craft Baselines" begin
    sc_a = make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0)
    sc_b = make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0)

    args_multi = build_config_multi(
        spacecraft=[sc_a, sc_b],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_a = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_b = build_config(
        spacecraft=make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df_multi = run_case_silent(args_multi)
    df_a = run_case_silent(args_a)
    df_b = run_case_silent(args_b)

    @test nrow(df_multi) > 10
    sample_idxs = round.(Int, range(1, nrow(df_multi), length=8))
    for idx in sample_idxs
        t = Float64(df_multi.time[idx])

        pa_m = SVector{3, Float64}(Float64(df_multi.sc1_pos_1[idx]), Float64(df_multi.sc1_pos_2[idx]), Float64(df_multi.sc1_pos_3[idx]))
        va_m = SVector{3, Float64}(Float64(df_multi.sc1_vel_1[idx]), Float64(df_multi.sc1_vel_2[idx]), Float64(df_multi.sc1_vel_3[idx]))
        pa_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_pos_1, t),
            interp_linear(df_a.time, df_a.sc1_pos_2, t),
            interp_linear(df_a.time, df_a.sc1_pos_3, t)
        )
        va_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_vel_1, t),
            interp_linear(df_a.time, df_a.sc1_vel_2, t),
            interp_linear(df_a.time, df_a.sc1_vel_3, t)
        )
        # Multi-satellite adaptive stepping can introduce modest trajectory differences vs single-body runs.
        @test norm(pa_m - pa_s) < 200.0
        @test norm(va_m - va_s) < 0.2

        pb_m = SVector{3, Float64}(Float64(df_multi.sc2_pos_1[idx]), Float64(df_multi.sc2_pos_2[idx]), Float64(df_multi.sc2_pos_3[idx]))
        vb_m = SVector{3, Float64}(Float64(df_multi.sc2_vel_1[idx]), Float64(df_multi.sc2_vel_2[idx]), Float64(df_multi.sc2_vel_3[idx]))
        pb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_pos_1, t),
            interp_linear(df_b.time, df_b.sc1_pos_2, t),
            interp_linear(df_b.time, df_b.sc1_pos_3, t)
        )
        vb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_vel_1, t),
            interp_linear(df_b.time, df_b.sc1_vel_2, t),
            interp_linear(df_b.time, df_b.sc1_vel_3, t)
        )
        @test norm(pb_m - pb_s) < 200.0
        @test norm(vb_m - vb_s) < 0.2
    end
end

@testset "Single-Link Drag Dissipates Specific Orbital Energy" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=180.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )

    df = run_case(args)
    eps = specific_energy(df, EARTH.μ)
    @test last(eps) < first(eps) - 1e5
end

@testset "Thruster Edge Cases (AI Second-Pass)" begin
    @testset "BaseThrusterModel Validation" begin
        model_ok = BaseThrusterModel(
            thrust=[1.0, 2.0],
            direction=[0.0, π],
            Δv=[10.0, 20.0],
            start_burn_time=[0.0, 0.0],
            stop_burn_time=[100.0, 120.0],
            Isp=[300.0, 300.0]
        )
        @test length(model_ok.thrust) == 2
        @test length(model_ok.direction) == 2

        @test_throws ArgumentError BaseThrusterModel(
            thrust=[1.0, 2.0],
            direction=[0.0],
            Δv=[10.0, 20.0],
            start_burn_time=[0.0, 0.0],
            stop_burn_time=[100.0, 120.0],
            Isp=[300.0, 300.0]
        )
    end

    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    p = ODEParams{1}(args=args)
    u = build_initial_conditions(args).sc[1]

    @testset "calcControlForceTorque" begin
        model = make_base_thruster_model(thrust=2.0, direction=0.0, start_burn_time=50.0, stop_burn_time=150.0)

        force_pre, torque_pre = calcControlForceTorque(model, u, p, 1, 49.9)
        @test force_pre == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_pre == SVector{3, Float64}(0.0, 0.0, 0.0)

        expected_dir = normalize(SVector{3, Float64}(u.vel))
        force_start, torque_start = calcControlForceTorque(model, u, p, 1, 50.0)
        @test isapprox(force_start, 2.0 .* expected_dir; atol=1e-12, rtol=0.0)
        @test torque_start == SVector{3, Float64}(0.0, 0.0, 0.0)

        force_stop, _ = calcControlForceTorque(model, u, p, 1, 150.0)
        @test isapprox(force_stop, 2.0 .* expected_dir; atol=1e-12, rtol=0.0)

        force_post, _ = calcControlForceTorque(model, u, p, 1, 150.1)
        @test force_post == SVector{3, Float64}(0.0, 0.0, 0.0)

        model.direction[1] = π
        force_retro, _ = calcControlForceTorque(model, u, p, 1, 100.0)
        @test dot(force_retro, expected_dir) < 0.0

        u_zero = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])
        force_zero, torque_zero = calcControlForceTorque(model, u_zero, p, 1, 100.0)
        @test force_zero == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_zero == SVector{3, Float64}(0.0, 0.0, 0.0)

        force_oob, torque_oob = calcControlForceTorque(model, u, p, 2, 100.0)
        @test force_oob == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_oob == SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    @testset "calcControlMassFlowRate" begin
        model = make_base_thruster_model(
            thrust=2.0,
            direction=0.0,
            start_burn_time=50.0,
            stop_burn_time=150.0,
            Isp=300.0
        )

        mdot_pre = calcControlMassFlowRate(model, u, p, 1, 49.9)
        @test mdot_pre == 0.0

        mdot_on = calcControlMassFlowRate(model, u, p, 1, 100.0)
        expected_mdot = -2.0 / (300.0 * 9.80665)
        @test isapprox(mdot_on, expected_mdot; atol=1e-14, rtol=0.0)

        mdot_post = calcControlMassFlowRate(model, u, p, 1, 150.1)
        @test mdot_post == 0.0

        model_bad_isp = make_base_thruster_model(
            thrust=2.0,
            direction=0.0,
            start_burn_time=50.0,
            stop_burn_time=150.0,
            Isp=0.0
        )
        @test calcControlMassFlowRate(model_bad_isp, u, p, 1, 100.0) == 0.0

        u_zero = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])
        @test calcControlMassFlowRate(model, u_zero, p, 1, 100.0) == 0.0
        @test calcControlMassFlowRate(model, u, p, 2, 100.0) == 0.0

        model_subtyped = TimedTangentialThrusterModel(1.0, +1.0, 10.0, 20.0)
        @test calcControlMassFlowRate(model_subtyped, u, p, 1, 15.0) == 0.0

        struct UntypedControlEffector end
        @test calcControlMassFlowRate(UntypedControlEffector(), u, p, 1, 100.0) == 0.0
    end

    @testset "calcControlEffect!" begin
        model = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        state = build_initial_conditions(args)
        calcControlEffect!(model, state, p, 100.0, 1)
        @test isfinite(model.start_burn_time[1])
        @test isfinite(model.stop_burn_time[1])
        @test model.start_burn_time[1] < model.stop_burn_time[1]
        expected_burn_duration = state.sc[1].mass * model.Δv[1] / model.thrust[1]
        @test isapprox(
            model.stop_burn_time[1] - model.start_burn_time[1],
            expected_burn_duration;
            atol=1e-9,
            rtol=0.0
        )

        # Pre-ignition tracking: a future scheduled window should be retimed.
        model_track = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=10_000.0, stop_burn_time=11_000.0)
        calcControlEffect!(model_track, state, p, 100.0, 1)
        @test model_track.start_burn_time[1] != 10_000.0
        @test model_track.stop_burn_time[1] != 11_000.0

        # Post-ignition lock: once burn start time has been reached, keep fixed.
        s_sched = model.start_burn_time[1]
        e_sched = model.stop_burn_time[1]
        calcControlEffect!(model, state, p, 120.0, 1)
        @test model.start_burn_time[1] == s_sched
        @test model.stop_burn_time[1] == e_sched

        sc_ineligible = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=210.0)
        args_ineligible = build_config(
            spacecraft=sc_ineligible,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_ineligible = ODEParams{1}(args=args_ineligible)
        state_ineligible = build_initial_conditions(args_ineligible)
        model_ineligible = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=11.0, stop_burn_time=22.0)
        calcControlEffect!(model_ineligible, state_ineligible, p_ineligible, 100.0, 1)
        @test model_ineligible.start_burn_time[1] == -1.0 || model_ineligible.start_burn_time[1] > 100.0
        @test model_ineligible.stop_burn_time[1] == -1.0 || model_ineligible.stop_burn_time[1] > model_ineligible.start_burn_time[1]

        model_zero_thrust = make_base_thruster_model(thrust=0.0, Δv=20.0, start_burn_time=33.0, stop_burn_time=44.0)
        calcControlEffect!(model_zero_thrust, state, p, 100.0, 1)
        @test model_zero_thrust.start_burn_time[1] == -1.0
        @test model_zero_thrust.stop_burn_time[1] == -1.0

        sc_edge_block = make_spacecraft(ra_alt_m=600e3, rp_alt_m=400e3, ν_deg=180.0)
        r_edge_block, _ = orbitalelemtorv(sc_edge_block.initial_condition, EARTH)
        ei_edge_block_km = (norm(SVector{3, Float64}(r_edge_block)) - EARTH.Rp_e) / 1e3
        args_edge_block = build_config(
            spacecraft=sc_edge_block,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=ei_edge_block_km,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_edge_block = ODEParams{1}(args=args_edge_block)
        state_edge_block = build_initial_conditions(args_edge_block)
        model_edge_block = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        calcControlEffect!(model_edge_block, state_edge_block, p_edge_block, 100.0, 1)
        @test model_edge_block.start_burn_time[1] == -1.0
        @test model_edge_block.stop_burn_time[1] == -1.0

        sc_edge_allow = make_spacecraft(ra_alt_m=600e3, rp_alt_m=400e3, ν_deg=170.0)
        r_edge_allow, _ = orbitalelemtorv(sc_edge_allow.initial_condition, EARTH)
        ei_edge_allow_km = (norm(SVector{3, Float64}(r_edge_allow)) - EARTH.Rp_e) / 1e3
        args_edge_allow = build_config(
            spacecraft=sc_edge_allow,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=ei_edge_allow_km,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_edge_allow = ODEParams{1}(args=args_edge_allow)
        state_edge_allow = build_initial_conditions(args_edge_allow)
        model_edge_allow = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        calcControlEffect!(model_edge_allow, state_edge_allow, p_edge_allow, 100.0, 1)
        @test model_edge_allow.start_burn_time[1] != -1.0
        @test model_edge_allow.stop_burn_time[1] != -1.0

        sc_circular = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=330.0)
        args_circular = build_config(
            spacecraft=sc_circular,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_circular = ODEParams{1}(args=args_circular)
        state_circular = build_initial_conditions(args_circular)
        model_circular = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        calcControlEffect!(model_circular, state_circular, p_circular, 100.0, 1)
        @test model_circular.start_burn_time[1] != -1.0
        @test model_circular.stop_burn_time[1] != -1.0

        state_hyperbolic = build_initial_conditions(args)
        rmag = norm(SVector{3, Float64}(state_hyperbolic.sc[1].pos))
        escape_speed = sqrt(2.0 * EARTH.μ / rmag)
        vhat = normalize(SVector{3, Float64}(state_hyperbolic.sc[1].vel))
        state_hyperbolic.sc[1].vel .= (1.2 * escape_speed) .* vhat
        model_hyperbolic = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=91.0, stop_burn_time=92.0)
        @test_nowarn calcControlEffect!(model_hyperbolic, state_hyperbolic, p, 100.0, 1)
        @test model_hyperbolic.start_burn_time[1] == -1.0
        @test model_hyperbolic.stop_burn_time[1] == -1.0

        state_near_parabolic = build_initial_conditions(args)
        rmag_parabolic = norm(SVector{3, Float64}(state_near_parabolic.sc[1].pos))
        vhat_parabolic = normalize(SVector{3, Float64}(state_near_parabolic.sc[1].vel))
        v_near_parabolic = 0.999999 * sqrt(2.0 * EARTH.μ / rmag_parabolic)
        state_near_parabolic.sc[1].vel .= v_near_parabolic .* vhat_parabolic
        model_near_parabolic = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=95.0, stop_burn_time=96.0)
        @test_nowarn calcControlEffect!(model_near_parabolic, state_near_parabolic, p, 100.0, 1)
        @test isfinite(model_near_parabolic.start_burn_time[1])
        @test isfinite(model_near_parabolic.stop_burn_time[1])

        state_singular = build_initial_conditions(args)
        state_singular.sc[1].pos .= 0.0
        state_singular.sc[1].vel .= 0.0
        model_singular = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=93.0, stop_burn_time=94.0)
        @test_nowarn calcControlEffect!(model_singular, state_singular, p, 100.0, 1)
        @test model_singular.start_burn_time[1] == -1.0
        @test model_singular.stop_burn_time[1] == -1.0

        state_tiny_a = build_initial_conditions(args)
        state_tiny_a.sc[1].pos .= SVector{3, Float64}(1.0, 0.0, 0.0)
        state_tiny_a.sc[1].vel .= SVector{3, Float64}(0.0, 0.1, 0.0)
        model_tiny_a = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=97.0, stop_burn_time=98.0)
        @test_nowarn calcControlEffect!(model_tiny_a, state_tiny_a, p, 100.0, 1)
        @test isfinite(model_tiny_a.start_burn_time[1])
        @test isfinite(model_tiny_a.stop_burn_time[1])

        sc_multi_1 = make_spacecraft(ra_alt_m=650e3, rp_alt_m=450e3, ν_deg=120.0)
        sc_multi_2 = make_spacecraft(ra_alt_m=700e3, rp_alt_m=500e3, ν_deg=150.0)
        args_multi = build_config_multi(
            spacecraft=[sc_multi_1, sc_multi_2],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_multi = ODEParams{2}(args=args_multi)
        state_multi = build_initial_conditions(args_multi)
        model_multi = BaseThrusterModel(
            thrust=[2.0, 3.0],
            direction=[0.0, π],
            Δv=[20.0, 10.0],
            start_burn_time=[-1.0, -2.0],
            stop_burn_time=[-3.0, -4.0],
            Isp=[300.0, 300.0]
        )
        calcControlEffect!(model_multi, state_multi, p_multi, 100.0, 1)
        @test model_multi.start_burn_time[1] != -1.0
        @test model_multi.stop_burn_time[1] != -3.0
        @test model_multi.start_burn_time[2] == -2.0
        @test model_multi.stop_burn_time[2] == -4.0

        s1 = model_multi.start_burn_time[1]
        e1 = model_multi.stop_burn_time[1]
        calcControlEffect!(model_multi, state_multi, p_multi, 100.0, 2)
        @test model_multi.start_burn_time[1] == s1
        @test model_multi.stop_burn_time[1] == e1
        @test model_multi.start_burn_time[2] != -2.0
        @test model_multi.stop_burn_time[2] != -4.0

        model_oob = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        @test_nowarn calcControlEffect!(model_oob, state, p, 100.0, 2)
        @test model_oob.start_burn_time[1] == -1.0
        @test model_oob.stop_burn_time[1] == -1.0

        throw_planet = ThrowingOrbitPlanet(EARTH.Rp_e)
        args_throw = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true,
            planet=throw_planet
        )
        p_throw = ODEParams{1}(args=args_throw)
        model_throw = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=123.0, stop_burn_time=124.0)
        withenv("SPACEAGORA_DEBUG_CONTROL" => "1", "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "0") do
            @test_logs (:warn, r"orbital-element conversion failed") calcControlEffect!(model_throw, state, p_throw, 100.0, 1)
        end
        @test model_throw.start_burn_time[1] == 123.0
        @test model_throw.stop_burn_time[1] == 124.0
        withenv("SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "1") do
            @test_throws ErrorException calcControlEffect!(model_throw, state, p_throw, 100.0, 1)
        end

        trace_key_helper = SimulationModel.ControlHooks._maneuver_trace_key
        trace_bool_helper = SimulationModel.ControlHooks._trace_bool_enabled
        trace_path_helper = SimulationModel.ControlHooks._maneuver_trace_path
        safe_orbit_counter_helper = SimulationModel.ControlHooks._safe_orbit_counter
        trace_event_helper = SimulationModel.ControlHooks._trace_maneuver_event!

        @test trace_bool_helper("yes")
        @test trace_bool_helper("ON")
        @test !trace_bool_helper("0")
        @test trace_key_helper(model_throw, 1) isa Tuple{UInt64, Int64}

        p_no_orbit_counter = (shared_buffers=(debug_control=Ref(false),),)
        @test safe_orbit_counter_helper(p_no_orbit_counter, 1) == -1

        mktempdir() do tmp
            trace_csv = joinpath(tmp, "maneuver_trace.csv")
            withenv(
                "SPACEAGORA_TRACE_MANEUVERS" => "1",
                "SPACEAGORA_MANEUVER_TRACE_CSV" => trace_csv
            ) do
                @test trace_path_helper() == trace_csv

                model_trace = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
                state_trace = build_initial_conditions(args)
                p_trace = ODEParams{1}(args=args)

                calcControlEffect!(model_trace, state_trace, p_trace, 100.0, 1)
                s_trace = model_trace.start_burn_time[1]
                e_trace = model_trace.stop_burn_time[1]
                @test s_trace < e_trace

                calcControlEffect!(model_trace, state_trace, p_trace, s_trace + 1e-6, 1)
                calcControlEffect!(model_trace, state_trace, p_trace, e_trace + 1.0, 1)
                @test !(model_trace.start_burn_time[1] == s_trace && model_trace.stop_burn_time[1] == e_trace)
                @test model_trace.start_burn_time[1] == -1.0 || model_trace.stop_burn_time[1] > model_trace.start_burn_time[1]

                @test isfile(trace_csv)
                trace_text = read(trace_csv, String)
                @test occursin("event,t_s,spacecraft_idx", trace_text)
                @test occursin("schedule_set", trace_text)
                @test occursin("burn_start", trace_text)
                @test occursin("burn_end", trace_text)
                @test occursin("schedule_clear", trace_text)

                trace_event_helper("manual_event", model_trace, p_no_orbit_counter, 1, 0.0)
                trace_text_after_manual = read(trace_csv, String)
                @test occursin("manual_event", trace_text_after_manual)
            end
        end

        helper = SimulationModel.ControlHooks._control_effector_exception_fallback
        p_dummy = (shared_buffers=(debug_control=Ref(false),),)
        err_eff = ErrorException("control-effector-fallback")
        @test helper(p_dummy, 1, err_eff, stacktrace()) === nothing
        withenv("SPACEAGORA_DEBUG_CONTROL" => "1", "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "0") do
            @test_logs (:warn, r"orbital-element conversion failed") helper(p_dummy, 1, err_eff, stacktrace())
        end
        withenv("SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "1") do
            @test_throws ErrorException helper(p_dummy, 1, err_eff, stacktrace())
        end
    end

    @testset "schmitt_trigger" begin
        @test schmitt_trigger(0.80, 0.75, 0.25) == 1.0
        @test schmitt_trigger(0.20, 0.75, 0.25) == 0.0
        @test schmitt_trigger(0.75, 0.75, 0.25) == 0.0
        @test schmitt_trigger(0.50, 0.75, 0.25) == 0.0

        # Sequential behavior check: current implementation is memoryless.
        state_high = schmitt_trigger(0.80, 0.75, 0.25)
        state_mid_after_high = schmitt_trigger(0.50, 0.75, 0.25)
        _ = schmitt_trigger(0.20, 0.75, 0.25)
        state_mid_after_low = schmitt_trigger(0.50, 0.75, 0.25)
        @test state_high == 1.0
        @test state_mid_after_high == 0.0
        @test state_mid_after_low == 0.0
        @test state_mid_after_high == state_mid_after_low
    end

    @testset "integrate_impulse!" begin
        link = Link{0}(root=true, attitude_control_rate=0.1)
        ω = 20.0

        thr_zero = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.0)
        impulse_zero = integrate_impulse!(link, thr_zero, 0.0, 0.0)
        @test impulse_zero == 0.0
        @test thr_zero.κ == 0.0

        thr_full = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.0)
        dt = link.attitude_control_rate
        expected_full = thr_full.max_thrust * (dt + (-1.0) / ω * (1 - exp(-ω * dt)))
        impulse_full = integrate_impulse!(link, thr_full, dt, 0.0)
        @test isapprox(impulse_full, expected_full; atol=1e-12, rtol=0.0)
        @test isapprox(thr_full.κ, 1 - exp(-ω * dt); atol=1e-12, rtol=0.0)

        thr_decay = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=1.0)
        expected_decay = thr_decay.max_thrust / ω * (1 - exp(-ω * dt))
        impulse_decay = integrate_impulse!(link, thr_decay, 0.0, 0.0)
        @test isapprox(impulse_decay, expected_decay; atol=1e-12, rtol=0.0)
        @test isapprox(thr_decay.κ, exp(-ω * dt); atol=1e-12, rtol=0.0)

        thr_small_ω = Thruster(max_thrust=10.0, cutoff_frequency=1e-12, κ=0.3)
        impulse_small_ω = integrate_impulse!(link, thr_small_ω, 0.05, 0.0)
        @test isfinite(impulse_small_ω)
        @test impulse_small_ω >= 0.0
        @test isfinite(thr_small_ω.κ)

        thr_ref_neg = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.4)
        thr_neg = deepcopy(thr_ref_neg)
        thr_zero_ref = deepcopy(thr_ref_neg)
        impulse_neg = integrate_impulse!(link, thr_neg, -0.01, 0.0)
        impulse_zero_ref = integrate_impulse!(link, thr_zero_ref, 0.0, 0.0)
        @test isapprox(impulse_neg, impulse_zero_ref; atol=1e-12, rtol=0.0)
        @test isapprox(thr_neg.κ, thr_zero_ref.κ; atol=1e-12, rtol=0.0)

        thr_ref_hi = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.4)
        thr_hi = deepcopy(thr_ref_hi)
        thr_dt_ref = deepcopy(thr_ref_hi)
        impulse_hi = integrate_impulse!(link, thr_hi, 10.0 * link.attitude_control_rate, 0.0)
        impulse_dt_ref = integrate_impulse!(link, thr_dt_ref, link.attitude_control_rate, 0.0)
        @test isapprox(impulse_hi, impulse_dt_ref; atol=1e-12, rtol=0.0)
        @test isapprox(thr_hi.κ, thr_dt_ref.κ; atol=1e-12, rtol=0.0)
    end

    @testset "thrust_calculation_schmitt_trigger!" begin
        mktempdir() do tmp
            cd(tmp) do
                link = Link{0}(root=true, attitude_control_rate=0.1)
                thr = Thruster(
                    max_thrust=10.0,
                    min_firing_time=0.05,
                    level_on=0.75,
                    level_off=0.25,
                    cutoff_frequency=100.0,
                    κ=0.0
                )
                debug_key = "SPACEAGORA_DEBUG_THRUSTER"
                old_debug = get(ENV, debug_key, nothing)
                try
                    if haskey(ENV, debug_key)
                        delete!(ENV, debug_key)
                    end
                    thrust_calculation_schmitt_trigger!(link, thr, 1.0, 0.0)
                    @test thr.thrust == 0.0
                    @test !isfile("thruster_debug.csv")

                    thrust_calculation_schmitt_trigger!(link, thr, 9.0, 0.1)
                    @test thr.thrust > 0.0
                    @test thr.thrust <= thr.max_thrust + 1e-9
                    @test !isfile("thruster_debug.csv")

                    ENV[debug_key] = "1"
                    thrust_calculation_schmitt_trigger!(link, thr, 9.0, 0.2)
                    @test isfile("thruster_debug.csv")
                finally
                    if old_debug === nothing
                        if haskey(ENV, debug_key)
                            delete!(ENV, debug_key)
                        end
                    else
                        ENV[debug_key] = old_debug
                    end
                end

                thr_zero_max = Thruster(max_thrust=0.0, cutoff_frequency=100.0, κ=0.0)
                @test_nowarn thrust_calculation_schmitt_trigger!(link, thr_zero_max, 1.0, 0.2)
                @test thr_zero_max.thrust == 0.0
            end
        end
    end

    @testset "update_thrusters!" begin
        mktempdir() do tmp
            cd(tmp) do
                link = Link{0}(root=true, attitude_control_rate=0.1)
                update_thrusters!(link, SVector{3, Float64}(0.0, 0.0, 0.0), 0.0)
                @test isempty(link.thrusters)
                @test size(link.J_thruster) == (3, 0)

                thr = Thruster(
                    max_thrust=1.0,
                    cutoff_frequency=100.0,
                    min_firing_time=0.0,
                    location=MVector{3, Float64}(0.0, 1.0, 0.0),
                    direction=MVector{3, Float64}(0.0, 0.0, 1.0)
                )
                push!(link.thrusters, thr)

                update_thrusters!(link, SVector{3, Float64}(-1.0, 0.0, 0.0), 0.1)
                @test link.thrusters[1].thrust >= 0.0

                link_full = Link{0}(root=true, attitude_control_rate=0.1)
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 1.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 0.0, 1.0), direction=MVector{3, Float64}(1.0, 0.0, 0.0)))
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(1.0, 0.0, 0.0), direction=MVector{3, Float64}(0.0, 1.0, 0.0)))
                τ_req_full = SVector{3, Float64}(1.0, 2.0, 3.0)
                update_thrusters!(link_full, τ_req_full, 0.2)
                @test rank(link_full.J_thruster) == 3
                @test all(isfinite, link_full.J_thruster)
                @test all(thr_i -> isfinite(thr_i.thrust) && thr_i.thrust >= 0.0, link_full.thrusters)
                @test any(thr_i -> thr_i.thrust > 0.0, link_full.thrusters)
                thrust_full = [thr_i.thrust for thr_i in link_full.thrusters]
                τ_ach_full = link_full.J_thruster * thrust_full
                @test norm(τ_ach_full - τ_req_full) / norm(τ_req_full) < 1e-6

                link_singular = Link{0}(root=true, attitude_control_rate=0.1)
                push!(link_singular.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 1.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                push!(link_singular.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 2.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                τ_req_singular = SVector{3, Float64}(0.0, 1.0, 0.0)
                update_thrusters!(link_singular, τ_req_singular, 0.3)
                @test rank(link_singular.J_thruster) == 1
                @test all(isfinite, link_singular.J_thruster)
                @test all(thr_i -> isfinite(thr_i.thrust) && thr_i.thrust >= 0.0, link_singular.thrusters)
                thrust_singular = [thr_i.thrust for thr_i in link_singular.thrusters]
                τ_ach_singular = link_singular.J_thruster * thrust_singular
                τ_lsq_singular = link_singular.J_thruster * (pinv(link_singular.J_thruster) * τ_req_singular)
                @test norm(τ_ach_singular - τ_lsq_singular) < 1e-6

                link_degenerate = Link{0}(root=true, attitude_control_rate=0.1)
                push!(
                    link_degenerate.thrusters,
                    Thruster(
                        max_thrust=50.0,
                        cutoff_frequency=100.0,
                        min_firing_time=0.0,
                        location=MVector{3, Float64}(0.5, 0.0, 0.0),
                        direction=MVector{3, Float64}(0.0, 0.0, 0.0)
                    )
                )
                @test_nowarn update_thrusters!(link_degenerate, SVector{3, Float64}(0.5, 0.0, 0.0), 0.4)
                @test all(isfinite, link_degenerate.J_thruster)
                @test link_degenerate.J_thruster[:, 1] == zeros(3)
                @test isfinite(link_degenerate.thrusters[1].thrust)
                @test link_degenerate.thrusters[1].thrust >= 0.0
            end
        end
    end
end

@testset "Control Burn Energy Sign (End-to-End)" begin
    sc0 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args0 = build_config(
        spacecraft=sc0,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df0 = run_case_silent(args0)
    eps0 = specific_energy(df0, EARTH.μ)
    Δeps0 = last(eps0) - first(eps0)
    @test abs(Δeps0) < 2e3

    scp = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_pro = TimedTangentialThrusterModel(800.0, +1.0, 100.0, 101.0)
    argsp = build_config(
        spacecraft=scp,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_pro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfp = run_case_silent(argsp)
    epsp = specific_energy(dfp, EARTH.μ)
    Δepsp = last(epsp) - first(epsp)
    @test Δepsp > 5e3

    scr = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_retro = TimedTangentialThrusterModel(800.0, -1.0, 100.0, 101.0)
    argsr = build_config(
        spacecraft=scr,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_retro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfr = run_case_silent(argsr)
    epsr = specific_energy(dfr, EARTH.μ)
    Δepsr = last(epsr) - first(epsr)
    @test Δepsr < -5e3
end

@testset "Control Mass Flow Coupling (End-to-End)" begin
    sc_no_control = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args_no_control = build_config(
        spacecraft=sc_no_control,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df_no_control = run_case_silent(args_no_control)
    @test "sc1_mass" in names(df_no_control)
    mass_no_control = Vector{Float64}(df_no_control.sc1_mass)
    @test maximum(abs.(mass_no_control .- first(mass_no_control))) < 1e-9

    function run_burn_case(isp::Float64)
        sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
        thruster = make_base_thruster_model(
            thrust=600.0,
            direction=0.0,
            Δv=0.0,
            start_burn_time=0.0,
            stop_burn_time=80.0,
            Isp=isp
        )
        args = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=180.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(thruster,),
            control_rates=[1.0],
            keplerian=true,
            tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
        )
        return run_case_silent(args)
    end

    df_low_isp = run_burn_case(200.0)
    df_high_isp = run_burn_case(400.0)

    mass_low = Vector{Float64}(df_low_isp.sc1_mass)
    mass_high = Vector{Float64}(df_high_isp.sc1_mass)
    Δm_low = first(mass_low) - last(mass_low)
    Δm_high = first(mass_high) - last(mass_high)

    @test Δm_low > 5.0
    @test Δm_high > 2.0
    @test all(diff(mass_low) .<= 1e-7)
    @test all(diff(mass_high) .<= 1e-7)
    @test isapprox(Δm_low / Δm_high, 2.0; atol=0.08, rtol=0.0)
end

@testset "Control Callback Multi-Spacecraft Mapping" begin
    sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    sc2 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    shared_thruster = BaseThrusterModel(
        thrust=[800.0, 800.0],
        direction=[0.0, π],
        Δv=[20.0, 20.0],
        start_burn_time=[-1.0, -1.0],
        stop_burn_time=[-1.0, -1.0],
        Isp=[300.0, 300.0]
    )
    args = build_config_multi(
        spacecraft=[sc1, sc2],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=1000.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(shared_thruster,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df = run_case_silent(args; isolate_state=false)

    for sat_idx in 1:2
        s = shared_thruster.start_burn_time[sat_idx]
        e = shared_thruster.stop_burn_time[sat_idx]
        @test (s == -1.0 && e == -1.0) || (isfinite(s) && isfinite(e) && e > s)
    end

    mass1 = Vector{Float64}(df.sc1_mass)
    mass2 = Vector{Float64}(df.sc2_mass)
    @test first(mass1) - last(mass1) > 0.1
    @test first(mass2) - last(mass2) > 0.1

    eps1 = 0.5 .* (df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    eps2 = 0.5 .* (df.sc2_vel_1.^2 .+ df.sc2_vel_2.^2 .+ df.sc2_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc2_pos_1.^2 .+ df.sc2_pos_2.^2 .+ df.sc2_pos_3.^2)
    Δeps1 = last(eps1) - first(eps1)
    Δeps2 = last(eps2) - first(eps2)
    @test Δeps1 > 2e3
    @test Δeps2 < -2e3
end

include(joinpath(REPO_ROOT, "test", "coverage_parallel_telemetry_probes.jl"))
include(joinpath(REPO_ROOT, "test", "coverage_targeted_90_probes.jl"))

@testset "JET Static Analysis" begin
    JET.@test_opt InitialCondition()
    JET.@test_opt Link()
    JET.@test_opt Joint()
    JET.@test_opt SpacecraftModel()

    sc = SpacecraftModel()
    JET.@test_opt rotate_to_inertial(sc, sc.root, 1)
end

@testset "Aqua Package Quality" begin
    Aqua.test_all(
        SimulationModel;
        ambiguities=false,
        stale_deps=false,
        deps_compat=false,
        project_extras=false,
        piracies=false,
        persistent_tasks=false,
        undocumented_names=false
    )
end
