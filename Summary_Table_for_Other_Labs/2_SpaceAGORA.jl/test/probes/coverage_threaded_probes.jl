using Test
using Dates
using DiffEqBase
using DiffEqCallbacks
using OrdinaryDiffEq
using Quaternions
using Serialization
using StaticArrays
using ComponentArrays
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel
include(joinpath(REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))

const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end
if !isdefined(@__MODULE__, :build_initial_conditions)
    const build_initial_conditions = SimulationEngine.build_initial_conditions
end

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

const HAS_GRAMSUITE = let
    vendored_gramsuite = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
    try
        if Base.find_package("GRAMSuite") === nothing && isdir(vendored_gramsuite)
            pushfirst!(LOAD_PATH, vendored_gramsuite)
        end
        @eval import GRAMSuite
        # Importing the package is not enough: the GRAM-backed probes construct
        # real models, which needs the native GRAM Suite root (Build/ + Julia/).
        # Dev machines often have the Julia wrapper but no native build; skip
        # cleanly there instead of erroring mid-testset. CI builds the root, so
        # this changes nothing in CI.
        @eval GRAMSuite._resolve_gram_root("", "")
        true
    catch err
        @info "Skipping GRAMSuite-backed threaded coverage probes" exception=(err, catch_backtrace())
        false
    end
end

if HAS_GRAMSUITE
    const EM = SimulationModel.EnvironmentModels
    const TEST_GRAM_LOCK = ReentrantLock()

    function EM.GRAMAtmosphereModel(; kwargs...)
        return EM.GRAMAtmosphereModel(GRAMSuite.GRAMAtmosphereModel(; kwargs...))
    end

    function Base.deepcopy_internal(model::EM.GRAMAtmosphereModel, stackdict::IdDict)
        haskey(stackdict, model) && return stackdict[model]
        copied = EM.GRAMAtmosphereModel(deepcopy(model.core))
        stackdict[model] = copied
        return copied
    end

    function EM._gram_core_density_state(
        core::GRAMSuite.GRAMAtmosphereModel,
        h::Float64,
        lat::Float64,
        lon::Float64,
        el_time::Float64,
        wind::Bool,
        lock_obj,
        vacuum_temperature::Float64
    )::Tuple{Float64, Float64, SVector{3, Float64}}
        return GRAMSuite.density_state(
            core,
            h,
            lat,
            lon,
            el_time,
            wind;
            lock_obj=lock_obj,
            vacuum_temperature=vacuum_temperature
        )
    end

    @inline function EM._gram_point_density(
        model::EM.GRAMAtmosphereModel,
        h::Float64,
        lat::Float64,
        lon::Float64,
        el_time::Float64,
        wind::Bool
    )::Tuple{Float64, Float64, SVector{3, Float64}}
        h_gram = max(h, -30.0)
        return GRAMSuite.point_density_state(model.core, h_gram, lat, lon, el_time, wind; lock_obj=TEST_GRAM_LOCK)
    end

    function EM.getDensity(
        model::EM.GRAMAtmosphereModel,
        h::Float64,
        lat::Float64,
        lon::Float64,
        el_time::Float64,
        wind::Bool,
        p::params
    )::Tuple{Float64, Float64, SVector{3, Float64}} where {params}
        EI = p.args.environment_model.EI * 1e3
        drag_state = h - EI <= 0.0

        if h > 2000.0e3
            rho = 0.0
            T = p.args.environment_model.planet.T_ref
            wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        elseif !drag_state && !p.args.mission_configuration.keplerian
            rho, T, wind_vec = EM.density_polyfit(h, p)
        else
            h_gram = max(h, -30.0)
            rho, T, wind_vec = GRAMSuite.density_state(
                model.core,
                h_gram,
                lat,
                lon,
                el_time,
                wind;
                lock_obj=TEST_GRAM_LOCK,
                vacuum_temperature=p.args.environment_model.planet.T_ref
            )
        end

        return rho, T, wind_vec
    end
end

struct ProbeDensityModel <: SimulationModel.AbstractDensityModel
end

struct ThrowDensityModel <: SimulationModel.AbstractDensityModel
end

mutable struct ProbeControlModel <: SimulationModel.AbstractControlEffectorModel
    hits::Vector{Int}
end

mutable struct ProbeGuidanceModel <: SimulationModel.AbstractTypes.AbstractGuidanceModel
    hits::Vector{Int}
end

struct ProbeImplicitForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
    torque::SVector{3, Float64}
end

function SimulationModel.calcControlEffect!(
    model::ProbeControlModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcGuidanceEffect!(
    model::ProbeGuidanceModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

SimulationModel.solver_partition(::ProbeImplicitForceModel) = :implicit

function SimulationModel.calcForceTorque(
    model::ProbeImplicitForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return model.force, model.torque
end

function SimulationModel.EnvironmentModels.getDensity(
    ::ThrowDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)::Tuple{Float64, Float64, SVector{3, Float64}}
    error("throw_density_probe")
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

function build_config(;
    spacecraft::SpacecraftModel,
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effecters::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    guidance_effecters::Tuple=(),
    guidance_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::InitialTime=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
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
        guidance_model=GuidanceModel(guidance_effecters, guidance_rates),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effecters, control_rates),
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
    control_effecters::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::InitialTime=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
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
        control_model=ControlModel(control_effecters, control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

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

@testset "Threaded Coverage Probes" begin
    @test Threads.nthreads() > 1

    callbacks = SimulationModel.SimulationCallbacks
    dyn = SimulationModel.DynamicEffectors

    # Drive callback helper branches that are gated behind multi-thread execution.
    thread_safe_args = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0),
            make_spacecraft(ra_alt_m=600e3, rp_alt_m=520e3, ν_deg=150.0),
            make_spacecraft(ra_alt_m=650e3, rp_alt_m=540e3, ν_deg=140.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),)
    )

    non_threadsafe_args = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=ProbeDensityModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),)
    )

    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off") do
        @test callbacks._density_callback_use_threads(thread_safe_args, 4) == false
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_ASSUME_THREADSAFE" => "0"
    ) do
        @test callbacks._density_callback_use_threads(non_threadsafe_args, 2) == false
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "8"
    ) do
        @test callbacks._density_callback_use_threads(thread_safe_args, 4) == true
    end
    # Pin the auto-mode budget floor (default 16 exceeds CI thread counts).
    # thread_safe_args uses NoAtmosphereModel, which routes through the
    # lock-free density-callback source (see the source= selection in
    # _density_callback_thread_decision), not the GRAM-locked one -- so the
    # lock-free-specific knob is what actually gates this decision.
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "3",
        "SPACEAGORA_DENSITY_CALLBACK_LOCKFREE_AUTO_THREAD_MIN_BUDGET" => "2"
    ) do
        @test callbacks._density_callback_use_threads(thread_safe_args, 4) == true
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test callbacks._density_callback_use_threads(thread_safe_args, 4) == false
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_LOCKFREE_AUTO_THREAD_MIN_BUDGET" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1"
    ) do
        @test callbacks._density_callback_use_threads(thread_safe_args, 4) == true
    end

    probe_control = ProbeControlModel(zeros(Int, 4))
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "0"
    ) do
        @test callbacks._control_callback_use_threads(probe_control, 4, false) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1"
    ) do
        @test callbacks._control_callback_use_threads(probe_control, 4, false) == true
    end
    # Pin the default auto-mode budget floor (4) at the 2-thread CI probe budget.
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "3",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1",
        "SPACEAGORA_AUTO_THREAD_MIN_BUDGET" => "2"
    ) do
        @test callbacks._control_callback_use_threads(probe_control, 4, false) == true
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test callbacks._control_callback_use_threads(probe_control, 4, false) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1",
        "SPACEAGORA_AUTO_THREAD_MIN_BUDGET" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1"
    ) do
        @test callbacks._control_callback_use_threads(probe_control, 4, false) == true
    end

    # Density callback threaded branch (line with Threads.@threads).
    p_density = ODEParams(n_sats=4, args=thread_safe_args)
    u_density = build_initial_conditions(thread_safe_args)
    integrator_density = MockCallbackIntegrator(
        p_density,
        u_density,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        density_cb = callbacks.get_density_callback(4, thread_safe_args)
        density_cb.affect!(integrator_density)
    end
    @test all(isfinite, p_density.shared_buffers.densities)

    if HAS_GRAMSUITE
        p_density_models = ODEParams(n_sats=4, args=thread_safe_args)
        u_density_models = build_initial_conditions(thread_safe_args)
        append!(p_density_models.shared_buffers.density_models, fill(GRAMAtmosphereModel(planet_name="earth"), 4))
        integrator_density_models = MockCallbackIntegrator(
            p_density_models,
            u_density_models,
            0.0,
            MockCallbackOpts(1.0, 1e-8, 1e-8),
            1,
            Inf
        )
        withenv(
            "SPACEAGORA_DENSITY_BATCH_PARALLEL" => "off",
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
            "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1"
        ) do
            density_cb = callbacks.get_density_callback(4, thread_safe_args)
            density_cb.affect!(integrator_density_models)
        end
        @test all(isfinite, p_density_models.shared_buffers.densities)
    end

    # Guidance invokelatest branch.
    probe_guidance = ProbeGuidanceModel([0])
    args_guidance = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        guidance_effecters=(probe_guidance,),
        guidance_rates=[1.0]
    )
    p_guidance = ODEParams(n_sats=1, args=args_guidance)
    u_guidance = build_initial_conditions(args_guidance)
    integrator_guidance = MockCallbackIntegrator(
        p_guidance,
        u_guidance,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv("SPACEAGORA_DEV_HOT_RELOAD" => "1") do
        guidance_cbs = callbacks.get_guidance_callbacks(1, args_guidance)
        guidance_cbs[1].affect!.affect!(integrator_guidance)
    end
    @test probe_guidance.hits == [1]

    # Control callback threaded branch (line with Threads.@threads).
    args_control = build_config_multi(
        spacecraft=thread_safe_args.dynamics_model.spacecraft,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effecters=(probe_control,),
        control_rates=[1.0]
    )
    p_control = ODEParams(n_sats=4, args=args_control)
    u_control = build_initial_conditions(args_control)
    integrator_control = MockCallbackIntegrator(
        p_control,
        u_control,
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
        control_cbs = callbacks.get_control_callbacks(4, args_control)
        control_cbs[1].affect!.affect!(integrator_control)
    end
    @test probe_control.hits == ones(Int, 4)

    # Aerodynamic helper branches and threaded accumulation branch.
    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test dyn._multibody_use_threads(8) == false
    end
    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY" => "1"
    ) do
        @test dyn._multibody_use_threads(8; heavy_work=false) == false
    end
    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_AUTO_THREAD_MIN_BUDGET" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY" => "0"
    ) do
        @test dyn._multibody_use_threads(8; heavy_work=true) == true
    end

    args_aero = build_config(
        spacecraft=make_spacecraft(
            ra_alt_m=520e3,
            rp_alt_m=420e3,
            ν_deg=165.0,
            orientation_state=(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), SVector{3, Float64}(0.0, 0.0, 0.0))
        ),
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM())
    )
    p_aero = ODEParams(n_sats=1, args=args_aero)
    u_aero = build_initial_conditions(args_aero)
    x_aero = u_aero.sc[1]
    p_aero.shared_buffers.densities[1] = 1e-6
    p_aero.shared_buffers.temperatures[1] = 220.0
    p_aero.shared_buffers.winds[1] = SVector{3, Float64}(0.0, 0.0, 0.0)

    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
        force_aero, torque_aero = calcForceTorque(AerodynamicCoefficientfM(), x_aero, p_aero, 1)
        @test length(force_aero) == 3
        @test length(torque_aero) == 3
    end

    p_aero.shared_buffers.densities[1] = NaN
    force_nan, torque_nan = calcForceTorque(AerodynamicCoefficientfM(), x_aero, p_aero, 1)
    @test force_nan == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test torque_nan == SVector{3, Float64}(0.0, 0.0, 0.0)

    # N-body perturbation threaded accumulation branch.
    nbody = NBodyGravityModel(["Sun", "Moon"], "Earth", SPICE_PATH)
    args_nbody = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=175.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),)
    )
    p_nbody = ODEParams(n_sats=1, args=args_nbody)
    x_nbody = build_initial_conditions(args_nbody).sc[1]
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
        force_nbody, torque_nbody = calcForceTorque(nbody, x_nbody, p_nbody, 1)
        @test all(isfinite, force_nbody)
        @test torque_nbody == SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    callbacks._gram_runtime_stats_reset!()

    # ParallelPolicy private-branch probes.
    policy = SimulationModel.ParallelPolicy
    req_stop = Channel{Any}(1)
    done_stop = Channel{Any}(1)
    worker_stop = Threads.@spawn policy._persistent_foreach_worker_loop(1, req_stop, done_stop)
    put!(req_stop, :stop)
    wait(worker_stop)
    @test istaskdone(worker_stop)

    req_err = Channel{Any}(1)
    done_err = Channel{Any}(1)
    worker_err = Threads.@spawn policy._persistent_foreach_worker_loop(1, req_err, done_err)
    put!(req_err, (num_items=1, active_workers=1, f=(idx -> error("persistent_worker_probe"))))
    @test take!(done_err) isa Base.CapturedException
    put!(req_err, :stop)
    wait(worker_err)

    req_dynamic = Channel{Any}(1)
    done_dynamic = Channel{Any}(1)
    worker_dynamic = Threads.@spawn policy._persistent_foreach_worker_loop(1, req_dynamic, done_dynamic)
    next_index = Base.Threads.Atomic{Int}(1)
    put!(
        req_dynamic,
        (
            num_items=5,
            active_workers=2,
            scheduler=:dynamic,
            chunk=2,
            next_index=next_index,
            f=identity,
        )
    )
    @test take!(done_dynamic) === nothing
    put!(req_dynamic, :stop)
    wait(worker_dynamic)

    pool_direct = policy._create_persistent_foreach_pool(2)
    @test pool_direct.workers >= 2
    policy._shutdown_persistent_foreach_pool!(pool_direct)

    scope_ctx = policy.PolicyContext()
    scope_id = policy._policy_scope_id(scope_ctx)
    lock(policy._persistent_foreach_lock) do
        policy._persistent_foreach_pools[(scope_id, :probe_scope)] = policy._create_persistent_foreach_pool(2)
    end
    policy._destroy_persistent_foreach_scope!(scope_id)
    @test !haskey(policy._persistent_foreach_pools, (scope_id, :probe_scope))

    withenv("SPACEAGORA_PARALLEL_POLICY_DELTA" => "oops") do
        @test_throws ArgumentError policy.adaptive_delta()
    end
    withenv("SPACEAGORA_PARALLEL_POLICY_RHO" => "1.0") do
        @test_throws ArgumentError policy.adaptive_rho()
    end
    withenv("SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA" => "oops") do
        @test_throws ArgumentError policy.adaptive_trim_quanta_budget()
    end
    withenv("SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "guided") do
        @test_throws ArgumentError policy.inner_scheduler_mode()
    end

    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(max(2, Threads.nthreads()))) do
        acc_foreach = Base.Threads.Atomic{Int}(0)
        policy.threaded_foreach(8, 8) do idx
            Base.Threads.atomic_add!(acc_foreach, idx)
        end
        @test acc_foreach[] == sum(1:8)

        worker_ids = zeros(Int, 4)
        policy.threaded_foreach_worker(4, 1) do worker_id, idx
            worker_ids[idx] = worker_id
        end
        @test all(==(1), worker_ids)

        reduced_mt = policy.threaded_reduce(
            12,
            8,
            () -> Ref(0),
            (acc, idx) -> begin
                acc[] += idx
                return nothing
            end,
            (dest, src) -> begin
                dest[] += src[]
                return nothing
            end
        )
        @test reduced_mt[] == sum(1:12)

        pool_err = policy._create_persistent_foreach_pool(2)
        persistent_err = try
            policy._threaded_foreach_persistent!(pool_err, 4, 2, idx -> error("persistent_dispatch_probe"))
            nothing
        catch err
            err
        end
        @test persistent_err !== nothing
        policy._shutdown_persistent_foreach_pool!(pool_err)
    end

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(max(2, Threads.nthreads())),
        "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "dynamic",
        "SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK" => "2"
    ) do
        acc_dynamic = Base.Threads.Atomic{Int}(0)
        policy.threaded_foreach(9, 9) do idx
            Base.Threads.atomic_add!(acc_dynamic, idx)
        end
        @test acc_dynamic[] == sum(1:9)

        reduced_dynamic = policy.threaded_reduce(
            9,
            9,
            () -> Ref(0),
            (acc, idx) -> begin
                acc[] += idx
                return nothing
            end,
            (dest, src) -> begin
                dest[] += src[]
                return nothing
            end
        )
        @test reduced_dynamic[] == sum(1:9)
    end

    args_partition = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(
            ProbeImplicitForceModel(SVector{3, Float64}(1.0, 0.0, 0.0), SVector{3, Float64}(0.0, 1.0, 0.0)),
            ProbeImplicitForceModel(SVector{3, Float64}(2.0, 0.0, 0.0), SVector{3, Float64}(0.0, 2.0, 0.0)),
        )
    )
    p_partition = ODEParams(n_sats=1, args=args_partition)
    u_partition = build_initial_conditions(args_partition)
    f_partition = MVector{3, Float64}(0.0, 0.0, 0.0)
    τ_partition = MVector{3, Float64}(0.0, 0.0, 0.0)
    SimulationEngine._accumulate_dynamic_effectors_partitioned!(
        f_partition,
        τ_partition,
        u_partition.sc[1],
        p_partition,
        1,
        0.0,
        args_partition.dynamics_model.dynamic_effectors,
        (use_threads=true, allotment=2, policy_applied=false, mode=:on),
        :implicit
    )
    @test f_partition == MVector{3, Float64}(3.0, 0.0, 0.0)
    @test τ_partition == MVector{3, Float64}(0.0, 3.0, 0.0)

    policy.reset_policy_telemetry!()
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
        "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "3",
        "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
        "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "2"
    ) do
        _ = policy.thread_policy_decision(4; mode=:auto, threshold=1, source=:probe_obs)
        policy.record_policy_observation!(
            :probe_obs;
            mode=:auto,
            num_items=1,
            use_threads=false,
            elapsed_ns=1
        )
        snap = policy.policy_telemetry_snapshot()
        @test snap.adaptation_updates_total == 0
    end

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION" => "1.5",
        "SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES" => "2"
    ) do
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            state.loaded = true
            state.dirty = false
            state.path = "inner_hint_probe_state.toml"
            empty!(state.history)
            state.history["profile=r5|machine=probe_machine|src=density_callback|items=3_4|thr=1|budget=2|outer=0|heavy_only=0|heavy=1"] = Dict(
                Int64(1) => policy.AdaptiveChoiceStats(samples=2, successes=2, failures=0, elapsed_sum_ns=200.0, elapsed_sq_sum_ns=20_000.0),
                Int64(2) => policy.AdaptiveChoiceStats(samples=2, successes=2, failures=0, elapsed_sum_ns=120.0, elapsed_sq_sum_ns=7_200.0)
            )
            state.history["profile=r5|machine=probe_machine|src=control_callback|items=3_4|thr=1|budget=2|outer=0|heavy_only=0|heavy=1"] = Dict(
                Int64(1) => policy.AdaptiveChoiceStats(samples=1, successes=1, failures=0, elapsed_sum_ns=50.0, elapsed_sq_sum_ns=2_500.0)
            )
        end

        layer_rows = policy.hint_layer_stats_snapshot(profile="R5", machine="probe_machine")
        @test length(layer_rows) == 2
        @test all(row -> row.profile == "r5", layer_rows)
        @test all(row -> row.machine == "probe_machine", layer_rows)
        density_row = only([row for row in layer_rows if row.layer == "density_callback"])
        @test density_row.samples_total == 4
        @test density_row.regret_mean_ns >= 0.0
        @test density_row.confidence_mean >= 0.0
        @test density_row.state_path == "inner_hint_probe_state.toml"
        @test isempty(policy.hint_layer_stats_snapshot(profile="R4", machine="probe_machine"))
    end
    lock(policy._persistent_hint_lock) do
        state = policy._persistent_hint_state[]
        state.loaded = false
        state.dirty = false
        state.path = ""
        empty!(state.history)
    end

    withenv("SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => joinpath("relative_hint_state", "policy.toml")) do
        rel_hint_path = policy._persistent_hint_path()
        @test isabspath(rel_hint_path)
        @test endswith(rel_hint_path, normpath(joinpath("relative_hint_state", "policy.toml")))
    end

    mktempdir() do tmp
        missing_path = joinpath(tmp, "missing_inner_policy.toml")
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => missing_path
        ) do
            lock(policy._persistent_hint_lock) do
                state = policy._persistent_hint_state[]
                state.loaded = false
                state.dirty = false
                state.path = ""
                empty!(state.history)
                policy._load_persistent_hint_state_locked!()
                @test state.loaded
                @test state.path == normpath(missing_path)
                @test isempty(state.history)
            end
        end

        broken_path = joinpath(tmp, "broken_inner_policy.toml")
        write(broken_path, "history = [")
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => broken_path
        ) do
            lock(policy._persistent_hint_lock) do
                state = policy._persistent_hint_state[]
                state.loaded = false
                state.dirty = false
                state.path = ""
                empty!(state.history)
                policy._load_persistent_hint_state_locked!()
                @test state.loaded
                @test isempty(state.history)
            end
        end

        good_path = joinpath(tmp, "good_inner_policy.toml")
        loaded_signature = "profile=r5|machine=probe_machine|src=density_callback|items=3_4|thr=1|budget=2|outer=0|heavy_only=0|heavy=1"
        payload = Dict(
            "schema_version" => 1,
            "history" => Any[
                Dict(
                    "signature" => loaded_signature,
                    "allotment" => 2,
                    "stats" => Dict(
                        "samples" => 3,
                        "successes" => 2,
                        "failures" => 1,
                        "elapsed_sum_ns" => 90.0,
                        "elapsed_sq_sum_ns" => 2_700.0
                    )
                ),
                Dict(
                    "signature" => "",
                    "allotment" => 1,
                    "stats" => Dict(
                        "samples" => 2,
                        "successes" => 1,
                        "failures" => 1,
                        "elapsed_sum_ns" => 10.0,
                        "elapsed_sq_sum_ns" => 100.0
                    )
                ),
                Dict(
                    "signature" => loaded_signature,
                    "allotment" => "bad",
                    "stats" => Dict(
                        "samples" => 1,
                        "successes" => 1,
                        "failures" => 0,
                        "elapsed_sum_ns" => 1.0,
                        "elapsed_sq_sum_ns" => 1.0
                    )
                ),
                Dict(
                    "signature" => "profile=r5|machine=probe_machine|src=control_callback|items=3_4|thr=1|budget=2|outer=0|heavy_only=0|heavy=1",
                    "allotment" => 1,
                    "stats" => "bad"
                )
            ]
        )
        open(good_path, "w") do io
            TOML.print(io, payload)
        end
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => good_path
        ) do
            lock(policy._persistent_hint_lock) do
                state = policy._persistent_hint_state[]
                state.loaded = false
                state.dirty = false
                state.path = ""
                empty!(state.history)
                policy._load_persistent_hint_state_locked!()
                @test haskey(state.history, loaded_signature)
                @test state.history[loaded_signature][Int64(2)].samples == 3
                @test policy._hint_entry_count(state) == 1
            end
        end

        withenv(
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => good_path,
            "SPACEAGORA_PARALLEL_POLICY_STATE_RESET" => "1"
        ) do
            lock(policy._persistent_hint_lock) do
                state = policy._persistent_hint_state[]
                state.loaded = false
                state.dirty = false
                state.path = ""
                empty!(state.history)
                policy._load_persistent_hint_state_locked!()
                @test state.loaded
                @test state.path == normpath(good_path)
                @test isempty(state.history)
            end
            @test policy.persistent_hints_state_reset_requested()
        end

        withenv(
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
            "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => missing_path
        ) do
            previous_registered = policy._persistent_hint_atexit_registered[]
            policy._persistent_hint_atexit_registered[] = false
            lock(policy._persistent_hint_lock) do
                state = policy._persistent_hint_state[]
                state.loaded = false
                state.dirty = false
                state.path = ""
                empty!(state.history)
            end
            policy._ensure_persistent_hint_state_loaded!()
            @test policy._persistent_hint_atexit_registered[]
            policy._persistent_hint_atexit_registered[] = previous_registered || policy._persistent_hint_atexit_registered[]
        end
    end

    lock(policy._persistent_hint_lock) do
        state = policy._persistent_hint_state[]
        state.loaded = true
        state.dirty = false
        state.path = "inner_hint_manual_reset.toml"
        state.history["manual_reset_sig"] = Dict(Int64(1) => policy.AdaptiveChoiceStats(samples=1, successes=1, failures=0, elapsed_sum_ns=1.0, elapsed_sq_sum_ns=1.0))
    end
    policy.reset_persistent_hint_state!()
    lock(policy._persistent_hint_lock) do
        state = policy._persistent_hint_state[]
        @test !state.loaded
        @test isempty(state.path)
        @test isempty(state.history)
    end

    @test policy._hint_bucket(12) == "9_16"
    @test policy._hint_bucket(24) == "17p"
    hint_candidates = policy._hint_candidate_allotments(12, 12)
    @test 4 in hint_candidates
    @test 8 in hint_candidates

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => "inner_hint_choose_probe.toml"
    ) do
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            state.loaded = true
            state.dirty = false
            state.path = "inner_hint_choose_probe.toml"
            empty!(state.history)
            state.history["sig_zero"] = Dict(Int64(1) => policy.AdaptiveChoiceStats())
            state.history["sig_choose"] = Dict(
                Int64(1) => policy.AdaptiveChoiceStats(samples=4, successes=4, failures=0, elapsed_sum_ns=360.0, elapsed_sq_sum_ns=32_400.0),
                Int64(2) => policy.AdaptiveChoiceStats(samples=4, successes=3, failures=1, elapsed_sum_ns=120.0, elapsed_sq_sum_ns=3_600.0)
            )
        end
        miss_choice = policy._hint_choose_allotment("sig_missing", Int64[1, 2])
        @test miss_choice.allotment == 1
        zero_choice = policy._hint_choose_allotment("sig_zero", Int64[1])
        @test zero_choice.allotment == 1
        chosen = policy._hint_choose_allotment("sig_choose", Int64[1, 2])
        @test chosen.allotment == 2
        @test chosen.confidence >= 0.0
        @test chosen.regret_ns >= 0.0

        policy._hint_record_observation!("sig_obs", Int64(0), Int64(10); success=true)
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            @test !haskey(state.history, "sig_obs")
        end
        policy._hint_record_observation!("sig_obs", Int64(2), Int64(10); success=false)
        policy._hint_record_observation!("sig_obs", Int64(2), Int64(20); success=true)
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            obs = state.history["sig_obs"][Int64(2)]
            @test obs.samples == 2
            @test obs.successes == 1
            @test obs.failures == 1
        end
    end

    withenv("SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0") do
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            state.loaded = true
            state.dirty = false
            state.path = "inner_hint_disabled_probe.toml"
            empty!(state.history)
        end
        policy._hint_record_observation!("sig_disabled", Int64(2), Int64(15); success=true)
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            @test isempty(state.history)
        end
    end

    @test policy._hint_signature_value("profile=|machine=probe", "profile") == ""
    @test isnothing(policy._hint_signature_value("profile=r5|machine=probe", "missing"))

    lock(policy._persistent_hint_lock) do
        state = policy._persistent_hint_state[]
        state.loaded = true
        state.dirty = false
        state.path = "inner_hint_snapshot_probe.toml"
        empty!(state.history)
    end
    @test isempty(policy.hint_layer_stats_snapshot())
    lock(policy._persistent_hint_lock) do
        state = policy._persistent_hint_state[]
        state.history["profile=r5|machine=machine_a|src=density_callback|items=1|thr=1|budget=1|outer=0|heavy_only=0|heavy=0"] = Dict(
            Int64(1) => policy.AdaptiveChoiceStats(samples=1, successes=1, failures=0, elapsed_sum_ns=10.0, elapsed_sq_sum_ns=100.0)
        )
    end
    filtered_rows = policy.hint_layer_stats_snapshot(profile="r5", machine="machine_b")
    @test filtered_rows isa Vector
    @test isempty(filtered_rows)

    withenv("SPACEAGORA_PARALLEL_POLICY_DELTA" => "1.2") do
        @test_throws ArgumentError policy.adaptive_delta()
    end

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
        "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "2",
        "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
        "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
        "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "1",
        "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "1",
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "4"
    ) do
        lock(policy._policy_telemetry_lock) do
            ctx = policy._active_policy_context()
            ctx.telemetry = policy.PolicyTelemetry()
            empty!(ctx.adaptive_state)
            empty!(ctx.decision_signature)
            empty!(ctx.decision_allotment)
        end
        signature_control = policy._hint_workload_signature(:control_callback, 4, 1, 4, false, false, true)
        lock(policy._persistent_hint_lock) do
            state = policy._persistent_hint_state[]
            state.loaded = true
            state.dirty = false
            state.path = "adaptive_hint_probe.toml"
            empty!(state.history)
            state.history[signature_control] = Dict(
                Int64(1) => policy.AdaptiveChoiceStats(samples=3, successes=3, failures=0, elapsed_sum_ns=300.0, elapsed_sq_sum_ns=30_000.0),
                Int64(2) => policy.AdaptiveChoiceStats(samples=3, successes=2, failures=1, elapsed_sum_ns=90.0, elapsed_sq_sum_ns=4_500.0)
            )
        end

        decision = policy.thread_policy_decision(
            4;
            mode=:auto,
            threshold=1,
            source=:control_callback,
            outer_active=false,
            allow_with_outer=true,
            heavy_only=false,
            heavy_work=true
        )
        @test decision.use_threads isa Bool
        @test decision.allotment >= 1
        if decision.use_threads
            @test decision.allotment >= 2
        else
            @test decision.allotment == 1
        end

        policy.record_policy_observation!(
            :control_callback;
            mode=:auto,
            num_items=4,
            use_threads=decision.use_threads,
            elapsed_ns=10
        )
        policy.record_policy_observation!(
            :control_callback;
            mode=:auto,
            num_items=4,
            use_threads=decision.use_threads,
            elapsed_ns=10
        )
        snap = policy.policy_telemetry_snapshot()
        @test snap.persistent_hints_updates >= 1
        @test snap.persistent_hints_entries >= 1
        @test snap.persistent_hints_loaded == true
        @test snap.persistent_hints_path == "adaptive_hint_probe.toml"
    end

    # Density_models helper probes.
    env_models = SimulationModel.EnvironmentModels
    args_density_helpers = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=false
    )
    p_density_helpers = ODEParams(n_sats=1, args=args_density_helpers)

    rho_no, T_no, wind_no = env_models.getDensity(NoAtmosphereModel(), 150e3, 0.0, 0.0, 0.0, true, p_density_helpers)
    @test rho_no == 0.0
    @test T_no == args_density_helpers.environment_model.planet.T_ref
    @test wind_no == SVector{3, Float64}(0.0, 0.0, 0.0)

    exp_model = env_models.ExponentialAtmosphereModel(1e-4, 120e3, 12e3)
    rho_exp, T_exp, wind_exp = env_models.getDensity(exp_model, 120e3, 0.0, 0.0, 0.0, false)
    @test isapprox(rho_exp, 1e-4; atol=0.0, rtol=1e-12)
    @test T_exp == 200.0
    @test wind_exp == SVector{3, Float64}(0.0, 0.0, 0.0)

    poly_model = env_models.PolynomialFitAtmosphereModel(
        [0.0, -1.0];
        valid_min_altitude_m=100e3,
        valid_max_altitude_m=200e3
    )
    @test_throws ArgumentError env_models.PolynomialFitAtmosphereModel(
        [0.0, -1.0];
        valid_min_altitude_m=200e3,
        valid_max_altitude_m=100e3
    )
    @test poly_model.valid_min_altitude_m == 100e3
    @test poly_model.valid_max_altitude_m == 200e3

    poly_planet = env_models.PolynomialFitAtmosphereModel(args_density_helpers.environment_model.planet)
    @test poly_planet.valid_min_altitude_m == 50e3
    @test poly_planet.valid_max_altitude_m == 2_000e3
    rho_scalar, T_scalar, wind_scalar = env_models.getDensity(poly_model, 120e3, 0.0, 0.0, 0.0, false, p_density_helpers)
    @test isfinite(rho_scalar)
    @test T_scalar == args_density_helpers.environment_model.planet.T_ref
    @test wind_scalar == SVector{3, Float64}(0.0, 0.0, 0.0)

    rho_vec, T_vec, wind_vec = env_models.getDensity(poly_model, [80e3, 120e3, 250e3], 0.0, 0.0, 0.0, false, p_density_helpers)
    @test rho_vec isa Vector{Float64}
    @test length(rho_vec) == 3
    @test all(isfinite, rho_vec)
    @test T_vec == args_density_helpers.environment_model.planet.T_ref
    @test wind_vec == SVector{3, Float64}(0.0, 0.0, 0.0)

    rho_at_min, _, _ = env_models.getDensity(poly_model, 100e3, 0.0, 0.0, 0.0, false, p_density_helpers)
    rho_below_min, _, _ = env_models.getDensity(poly_model, 80e3, 0.0, 0.0, 0.0, false, p_density_helpers)
    rho_at_max, _, _ = env_models.getDensity(poly_model, 200e3, 0.0, 0.0, 0.0, false, p_density_helpers)
    rho_above_max, _, _ = env_models.getDensity(poly_model, 250e3, 0.0, 0.0, 0.0, false, p_density_helpers)
    @test rho_below_min == rho_at_min
    @test rho_above_max == rho_at_max

    poly_overflow = env_models.PolynomialFitAtmosphereModel(
        [1.0e6];
        valid_min_altitude_m=100e3,
        valid_max_altitude_m=200e3
    )
    rho_overflow, _, _ = env_models.getDensity(poly_overflow, 150e3, 0.0, 0.0, 0.0, false, p_density_helpers)
    @test isfinite(rho_overflow)
    @test rho_overflow > 0.0

    hs_batch = [120e3, 130e3, 140e3]
    lats_batch = [0.0, 0.05, -0.02]
    lons_batch = [0.0, 0.1, -0.2]
    ts_batch = [0.0, 20.0, 40.0]

    rhos_batch = zeros(Float64, length(hs_batch))
    Ts_batch = zeros(Float64, length(hs_batch))
    winds_batch = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in eachindex(hs_batch)]
    env_models.getDensityBatch!(
        rhos_batch,
        Ts_batch,
        winds_batch,
        NoAtmosphereModel(),
        hs_batch,
        lats_batch,
        lons_batch,
        ts_batch,
        true,
        p_density_helpers
    )
    @test all(==(0.0), rhos_batch)
    @test all(x -> x == args_density_helpers.environment_model.planet.T_ref, Ts_batch)
    @test all(==(SVector{3, Float64}(0.0, 0.0, 0.0)), winds_batch)

    rhos_poly_batch = zeros(Float64, length(hs_batch))
    Ts_poly_batch = zeros(Float64, length(hs_batch))
    winds_poly_batch = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in eachindex(hs_batch)]
    env_models.getDensityBatch!(
        rhos_poly_batch,
        Ts_poly_batch,
        winds_poly_batch,
        poly_model,
        hs_batch,
        lats_batch,
        lons_batch,
        0.0,
        false,
        p_density_helpers
    )
    @test all(isfinite, rhos_poly_batch)
    @test all(x -> x == args_density_helpers.environment_model.planet.T_ref, Ts_poly_batch)
    @test all(==(SVector{3, Float64}(0.0, 0.0, 0.0)), winds_poly_batch)
    @inbounds for i in eachindex(hs_batch)
        rho_i, T_i, wind_i = env_models.getDensity(poly_model, hs_batch[i], lats_batch[i], lons_batch[i], 0.0, false, p_density_helpers)
        @test isapprox(rhos_poly_batch[i], rho_i; atol=0.0, rtol=1e-12)
        @test Ts_poly_batch[i] == T_i
        @test winds_poly_batch[i] == wind_i
    end

    @test_throws ArgumentError env_models.getDensityBatch!(
        zeros(Float64, 2),
        zeros(Float64, 3),
        [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:3],
        poly_model,
        hs_batch,
        lats_batch,
        lons_batch,
        0.0,
        false,
        p_density_helpers
    )

    @test isfinite(env_models.interp(5.0, 355.0, 0.25))
    @test isfinite(env_models.interp(355.0, 5.0, 0.25))
    @test env_models.temperature_linear(10.0, (T_ref=123.0,)) == 123.0
    @test env_models._gram_use_global_lock() isa Bool
    @test_throws MethodError env_models._gram_point_density(:bad_model, 0.0, 0.0, 0.0, 0.0, false)

    gram_model = nothing
    if HAS_GRAMSUITE
        gram_model = env_models.GRAMAtmosphereModel(
            planet_name="earth",
            initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
        )
        @test gram_model.core === gram_model.core
        @test gram_model.planet_name == "earth"
        copied_gram = Base.deepcopy_internal(gram_model, IdDict())
        @test copied_gram !== gram_model

        rho_hi, T_hi, wind_hi = env_models.getDensity(gram_model, 2_500e3, 0.0, 0.0, 0.0, true, p_density_helpers)
        @test rho_hi == 0.0
        @test T_hi == args_density_helpers.environment_model.planet.T_ref
        @test wind_hi == SVector{3, Float64}(0.0, 0.0, 0.0)

        rho_mid, T_mid, wind_mid = env_models.getDensity(gram_model, 150e3, 0.0, 0.0, 0.0, true, p_density_helpers)
        @test isfinite(rho_mid)
        @test T_mid == args_density_helpers.environment_model.planet.T_ref
        @test wind_mid == SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    rho_polyfit, T_polyfit, wind_polyfit = env_models.density_polyfit(150e3, p_density_helpers)
    @test isfinite(rho_polyfit)
    @test T_polyfit == args_density_helpers.environment_model.planet.T_ref
    @test wind_polyfit == SVector{3, Float64}(0.0, 0.0, 0.0)

    # Deep callback branch probes around entry targeting and cache refresh.
    withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "off") do
        @test callbacks._gram_entry_target_mode() == :off
    end
    withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "allen_eggers") do
        @test callbacks._gram_entry_target_mode() == :allen_eggers
    end
    withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "invalid_mode") do
        @test_throws ArgumentError callbacks._gram_entry_target_mode()
    end
    withenv("SPACEAGORA_GRAM_ENTRY_TARGET_CD" => "0.01") do
        @test callbacks._gram_entry_target_cd() == 0.05
    end
    withenv("SPACEAGORA_GRAM_ENTRY_TARGET_DT_S" => "0.01") do
        @test callbacks._gram_entry_target_dt() == 0.05
    end
    withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MAX_STEPS" => "4") do
        @test callbacks._gram_entry_target_max_steps() == 8
    end
    withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MAX_STEPS" => "oops") do
        @test_throws ArgumentError callbacks._gram_entry_target_max_steps()
    end
    withenv("SPACEAGORA_CB_TEST_FLOAT_OPT_VALID" => "1.25") do
        @test callbacks._parse_float_env_optional("SPACEAGORA_CB_TEST_FLOAT_OPT_VALID") == 1.25
    end
    withenv("SPACEAGORA_CB_TEST_INT_OPT_VALID" => "7") do
        @test callbacks._parse_int_env_optional("SPACEAGORA_CB_TEST_INT_OPT_VALID") == 7
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "1") do
        @test callbacks._gram_track_cache_max_npos() == 2
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "oops") do
        @test_throws ArgumentError callbacks._gram_track_cache_max_npos()
    end
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "invalid_mode") do
        @test_throws ArgumentError callbacks._gram_isolated_pool_mode()
    end
    withenv(
        "SPACEAGORA_GRAM_ISOLATED_POOL" => "auto",
        "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "2"
    ) do
        @test callbacks._gram_isolated_pool_enabled(2) == (Threads.nthreads() > 1)
    end
    if HAS_GRAMSUITE
        withenv(
            "SPACEAGORA_GRAM_ISOLATED_POOL" => "off",
            "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => "2"
        ) do
            rho_pool = zeros(Float64, length(hs_batch))
            T_pool = zeros(Float64, length(hs_batch))
            wind_pool = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in eachindex(hs_batch)]
            @test !callbacks._gram_isolated_pool_batch_eval!(
                rho_pool,
                T_pool,
                wind_pool,
                gram_model,
                hs_batch,
                lats_batch,
                lons_batch,
                ts_batch,
                true,
                p_density_helpers
            )
        end
        withenv(
            "SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
            "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => "2"
        ) do
            rho_pool = zeros(Float64, length(hs_batch))
            T_pool = zeros(Float64, length(hs_batch))
            wind_pool = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in eachindex(hs_batch)]
            pooled = callbacks._gram_isolated_pool_batch_eval!(
                rho_pool,
                T_pool,
                wind_pool,
                gram_model,
                hs_batch,
                lats_batch,
                lons_batch,
                ts_batch,
                true,
                p_density_helpers;
                allotment_hint=2
            )
            @test pooled || Threads.nthreads() == 1
            @test all(isfinite, rho_pool)
            @test all(isfinite, T_pool)
            @test all(isfinite, getindex.(wind_pool, 1))
            if Threads.nthreads() > 1
                @test length(p_density_helpers.shared_buffers.gram_isolated_pool_models) >= 2
                @test length(p_density_helpers.shared_buffers.gram_isolated_pool_locks) >= 2
            end
        end
    end
    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on") do
        rho_pool = zeros(Float64, length(hs_batch))
        T_pool = zeros(Float64, length(hs_batch))
        wind_pool = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in eachindex(hs_batch)]
        @test !callbacks._gram_isolated_pool_batch_eval!(
            rho_pool,
            T_pool,
            wind_pool,
            NoAtmosphereModel(),
            hs_batch,
            lats_batch,
            lons_batch,
            ts_batch,
            true,
            p_density_helpers
        )
    end
    @test callbacks._gram_track_cache_target_spacing_m(500.0, deg2rad(1.0), 6.5e6) >= 1.0

    cache_unset = callbacks.GramTrackCache()
    @test callbacks._gram_track_cache_segment(cache_unset, 0.0) === nothing

    cache_seek = callbacks.GramTrackCache()
    cache_seek.valid = true
    cache_seek.t0 = 0.0
    cache_seek.t1 = 2.0
    cache_seek.index_hint = 1
    cache_seek.times = [0.0, 1.0, 2.0]
    cache_seek.alts = [100.0, 200.0, 300.0]
    cache_seek.lats = [0.0, 0.1, 0.2]
    cache_seek.lons = [0.0, 0.1, 0.2]
    cache_seek.rhos = [1.0, 2.0, 3.0]
    cache_seek.Ts = [200.0, 220.0, 240.0]
    cache_seek.winds = [SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)]
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "1") do
        @test callbacks._gram_track_cache_segment(cache_seek, 1.8) == (2, 0.8)
    end

    cache_nonmono = callbacks.GramTrackCache()
    cache_nonmono.valid = true
    cache_nonmono.t0 = 0.0
    cache_nonmono.t1 = 1.0
    cache_nonmono.index_hint = 2
    cache_nonmono.times = [0.0, 1.0, 0.5]
    cache_nonmono.alts = [100.0, 200.0, 300.0]
    cache_nonmono.lats = [0.0, 0.1, 0.2]
    cache_nonmono.lons = [0.0, 0.1, 0.2]
    cache_nonmono.rhos = [1.0, 2.0, 3.0]
    cache_nonmono.Ts = [200.0, 220.0, 240.0]
    cache_nonmono.winds = [SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)]
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "1") do
        seg_nonmono = callbacks._gram_track_cache_segment(cache_nonmono, 0.75)
        @test seg_nonmono === nothing || seg_nonmono isa Tuple{Int, Float64}
    end

    @test callbacks._gram_entry_reference_area_m2(p_density_helpers, 1) > 0.0
    @test callbacks._gram_entry_reference_area_m2(p_density_helpers, 999) == 1.0
    @test callbacks._gram_entry_mass_kg(p_density_helpers, 1, 42.0) == 42.0
    @test callbacks._gram_entry_mass_kg(p_density_helpers, 1, -1.0) > 0.0
    @test callbacks._gram_entry_mass_kg(p_density_helpers, 999, -1.0) == 100.0

    sc_zero_area = make_spacecraft(ra_alt_m=450e3, rp_alt_m=420e3, ν_deg=175.0)
    @inbounds for link in sc_zero_area.links
        link.ref_area = 0.0
    end
    args_zero_area = build_config(
        spacecraft=sc_zero_area,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=false
    )
    p_zero_area = ODEParams(n_sats=1, args=args_zero_area)
    @test callbacks._gram_entry_reference_area_m2(p_zero_area, 1) == 1.0

    u_refresh = build_initial_conditions(args_density_helpers)
    planet_refresh = p_density_helpers.args.environment_model.planet
    planet_refresh.L_PI .= [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    pos_refresh = SVector{3, Float64}(planet_refresh.Rp_e + 250e3, 0.0, 0.0)
    vel_refresh = SVector{3, Float64}(0.0, 7.6e3, 0.0)
    rp_refresh, _ = r_intor_p!(pos_refresh, vel_refresh, p_density_helpers.args.environment_model.planet)
    alt_refresh, lat_refresh, lon_refresh = rtolatlong(rp_refresh, p_density_helpers.args.environment_model.planet)
    cache_refresh = callbacks.GramTrackCache()
    @test callbacks._gram_entry_reference_density(planet_refresh, 120e3) >= 0.0
    @test callbacks._gram_entry_target_allen_eggers(
        pos_refresh,
        SVector{3, Float64}(0.0, 7.6e3, 0.0),
        planet_refresh,
        5.0,
        500.0,
        2.0
    ) !== nothing
    @test callbacks._gram_entry_target_allen_eggers(
        pos_refresh,
        SVector{3, Float64}(0.0, 7.6e3, 0.0),
        planet_refresh,
        0.0,
        500.0,
        2.0
    ) === nothing

    withenv(
        "SPACEAGORA_GRAM_PROFILE" => "1",
        "SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => "0",
        "SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "64"
    ) do
        rho_ok, T_ok, _ = callbacks._gram_track_cache_refresh!(
            callbacks.GramTrackCache(),
            NoAtmosphereModel(),
            p_density_helpers,
            pos_refresh,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            alt_refresh,
            lat_refresh,
            lon_refresh,
            0.0,
            10.0,
            6,
            500.0,
            deg2rad(2.0),
            0.0,
            10.0,
            false,
            1,
            Float64(u_refresh.sc[1].mass)
        )
        @test isfinite(rho_ok)
        @test isfinite(T_ok)

        rho_a, T_a, _ = callbacks._gram_track_cache_refresh!(
            cache_refresh,
            NoAtmosphereModel(),
            p_density_helpers,
            pos_refresh,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            alt_refresh,
            lat_refresh,
            lon_refresh,
            0.0,
            30.0,
            8,
            500.0,
            deg2rad(2.0),
            0.0,
            20.0,
            false,
            1,
            Float64(u_refresh.sc[1].mass)
        )
        @test isfinite(rho_a)
        @test isfinite(T_a)

        rho_b, T_b, _ = callbacks._gram_track_cache_refresh!(
            cache_refresh,
            NoAtmosphereModel(),
            p_density_helpers,
            pos_refresh,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            alt_refresh,
            lat_refresh,
            lon_refresh,
            0.0,
            30.0,
            8,
            500.0,
            deg2rad(2.0),
            0.0,
            NaN,
            false,
            1,
            Float64(u_refresh.sc[1].mass)
        )
        @test isfinite(rho_b)
        @test isfinite(T_b)
    end

    mission_orbits = MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=false,
        number_of_orbits=1,
        mission_time=120.0,
        orientation_sim=false,
        num_steps_to_save=1000
    )
    args_orbit_refresh = SimulationConfiguration(
        simulation_settings=args_density_helpers.simulation_settings,
        mission_configuration=mission_orbits,
        environment_model=args_density_helpers.environment_model,
        dynamics_model=args_density_helpers.dynamics_model,
        guidance_model=args_density_helpers.guidance_model,
        navigation_model=args_density_helpers.navigation_model,
        control_model=args_density_helpers.control_model,
        initial_time=args_density_helpers.initial_time,
        integration_tolerances=args_density_helpers.integration_tolerances
    )
    p_orbit_refresh = ODEParams(n_sats=1, args=args_orbit_refresh)
    u_orbit_refresh = build_initial_conditions(args_orbit_refresh)
    pos_orbit = SVector{3, Float64}(u_orbit_refresh.sc[1].pos)
    vel_orbit = SVector{3, Float64}(u_orbit_refresh.sc[1].vel)
    rp_orbit, _ = r_intor_p!(pos_orbit, vel_orbit, p_orbit_refresh.args.environment_model.planet)
    alt_orbit, lat_orbit, lon_orbit = rtolatlong(rp_orbit, p_orbit_refresh.args.environment_model.planet)
    withenv("SPACEAGORA_GRAM_PROFILE" => "1", "SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => "0") do
        rho_orbit, T_orbit, _ = callbacks._gram_track_cache_refresh!(
            callbacks.GramTrackCache(),
            NoAtmosphereModel(),
            p_orbit_refresh,
            pos_orbit,
            vel_orbit,
            alt_orbit,
            lat_orbit,
            lon_orbit,
            0.0,
            40.0,
            8,
            500.0,
            deg2rad(2.0),
            0.0,
            NaN,
            false,
            1,
            Float64(u_orbit_refresh.sc[1].mass)
        )
        @test isfinite(rho_orbit)
        @test isfinite(T_orbit)
    end

    withenv("SPACEAGORA_GRAM_PROFILE" => "1", "SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => "1") do
        rho_entry_peri, T_entry_peri, _ = callbacks._gram_track_cache_refresh!(
            callbacks.GramTrackCache(),
            NoAtmosphereModel(),
            p_density_helpers,
            pos_refresh,
            vel_refresh,
            80e3,
            lat_refresh,
            lon_refresh,
            0.0,
            20.0,
            8,
            500.0,
            deg2rad(2.0),
            40e3,
            NaN,
            false,
            1,
            Float64(u_refresh.sc[1].mass)
        )
        @test isfinite(rho_entry_peri)
        @test isfinite(T_entry_peri)
    end

    withenv(
        "SPACEAGORA_GRAM_PROFILE" => "1",
        "SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => "1",
        "SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "off"
    ) do
        rho_entry_off, T_entry_off, _ = callbacks._gram_track_cache_refresh!(
            callbacks.GramTrackCache(),
            NoAtmosphereModel(),
            p_density_helpers,
            pos_refresh,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            80e3,
            lat_refresh,
            lon_refresh,
            0.0,
            20.0,
            8,
            500.0,
            deg2rad(2.0),
            40e3,
            15.0,
            false,
            1,
            Float64(u_refresh.sc[1].mass)
        )
        @test isfinite(rho_entry_off)
        @test isfinite(T_entry_off)
    end

    withenv(
        "SPACEAGORA_GRAM_PROFILE" => "1",
        "SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => "1",
        "SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "allen_eggers"
    ) do
        rho_entry_ae, T_entry_ae, _ = callbacks._gram_track_cache_refresh!(
            callbacks.GramTrackCache(),
            NoAtmosphereModel(),
            p_density_helpers,
            pos_refresh,
            SVector{3, Float64}(0.0, 0.0, 0.0),
            80e3,
            lat_refresh,
            lon_refresh,
            0.0,
            20.0,
            8,
            500.0,
            deg2rad(2.0),
            40e3,
            NaN,
            false,
            1,
            Float64(u_refresh.sc[1].mass)
        )
        @test isfinite(rho_entry_ae)
        @test isfinite(T_entry_ae)
    end

    callbacks._gram_track_cache_warning_emitted[] = false
    withenv("SPACEAGORA_GRAM_PROFILE" => "1", "SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => "0") do
        @test_throws ErrorException callbacks._gram_track_cache_refresh!(
            callbacks.GramTrackCache(),
            ThrowDensityModel(),
            p_density_helpers,
            pos_refresh,
            vel_refresh,
            alt_refresh,
            lat_refresh,
            lon_refresh,
            0.0,
            5.0,
            4,
            500.0,
            deg2rad(2.0),
            0.0,
            5.0,
            false,
            1,
            Float64(u_refresh.sc[1].mass)
        )
    end

    args_thermal_threaded = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=530e3, rp_alt_m=460e3, ν_deg=165.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),)
    )
    p_thermal_threaded = ODEParams(n_sats=2, args=args_thermal_threaded)
    u_thermal_threaded = build_initial_conditions(args_thermal_threaded)
    p_thermal_threaded.shared_buffers.densities .= 1e-6
    p_thermal_threaded.shared_buffers.temperatures .= 250.0
    withenv(
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test callbacks._thermal_callback_thread_decision(2).use_threads == false
    end
    # Pin the auto-mode budget floor (default 16 exceeds CI thread counts).
    withenv(
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_THERMAL_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1"
    ) do
        @test callbacks._thermal_callback_thread_decision(2).use_threads == true
    end
    withenv(
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        thermal_cb_threaded = callbacks.get_thermal_callback(2, args_thermal_threaded)
        thermal_cb_threaded.affect!((p=p_thermal_threaded, u=u_thermal_threaded, t=0.0))
    end

    args_quat = build_config_multi(
        spacecraft=[
            make_spacecraft(
                ra_alt_m=500e3,
                rp_alt_m=450e3,
                ν_deg=170.0,
                orientation_state=(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), SVector{3, Float64}(0.0, 0.0, 0.0))
            ),
            make_spacecraft(
                ra_alt_m=530e3,
                rp_alt_m=460e3,
                ν_deg=165.0,
                orientation_state=(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), SVector{3, Float64}(0.0, 0.0, 0.0))
            )
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),)
    )
    p_quat = ODEParams(n_sats=2, args=args_quat)
    u_quat = build_initial_conditions(args_quat)
    p_quat.is_active[1] = false
    u_quat.sc[2].q .= (NaN, NaN, NaN, NaN)
    quat_cb = callbacks.get_quaternion_projection_callback(2, args_quat)
    @test quat_cb.condition(u_quat, 0.0, (p=p_quat, u=u_quat, t=0.0))
    quat_cb.affect!((p=p_quat, u=u_quat, t=0.0))
    @test u_quat.sc[2].q == SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

    withenv("SPACEAGORA_ENTRY_TARGET_COUNT" => "1") do
        entry_cb = callbacks.get_entry_end_callback(1, args_density_helpers)
        @test entry_cb isa DiffEqBase.VectorContinuousCallback
    end
end

println("coverage_threaded_probes_ok")
