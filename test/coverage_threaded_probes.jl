using Test
using StaticArrays
using ComponentArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

const quat_mult = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

struct ProbeDensityModel <: SimulationModel.AbstractDensityModel
end

mutable struct ProbeControlModel <: SimulationModel.AbstractControlEffectorModel
    hits::Vector{Int}
end

mutable struct ProbeGuidanceModel <: SimulationModel.AbstractTypes.AbstractGuidanceModel
    hits::Vector{Int}
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
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "3"
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
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "3",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1"
    ) do
        @test callbacks._control_callback_use_threads(probe_control, 4, false) == true
    end

    # Density callback threaded branch (line with Threads.@threads).
    p_density = ODEParams{4}(args=thread_safe_args)
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
    p_control = ODEParams{4}(args=args_control)
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
    p_aero = ODEParams{1}(args=args_aero)
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
    p_nbody = ODEParams{1}(args=args_nbody)
    x_nbody = build_initial_conditions(args_nbody).sc[1]
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
        force_nbody, torque_nbody = calcForceTorque(nbody, x_nbody, p_nbody, 1)
        @test all(isfinite, force_nbody)
        @test torque_nbody == SVector{3, Float64}(0.0, 0.0, 0.0)
    end
end

println("coverage_threaded_probes_ok")
