module EnvironmentWindFlagTests
# `EnvironmentModel.wind = false` means the simulation uses no atmospheric wind:
# density queries are made with `wind=false`, and the wind that reaches
# aerodynamics, the density buffers and the saved `wind` field is zero, for every
# density model. Native GRAM needs the masking because its own `wind=false`
# query still returns nominal winds.
#
# The same file pins the wind-history rule shared by the two native-GRAM pools
# (in-process isolated pool and process-backed density service): under `auto` a
# history-dependent wind request stays on the locked path, `on` accepts separate
# histories, and density-only runs stay eligible.
#
# The native-backed part is skipped when GRAM is not built on the host, unless
# SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES=1, in which case a missing GRAM fails.
using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using StaticArrays

include(joinpath(@__DIR__, "..", "..", "helpers", "native_probe_reporting.jl"))

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const EM = SM.EnvironmentModels
const CB = SM.SimulationCallbacks

const WF_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WF_GRAM_ROOT = joinpath(WF_REPO, "data", "GRAMSuite.jl")
const WF_NATIVE_ROOT = get(ENV, "GRAM_ROOT", joinpath(WF_GRAM_ROOT, "GRAM Suite 2.0"))
const WF_SPICE_PATH = joinpath(WF_NATIVE_ROOT, "SPICE")
const WF_GRAM_LIB = joinpath(WF_NATIVE_ROOT, "Build", "lib",
    Sys.isapple() ? "libGRAM.dylib" : (Sys.iswindows() ? "libGRAM.dll" : "libGRAM.so"))
const WF_GRAM_READY = isfile(WF_GRAM_LIB) && isdir(WF_SPICE_PATH)

if WF_GRAM_READY
    if Base.find_package("GRAMSuite") === nothing
        pushfirst!(LOAD_PATH, WF_GRAM_ROOT)
    end
    @eval import GRAMSuite
end

const WF_ZERO = SVector{3, Float64}(0.0, 0.0, 0.0)
const WF_WIND = SVector{3, Float64}(40.0, -25.0, 3.0)

# A non-GRAM model that, like native GRAM, returns a wind whatever the flag
# says, and counts how often it was asked for one.
struct WindyDensityModel <: SM.AbstractDensityModel
    wind::SVector{3, Float64}
    wind_requests::Base.RefValue{Int}
    queries::Base.RefValue{Int}
end
WindyDensityModel(w) = WindyDensityModel(w, Ref(0), Ref(0))
function EM.getDensity(m::WindyDensityModel, h::Float64, lat::Float64, lon::Float64,
        t::Float64, wind::Bool)
    m.queries[] += 1
    wind && (m.wind_requests[] += 1)
    return 1.0e-9 * exp(-(h - 180.0e3) / 30.0e3), 700.0, m.wind
end
EM.getDensity(m::WindyDensityModel, h::Float64, lat::Float64, lon::Float64,
    t::Float64, wind::Bool, p) = EM.getDensity(m, h, lat, lon, t, wind)

function wf_config(density_model; wind::Bool, n_sats::Int=2, mission_s::Float64=20.0,
        planet=SM.Earth(), ephemerides_model=SimpleEphemeridesModel())
    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link(root=true, m=500.0, ref_area=12.0)
        ic = InitialCondition(
            ra=planet.Rp_e + 400.0e3, rp=planet.Rp_e + 170.0e3 + 5.0e3 * i,
            i=53.0, ω=0.0, Ω=10.0 + 30.0 * i, ν=0.0,
        )
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0,
            root.inertia, 0, 0, ic, i))
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=MissionConfiguration(MissionTime, true, 1, mission_s, false, 20, 2.0),
        environment_model=EnvironmentModel(
            planet=planet, EI=600.0, density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=wind, ephemerides_model=ephemerides_model,
        ),
        dynamics_model=DynamicsModel(spacecraft,
            (InverseSquaredGravityModel(), AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-10,
            dt_max_orbit=1.0),
    )
end

const WF_RUN_ENV = (
    "SPACEAGORA_RHS_CALIBRATE" => "off",
    "SPACEAGORA_RHS_IDENTIFY" => "0",
    "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
    "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
    "SPACEAGORA_GRAM_ISOLATED_POOL" => "off",
    "SPACEAGORA_GRAM_PROCESS_POOL" => "off",
)

# Solve and return the recorded trajectory and the saved winds and drag.
function wf_run(args)
    recorder = TrajectoryRecorder(args; capacity=1)
    result = withenv(WF_RUN_ENV...) do
        run_simulation(args; isolate_state=false, return_solver_metadata=true,
            visualization=false, extra_callbacks=(get_trajectory_recorder_callback(recorder),))
    end
    @test result.retcode == "Success"
    rows = trajectory_save_data(recorder)
    return (
        times=collect(trajectory_times(recorder)),
        pos=copy(trajectory_positions(recorder)),
        vel=copy(trajectory_velocities(recorder)),
        wind=[collect(row[:wind]) for row in rows],
        drag=[collect(row[:drag]) for row in rows],
    )
end

wf_flat(x::Real) = Float64[x]
wf_flat(x) = reduce(vcat, map(wf_flat, collect(x)); init=Float64[])
wf_bits(a) = reinterpret(UInt64, wf_flat(a))
wf_all_zero_wind(run) = all(ws -> all(w -> w == WF_ZERO, ws), run.wind)

# One density-callback firing and one staged environment sample on a fresh
# ODEParams, the two places aerodynamics reads its wind from.
function wf_stage_samples(args; batch::String="off")
    n = length(args.dynamics_model.spacecraft)
    p = SM.ODEParams(n_sats=n, args=args)
    p.shared_buffers.et_start[] = SM.ephemerides_time_seconds(args.initial_time,
        args.environment_model.ephemerides_model)
    withenv("SPACEAGORA_DENSITY_BATCH_PARALLEL" => batch, WF_RUN_ENV...) do
        p.shared_buffers.callback_env_config[] = CB._snapshot_callback_env_config()
        u = SE.build_initial_conditions(args)
        cb = CB.get_density_callback(n, args.dynamics_model.dynamic_effectors, args)
        cb.affect!((p=p, u=u, t=0.0))
        buffered = copy(p.shared_buffers.winds[1:n])
        staged = [CB._stage_environment_state(u.sc[i], p, i, 0.5; write_buffers=false).wind
                  for i in 1:n]
        sampled = [SE.sample_atmosphere(u.sc[i], p, i, 0.75; write_buffers=false).wind_pp
                   for i in 1:n]
        return (buffered=buffered, staged=staged, sampled=sampled, p=p)
    end
end

@testset "environment wind helpers" begin
    on = (args=(environment_model=(wind=true,),),)
    off = (args=(environment_model=(wind=false,),),)
    legacy = (args=(environment_model=(EI=100.0,),),)
    @test EM._environment_wind_enabled(on)
    @test !EM._environment_wind_enabled(off)
    @test EM._environment_wind_enabled(legacy)
    @test EM._environment_wind(on, WF_WIND) == WF_WIND
    @test EM._environment_wind(off, WF_WIND) == WF_ZERO
    @test EM._environment_wind(false, (1.0, 2.0, 3.0)) == WF_ZERO
    buf = fill(WF_WIND, 3)
    EM._zero_environment_winds!(on, buf)
    @test buf == fill(WF_WIND, 3)
    EM._zero_environment_winds!(off, buf)
    @test buf == fill(WF_ZERO, 3)
end

@testset "wind=false zeroes non-GRAM winds in aerodynamics and output" begin
    for batch in ("off", "on")
        windy = WindyDensityModel(WF_WIND)
        s = wf_stage_samples(wf_config(windy; wind=false); batch=batch)
        @test all(==(WF_ZERO), s.buffered)
        @test all(==(WF_ZERO), s.staged)
        @test all(==(WF_ZERO), s.sampled)
        @test windy.queries[] > 0
        @test windy.wind_requests[] == 0

        # Enabled winds are passed through unchanged.
        windy_on = WindyDensityModel(WF_WIND)
        s_on = wf_stage_samples(wf_config(windy_on; wind=true); batch=batch)
        @test all(==(WF_WIND), s_on.buffered)
        @test all(==(WF_WIND), s_on.staged)
        @test all(==(WF_WIND), s_on.sampled)
        @test windy_on.wind_requests[] == windy_on.queries[] > 0
    end

    # Uniform-light batch prefill (flat route), called directly.
    windy = WindyDensityModel(WF_WIND)
    args = wf_config(windy; wind=false, n_sats=3)
    p = SM.ODEParams(n_sats=3, args=args)
    fill!(p.shared_buffers.winds, WF_WIND)
    SE._fill_uniform_light_atmosphere!(p, 0.0, 3, windy,
        [180.0e3, 190.0e3, 200.0e3], zeros(3), zeros(3))
    @test all(==(WF_ZERO), p.shared_buffers.winds[1:3])
    @test windy.wind_requests[] == 0

    # Per-link aerodynamic atmosphere query.
    rho, T, w = SM.DynamicEffectors.AerodynamicEffectors._aero_link_atmosphere_query(
        p, 1, 0.0, SVector{3, Float64}(args.environment_model.planet.Rp_e + 180.0e3, 0.0, 0.0),
        args.environment_model.planet)
    @test w == WF_ZERO
    @test rho > 0.0
end

@testset "wind=false equals a run with nominally zero winds" begin
    # Wind disabled on a model that returns winds is the same solve, bit for
    # bit, as winds enabled on the same model with zero winds.
    disabled = wf_run(wf_config(WindyDensityModel(WF_WIND); wind=false))
    calm = wf_run(wf_config(WindyDensityModel(WF_ZERO); wind=true))
    windy = wf_run(wf_config(WindyDensityModel(WF_WIND); wind=true))
    @test length(disabled.times) > 2
    @test disabled.times == calm.times
    @test wf_bits(disabled.pos) == wf_bits(calm.pos)
    @test wf_bits(disabled.vel) == wf_bits(calm.vel)
    @test wf_bits(disabled.drag) == wf_bits(calm.drag)
    @test wf_all_zero_wind(disabled)
    # The flag is what removed the wind: enabled, the same model changes both
    # the saved wind and the drag.
    @test all(ws -> all(==(WF_WIND), ws), windy.wind)
    @test wf_bits(windy.drag) != wf_bits(disabled.drag)
    @test any(!iszero, wf_flat(disabled.drag))

    # Analytic models already return zero wind; the flag must not change them.
    planet = SM.Earth()
    expo = ExponentialAtmosphereModel(1.0e-9, 180.0e3, 30.0e3; temperature_k=700.0)
    e_off = wf_run(wf_config(expo; wind=false, planet=planet))
    e_on = wf_run(wf_config(expo; wind=true, planet=planet))
    @test wf_bits(e_off.pos) == wf_bits(e_on.pos)
    @test wf_bits(e_off.drag) == wf_bits(e_on.drag)
    @test wf_all_zero_wind(e_off)
    @test wf_all_zero_wind(e_on)
end

# A raw GRAM core without native GRAM: exercises both pools' routing guards.
mutable struct WindHistoryCore
    history_dependent::Bool
    calls::Int
end
EM._gram_core_wind_is_history_dependent(core::WindHistoryCore) = core.history_dependent
function EM._gram_core_density_state(core::WindHistoryCore, h::Float64, lat::Float64,
        lon::Float64, t::Float64, wind::Bool, lk, temperature::Float64)
    lock(lk) do
        core.calls += 1
        return 1.0e-9, 700.0, WF_WIND
    end
end

@testset "pool wind-history rule" begin
    for history in (false, true), wind in (false, true)
        model = EM.GRAMAtmosphereModel(WindHistoryCore(history, 0))
        @test !CB._gram_pool_declines_history_dependent_winds(:off, wind, model)
        @test !CB._gram_pool_declines_history_dependent_winds(:on, wind, model)
        @test CB._gram_pool_declines_history_dependent_winds(:auto, wind, model) ==
              (history && wind)
    end
end

@testset "process pool applies the wind-history rule per mode" begin
    function candidate(mode; wind::Bool, history::Bool)
        model = EM.GRAMAtmosphereModel(WindHistoryCore(history, 0), ReentrantLock(),
            Dict{Symbol, Any}(:planet_name => "earth"))
        args = wf_config(model; wind=wind, n_sats=4)
        p = SM.ODEParams(n_sats=4, args=args)
        withenv("SPACEAGORA_GRAM_PROCESS_POOL" => mode,
                "SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => "1",
                "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
                "SPACEAGORA_VACUUM_GRAM_CACHE" => nothing,
                "SPACEAGORA_GRAM_TRACK_CACHE" => "off") do
            p.shared_buffers.callback_env_config[] = CB._snapshot_callback_env_config()
            return CB._rhs_density_service_candidate(p, 4)
        end
    end
    for wind in (false, true), history in (false, true)
        @test !candidate("off"; wind=wind, history=history)
        @test candidate("on"; wind=wind, history=history)
        @test candidate("auto"; wind=wind, history=history) == !(wind && history)
    end

    # The batch evaluator declines under auto before it asks for workers or
    # writes an output.
    withenv("SPACEAGORA_GRAM_PROCESS_POOL" => "auto",
            "SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => "1",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        model = EM.GRAMAtmosphereModel(WindHistoryCore(true, 0), ReentrantLock(),
            Dict{Symbol, Any}(:planet_name => "earth"))
        p = SM.ODEParams(n_sats=4, args=wf_config(model; wind=true, n_sats=4))
        n = 4
        rhos, Ts, ws = fill(-1.0, n), fill(-2.0, n), fill(WF_WIND, n)
        @test !CB._gram_process_pool_batch_eval!(rhos, Ts, ws, model,
            fill(150.0e3, n), zeros(n), zeros(n), 0.0, true, p)
        @test rhos == fill(-1.0, n)
        @test Ts == fill(-2.0, n)
        @test ws == fill(WF_WIND, n)
        @test model.core.calls == 0
        @test !CB.ParallelProcess.density_service_failed(model.constructor_kwargs)
    end
end

# The two fixes compose: a history-dependent GRAM model in a run with
# `wind = false` is not history-dependent for that run, so the run-scoped
# snapshot (the one `setup.jl` installs) leaves the explicit `auto` freeze off,
# and both pools stay eligible under `auto`. Winds on engage both pool guards;
# freezing additionally requires an explicit sampling opt-in.
@testset "wind=false lifts the freeze and the pool guard for a history-dependent model" begin
    n = 4
    for wind in (false, true), freeze in (nothing, "auto")
        model = EM.GRAMAtmosphereModel(WindHistoryCore(true, 0), ReentrantLock(),
            Dict{Symbol, Any}(:planet_name => "earth"))
        args = wf_config(model; wind=wind, n_sats=n)
        p = SM.ODEParams(n_sats=n, args=args)
        withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => freeze,
                "SPACEAGORA_GRAM_PROCESS_POOL" => "auto",
                "SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD" => "1",
                "SPACEAGORA_GRAM_ISOLATED_POOL" => "auto",
                "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
                "SPACEAGORA_VACUUM_GRAM_CACHE" => nothing,
                "SPACEAGORA_GRAM_TRACK_CACHE" => "off") do
            cfg = CB._snapshot_callback_env_config(args)
            p.shared_buffers.callback_env_config[] = cfg
            @test EM.density_model_history_dependent(model)
            @test cfg.density_history_dependent == wind
            @test cfg.density_freeze_per_step == (wind && freeze == "auto")
            @test CB._rhs_density_service_candidate(p, n) == !wind

            # The isolated pool's callers pass the run's wind flag to the
            # shared guard (`gram_isolated_pool_tests.jl` covers the batch
            # evaluator itself; its persistent workers cannot see a core type
            # defined after they started, so it is not driven from here).
            wind_requested = EM._environment_wind_enabled(p)
            @test wind_requested == wind
            @test cfg.gram_isolated_pool_mode === :auto
            @test CB._gram_pool_declines_history_dependent_winds(
                cfg.gram_isolated_pool_mode, wind_requested, model) == wind
        end
    end
end

if !WF_GRAM_READY
    if NativeProbeReporting.native_probe_required()
        error("Native GRAM wind-flag tests are required but no libGRAM was found at $(WF_GRAM_LIB)")
    end
    @info "Skipping native GRAM wind-flag tests: no libGRAM on this host." lib = WF_GRAM_LIB
else
    @testset "wind=false zeroes native GRAM winds" begin
        planet = SM.Earth("", WF_SPICE_PATH)
        gram = EM.GRAMAtmosphereModel(planet_name="earth")
        function gram_args(wind)
            return wf_config(deepcopy(gram); wind=wind, planet=planet, mission_s=10.0)
        end
        # Nominal winds are deterministic, so a run with winds enabled is the
        # reference that shows GRAM does return a wind here.
        withenv("SPACEAGORA_GRAM_WIND_MODE" => "nominal") do
            s_on = wf_stage_samples(gram_args(true))
            @test any(w -> w != WF_ZERO, s_on.buffered)
            @test any(w -> w != WF_ZERO, s_on.staged)
            @test any(w -> w != WF_ZERO, s_on.sampled)
            s_off = wf_stage_samples(gram_args(false))
            @test all(==(WF_ZERO), s_off.buffered)
            @test all(==(WF_ZERO), s_off.staged)
            @test all(==(WF_ZERO), s_off.sampled)
            # Density does not depend on the wind request.
            @test s_off.p.shared_buffers.densities == s_on.p.shared_buffers.densities
        end
        # The default (perturbed) mode is masked the same way.
        for mode in (nothing, "nominal")
            withenv("SPACEAGORA_GRAM_WIND_MODE" => mode) do
                off = wf_run(gram_args(false))
                @test length(off.times) > 2
                @test wf_all_zero_wind(off)
                @test any(!iszero, wf_flat(off.drag))
            end
        end
        withenv("SPACEAGORA_GRAM_WIND_MODE" => "nominal") do
            on = wf_run(gram_args(true))
            @test !wf_all_zero_wind(on)
        end
    end
end

end # module
