using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel

# SolverConfig(parallel=true): the single user-facing switch for parallel
# execution. It must (1) change nothing when off, (2) run a single simulation
# and a campaign under the R7 profile scoped to the call, (3) calibrate the
# machine once when its cost constants are missing, and (4) leave nested runs
# inside an enclosing outer split alone.

const PF_SM = SpaceAGORA.SimulationModel
const PF_SE = SpaceAGORA.SimulationEngine
const PF_SC = SpaceAGORA.SimulationCampaigns
const PF_PR = SpaceAGORA.ParallelProfiles
const PF_PC = SpaceAGORA.SimulationModel.ParallelCost
const PF_DiscreteCallback = SpaceAGORA.SimulationEngine.DiffEqBase.DiscreteCallback

# Everything a parallel profile writes, plus the two outer-split variables.
const PF_WATCHED_KEYS = sort!(unique!(vcat(
    first.(PF_PR.profile_env_pairs(PF_PR.R7; preserve_existing = false)),
    ["SPACEAGORA_INNER_THREAD_BUDGET", "SPACEAGORA_COST_CONSTANTS_PATH"],
)))

pf_env_snapshot() = Dict(k => get(ENV, k, nothing) for k in PF_WATCHED_KEYS)

function pf_config(; solver_config = nothing, n_sats::Int = 1, mission_s::Float64 = 30.0)
    planet = PF_SM.Earth()
    spacecraft = map(1:n_sats) do k
        root = Link(root = true, m = 500.0, ref_area = 12.0)
        ic = InitialCondition(ra = planet.Rp_e + 550_000.0 + 1_000.0 * k,
                              rp = planet.Rp_e + 550_000.0 + 1_000.0 * k,
                              i = 53.0, ω = 0.0, Ω = 10.0 * k, ν = 0.0)
        SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, k)
    end
    return SimulationConfiguration(
        simulation_settings = SimulationSettings(
            results = false, verbose = false, generate_plots = false, normalize = false,
            save_csv = false, results_directory = mktempdir(),
        ),
        mission_configuration = MissionConfiguration(
            mission_type = MissionTime, keplerian = true, number_of_orbits = 1,
            mission_time = mission_s, orientation_sim = false, num_steps_to_save = 20,
            data_rate = 5.0,
        ),
        environment_model = EnvironmentModel(
            planet = planet, EI = 300.0,
            density_model = NoAtmosphereModel(),
            thermal_model = MaxwellianHeat(thermal_accomodation_factor = 1.0, planet = planet),
            topography = false, wind = false,
            ephemerides_model = SimpleEphemeridesModel(),
        ),
        dynamics_model = DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model = GuidanceModel(guidance_effectors = (), guidance_rates = Float64[]),
        navigation_model = NavigationModel(navigation_effectors = (), navigation_rates = Float64[]),
        control_model = ControlModel(control_effectors = (), control_rates = Float64[]),
        initial_time = InitialTime(year = 2020, month = 1, day = 1, hour = 0, minute = 0, second = 0.0),
        integration_tolerances = IntegrationTolerances(reltol_orbit = 1e-9, abstol_orbit = 1e-9, dt_max_orbit = 2.0),
        solver_config = solver_config,
    )
end

# A callback that records the environment the solve actually ran under, on its
# first step. `sink` collects one snapshot per solve (a campaign has several).
function pf_env_probe(sink::Vector, lock_ = ReentrantLock())
    seen = Ref(false)
    condition = (u, t, integrator) -> !seen[]
    affect! = integrator -> begin
        seen[] = true
        lock(lock_) do
            push!(sink, pf_env_snapshot())
        end
        return nothing
    end
    return PF_DiscreteCallback(condition, affect!; save_positions = (false, false))
end

pf_throwing_probe() = PF_DiscreteCallback((u, t, i) -> true,
                                          i -> error("parallel flag probe failure");
                                          save_positions = (false, false))

# Every test below runs against a scratch constants path, so the real store
# under output/parallel_policy_state is never read or written.
const PF_TMP = mktempdir()
const PF_CONSTANTS_PATH = joinpath(PF_TMP, "cost_constants_test.toml")
# A cheap calibration: the real kernels with the fewest repeats. The tests
# assert when it runs, not what it measures.
const PF_CALIBRATIONS = Ref(0)
pf_calibrate() = (PF_CALIBRATIONS[] += 1; PF_PC.calibrate_machine(k = 1))

# The flag turns on persisted policy state (hints, RHS calibration, route
# state, planner corrections); point all of it at the scratch directory too.
const PF_ISOLATION = [
    "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => joinpath(PF_TMP, "policy_state.toml"),
    "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(PF_TMP, "rhs_calibration.toml"),
    "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => joinpath(PF_TMP, "outer_route_state.toml"),
    "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => joinpath(PF_TMP, "campaign_corrections.toml"),
]

withenv(PF_ISOLATION...) do

@testset "SolverConfig carries the flag, off by default" begin
    @test PF_SM.SolverConfig().parallel === false
    @test PF_SM.SolverConfig(parallel = true).parallel === true
    @test SpaceAGORA.SolverConfig === PF_SM.SolverConfig
    # The flag's profile is R7, applied with preserve_existing=false.
    @test PF_PR.PARALLEL_FLAG_PROFILE === PF_PR.R7
    withenv("SPACEAGORA_CAMPAIGN_PLANNER" => "bandit", "SPACEAGORA_PARALLEL_POLICY_V2" => "0") do
        @test Dict(PF_PR.parallel_flag_env_pairs()) ==
              Dict(PF_PR.profile_env_pairs(PF_PR.R7; preserve_existing = false))
    end
end

@testset "machine constants are calibrated once, and never when present" begin
    path_a = joinpath(PF_TMP, "a", "constants.toml")
    PF_CALIBRATIONS[] = 0
    @test PF_PC.ensure_machine_constants!(path = path_a, calibrate = pf_calibrate) === :calibrated
    @test PF_CALIBRATIONS[] == 1
    @test isfile(path_a)
    @test PF_PC.load_machine_constants(path_a) !== nothing
    # Settled for this process: no second measurement, no second file read.
    @test PF_PC.ensure_machine_constants!(path = path_a, calibrate = pf_calibrate) === :checked
    @test PF_CALIBRATIONS[] == 1

    # A current file written elsewhere is used as is.
    path_b = joinpath(PF_TMP, "b", "constants.toml")
    mkpath(dirname(path_b))
    cp(path_a, path_b)
    @test PF_PC.ensure_machine_constants!(path = path_b, calibrate = pf_calibrate) === :present
    @test PF_CALIBRATIONS[] == 1

    # A file of another schema is not current: it is recalibrated.
    path_c = joinpath(PF_TMP, "c", "constants.toml")
    mkpath(dirname(path_c))
    write(path_c, "schema_version = 1\n")
    @test PF_PC.ensure_machine_constants!(path = path_c, calibrate = pf_calibrate) === :calibrated
    @test PF_CALIBRATIONS[] == 2

    # A failing calibration does not fail the caller and is not retried.
    path_d = joinpath(PF_TMP, "d", "constants.toml")
    failing = () -> error("calibration unavailable")
    @test (@test_logs (:warn,) match_mode = :any PF_PC.ensure_machine_constants!(path = path_d, calibrate = failing)) === :failed
    @test !isfile(path_d)
    @test PF_PC.ensure_machine_constants!(path = path_d, calibrate = failing) === :checked

    # The planner's cached constants follow a newly written file.
    generation = PF_PC.machine_constants_generation()
    PF_PC.save_machine_constants(PF_PC.load_machine_constants(path_a), joinpath(PF_TMP, "e.toml"))
    @test PF_PC.machine_constants_generation() == generation + 1
end

@testset "parallel=false is the default path: same environment, bit-identical trajectory" begin
    withenv("SPACEAGORA_COST_CONSTANTS_PATH" => joinpath(PF_TMP, "never", "constants.toml")) do
        before = pf_env_snapshot()
        seen_nothing = Any[]
        seen_false = Any[]
        sol_nothing = run_simulation(pf_config(); return_solution = true,
                                     extra_callbacks = (pf_env_probe(seen_nothing),))
        sol_false = run_simulation(pf_config(solver_config = PF_SM.SolverConfig(parallel = false));
                                   return_solution = true,
                                   extra_callbacks = (pf_env_probe(seen_false),))
        @test sol_nothing.t == sol_false.t
        @test sol_nothing.u == sol_false.u
        @test only(seen_nothing) == only(seen_false) == before
        @test pf_env_snapshot() == before
        # parallel=false never calibrates.
        @test !isfile(joinpath(PF_TMP, "never", "constants.toml"))
    end
end

@testset "parallel=true runs one simulation under R7 and restores the environment" begin
    path = joinpath(PF_TMP, "single", "constants.toml")
    withenv("SPACEAGORA_COST_CONSTANTS_PATH" => path,
            # Stale shell values the flag must override for the run and restore after.
            "SPACEAGORA_CAMPAIGN_PLANNER" => "bandit",
            "SPACEAGORA_PARALLEL_POLICY_V2" => "0",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        before = pf_env_snapshot()
        seen = Any[]
        sol_off = run_simulation(pf_config(); return_solution = true)
        sol_on = run_simulation(pf_config(solver_config = PF_SM.SolverConfig(parallel = true));
                                return_solution = true, extra_callbacks = (pf_env_probe(seen),))
        inside = only(seen)
        expected = Dict(PF_PR.profile_env_pairs(PF_PR.R7; preserve_existing = false))
        for (k, v) in expected
            @test inside[k] == v
        end
        @test inside["SPACEAGORA_CAMPAIGN_PLANNER"] == "predictive"
        @test inside["SPACEAGORA_PARALLEL_POLICY_V2"] == "1"
        @test pf_env_snapshot() == before
        # First parallel run on a "machine" without constants: calibrated once.
        @test isfile(path)
        @test PF_PC.ensure_machine_constants!(path = path) === :checked
        # Same physics; the flag changes scheduling only (a one-satellite,
        # vacuum run threads nothing, so the answer is the same bits).
        @test sol_on.u[end] ≈ sol_off.u[end] rtol = 1e-9

        # An exception inside the run still restores everything.
        @test_throws Exception run_simulation(
            pf_config(solver_config = PF_SM.SolverConfig(parallel = true));
            extra_callbacks = (pf_throwing_probe(),))
        @test pf_env_snapshot() == before
        @test_throws ErrorException PF_SE._with_parallel_flag(true) do
            @test ENV["SPACEAGORA_CAMPAIGN_PLANNER"] == "predictive"
            error("inside the flag scope")
        end
        @test pf_env_snapshot() == before
    end
end

@testset "parallel=true through a SimulationEngineConfig" begin
    withenv("SPACEAGORA_COST_CONSTANTS_PATH" => PF_CONSTANTS_PATH,
            "SPACEAGORA_CAMPAIGN_PLANNER" => nothing) do
        PF_PC.ensure_machine_constants!(path = PF_CONSTANTS_PATH, calibrate = pf_calibrate)
        before = pf_env_snapshot()
        engine_on = PF_SE.SimulationEngineConfig(solver = PF_SM.SolverConfig(parallel = true))

        seen = Any[]
        run_simulation(engine_on, pf_config(); extra_callbacks = (pf_env_probe(seen),))
        @test only(seen)["SPACEAGORA_CAMPAIGN_PLANNER"] == "predictive"
        @test only(seen)["SPACEAGORA_PARALLEL_POLICY_V2"] == "1"
        @test pf_env_snapshot() == before

        # An explicit args.solver_config wins over the engine config, as it does
        # for every other solver field.
        seen_off = Any[]
        run_simulation(engine_on, pf_config(solver_config = PF_SM.SolverConfig(parallel = false));
                       extra_callbacks = (pf_env_probe(seen_off),))
        @test only(seen_off)["SPACEAGORA_CAMPAIGN_PLANNER"] === nothing
        seen_on = Any[]
        run_simulation(PF_SE.SimulationEngineConfig(),
                       pf_config(solver_config = PF_SM.SolverConfig(parallel = true));
                       extra_callbacks = (pf_env_probe(seen_on),))
        @test only(seen_on)["SPACEAGORA_CAMPAIGN_PLANNER"] == "predictive"
        @test pf_env_snapshot() == before

        # The flag and a different explicit profile contradict each other.
        conflicting = PF_SE.SimulationEngineConfig(parallel = PF_SE.ParallelConfig(profile = "R2"),
                                                   solver = PF_SM.SolverConfig(parallel = true))
        @test_throws ArgumentError run_simulation(conflicting, pf_config())
        @test pf_env_snapshot() == before
    end
end

@testset "a nested parallel run keeps the enclosing split's environment" begin
    path = joinpath(PF_TMP, "nested", "constants.toml")
    withenv("SPACEAGORA_COST_CONSTANTS_PATH" => path,
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
            "SPACEAGORA_INNER_THREAD_BUDGET" => "2",
            "SPACEAGORA_CAMPAIGN_PLANNER" => nothing,
            "SPACEAGORA_PARALLEL_POLICY_V2" => nothing) do
        before = pf_env_snapshot()
        seen = Any[]
        run_simulation(pf_config(solver_config = PF_SM.SolverConfig(parallel = true));
                       extra_callbacks = (pf_env_probe(seen),))
        @test only(seen) == before
        @test !PF_SE._parallel_flag_applies(true)
        # A nested run never starts a calibration of its own.
        @test !isfile(path)
    end
end

@testset "ParallelConfig(profile=...) expands the whole profile" begin
    cfg_r7 = PF_SE.SimulationEngineConfig(parallel = PF_SE.ParallelConfig(profile = "R7"))
    overrides = PF_SE._engine_env_overrides(cfg_r7)
    for (k, v) in PF_PR.profile_env_pairs(PF_PR.R7; preserve_existing = false)
        @test overrides[k] == v
    end
    @test overrides["SPACEAGORA_PARALLEL_POLICY_ADAPTIVE"] == "1"
    @test overrides["SPACEAGORA_CAMPAIGN_PLANNER"] == "predictive"

    # An explicit non-default field still overrides its profile value.
    cfg_override = PF_SE.SimulationEngineConfig(parallel = PF_SE.ParallelConfig(
        profile = "R7", density_callback_parallel_mode = "off", outer_parallel_active = true))
    overrides = PF_SE._engine_env_overrides(cfg_override)
    @test overrides["SPACEAGORA_DENSITY_CALLBACK_PARALLEL"] == "off"
    @test overrides["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "1"
    @test overrides["SPACEAGORA_CONTROL_CALLBACK_PARALLEL"] == PF_PR.profile_config(PF_PR.R7).control_mode
    @test overrides["SPACEAGORA_PARALLEL_POLICY_V2"] == "1"

    # No profile: exactly the fixed set it always wrote.
    overrides = PF_SE._engine_env_overrides(PF_SE.SimulationEngineConfig())
    @test !haskey(overrides, "SPACEAGORA_PARALLEL_PROFILE")
    @test !haskey(overrides, "SPACEAGORA_CAMPAIGN_PLANNER")
    @test overrides["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] == "0"
    @test overrides["SPACEAGORA_PARALLEL_POLICY_ADAPTIVE"] == "0"
    for key in ("SPACEAGORA_EFFECTOR_PARALLEL", "SPACEAGORA_RHS_BATCH_PARALLEL",
                "SPACEAGORA_DENSITY_CALLBACK_PARALLEL", "SPACEAGORA_CONTROL_CALLBACK_PARALLEL",
                "SPACEAGORA_THERMAL_CALLBACK_PARALLEL")
        @test overrides[key] == "auto"
    end

    # An unknown profile is an error rather than a silently ignored string.
    @test_throws ArgumentError PF_SE._engine_env_overrides(
        PF_SE.SimulationEngineConfig(parallel = PF_SE.ParallelConfig(profile = "R9")))
end

@testset "campaigns take the flag" begin
    no_pool = PF_PR.OuterRouteTuning(process_max_workers = 1)
    withenv("SPACEAGORA_COST_CONSTANTS_PATH" => PF_CONSTANTS_PATH,
            "SPACEAGORA_CAMPAIGN_PLANNER" => nothing,
            "SPACEAGORA_PARALLEL_POLICY_V2" => nothing,
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
            # Keep the campaign's route state out of the persisted store.
            "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0") do
        PF_PC.ensure_machine_constants!(path = PF_CONSTANTS_PATH, calibrate = pf_calibrate)
        before = pf_env_snapshot()

        # run_monte_carlo: the keyword.
        result = run_monte_carlo(1:3; parallel = true, route_tuning = no_pool) do seed
            (seed, PF_SC.campaign_planner_mode(), get(ENV, "SPACEAGORA_PARALLEL_POLICY_V2", ""))
        end
        @test length(result.successful) == 3
        @test all(s -> s.value[2] === :predictive && s.value[3] == "1", result.successful)
        @test pf_env_snapshot() == before
        # Without it, nothing changes: serial, default planner.
        result = run_monte_carlo(1:2) do seed
            PF_SC.campaign_planner_mode()
        end
        @test result.threads == 1
        @test all(s -> s.value === :bandit, result.successful)
        @test_throws ArgumentError run_monte_carlo(identity, 1:2; parallel = true, threads = 2)
        @test_throws ArgumentError run_monte_carlo(identity, 1:2; parallel = true, threads = 1)
        @test run_monte_carlo(identity, 1:2; threads = 1).threads == 1

        # run_constellation_ensemble: inferred from the members' configuration.
        args_on = pf_config(solver_config = PF_SM.SolverConfig(parallel = true), n_sats = 2)
        seen = Any[]
        probe = pf_env_probe(seen)
        result = run_constellation_ensemble(args_on; route_tuning = no_pool,
                                            extra_callbacks = (probe,))
        @test length(result.successful) == 2
        @test !isempty(seen)
        @test all(s -> s["SPACEAGORA_CAMPAIGN_PLANNER"] == "predictive", seen)
        @test pf_env_snapshot() == before
        @test_throws ArgumentError run_constellation_ensemble(args_on; threads = 2)

        args_off = pf_config(n_sats = 2)
        seen_off = Any[]
        result = run_constellation_ensemble(args_off; extra_callbacks = (pf_env_probe(seen_off),))
        @test result.threads == 1
        @test all(s -> s["SPACEAGORA_CAMPAIGN_PLANNER"] === nothing, seen_off)

        # Members that disagree are an error.
        @test PF_SC._campaign_parallel_flag([args_on, args_on]) === true
        @test PF_SC._campaign_parallel_flag([args_off, pf_config(solver_config = PF_SM.SolverConfig())]) === false
        @test_throws ArgumentError PF_SC._campaign_parallel_flag([args_on, args_off])
    end
end

end # withenv(PF_ISOLATION...)
