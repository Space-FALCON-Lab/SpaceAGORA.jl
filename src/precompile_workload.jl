function _spaceagora_precompile_args()
    SM = SimulationModel
    planet = SM.Mars()
    spacecraft = TelemetryVerification.make_three_body_spacecraft(
        bus_dims=(1.2, 1.1, 0.9),
        panel_dims=(0.01, 0.8, 0.4),
        bus_mass=150.0,
        panel_mass_each=2.0,
        panel_offset_y=0.7,
        ic=SM.InitialCondition(
            ra=planet.Rp_e + 220e3,
            rp=planet.Rp_e + 150e3,
            i=28.0,
            ω=10.0,
            Ω=15.0,
            ν=165.0
        ),
        prop_mass=15.0,
        id=1
    )

    return TelemetryVerification.make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=5.0,
        initial_time=SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredGravityModel(),),
        density_model=SM.ExponentialAtmosphereModel(planet),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=true,
        EI_km=140.0,
        verbose=false,
        results=false,
        results_directory=joinpath(tempdir(), "spaceagora_precompile")
    )
end

const _SPACEAGORA_PRECOMPILE_ENV = Dict("SPACEAGORA_PARALLEL_PROFILE" => "R2", "SPACEAGORA_SAVE_BUNDLE" => "0", "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0")

function _run_spaceagora_precompile_workload(; workspace::AbstractString=tempdir())
    ParallelProfiles.parse_parallel_profile("R2")
    engine_config = simulation_engine_config_from_env(_SPACEAGORA_PRECOMPILE_ENV)
    args = _spaceagora_precompile_args()
    mktempdir(workspace) do tmp
        cd(tmp) do
            run_simulation(engine_config, args; return_solution=true)
        end
    end
end

# `_warm_campaign_dispatchers` (SimulationCampaigns, monte_carlo.jl) exercises
# the low-level Monte Carlo dispatchers directly -- the serial loop, the mixed
# dispatcher, the threaded dispatcher -- but never the campaign-level PLANNING
# layer above them (`_campaign_route_plan`/`_run_campaign_with_route_env` for
# the default bandit route, `predictive_plan`/`_run_campaign_predictive` for
# R7), because it calls those dispatchers straight, not through
# `run_monte_carlo`. Both warmups below go through that public entry instead,
# with two trivial samples so a route plan resolves without ever reaching a
# process pool (`_mc_process_worth_exploring` requires >= 16 samples or a
# 3600 s mission; two samples with the default zero mission time clears
# neither, on any tuning). Route persistence is forced off for the same
# reason `_reset_furnished_kernels!` exists below: `ensure_campaign_route_state_loaded!`
# and `predictive_machine_constants` both cache into module-level `const`s the
# first time a campaign reads them, and an untouched cache would otherwise be
# serialised into the pkgimage -- an "already loaded" flag with nothing
# actually loaded, or this precompiling machine's calibration file (or absence
# of one) served to every process that later loads the same pkgimage. Each
# warmup uses its own throwaway `OuterRouteState()` rather than the process
# default, so a trivial `seed -> seed * 2` timing never pollutes the route
# bandit's history for a real campaign that runs later in the same process.
const _SPACEAGORA_PRECOMPILE_ROUTE_ENV = Dict(
    "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
    "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
)

function _warm_predictive_campaign()::Nothing
    sample = seed -> seed * 2
    seeds = collect(1:2)
    withenv(_SPACEAGORA_PRECOMPILE_ROUTE_ENV..., "SPACEAGORA_CAMPAIGN_PLANNER" => "predictive") do
        run_monte_carlo(
            sample, seeds;
            threads=:auto,
            route_state=ParallelProfiles.OuterRouteState(),
            route_tuning=ParallelProfiles.OuterRouteTuning(process_max_workers=1),
        )
    end
    return nothing
end

function _warm_mixed_dispatch_campaign()::Nothing
    sample = seed -> seed * 2
    seeds = collect(1:2)
    withenv(_SPACEAGORA_PRECOMPILE_ROUTE_ENV...) do
        run_monte_carlo(
            sample, seeds;
            threads=:auto,
            route_state=ParallelProfiles.OuterRouteState(),
            route_tuning=ParallelProfiles.OuterRouteTuning(process_max_workers=1),
        )
    end
    return nothing
end

"""
    _reset_process_local_state!()

Clear every module-level cache that records a fact about the process or host
it was filled in: the machine topology (cores, affinity, cgroup quota and
memory limit, and the `SPACEAGORA_CORE_BUDGET` / `SPACEAGORA_PHYSICAL_CORES`
overrides read when it was filled), the RHS calibration machine label and its
loaded store, the measured hint-consultation overhead, the planner's machine
constants, corrections and warm-pool record, the "already loaded" and
"exit hook registered" flags of the persisted route and hint state, and the
native-lock counters with the time their window started.

A `const` cache filled during precompilation is serialized into the pkgimage
and served, unchanged, to every process that loads it: a container image built
on one host then plans for that host's cores and memory, the overrides stop
taking effect at run time, and an "exit hook registered" flag is set in a
process that never registered one. So this runs at the end of the precompile
workload, where nothing it cleared can reach the image, and again from
`__init__`, so a cache filled by any future workload that forgets to clear it
is still refilled from the process that loads the image. Everything here is
recomputed or reloaded on first use.
"""
function _reset_process_local_state!()::Nothing
    lock(ParallelProfiles._TOPOLOGY_LOCK) do
        ParallelProfiles._TOPOLOGY_CACHE[] = nothing
    end

    SE = SimulationEngine
    lock(SE._rhs_calib_lock) do
        SE._CALIB_MACHINE_LABEL[] = ""
        empty!(SE._rhs_calib_cache)
        SE._rhs_calib_loaded[] = false
        SE._rhs_calib_loaded_path[] = ""
        empty!(SE._rhs_calib_solve_start)
        empty!(SE._rhs_calib_solve_honoured)
    end
    SE._PARALLEL_FLAG_ONE_THREAD_NOTED[] = false

    PPol = SimulationModel.ParallelPolicy
    PPol._HINT_OVERHEAD_NS[] = -1.0
    lock(PPol._persistent_hint_lock) do
        PPol._persistent_hint_state[] = PPol._PersistentHintState()
        PPol._persistent_hint_atexit_registered[] = false
    end
    PPol._global_policy_context[] = PPol.PolicyContext()

    PC = SimulationModel.ParallelCost
    lock(PC._ENSURE_CONSTANTS_LOCK) do
        empty!(PC._ENSURED_CONSTANTS_PATHS)
    end

    SC = SimulationCampaigns
    SC.reset_predictive_machine_constants!()
    lock(SC._CAMPAIGN_ROUTE_STATE_LOCK) do
        SC._CAMPAIGN_ROUTE_STATE_LOADED[] = false
        SC._CAMPAIGN_ROUTE_STATE_ATEXIT[] = false
    end
    ParallelProfiles.reset_outer_route_state!(SC._CAMPAIGN_OUTER_ROUTE_STATE)
    lock(SC._CAMPAIGN_CORRECTIONS_LOCK) do
        SC._CAMPAIGN_CORRECTIONS[] = nothing
        SC._CAMPAIGN_CORRECTIONS_PATH[] = ""
    end
    lock(SC._PREDICTIVE_WARM_LOCK) do
        empty!(SC._PREDICTIVE_WARM_WORKERS)
    end
    SC._GC_DEBT[] = false
    SC._POOL_PROBE_DONE[] = false

    RuntimeServices.reset_native_lock_stats!()
    return nothing
end

@setup_workload begin
    @compile_workload begin
        _run_spaceagora_precompile_workload()
        _warm_predictive_campaign()
        _warm_mixed_dispatch_campaign()
        # The Monte Carlo dispatchers compile on their first campaign in a
        # process: the job channel, the feeders and local consumers of the mixed
        # dispatcher, the sample wrapper, the steady-cost estimator. Measured on
        # the paper harness (L12, independent_1sat_1hr, 64 samples), the first
        # pool campaign cost 3.1-3.2 s against 1.8-2.2 s for the static pool
        # path's own cold start on both machines, and 0.2-0.6 s warm. Run on a
        # trivial sample; the user's sample closure still specializes on first call.
        SimulationCampaigns._warm_campaign_dispatchers()
    end
    # MANDATORY whenever a workload above touches a SPICE-backed planet, and
    # cheap insurance when none does. `_FURNISHED_KERNELS` and the planet
    # instance caches are module-level `const`s, so anything furnished HERE is
    # serialised into the pkgimage; at run time `_furnsh_once` would then skip
    # the furnish for a kernel CSPICE never actually loaded, and every lookup
    # needing it fails against an empty pool (`utc2et` losing the leapseconds
    # kernel is the first symptom). Same hazard `_reset_furnished_kernels!`
    # documents for `kclear()`, reached by a different route. Verified by
    # observation, not theory: a workload constructing `Earth(...)` here left
    # every subsequent process unable to resolve a UTC epoch until the pkgimage
    # was rebuilt.
    SimulationModel.Planets._reset_furnished_kernels!()
    # Same hazard, for every cache the workloads above fill with a fact about
    # this precompiling process or host -- the campaign warmups reach
    # `machine_topology()` through the route planners, for one. See
    # `_reset_process_local_state!`.
    _reset_process_local_state!()
end
