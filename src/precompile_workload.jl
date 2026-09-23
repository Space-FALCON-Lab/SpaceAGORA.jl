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
    parse_parallel_profile("R2")
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
            route_state=OuterRouteState(),
            route_tuning=OuterRouteTuning(process_max_workers=1),
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
            route_state=OuterRouteState(),
            route_tuning=OuterRouteTuning(process_max_workers=1),
        )
    end
    return nothing
end

@setup_workload begin
    @compile_workload begin
        _run_spaceagora_precompile_workload()
        _warm_predictive_campaign()
        _warm_mixed_dispatch_campaign()
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
    # Same hazard, for the two campaign warmups just above: reset the R7
    # machine-constants cache and the route-state "already loaded" flag they
    # populate, even though route persistence was held off above, as cheap
    # insurance against the same class of bug.
    SimulationCampaigns.reset_predictive_machine_constants!()
    SimulationCampaigns.reset_campaign_route_state_persistence!()
end
