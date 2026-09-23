# Shared configuration for the WS11b GRAM thread-scaling study.
#
# `run_scaling.jl` (timing) and `dump_states.jl` (bit-identity) must build the
# same configurations from the same env, or the dumps would prove identity of
# something the timings never ran. Everything both need lives here; neither file
# defines a case or an env pair of its own.
#
# Include AFTER parallelization_performance/{cli,modes,cases}.jl: the builders
# below call `ppc_build_config`, `ppc_spacecraft` and `ppc_gram_atmosphere_model`.

# Entry interface for the pool-engaged configuration, in km. DERIVED: the
# constellation below tops out at 480 km apoapsis altitude, and `in_atmosphere`
# -- which gates the vacuum-predicted look-ahead cache -- is set from
# `altitude <= EI`. 600 km is the round number above that ceiling.
const GTS_ENGAGED_EI_KM = 600.0

"""
    gts_density_env(path, mission_s)

The env pairs for one density path. Values copied from
`_ppc_p6_gram_density_env!` in parallelization_performance/cases.jl, which is
where the S2 scenario's settings live; keeping them identical is what makes this
study's rows comparable with the P6 and S2 traces rather than a separate
configuration that happens to have the same name.
"""
function gts_density_env(path::String, mission_s::Float64)
    if path == "lookahead"
        return [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
            "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(mission_s + 500.0),
            "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
        ]
    elseif path == "freeze"
        return [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
        ]
    end
    error("Unknown density path '$path'. Use lookahead or freeze.")
end

function gts_pool_env(workers::Int)
    workers <= 0 && return ["SPACEAGORA_GRAM_ISOLATED_POOL" => "off"]
    return [
        "SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
        "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => string(workers),
        "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "1",
    ]
end

"""
    GTS_WIDTH_ENV

Both arms are run with the density callback's automatic minimum thread budget
lowered from its default of 16 to 1.

Without it there is nothing to measure below 16 threads. The pool's width is
`_density_callback_thread_decision(...; heavy_work=true).allotment`
(density_callbacks/runtime.jl), and `auto_thread_min_budget(:density_callback)`
(parallel/policy/env_config.jl) pins that decision to a single thread on any
process with fewer than 16 threads -- so `_gram_isolated_pool_batch_eval!`
computes `workers = 1`, fails its own `workers > 1` guard and returns false,
and the run silently takes the locked path however `SPACEAGORA_GRAM_ISOLATED_POOL`
is set. Verified directly on this workstation at 4 threads and 64 spacecraft:
the decision returns `(use_threads = false, allotment = 1)`, and the locked and
pooled arms record the same `gram_density` acquisition count to the unit.

That floor exists, by its own comment, because native GRAM is serialized behind
a process-wide lock and oversubscribing it wastes cycles fighting for that lock.
Which makes it circular for this measurement: it is the lock that motivates the
gate, and the pool is the thing that removes the lock. Setting it to 1 in BOTH
arms keeps the comparison single-variable -- it is the same number in the locked
run and the pooled run, and in the locked run it can only affect the kinematics
pre-fill, never the GRAM evaluation, which `getDensityBatch!` performs serially
whatever the width.
"""
const GTS_WIDTH_ENV = ["SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "1"]

"""
    gts_constellation(planet, n)

A low-Earth constellation whose every member samples native GRAM for the whole
mission. It is NOT `ppc_constellation`, and the difference is the point:

  * `ppc_constellation` places member i at 540 + 2(i-1) km apoapsis, so at
    n = 1024 its highest members sit above 2000 km, where
    `getDensity(::GRAMAtmosphereModel, ...)` returns zero without calling GRAM
    at all. Measuring a GRAM concurrency change on a constellation that skips
    GRAM for part of its membership would understate the effect by an amount
    that depends on n, which is the one axis being swept.
  * every member of `ppc_constellation` also starts far above the 120 km entry
    interface `ppc_build_config` sets, so `in_atmosphere` is false for all of
    them and the vacuum-predicted look-ahead cache -- one of the two density
    paths this study is about -- never activates.

So: altitudes in a 300-480 km band (ASSUMED; chosen as a drag-relevant LEO band
that is entirely below the 2000 km GRAM cut-off and entirely below the entry
interface set by `--ei-km`), spread in true anomaly and right ascension so the
members do not share a ground track. Everything else -- the L50 harmonics, the
`AerodynamicCoefficientfM` aero effector, the 5 s step cap, the tolerances --
is `ppc_build_config`'s and matches the P6 aero traces.
"""
function gts_constellation(planet, n::Int)
    sats = SpacecraftModel[]
    for i in 1:n
        k = mod(i - 1, 5)
        push!(sats, ppc_spacecraft(
            planet;
            id=i,
            ra_alt_m=400e3 + 20e3 * k,
            rp_alt_m=300e3 + 10e3 * k,
            raan_deg=360.0 * (i - 1) / n,
            nu_deg=120.0 + 240.0 * (i - 1) / max(1, n),
        ))
    end
    return sats
end

# L50 harmonics plus aero, live native GRAM, dt_max_orbit 5 s: the P6 aero
# trace's force model on the constellation above.
function gts_build_config(n::Int, mission_s::Float64, ei_km::Float64)
    planet = Earth("", PPC_SPICE_PATH)
    args = ppc_build_config(
        planet=planet,
        spacecraft=gts_constellation(planet, n),
        mission_time_s=mission_s,
        orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
        density_model=ppc_gram_atmosphere_model("earth"),
        dt_max_orbit=5.0
    )
    # The entry interface decides `in_atmosphere`, which is what gates the
    # vacuum-predicted look-ahead cache. Put it above the constellation so every
    # member is inside it from the first step.
    return SimulationConfiguration(
        simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration,
        environment_model=EnvironmentModel(
            planet=args.environment_model.planet,
            EI=ei_km,
            density_model=args.environment_model.density_model,
            thermal_model=args.environment_model.thermal_model,
            topography=false,
            wind=false
        ),
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

