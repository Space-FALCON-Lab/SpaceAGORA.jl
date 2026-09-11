# Single measured subprocess for the paper_scenarios suite. One invocation = one
# benchmark point: it builds the requested Earth workload, runs PS_WARMUP warmup
# solves + PS_REPEATS timed solves, and prints exactly one machine-parseable line:
#
#   PS_RESULT ok=<bool> median_s=... min_s=... max_s=... times_s=a|b|c \
#             maxrss_mb=... workers_rss_mb=... detail=...
#
# All configuration arrives via environment variables (set by the s*_*.jl
# controllers; SPACEAGORA_* routing/GRAM knobs are set as real process env so
# Distributed pool workers inherit them):
#
#   PS_WORKLOAD       constellation | mc
#   PS_N_SATS         satellites per propagation (constellation size, or sats/sample)
#   PS_MISSION_S      mission length [s]
#   PS_ALT_KM         circular orbit altitude [km] (default 150 with GRAM, 550 without)
#   PS_GRAVITY        invsq | l20 | l50
#   PS_DENSITY        none | gram_standard | gram_lookahead | gram_surrogate
#   PS_PROFILE        "" or R0|R1_a|R2|R3|R4|R5 (wraps the run in with_parallel_profile)
#   PS_MC_SAMPLES     Monte Carlo sample count            (mc only)
#   PS_MC_BACKEND     serial | threads | process          (mc only)
#   PS_OUTER_WORKERS  outer campaign workers              (mc only)
#   PS_MC_MEMBERS     1 = seeds are constellation member indices, each sample
#                     propagates one member (process-ensemble mode); 0 = each
#                     sample is a full PS_N_SATS constellation with per-seed jitter
#   PS_WARMUP / PS_REPEATS
#
# Julia's --threads is set by the controller and is part of the point definition.

const PS_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PS_REPO_ROOT, "examples", "common.jl"))

using Statistics
using Distributed

const PS_WORKLOAD = get(ENV, "PS_WORKLOAD", "constellation")
const PS_N_SATS = parse(Int, get(ENV, "PS_N_SATS", "1"))
const PS_MISSION_S = parse(Float64, get(ENV, "PS_MISSION_S", "600.0"))
const PS_GRAVITY = get(ENV, "PS_GRAVITY", "l20")
const PS_DENSITY = get(ENV, "PS_DENSITY", "none")
const PS_PROFILE = get(ENV, "PS_PROFILE", "")
const PS_MC_SAMPLES = parse(Int, get(ENV, "PS_MC_SAMPLES", "8"))
const PS_MC_BACKEND = get(ENV, "PS_MC_BACKEND", "serial")
const PS_OUTER_WORKERS = parse(Int, get(ENV, "PS_OUTER_WORKERS", "1"))
const PS_MC_MEMBERS = get(ENV, "PS_MC_MEMBERS", "0") == "1"
const PS_WARMUP = parse(Int, get(ENV, "PS_WARMUP", "1"))
const PS_REPEATS = parse(Int, get(ENV, "PS_REPEATS", "3"))
const PS_ALT_KM = parse(Float64, get(ENV, "PS_ALT_KM", PS_DENSITY == "none" ? "550.0" : "150.0"))

const _NEEDS_GRAM = startswith(PS_DENSITY, "gram")
_NEEDS_GRAM && ensure_gramsuite_loaded!()

function _ps_density_model()
    PS_DENSITY == "none" && return NoAtmosphereModel()
    PS_DENSITY == "gram_surrogate" && return SimulationModel.GRAMAtmosphereModelSurrogate(planet_name="earth")
    # gram_standard and gram_lookahead share the native model; the calling method
    # (direct vs. vacuum-predicted cache) is selected purely via SPACEAGORA_* env.
    return SimulationModel.GRAMAtmosphereModel(planet_name="earth")
end

function _ps_gravity_effector(planet)
    PS_GRAVITY == "invsq" && return InverseSquaredGravityModel()
    lm = PS_GRAVITY == "l50" ? 50 : 20
    return GravitationalHarmonicsModel(
        lm, lm, joinpath(PS_REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )
end

# One Earth LEO constellation configuration. `seed` jitters altitude/RAAN so Monte
# Carlo samples are distinct constellations; `only_member` restricts the build to a
# single satellite of the constellation (ensemble-member mode), keeping the exact
# same orbit geometry formula so member and monolithic runs are comparable.
function ps_build_config(; n_sats::Int, seed::Int=0, only_member::Int=0)
    planet = Earth("", SPICE_PATH)
    gravity = _ps_gravity_effector(planet)
    effectors = PS_DENSITY == "none" ? (gravity,) : (gravity, AerodynamicCoefficientfM())
    alt_m = PS_ALT_KM * 1e3
    alt_jitter_m = 200.0 * seed
    raan_jitter_deg = 0.75 * seed

    member_range = only_member > 0 ? (only_member:only_member) : (1:n_sats)
    spacecraft = SpacecraftModel[]
    for i in member_range
        root = Link{0}(root=true, m=500.0, ref_area=12.0)
        phase = 50.0 * (i - 1) / max(n_sats, 1)
        ic = InitialCondition(
            ra=planet.Rp_e + alt_m + phase + alt_jitter_m,
            rp=planet.Rp_e + alt_m + phase + alt_jitter_m,
            i=53.0,
            ω=0.0,
            Ω=10.0 + raan_jitter_deg,
            ν=360.0 * (i - 1) / max(n_sats, 1)
        )
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, i))
    end

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=PS_MISSION_S,
            orientation_sim=false,
            num_steps_to_save=20
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=300.0,
            density_model=_ps_density_model(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

# ── Workload closures ─────────────────────────────────────────────────────────────
# Each returns (run_once::Function, describe::String). run_once returns true on success.

function ps_constellation_workload()
    args = ps_build_config(n_sats=PS_N_SATS)
    solve = () -> begin
        result = SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=true, return_solver_metadata=true)
        string(result.solution.retcode) == "Success"
    end
    run_once = if isempty(PS_PROFILE)
        solve
    else
        () -> SpaceAGORA.with_parallel_profile(solve, PS_PROFILE)
    end
    return run_once, "constellation n_sats=$(PS_N_SATS) gravity=$(PS_GRAVITY) density=$(PS_DENSITY) profile=$(PS_PROFILE)"
end

# Expression eval'd on each Distributed pool worker to define Main._ps_mc_sample.
# Interpolates every scalar it needs so the worker never depends on this script's
# globals; body only references SpaceAGORA/Base names (pool bootstrap runs
# `using SpaceAGORA` + GRAM warmup on every worker).
function _ps_mc_sample_defn_expr()::Expr
    density = PS_DENSITY
    gravity = PS_GRAVITY
    return quote
        function _ps_mc_sample(seed::Int)::Bool
            planet = Earth("", $(SPICE_PATH))
            gravity_eff = if $(gravity) == "invsq"
                InverseSquaredGravityModel()
            else
                lm = $(gravity) == "l50" ? 50 : 20
                GravitationalHarmonicsModel(
                    lm, lm,
                    joinpath($(PS_REPO_ROOT), "data", "Gravity_harmonics_data", "EarthGGM05C.csv"),
                    planet,
                )
            end
            density_model = if $(density) == "none"
                NoAtmosphereModel()
            elseif $(density) == "gram_surrogate"
                SpaceAGORA.SimulationModel.GRAMAtmosphereModelSurrogate(planet_name="earth")
            else
                SpaceAGORA.SimulationModel.GRAMAtmosphereModel(planet_name="earth")
            end
            effectors = $(density) == "none" ? (gravity_eff,) : (gravity_eff, AerodynamicCoefficientfM())
            alt_m = $(PS_ALT_KM) * 1e3
            n_sats = $(PS_N_SATS)
            members = $(PS_MC_MEMBERS)
            # Member mode: seed selects one satellite of the fixed constellation.
            # Sample mode: seed jitters a whole (possibly multi-sat) constellation.
            member_range = members ? (seed:seed) : (1:n_sats)
            alt_jitter_m = members ? 0.0 : 200.0 * seed
            raan_jitter_deg = members ? 0.0 : 0.75 * seed
            spacecraft = SpacecraftModel[]
            for i in member_range
                root = Link{0}(root=true, m=500.0, ref_area=12.0)
                phase = 50.0 * (i - 1) / max(n_sats, 1)
                ic = InitialCondition(
                    ra=planet.Rp_e + alt_m + phase + alt_jitter_m,
                    rp=planet.Rp_e + alt_m + phase + alt_jitter_m,
                    i=53.0,
                    ω=0.0,
                    Ω=10.0 + raan_jitter_deg,
                    ν=360.0 * (i - 1) / max(n_sats, 1)
                )
                push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, i))
            end
            cfg = SimulationConfiguration(
                simulation_settings=SimulationSettings(
                    results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
                ),
                mission_configuration=MissionConfiguration(
                    mission_type=MissionTime,
                    keplerian=true,
                    number_of_orbits=1,
                    mission_time=$(PS_MISSION_S),
                    orientation_sim=false,
                    num_steps_to_save=20
                ),
                environment_model=EnvironmentModel(
                    planet=planet,
                    EI=300.0,
                    density_model=density_model,
                    thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
                    topography=false,
                    wind=false
                ),
                dynamics_model=DynamicsModel(spacecraft, effectors),
                guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
                navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
                control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
                initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
                integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
            )
            result = SpaceAGORA.run_simulation(cfg; isolate_state=false, return_solution=true, return_solver_metadata=true)
            return string(result.solution.retcode) == "Success"
        end
        nothing
    end
end

const _PS_POOL_WORKER_IDS = Int[]

function ps_mc_workload()
    n_samples = PS_MC_MEMBERS ? PS_N_SATS : PS_MC_SAMPLES
    seeds = collect(1:n_samples)
    describe = "mc backend=$(PS_MC_BACKEND) samples=$(n_samples) sats_per_sample=$(PS_MC_MEMBERS ? 1 : PS_N_SATS) " *
               "outer=$(PS_OUTER_WORKERS) density=$(PS_DENSITY) members=$(PS_MC_MEMBERS)"

    if PS_MC_BACKEND == "process"
        pool = SpaceAGORA.campaign_process_pool()
        worker_ids = SpaceAGORA.ensure_process_workers!(pool, PS_OUTER_WORKERS)
        active = worker_ids[1:min(PS_OUTER_WORKERS, length(worker_ids))]
        append!(empty!(_PS_POOL_WORKER_IDS), active)
        # Pool bootstrap only runs `using SpaceAGORA`; the bare model names the
        # sample builder uses (Earth, SpacecraftModel, ...) are SimulationModel
        # exports, so bring those into each worker's Main first.
        Distributed.remotecall_eval(Main, active, :(using SpaceAGORA.SimulationModel))
        Distributed.remotecall_eval(Main, active, _ps_mc_sample_defn_expr())
        # Pay each worker's one-time run_simulation JIT cost here, outside the timed
        # region (see ensure_process_workers!'s warmup_fn rationale) — the pool's own
        # bootstrap only warms GRAM construction, not the full solve surface.
        @sync for w in active
            @async begin
                ok = remotecall_fetch(s -> Base.invokelatest(getfield(Main, :_ps_mc_sample), s), w, 1)
                ok || @warn "process warmup sample failed on worker $(w)"
            end
        end
        f = seed -> Base.invokelatest(getfield(Main, :_ps_mc_sample), seed)
        spec = SpaceAGORA.MonteCarloSpec(seeds=seeds, threads=1)
        run_once = () -> begin
            samples = SpaceAGORA.SimulationCampaigns._run_monte_carlo_process(f, seeds, spec, active)
            count(s -> s.success, samples) == n_samples
        end
        return run_once, describe
    end

    # serial/threads backends: sample function lives in this process.
    Core.eval(Main, _ps_mc_sample_defn_expr())
    f = seed -> Base.invokelatest(getfield(Main, :_ps_mc_sample), seed)
    mc_threads = PS_MC_BACKEND == "threads" ? PS_OUTER_WORKERS : 1
    run_once = () -> begin
        result = SpaceAGORA.run_monte_carlo(f, seeds; threads=mc_threads)
        length(result.successful) == n_samples
    end
    return run_once, describe
end

# ── Measurement ────────────────────────────────────────────────────────────────────

function main()
    run_once, describe = PS_WORKLOAD == "mc" ? ps_mc_workload() : ps_constellation_workload()
    println("[worker] $(describe) julia_threads=$(Threads.nthreads()) mission_s=$(PS_MISSION_S) " *
            "warmup=$(PS_WARMUP) repeats=$(PS_REPEATS)")
    flush(stdout)

    ok = true
    for w in 1:PS_WARMUP
        warm_ok = run_once()
        println("[worker] warmup $(w): ok=$(warm_ok)")
        flush(stdout)
        ok &= warm_ok
    end

    times = Float64[]
    for r in 1:PS_REPEATS
        GC.gc()
        t = @elapsed begin
            repeat_ok = run_once()
            ok &= repeat_ok
        end
        push!(times, t)
        println("[worker] repeat $(r): $(round(t; digits=3)) s")
        flush(stdout)
    end

    # Memory: this process's peak RSS, plus (process backend only) the summed peak
    # RSS of the Distributed pool workers — the paper's memory-vs-speed column.
    maxrss_mb = Sys.maxrss() / 1e6
    workers_rss_mb = 0.0
    if !isempty(_PS_POOL_WORKER_IDS)
        for w in _PS_POOL_WORKER_IDS
            workers_rss_mb += remotecall_fetch(() -> Sys.maxrss() / 1e6, w)
        end
    end

    med = median(times)
    println("PS_RESULT ok=$(ok) median_s=$(round(med; digits=4)) min_s=$(round(minimum(times); digits=4)) " *
            "max_s=$(round(maximum(times); digits=4)) times_s=$(join(round.(times; digits=4), '|')) " *
            "maxrss_mb=$(round(maxrss_mb; digits=1)) workers_rss_mb=$(round(workers_rss_mb; digits=1))")
    flush(stdout)

    # Shut the pool down explicitly: orphaned pool workers holding this process's
    # stdout pipe open were observed to hang controller-side wait() in prior studies.
    isempty(_PS_POOL_WORKER_IDS) || SpaceAGORA.shutdown_process_pool!(SpaceAGORA.campaign_process_pool())
    return nothing
end

main()
