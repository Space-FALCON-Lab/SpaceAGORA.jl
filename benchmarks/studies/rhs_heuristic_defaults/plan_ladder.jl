# WS11a: what every candidate RHS plan costs, size by size, against the plan the
# heuristic picks.
#
# The pinned parallel routes run whatever `_rhs_execution_plan_uncached`
# (src/simulation/engine/setup.jl) decides; the adaptive routes run whatever the
# calibration sweep measured. On the archived TRX50 P1 ladder those two answers
# differ by up to 10x in the middle of the constellation-size range, so the
# default is what needs fixing, not the sweep. This ladder measures the whole
# candidate set the sweep would consider (`_rhs_plan_candidates`,
# src/simulation/engine/rhs_calibration.jl) plus the heuristic's own plan, at
# one thread budget, in one process, back to back.
#
# How a plan is pinned: through the production calibration cache, not a private
# hook. Each plan gets its own one-row TOML store, the path is handed to the
# engine with SPACEAGORA_RHS_CALIBRATION_PATH, and
# SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S is raised past any solve length here so
# the cached verdict is honored rather than re-swept (`_rhs_calib_cached_verdict`).
# The `heuristic` row instead runs SPACEAGORA_RHS_CALIBRATE=off, which is exactly
# what the pinned routes in the paper harness run.
#
# Case: the P1 iso-work L50 vacuum ladder (`gravity_<N>sat_l50_vacuum_<S>s` in
# benchmarks/studies/parallelization_performance/cases.jl), built here from the
# same helpers so the shape is identical.
#
# Usage:
#   julia --project=. --threads=8 benchmarks/studies/rhs_heuristic_defaults/plan_ladder.jl \
#       --sizes=16,32,64,128,256,512,1024 --repeats=3 [--out=results/<file>.csv]
#
# Only ratios within one invocation are meaningful; the machine is shared.

using Printf
using Dates
using TOML
using ComponentArrays: getdata

const RHD_DIR = @__DIR__
const RHD_REPO_ROOT = normpath(joinpath(RHD_DIR, "..", "..", ".."))
const RHD_PPC_DIR = joinpath(RHD_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(RHD_PPC_DIR, "cli.jl"))
include(joinpath(RHD_PPC_DIR, "modes.jl"))
include(joinpath(RHD_PPC_DIR, "cases.jl"))

rhd_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const RHD_SIZES = parse.(Int, split(rhd_arg(ARGS, "sizes", "16,32,64,128,256,512,1024"), ","))
const RHD_REPEATS = parse(Int, rhd_arg(ARGS, "repeats", "3"))
const RHD_WARMUP_S = parse(Float64, rhd_arg(ARGS, "warmup-mission", "60"))
# Functional smoke only: overrides every rung's mission length so the whole
# matrix runs in seconds. Never use for a measured row.
const RHD_MISSION_OVERRIDE = parse(Float64, rhd_arg(ARGS, "mission", "0"))
# Which shape the ladder runs. `vacuum` is the P1 L50 vacuum constellation;
# `aero` is the L50 harmonics + aerodynamic-coefficient stack over an analytic
# exponential atmosphere, i.e. the two-effector reference case.
const RHD_CASE = rhd_arg(ARGS, "case", "vacuum")
# When set, every timed solve also writes its full state history (every saved
# time and every component of every spacecraft, raw Float64) to
# <dir>/<case>_<N>_<plan>_rep<k>.bin, for a byte-for-byte `cmp`.
const RHD_DUMP_DIR = rhd_arg(ARGS, "dump-dir", "")
const RHD_OUT = rhd_arg(ARGS, "out", joinpath(RHD_DIR, "results",
    "plan_ladder_$(rhd_arg(ARGS, "case", "vacuum"))_$(gethostname())_t$(Threads.nthreads()).csv"))
const RHD_PINDIR = mktempdir(; prefix="rhd_pins_")

# Simulated mission seconds per constellation size.
#
# SOURCED for 16/64/256/1024/4096 from PPC_L50_ISO_MISSION_S
# (benchmarks/studies/parallelization_performance/cases.jl), so those rungs are
# the same workload the archived P1 runs measured. DERIVED for 32/128/512 as the
# geometric mean of the two neighboring sourced rungs, which keeps the ladder's
# own shape (its duration falls with N) without inventing a new one.
const RHD_MISSION_S = Dict{Int, Float64}(
    1 => 4150000.0,
    16 => 514000.0,
    32 => 462000.0,     # DERIVED: sqrt(514000 * 415000)
    64 => 415000.0,
    128 => 227000.0,    # DERIVED: sqrt(415000 * 124000)
    256 => 124000.0,
    512 => 55000.0,     # DERIVED: sqrt(124000 * 24600)
    1024 => 24600.0,
    4096 => 5800.0,
)

const RHD_PLANET = Earth("", PPC_SPICE_PATH)

function rhd_effectors()
    RHD_CASE == "vacuum" && return (ppc_harmonics_model(RHD_PLANET, 50),)
    RHD_CASE == "aero" &&
        return (ppc_harmonics_model(RHD_PLANET, 50), AerodynamicCoefficientfM())
    error("unknown --case=$(RHD_CASE); expected vacuum or aero")
end

rhd_density_model() =
    RHD_CASE == "aero" ? ExponentialAtmosphereModel(RHD_PLANET) : NoAtmosphereModel()

function rhd_config(n::Int, mission_s::Float64)
    return ppc_build_config(
        planet=RHD_PLANET,
        spacecraft=ppc_constellation(RHD_PLANET, n),
        mission_time_s=mission_s,
        orientation_sim=false,
        dynamic_effectors=rhd_effectors(),
        density_model=rhd_density_model(),
        dt_max_orbit=20.0,
    )
end

# The plan set the sweep would consider at this budget and size, mirroring
# `_rhs_plan_candidates`, plus the heuristic. Kept as a local mirror rather than
# a call into the private function because that one needs a live ODEParams.
function rhd_plans(n::Int, budget::Int)
    plans = Any[(label="heuristic", mode="", allotment=0, scheduler="")]
    push!(plans, (label="batch@1", mode="satellite_batch", allotment=1, scheduler="static"))
    if budget >= 4
        for w in (budget ÷ 2, budget ÷ 4)
            w >= 2 && push!(plans,
                (label="batch@$(w)", mode="satellite_batch", allotment=w, scheduler="static"))
        end
    end
    min_sats_floor = SimulationModel.ParallelPolicy.harmonics_batch_spin_barrier_enabled() ? 1 : 4
    viable_workers = fld(n, max(1, min_sats_floor))
    if viable_workers >= 2
        max_workers = max(1, min(budget, viable_workers))
        allotments = Int[1]
        a = 2
        while a < max_workers
            push!(allotments, a)
            a *= 2
        end
        push!(allotments, max_workers)
        sort!(unique!(allotments))
        for scheduler in ("static", "dynamic"), al in allotments
            al >= 2 || continue
            push!(plans, (label="flat@$(al)/$(scheduler)",
                          mode="flat_constellation_effector_queue",
                          allotment=al, scheduler=scheduler))
        end
        1 in allotments && push!(plans,
            (label="flat@1/static", mode="flat_constellation_effector_queue",
             allotment=1, scheduler="static"))
    end
    return plans
end

# The calibration signature this shape will produce, assembled from the engine's
# own pieces so a pinned row cannot miss.
function rhd_signature(n::Int, effectors, density_model)
    return join([
        "v6",
        "machine=$(SimulationEngine._calib_machine_label())",
        "budget=$(SimulationModel.ParallelPolicy.effective_inner_thread_budget())",
        "sats=$(SimulationEngine._calib_sat_bucket(n))",
        "effs=$(length(effectors))",
        "harm=$(length(effectors) == 1 && effectors[1] isa SimulationModel.GravitationalHarmonicsModel ? "1" : "0")",
        "eff=$(SimulationEngine._rhs_calib_effector_token(effectors))",
        "dens=$(string(nameof(typeof(density_model))))",
        "outer=$(SimulationModel.ParallelPolicy.outer_parallel_active() ? "1" : "0")",
    ], "|")
end

function rhd_write_pin(path::String, sig::String, plan)
    payload = Dict{String, Any}(
        "schema_version" => 1,
        "calibrations" => [Dict{String, Any}(
            "signature" => sig,
            "mode" => plan.mode,
            "allotment" => plan.allotment,
            "scheduler" => plan.scheduler,
            "elapsed_mean_ns" => 1.0,
            # Non-zero so the solve-length gate reads "measured"; tiny so the
            # verdict is honored rather than re-swept. See the header.
            "solve_ns" => 1.0,
            "heuristic_votes" => 0,
            "sweep_ns" => 1.0,
            "honoured_ns" => 0.0,
            "plan_votes" => 9,
        )],
    )
    open(path, "w") do io
        TOML.print(io, payload)
    end
    return path
end

function rhd_env_pairs(plan, sig::String, n::Int)
    mode = ppc_mode_specs()["inner_only"]
    cfg = PPCConfig(profile="full", solver_mode="auto_stiff")
    pairs = copy(ppc_mode_env_pairs(mode, cfg; outer_tasks=1))
    if plan.label == "heuristic"
        push!(pairs, "SPACEAGORA_RHS_CALIBRATE" => "off")
    else
        path = joinpath(RHD_PINDIR, "pin_$(n)_$(replace(plan.label, "/" => "_", "@" => "at")).toml")
        rhd_write_pin(path, sig, plan)
        push!(pairs, "SPACEAGORA_RHS_CALIBRATE" => "auto")
        push!(pairs, "SPACEAGORA_RHS_CALIBRATION_PATH" => path)
        # Every solve here is shorter than this, so the cached verdict is always
        # honored and the sweep never runs.
        push!(pairs, "SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S" => "1.0e9")
    end
    return pairs
end

function rhd_solve(args)
    timed = @timed SimulationEngine.run_simulation(
        args; isolate_state=false, return_solution=true, return_solver_metadata=true
    )
    r = timed.value
    sol = r.solution
    # From the solve's own return value: the engine records every routing
    # decision into its own scoped PolicyContext, so a snapshot taken out here
    # would read the global context and report zeros (see
    # benchmarks/studies/parallelization_performance/execution.jl).
    telem = get(r, :parallel_policy, nothing)
    tget(k, default) = telem === nothing ? default :
        (hasproperty(telem, k) ? getproperty(telem, k) : default)
    geti(name) = try Int(getproperty(sol.stats, name)) catch; -1 end
    sc = sol.u[end].sc[1]
    return (
        wall_s=Float64(timed.time), gc_s=Float64(timed.gctime), bytes=Int(timed.bytes),
        nf=geti(:nf), naccept=geti(:naccept), nreject=geti(:nreject),
        steps=length(sol.t), retcode=string(sol.retcode),
        pos=Float64.(sc.pos), vel=Float64.(sc.vel),
        applied_source=string(tget(:rhs_plan_source, "none")),
        applied_mode=string(tget(:rhs_plan_mode, "none")),
        applied_allotment=Int(tget(:rhs_plan_allotment, 0)),
        applied_scheduler=string(tget(:rhs_plan_scheduler, "none")),
        sol=sol,
    )
end

# Full state history, raw Float64, no formatting: every saved time followed by
# every component of every spacecraft at that time.
function rhd_dump_states(path::String, sol)
    open(path, "w") do io
        write(io, Float64.(sol.t))
        for u in sol.u
            write(io, Float64.(vec(getdata(u))))
        end
    end
    return path
end

function main()
    budget = SimulationModel.ParallelPolicy.effective_inner_thread_budget()
    @printf("host=%s julia_threads=%d inner_budget=%d case=%s sizes=%s repeats=%d\n",
            gethostname(), Threads.nthreads(), budget, RHD_CASE,
            join(RHD_SIZES, ","), RHD_REPEATS)
    println("load-at-start: ", strip(read(`uptime`, String)))
    println("pin store dir: ", RHD_PINDIR)
    flush(stdout)

    mkpath(dirname(RHD_OUT))
    io = open(RHD_OUT, "w")
    println(io, join([
        "timestamp_utc", "host", "julia_threads", "inner_budget", "case", "n_sats", "mission_s",
        "plan_label", "plan_mode", "plan_allotment", "plan_scheduler", "rep",
        "wall_s", "gc_s", "alloc_bytes", "nf", "naccept", "nreject", "steps", "retcode",
        "applied_source", "applied_mode", "applied_allotment", "applied_scheduler",
        "final_pos_x", "final_pos_y", "final_pos_z",
        "final_vel_x", "final_vel_y", "final_vel_z",
    ], ","))
    flush(io)

    for n in RHD_SIZES
        mission_s = RHD_MISSION_OVERRIDE > 0 ? RHD_MISSION_OVERRIDE : get(RHD_MISSION_S, n) do
            error("no mission length registered for N=$(n); add one to RHD_MISSION_S")
        end
        plans = rhd_plans(n, budget)
        args_probe = rhd_config(n, 10.0)
        sig = rhd_signature(n, args_probe.dynamics_model.dynamic_effectors,
                            args_probe.environment_model.density_model)
        @printf("\n== N=%d mission=%.0fs plans=%d\n   sig=%s\n",
                n, mission_s, length(plans), sig)
        flush(stdout)

        # Warm every plan on a short mission first, so no timed row pays
        # compilation for its own route.
        for plan in plans
            withenv(rhd_env_pairs(plan, sig, n)...) do
                rhd_solve(rhd_config(n, RHD_WARMUP_S))
            end
        end
        println("   warmed"); flush(stdout)

        for rep in 1:RHD_REPEATS, plan in plans
            r = withenv(rhd_env_pairs(plan, sig, n)...) do
                rhd_solve(rhd_config(n, mission_s))
            end
            println(io, join([
                string(Dates.now(Dates.UTC)), gethostname(), Threads.nthreads(), budget,
                RHD_CASE, n, mission_s, plan.label, plan.mode, plan.allotment, plan.scheduler, rep,
                @sprintf("%.6f", r.wall_s), @sprintf("%.6f", r.gc_s), r.bytes,
                r.nf, r.naccept, r.nreject, r.steps, r.retcode,
                r.applied_source, r.applied_mode, r.applied_allotment, r.applied_scheduler,
                (@sprintf("%.17g", x) for x in r.pos)...,
                (@sprintf("%.17g", x) for x in r.vel)...,
            ], ","))
            flush(io)
            if !isempty(RHD_DUMP_DIR)
                mkpath(RHD_DUMP_DIR)
                rhd_dump_states(joinpath(RHD_DUMP_DIR, string(
                    RHD_CASE, "_", n, "_",
                    replace(plan.label, "/" => "-", "@" => "at"), "_rep", rep, ".bin")), r.sol)
            end
            @printf("   rep%d %-18s wall=%8.3fs gc=%6.3fs nf=%8d acc=%6d rej=%5d %s applied=%s/%s@%d/%s\n",
                    rep, plan.label, r.wall_s, r.gc_s, r.nf, r.naccept, r.nreject, r.retcode,
                    r.applied_source, r.applied_mode, r.applied_allotment, r.applied_scheduler)
            flush(stdout)
        end
    end
    close(io)
    println("\nwrote ", RHD_OUT)
    println("load-at-end: ", strip(read(`uptime`, String)))
end

main()
