# WS10b variant ladder: which change actually makes the SRP + third-body
# constellation expensive.
#
# Every rung keeps the P6 trace's constellation, tolerances, step cap and
# mission and changes exactly one thing, so the ladder reads as an attribution
# rather than a benchmark. All rungs are solved back to back in one process at
# one thread; only ratios between rungs are meaningful.
#
#   vacuum        harmonics L50 only                      (the P6 trace 2 rung)
#   nbody         + Sun/Moon third-body gravity
#   srp           + solar radiation pressure
#   srp_nbody     + both                                  (the P6 trace 3 rung)
#
# --rhs=satellite forces every rung onto the per-satellite RHS route, and
# --solver=tsit5 forces every rung onto the explicit solver, which is how the
# route effect and the solver effect are separated from each other.
#
# Usage:
#   julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl \
#       --n=256 --mission=5800 [--variants=vacuum,nbody,srp,srp_nbody]
#       [--solver=auto_stiff|tsit5] [--rhs=auto|satellite|flat] [--repeats=1]

using Printf
using Profile

const TBC_STUDY_DIR = @__DIR__
const TBC_REPO_ROOT = normpath(joinpath(TBC_STUDY_DIR, "..", "..", ".."))
const TBC_PPC_DIR = joinpath(TBC_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(TBC_PPC_DIR, "cli.jl"))
include(joinpath(TBC_PPC_DIR, "modes.jl"))
include(joinpath(TBC_PPC_DIR, "cases.jl"))
using ComponentArrays: getdata

tbc_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const TBC_N       = parse(Int, tbc_arg(ARGS, "n", "256"))
const TBC_MISSION = parse(Float64, tbc_arg(ARGS, "mission", "5800"))
const TBC_REPEATS = parse(Int, tbc_arg(ARGS, "repeats", "1"))
const TBC_SOLVER  = tbc_arg(ARGS, "solver", "auto_stiff")
const TBC_RHS     = tbc_arg(ARGS, "rhs", "auto")
const TBC_MODE = tbc_arg(ARGS, "mode", "serial")
const TBC_DUMP = tbc_arg(ARGS, "dump", "")
const TBC_PROFILE = any(a -> a == "--profile", ARGS)
const TBC_VARIANTS = String.(split(tbc_arg(ARGS, "variants", "vacuum,nbody,srp,srp_nbody"), ","))

const TBC_PLANET = Earth("", PPC_SPICE_PATH)

function tbc_effectors(variant::String)
    harmonics = ppc_harmonics_model(TBC_PLANET, 50)
    srp = SolarRadiationPressureModel(1.2, 12.0)
    nbody = NBodyGravityModel(
        body_names=("Sun", "Moon"), primary_body_name="Earth", planet=TBC_PLANET
    )
    variant == "vacuum"    && return (harmonics,)
    variant == "nbody"     && return (harmonics, nbody)
    variant == "srp"       && return (harmonics, srp)
    variant == "srp_nbody" && return (harmonics, srp, nbody)
    error("unknown variant $variant")
end

function tbc_config(variant::String, mission_s::Float64)
    return ppc_build_config(
        planet=TBC_PLANET,
        spacecraft=ppc_constellation(TBC_PLANET, TBC_N),
        mission_time_s=mission_s,
        orientation_sim=false,
        dynamic_effectors=tbc_effectors(variant),
        density_model=NoAtmosphereModel(),
        dt_max_orbit=20.0
    )
end

function tbc_counters(counters)
    counters === nothing && return (nbody=0, srp=0, pxform=0)
    geti(k) = (v = getproperty(counters, k); v isa Base.Threads.Atomic ? v[] : Int(v))
    return (nbody=geti(:nbody_spkpos_runtime_calls),
            srp=geti(:srp_spkpos_runtime_calls),
            pxform=geti(:planet_pxform_runtime_calls))
end

# solver_trace is whatever the engine returns (a vector of per-phase entries
# in the shipped path); reduce it to the solver label(s) actually used.
function tbc_solver_label(trace)
    trace === nothing && return "?"
    entries = trace isa AbstractVector ? trace : [trace]
    labels = String[]
    for e in entries
        lbl = try
            string(getproperty(e, :solver))
        catch
            string(e)
        end
        fb = try
            getproperty(e, :fallback_used) === true ? "+switched" : ""
        catch
            ""
        end
        push!(labels, lbl * fb)
    end
    return join(labels, "|")
end

function tbc_solve(variant::String, mission_s::Float64)
    args = tbc_config(variant, mission_s)
    timed = @timed SimulationEngine.run_simulation(
        args; isolate_state=false, return_solution=true, return_solver_metadata=true
    )
    r = timed.value
    sol = r.solution
    stats = sol.stats
    geti(name) = try Int(getproperty(stats, name)) catch; -1 end
    sc = sol.u[end].sc[1]
    return (
        wall_s=Float64(timed.time), gc_s=Float64(timed.gctime), bytes=Int(timed.bytes),
        nf=geti(:nf), njacs=geti(:njacs), nw=geti(:nw), nsolve=geti(:nsolve),
        naccept=geti(:naccept), nreject=geti(:nreject),
        steps=length(sol.t), retcode=string(sol.retcode),
        solver=tbc_solver_label(r.solver_trace),
        counters=tbc_counters(get(r, :spice_counters, nothing)),
        pos=Float64.(sc.pos), vel=Float64.(sc.vel),
        sol=sol,
    )
end

function main()
    mode = ppc_mode_specs()[TBC_MODE]
    cfg = PPCConfig(profile="full", solver_mode=TBC_SOLVER)
    envpairs = copy(ppc_mode_env_pairs(mode, cfg; outer_tasks=1))
    if TBC_RHS != "auto"
        push!(envpairs, "SPACEAGORA_RHS_EXECUTION_MODE" => TBC_RHS)
    end

    @printf("host=%s threads=%d N=%d mission=%.0fs solver=%s rhs_mode=%s parallel_mode=%s\n",
            gethostname(), Threads.nthreads(), TBC_N, TBC_MISSION, TBC_SOLVER, TBC_RHS, TBC_MODE)
    println("load-at-start: ", strip(read(`uptime`, String)))
    flush(stdout)

    for v in TBC_VARIANTS
        withenv(envpairs...) do
            tbc_solve(v, 10.0)
        end
    end
    println("warmed all variants"); flush(stdout)

    rows = Dict{String, Vector{NamedTuple}}()
    for rep in 1:TBC_REPEATS, v in TBC_VARIANTS
        r = withenv(envpairs...) do
            tbc_solve(v, TBC_MISSION)
        end
        if !isempty(TBC_DUMP)
            # Full state history, raw Float64, for a byte-for-byte before/after
            # comparison of the trajectory (no tolerance).
            open("$(TBC_DUMP)_$(v).bin", "w") do io
                write(io, Float64.(r.sol.t))
                for u in r.sol.u
                    write(io, Float64.(vec(getdata(u))))
                end
            end
        end
        push!(get!(rows, v, NamedTuple[]), r)
        @printf("rep%d %-10s wall=%8.3fs gc=%6.3fs alloc=%7.2fGiB nf=%7d njacs=%5d nw=%5d nsolve=%6d acc=%5d rej=%4d solver=%-18s spkpos(nb=%d srp=%d) pxform=%d\n",
                rep, v, r.wall_s, r.gc_s, r.bytes/2^30, r.nf, r.njacs, r.nw, r.nsolve,
                r.naccept, r.nreject, string(r.solver), r.counters.nbody, r.counters.srp, r.counters.pxform)
        flush(stdout)
    end

    if haskey(rows, "vacuum")
        base = minimum(r.wall_s for r in rows["vacuum"])
        base_nf = rows["vacuum"][1].nf
        println("\n-- relative to vacuum (same process, back to back) --")
        for v in TBC_VARIANTS
            haskey(rows, v) || continue
            w = minimum(r.wall_s for r in rows[v])
            nf = rows[v][1].nf
            @printf("%-10s wall x%.2f   nf x%.2f   per-RHS-eval x%.2f\n",
                    v, w/base, nf/base_nf, (w/nf)/(base/base_nf))
        end
    end
    println("\n-- terminal state of spacecraft 1 (for cross-run comparison) --")
    for v in TBC_VARIANTS
        haskey(rows, v) || continue
        r = rows[v][1]
        @printf("%-10s pos=[%.17g, %.17g, %.17g] vel=[%.17g, %.17g, %.17g]\n",
                v, r.pos[1], r.pos[2], r.pos[3], r.vel[1], r.vel[2], r.vel[3])
    end
    if TBC_PROFILE
        for v in TBC_VARIANTS
            Profile.clear()
            Profile.init(n = 20_000_000, delay = 0.0005)
            withenv(envpairs...) do
                @profile tbc_solve(v, TBC_MISSION)
            end
            tag = "$(v)_$(TBC_N)sat_$(round(Int, TBC_MISSION))s_$(TBC_SOLVER)_$(TBC_RHS)"
            open(joinpath(TBC_STUDY_DIR, "profile_flat_$(tag).txt"), "w") do io
                Profile.print(IOContext(io, :displaysize => (24, 2000)); format=:flat, C=true, sortedby=:count, mincount=20)
            end
            println("profile written: ", tag); flush(stdout)
        end
    end
    println("load-at-end: ", strip(read(`uptime`, String)))
end

main()
