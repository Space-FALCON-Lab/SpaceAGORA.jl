# Byte-for-byte dump of what one derivative evaluation writes besides the
# derivative: the aerodynamic save caches (drag, lift, cross) that the trajectory
# outputs record, and the density/temperature/wind/sample-time buffers the RHS
# atmosphere pre-sample fills. `variants.jl --dump` compares the state history;
# these outputs are not part of the state, so a change could alter them without
# moving a single state bit. This covers that gap.
#
# Each case is evaluated at its initial state and at a handful of later times
# with the state advanced by a crude, fixed Euler step (the step only has to be
# identical before and after, not accurate), so the evaluations see several
# distinct altitudes and densities.
#
# Usage (one Julia process at a time; check `uptime` first):
#   julia --project=. --threads=8 benchmarks/studies/aero_batch/rhs_dump.jl \
#       --cases=aero_1024sat_l50_expatm_100s,atmo256_exponential_10min \
#       --rhs=flat --dump=/path/prefix

using Printf

const AB_STUDY_DIR = @__DIR__
const AB_REPO_ROOT = normpath(joinpath(AB_STUDY_DIR, "..", "..", ".."))
const AB_PPC_DIR = joinpath(AB_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(AB_PPC_DIR, "cli.jl"))
include(joinpath(AB_PPC_DIR, "modes.jl"))
include(joinpath(AB_PPC_DIR, "cases.jl"))
using ComponentArrays: getdata

ab_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const AB_CASES = String.(split(ab_arg(ARGS, "cases",
    "aero_64sat_l50_expatm_100s,atmo256_exponential_10min,aero_1024sat_l50_expatm_100s,multi_64_high_fidelity"), ","))
const AB_RHS   = ab_arg(ARGS, "rhs", "flat")
const AB_DUMP  = ab_arg(ARGS, "dump", "")
const AB_EVALS = parse(Int, ab_arg(ARGS, "evals", "6"))
const AB_STEP  = parse(Float64, ab_arg(ARGS, "step", "20.0"))

const SE_AB = SimulationEngine

function ab_write_vec3s(io, v, n)
    for i in 1:n
        write(io, Float64(v[i][1]), Float64(v[i][2]), Float64(v[i][3]))
    end
end

function main()
    isempty(AB_DUMP) && error("--dump=<prefix> is required")
    @printf("host=%s threads=%d rhs_mode=%s cases=%s\n",
            gethostname(), Threads.nthreads(), AB_RHS, join(AB_CASES, ","))
    for case in AB_CASES
        args = ppc_single_config(case, PPCConfig(profile="full"))
        n = length(args.dynamics_model.spacecraft)
        u = SE_AB.build_initial_conditions(args)
        du = zero(u)
        p = ODEParams(n_sats=n, args=args)
        SE_AB._initialize_save_cache_buffers!(p)
        SE_AB._initialize_heat_rate_buffers!(p)
        path = "$(AB_DUMP)_$(case).bin"
        mode = ""
        withenv("SPACEAGORA_RHS_EXECUTION_MODE" => AB_RHS,
                "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())) do
            SE_AB._initialize_runtime_env_config!(p)
            mode = string(SE_AB._rhs_execution_plan(args, p, args.dynamics_model.dynamic_effectors, n).mode)
            open(path, "w") do io
                t = 0.0
                for k in 1:AB_EVALS
                    SE_AB.spacecraft_dynamics!(du, u, p, t)
                    sb = p.shared_buffers
                    write(io, t)
                    write(io, Float64.(vec(getdata(du))))
                    ab_write_vec3s(io, p.save_cache.drag_cache, n)
                    ab_write_vec3s(io, p.save_cache.lift_cache, n)
                    ab_write_vec3s(io, p.save_cache.cross_cache, n)
                    write(io, Float64.(sb.densities[1:n]))
                    write(io, Float64.(sb.temperatures[1:n]))
                    ab_write_vec3s(io, sb.winds, n)
                    write(io, Float64.(sb.density_sample_t[1:n]))
                    # Advance: identical arithmetic before and after, so the
                    # next evaluation's inputs agree exactly iff this one's
                    # derivative did.
                    getdata(u) .+= AB_STEP .* getdata(du)
                    t += AB_STEP
                end
            end
        end
        @printf("%-34s N=%5d route=%s -> %s (%d bytes)\n", case, n, mode, path, filesize(path))
    end
end

main()
