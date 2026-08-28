#!/usr/bin/env julia
#
# Paired in-process A/B between two routing profiles on a MONTE CARLO CAMPAIGN.
#
#   julia --project=. --threads=12 scripts/paired_campaign_probe.jl \
#       --cases=montecarlo_mars_aerobraking --a=full_smart --b=outer_inner_static \
#       --mc-samples=64 --pairs=15
#
# scripts/paired_profile_probe.jl does this for a single simulation and settled
# the single-simulation questions to about one percent. It cannot run a campaign:
# it builds one SimulationConfiguration and calls run_simulation, so the outer
# axis -- which is the whole of the Monte Carlo question -- is never exercised.
# The benchmark harness does run campaigns, but it runs modes in BLOCKS with a
# fresh process per point, which is the ordering this project's own methodology
# notes say fakes effects, and on short campaigns its floor is around +/-15%.
#
# That gap is why the Monte Carlo residuals could not be read. Across fourteen
# measured points the adaptive profiles scatter roughly +/-9% around zero, R4 and
# R5 disagree in SIGN between adjacent thread counts on the same case, and
# per-sample cost does not predict the regret -- all signatures of noise rather
# than mechanism. Deciding whether Monte Carlo warrants its own policy needs an
# instrument that can tell a real few-percent effect from that scatter.
#
# The two changes that make it one, both borrowed from the single-simulation
# probe:
#
#   INTERLEAVING. Profiles alternate A,B,A,B and the reduction is over per-pair
#   ratios, so anything that drifts slowly -- thermal state, frequency, a
#   background process, the allocator settling -- scales both halves of a pair
#   together and cancels. The order within each pair alternates too, so a
#   first-position or second-position bias cancels across pairs instead of
#   accumulating.
#
#   A SIGNIFICANCE VERDICT. A two-sided sign test on pair wins, not a difference
#   of medians. It needs no distributional assumption, is immune to a single wild
#   pair, and can return "too close to call" -- which is the answer the Monte
#   Carlo residuals may well deserve, and which a median cannot express.
#
# The campaign itself is run by ppc_run_sample_batch from the shipped harness
# rather than reimplemented here. That function already resolves the outer route,
# applies the profile's environment, provisions and warms the Distributed pool
# outside the timed region, and returns the batch wall time. Duplicating any of
# it would drift, which is exactly what the single-simulation probe's header
# records happening to a copied env table within an hour.

const PPC_DIR = normpath(joinpath(@__DIR__, "..", "benchmarks", "studies",
                                  "parallelization_performance"))
include(joinpath(PPC_DIR, "cli.jl"))
include(joinpath(PPC_DIR, "modes.jl"))
include(joinpath(PPC_DIR, "cases.jl"))
include(joinpath(PPC_DIR, "trajectory_parity.jl"))
include(joinpath(PPC_DIR, "reporting.jl"))
include(joinpath(PPC_DIR, "execution.jl"))

using Printf
using Statistics

# ── CLI ───────────────────────────────────────────────────────────────────────
cases      = ["montecarlo_mars_aerobraking"]
a_name     = "full_smart"
b_name     = "outer_inner_static"
pairs      = 15
mc_samples = 64
seed       = 20260615
procs      = 2
cold_calib = false
isolate    = false
for arg in ARGS
    startswith(arg, "--cases=")      && (global cases = String.(split(split(arg, "=", limit=2)[2], ",")))
    startswith(arg, "--a=")          && (global a_name = String(split(arg, "=", limit=2)[2]))
    startswith(arg, "--b=")          && (global b_name = String(split(arg, "=", limit=2)[2]))
    startswith(arg, "--pairs=")      && (global pairs = parse(Int, split(arg, "=", limit=2)[2]))
    startswith(arg, "--mc-samples=") && (global mc_samples = parse(Int, split(arg, "=", limit=2)[2]))
    startswith(arg, "--seed=")       && (global seed = parse(Int, split(arg, "=", limit=2)[2]))
    startswith(arg, "--procs=")      && (global procs = parse(Int, split(arg, "=", limit=2)[2]))
    arg == "--cold-calib"            && (global cold_calib = true)
    arg == "--isolate-calib"         && (global isolate = true)
end

# A two-sided sign test on n pairs cannot return below 2/2^n, so fewer than six
# pairs can never reach 0.05 however lopsided the result: five wins out of five
# gives p = 0.0625. Refusing is better than reporting "not distinguishable" for a
# comparison the test was never able to resolve.
if pairs < 6
    error("--pairs must be at least 6; a sign test on $(pairs) pairs bottoms out at " *
          "p = $(round(2.0^(1 - pairs), digits=4)) and can never clear 0.05.")
end

const SPECS   = ppc_mode_specs()
const CATALOG = ppc_case_catalog()
haskey(SPECS, a_name) || error("unknown mode '$(a_name)'")
haskey(SPECS, b_name) || error("unknown mode '$(b_name)'")

# Per-arm calibration stores, so the two arms cannot feed each other cache hits.
# The calibration TOML is shared mutable state that both arms read and write, and
# the sweep's verdict is not stable run to run.
function _arm_calib_path(tag::String)
    dir = mktempdir(; prefix="spaceagora_campaign_probe_")
    return joinpath(dir, "calib_$(tag).toml")
end

# Drop the in-process RHS-calibration memo so the next campaign re-enters
# calibration as a fresh process would. The on-disk store is left alone, so a
# signature that HAS a stored plan reloads it cheaply while one that has none
# re-sweeps -- which is what a second `julia` invocation actually does.
function _drop_calibration_memo!()
    SE = SpaceAGORA.SimulationEngine
    empty!(SE._rhs_calib_cache)
    SE._rhs_calib_loaded[] = false
    return nothing
end

function make_cfg(mode_name::String, calib_path::Union{Nothing, String})
    over = Dict{String, String}()
    calib_path === nothing || (over["SPACEAGORA_RHS_CALIBRATION_PATH"] = calib_path)
    return (_ppc_with(PPCConfig();
                      worker_case="", worker_mode=mode_name, worker_seed=seed,
                      worker_threads=Threads.nthreads(), process_workers=procs,
                      warmup=0, mc_samples=[mc_samples], mc_samples_explicit=true),
            over)
end

# One timed campaign. Everything outside the batch -- pool provisioning, warm-up,
# route resolution -- is handled by ppc_run_sample_batch itself and, for the
# process pool, deliberately kept outside its timed region.
function run_campaign(case_name::String, mode_name::String,
                      calib_path::Union{Nothing, String})::Float64
    cfg, over = make_cfg(mode_name, calib_path)
    cfg = _ppc_with(cfg; worker_case=case_name)
    cold_calib && _drop_calibration_memo!()
    return withenv((k => v for (k, v) in over)...) do
        batch = ppc_run_sample_batch(CATALOG[case_name], cfg, SPECS[mode_name], mc_samples)
        all(r -> r.success, batch.results) ||
            error("campaign failed on $(case_name)/$(mode_name)")
        return Float64(batch.batch_wall_time_s)
    end
end

@printf("paired campaign probe: %s (A) vs %s (B)\n", a_name, b_name)
@printf("  %d pairs, %d samples/campaign, %d threads, %d process workers%s%s\n\n",
        pairs, mc_samples, Threads.nthreads(), procs,
        cold_calib ? ", cold calibration memo" : "",
        isolate ? ", isolated calibration stores" : "")

for case in cases
    calib_a = isolate ? _arm_calib_path("a_" * string(getpid())) : nothing
    calib_b = isolate ? _arm_calib_path("b_" * string(getpid())) : nothing

    # Warm both arms once. The first campaign in a process pays cold JIT of the
    # whole RHS/solver stack and, on the process route, worker startup; neither
    # belongs in a comparison of steady-state routing.
    run_campaign(case, a_name, calib_a)
    run_campaign(case, b_name, calib_b)

    ratios = Float64[]; wins_a = 0; wins_b = 0
    for i in 1:pairs
        ta, tb = if isodd(i)
            x = run_campaign(case, a_name, calib_a); y = run_campaign(case, b_name, calib_b); (x, y)
        else
            y = run_campaign(case, b_name, calib_b); x = run_campaign(case, a_name, calib_a); (x, y)
        end
        push!(ratios, ta / tb)
        ta < tb ? (wins_a += 1) : (tb < ta ? (wins_b += 1) : nothing)
    end

    sort!(ratios)
    med     = ratios[cld(length(ratios), 2)]
    decided = wins_a + wins_b
    p       = SpaceAGORA.SimulationModel.ParallelCost._sign_test_two_sided(
                  max(wins_a, wins_b), decided)
    verdict = p <= 0.05 ? (wins_a > wins_b ? "A FASTER" : "B FASTER") : "NOT DISTINGUISHABLE"
    @printf("%-30s  A wins %2d / B wins %2d   median A/B = %.3f (%+.1f%%)  p = %.4f  -> %s\n",
            case, wins_a, wins_b, med, 100 * (med - 1), p, verdict)
end
