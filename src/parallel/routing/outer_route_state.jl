"""
    OuterRouteFeatures

Feature vector describing a simulation workload for outer parallel-route
selection.
"""
Base.@kwdef struct OuterRouteFeatures
    category::String = "deterministic"
    n_sats::Int = 1
    n_links::Int = 1
    max_links_per_sat::Int = 1
    mission_time_s::Float64 = 0.0
    has_nbody::Bool = false
    has_srp::Bool = false
    harmonics_degree::Int = 0
    has_control::Bool = false
    orientation_on::Bool = false
    density_family::String = "unknown"
    solver_mode::String = "auto"
    dt_max_orbit_s::Float64 = 0.0
    control_rate_s::Float64 = 0.0
    guidance_rate_s::Float64 = 0.0
    navigation_rate_s::Float64 = 0.0
    gram_surrogate_enabled::Bool = false
    gram_static_grid_enabled::Bool = false
    control_effector_count::Int = 0
    thermal_enabled::Bool = false
    dynamic_effector_count::Int = 0
    effector_cost_class::String = "unknown"
    montecarlo_samples::Int = 0
end

"""
    OuterRouteTuning

Tunable thresholds and adaptive-selection parameters for outer parallel-route
policy.
"""
Base.@kwdef struct OuterRouteTuning
    inner_sat_threshold::Int = 8
    inner_link_threshold::Int = 12
    outer_light_sat_threshold::Int = 2
    outer_light_link_threshold::Int = 4
    outer_light_mission_threshold_s::Float64 = 14_400.0
    spice_constellation_process_enabled::Bool = true
    spice_constellation_min_sats::Int = 4
    adaptive_enabled::Bool = true
    adaptive_min_samples::Int = 2
    adaptive_exploration_c::Float64 = 1.25
    failure_penalty_s::Float64 = 120.0
    mc_process_min_samples::Int = 16
    mc_process_min_mission_s::Float64 = 3600.0
    # Process workers run --threads=1 each and don't share the coordinator's
    # thread pool, so they are capped by physical parallelism rather than by
    # Threads.nthreads() like the thread route.
    #
    # `usable_core_budget()`, not `Sys.CPU_THREADS`. The latter counts SMT
    # siblings, ignores an affinity mask and ignores a cgroup quota, so on the
    # 12-physical/24-logical reference box this default was twice the useful
    # width -- and every workload measured there regressed from 16 to 24
    # threads. Under `taskset -c 0-3` it claimed 24 where four cores were
    # available, and in a container it claims the host's count while the process
    # is throttled to a fraction of it.
    #
    # This CHANGES ROUTING DECISIONS on any SMT or constrained machine, which is
    # the point of the change rather than a side effect of it: numbers measured
    # against the previous default describe a configuration that was asking for
    # more parallelism than the machine could deliver.
    process_max_workers::Int = _outer_process_worker_cap()
    # Stop forced exploration once ANY candidate is proven best, rather than
    # only when the default is. See _any_candidate_proven. Follows
    # SPACEAGORA_PARALLEL_POLICY_V2 so the shipped selector stays reproducible.
    explore_until_any_proven::Bool = outer_route_policy_v2()
    # Threads the campaign's thread route can split across. The Monte Carlo
    # default compares it against process_max_workers (see
    # _priority_outer_route_montecarlo). Fixed for the life of the process.
    outer_thread_budget::Int = Base.Threads.nthreads()
    # Route a Monte Carlo campaign to whichever outer axis has more cores.
    # Follows SPACEAGORA_PARALLEL_POLICY_V2; off keeps the shipped threads
    # default.
    mc_route_by_core_budget::Bool = outer_route_policy_v2()
    # Race the split widths inside the first campaign of an unseen signature
    # instead of trying one width per campaign. Follows
    # SPACEAGORA_PARALLEL_POLICY_V2. See SimulationCampaigns._run_campaign_split_race.
    split_race::Bool = outer_route_policy_v2()
    trace::Bool = false
end

# Off restores the previous behaviour, in which an arm's cold first timing was
# recorded as if it were its steady-state cost. Exists so the isolating A/B has
# a B side; there is no reason to run with it off.
# Whether a default arm must be measured against EVERY enumerated alternative
# before it counts as proven, rather than against any one of them.
#
# DEFAULT OFF, and that is a measured reversal of the obvious guess. Forcing the
# selector through the whole ladder looked like the fix for a width optimum it
# never found, and it is not: with the ladder trimmed to plausible widths (see
# outer_split_candidates) the ORDINARY guard finds the optimum on every campaign
# measured, and forcing full exploration is worse on all of them. Steady state
# over runs 11-20, 64 samples at 12 threads, off against on:
#
#     montecarlo_multi_sat        w8  0.0155 s   vs   w4  0.0220 s
#     montecarlo_high_accuracy    w8  0.0230 s   vs   w4  0.0345 s
#     montecarlo_mars_aerobraking w12 0.0435 s   vs   w12 0.0445 s
#
# Three noisy arms can be ranked; five cannot, and the extra exploration buys a
# worse answer as well as costing more to reach. The real defect was the useless
# narrow arms satisfying the guard's "tested against an alternative" without
# being informative, which trimming the ladder removes at the source.
#
# Kept as a knob because it is the natural hypothesis and someone will have it
# again; the numbers above are why it is not the default.
@inline function proven_requires_all_candidates()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_OUTER_ROUTE_PROVEN_ALL", "0")))
    return raw in ("1", "true", "yes", "on")
end

# Routing-layer read of SPACEAGORA_PARALLEL_POLICY_V2. ParallelProfiles is
# assembled before ParallelPolicy, so it cannot borrow policy_v2_enabled; the
# parse is the same one proven_requires_all_candidates uses. Read once per
# OuterRouteTuning construction, i.e. once per campaign, never on a hot path.
@inline function outer_route_policy_v2()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PARALLEL_POLICY_V2", "0")))
    return raw in ("1", "true", "yes", "on")
end

# Process workers the campaign pool may use.
#
# `usable_core_budget()` as shipped. Under SPACEAGORA_PARALLEL_POLICY_V2 also
# capped by SPACEAGORA_PERF_PROCS when it is set: that variable is how the
# benchmark harness (and the HPC launch scripts) state how many worker
# processes a run was given, and a router that ignores it sizes the pool from
# the whole machine. On the B13 budget grid the runtime would have spawned
# twelve workers at a point whose budget was two, and the core-budget default
# below would then have compared twelve against six threads and chosen wrong.
@inline function _outer_process_worker_cap()::Int
    cap = usable_core_budget()
    outer_route_policy_v2() || return cap
    raw = strip(get(ENV, "SPACEAGORA_PERF_PROCS", ""))
    isempty(raw) && return cap
    v = tryparse(Int, raw)
    (v === nothing || v <= 0) && return cap
    return max(1, min(cap, v))
end

@inline function outer_route_discard_cold_observation()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_OUTER_ROUTE_DISCARD_COLD", "1")))
    return raw in ("1", "true", "yes", "on")
end

Base.@kwdef mutable struct OuterRouteStats
    samples::Int = 0
    successes::Int = 0
    failures::Int = 0
    elapsed_sum_s::Float64 = 0.0
    elapsed_sq_sum_s::Float64 = 0.0
    # How many observations this arm has taken, used to evict its first, cold
    # timing once a warm one exists.
    #
    # The first campaign to run on an arm pays costs that will never recur: JIT
    # of the threaded or distributed dispatch path, thread-pool spin-up, worker
    # provisioning. Recording that timing as if it were the arm's steady-state
    # cost is how a good route gets written off permanently after one look, and
    # it is exactly what happened. Measured on montecarlo_multi_sat, 64 samples,
    # 12 threads, one persisted state across repeated campaigns: the threads
    # arm's first observation came in at 0.590 s against serial's 0.060 s, the
    # selector stopped choosing threads, and every subsequent campaign ran
    # serially -- at 0.060 s where a warm threads run costs 0.024 s and the best
    # width costs 0.017 s. One cold sample cost a factor of 3.4, permanently.
    #
    # Every other timed path in this codebase already discards a warm-up: the
    # RHS calibration sweep warms before it measures, and the benchmark harness
    # takes --warmup=1. This is the same discipline for the route bandit.
    #
    # The cold observation is EVICTED rather than skipped. Skipping it entirely
    # would also discard whether the route succeeded, so a route that failed on
    # its only attempt would be retried forever, and it would leave a campaign
    # that ran with no record of having run. Recording it and then replacing it
    # wholesale on the second observation keeps the success and failure counts
    # from the first attempt available to the selector while ensuring no cold
    # timing survives into the average.
    observations::Int = 0
    # How many CAMPAIGNS this arm has been observed over, as distinct from how
    # many Monte Carlo samples those campaigns contained.
    #
    # The two were conflated, and the exploration gate read the wrong one. One
    # campaign of 64 samples credited the arm with 64 `samples`, which satisfies
    # any reasonable min-samples guarantee immediately -- so a single campaign
    # was enough to mark a route "sufficiently observed" and stop the selector
    # ever trying it again. Those 64 numbers are not 64 observations of the
    # route; they are one, and the route bandit's whole exploration budget was
    # being spent by it.
    #
    # `samples` still normalises the mean, which is a per-sample time and wants
    # the sample count. `campaigns` gates exploration, which wants the number of
    # independent timings.
    campaigns::Int = 0
end

"""
    OuterRouteState

Mutable state used to accumulate historical routing outcomes for adaptive outer
parallel-route selection.
"""
Base.@kwdef mutable struct OuterRouteState
    lock::ReentrantLock = ReentrantLock()
    history::Dict{String, Dict{Symbol, OuterRouteStats}} = Dict{String, Dict{Symbol, OuterRouteStats}}()
end

"""
    reset_outer_route_state!(state)

Clear all accumulated adaptive outer-route history from `state`.
"""
function reset_outer_route_state!(state::OuterRouteState)
    lock(state.lock) do
        empty!(state.history)
    end
    return nothing
end

# `observations` and `campaigns` are persisted, and leaving them out was not a
# cosmetic omission -- it defeated persistence twice over.
#
# Both fields gate behaviour, and both restored as zero:
#
#   campaigns == 0     `_route_is_proven` and `_under_sampled_candidate` read
#                      campaigns, not samples, so every restored arm looked
#                      never-observed and the selector was obliged to spend a
#                      full exploratory campaign on each one -- which is exactly
#                      the cost persistence exists to remove. The docstring on
#                      `ensure_campaign_route_state_loaded!` claims a restored
#                      signature "should exploit immediately, not re-run the
#                      round-robin its predecessor already paid for"; it could
#                      not.
#
#   observations == 0  worse, and silent. The cold-eviction branch in
#                      `record_outer_route_feedback!` fires when observations
#                      reaches 2 and REPLACES the arm's statistics wholesale. So
#                      the second live observation after a load discarded the
#                      entire restored history, however many campaigns it
#                      represented, and left the arm holding one fresh timing.
#
# Test `outer_route_persistence_tests.jl` asserted the intended behaviour and
# was failing on the first of these.
#
# RESIDUAL, deliberately not addressed here: eviction guards a PER-PROCESS cost
# (JIT of the dispatch path, thread-pool spin-up, worker provisioning), but the
# counter is now restored across processes -- so a restored arm with
# observations >= 2 will average this process's own cold first timing into its
# history rather than evicting it. It is diluted by the restored weight and is
# strictly better than either previous behaviour, but the correct fix is a
# separate per-process observation counter that resets on load, which is a
# semantics change rather than a serialisation fix.
@inline function _route_stats_payload(stats::OuterRouteStats)::Dict{String, Any}
    return Dict{String, Any}(
        "samples" => max(0, Int(stats.samples)),
        "successes" => max(0, Int(stats.successes)),
        "failures" => max(0, Int(stats.failures)),
        "elapsed_sum_s" => max(0.0, Float64(stats.elapsed_sum_s)),
        "elapsed_sq_sum_s" => max(0.0, Float64(stats.elapsed_sq_sum_s)),
        "observations" => max(0, Int(stats.observations)),
        "campaigns" => max(0, Int(stats.campaigns))
    )
end

@inline function _route_payload_stats(payload)::Union{Nothing, OuterRouteStats}
    payload isa AbstractDict || return nothing
    samples = try
        max(0, Int(get(payload, "samples", 0)))
    catch
        0
    end
    successes = try
        max(0, Int(get(payload, "successes", 0)))
    catch
        0
    end
    failures = try
        max(0, Int(get(payload, "failures", 0)))
    catch
        0
    end
    elapsed_sum_s = try
        max(0.0, Float64(get(payload, "elapsed_sum_s", 0.0)))
    catch
        0.0
    end
    elapsed_sq_sum_s = try
        max(0.0, Float64(get(payload, "elapsed_sq_sum_s", NaN)))
    catch
        NaN
    end
    samples <= 0 && return nothing
    # Schema-2 files carry neither counter. One is the honest default for both:
    # a row with samples > 0 is evidence the arm ran at least once, and nothing
    # in the file says how many times.
    #
    # One rather than zero for `campaigns`, so a legacy file is not treated as
    # never-observed. One rather than two for `observations`, so the next live
    # observation still evicts a restored timing whose warmth is unknowable --
    # the conservative reading, and the same discipline the eviction applies to
    # a first in-process timing.
    default_count = 1
    observations = try
        max(0, Int(get(payload, "observations", default_count)))
    catch
        default_count
    end
    campaigns = try
        max(0, Int(get(payload, "campaigns", default_count)))
    catch
        default_count
    end
    successes = min(samples, successes)
    failures = min(samples - successes, failures)
    elapsed_sum_s = max(0.0, elapsed_sum_s)
    if !isfinite(elapsed_sq_sum_s)
        # Backward-compatibility path for legacy state files without second-moment data.
        elapsed_sq_sum_s = (elapsed_sum_s^2) / samples
    end
    elapsed_sq_sum_s = max(elapsed_sq_sum_s, (elapsed_sum_s^2) / samples)
    return OuterRouteStats(
        samples=samples,
        successes=successes,
        failures=failures,
        elapsed_sum_s=elapsed_sum_s,
        elapsed_sq_sum_s=elapsed_sq_sum_s,
        observations=observations,
        campaigns=campaigns
    )
end

function save_outer_route_state(
    state::OuterRouteState,
    path::AbstractString;
    metadata::AbstractDict=Dict{String, Any}()
)::NamedTuple{(:path, :signatures, :rows), Tuple{String, Int, Int}}
    path_s = String(path)
    rows = Dict{String, Any}[]
    signatures = 0
    lock(state.lock) do
        for signature in sort!(collect(keys(state.history)))
            bucket = state.history[signature]
            isempty(bucket) && continue
            signature_rows = 0
            # Every arm in the bucket, not a fixed three. The split selector adds
            # `split_<route>_w<N>` arms alongside the route arms, and a fixed
            # enumeration would silently drop them on save -- so a restored state
            # would carry route history but no split history, and the split
            # bandit would re-explore from cold on every run while the route
            # bandit did not.
            for route in sort!(collect(keys(bucket)); by = String)
                stats = get(bucket, route, nothing)
                stats isa OuterRouteStats || continue
                stats.samples > 0 || continue
                push!(rows, Dict{String, Any}(
                    "signature" => signature,
                    "route" => String(route),
                    "stats" => _route_stats_payload(stats)
                ))
                signature_rows += 1
            end
            signature_rows > 0 && (signatures += 1)
        end
    end

    metadata_out = Dict{String, Any}()
    for (k, v) in metadata
        key = String(k)
        metadata_out[key] = if v isa Number || v isa Bool
            v
        else
            String(v)
        end
    end

    payload = Dict{String, Any}(
        # 3: adds the `observations` and `campaigns` counters. Not gated on when
        # loading -- a schema-2 file is still read, with the defaults documented
        # in `_route_payload_stats` -- so the bump records what a file contains
        # rather than rejecting older ones.
        "schema_version" => 3,
        "updated_utc" => string(now(UTC)),
        "metadata" => metadata_out,
        "history" => rows
    )

    mkpath(dirname(path_s))
    open(path_s, "w") do io
        TOML.print(io, payload)
    end
    return (path=path_s, signatures=signatures, rows=length(rows))
end

function load_outer_route_state!(
    state::OuterRouteState,
    path::AbstractString;
    replace::Bool=true
)::NamedTuple{(:path, :signatures, :rows), Tuple{String, Int, Int}}
    path_s = String(path)
    isfile(path_s) || return (path=path_s, signatures=0, rows=0)
    parsed = TOML.parsefile(path_s)
    rows_in = get(parsed, "history", Any[])
    rows_in isa AbstractVector || return (path=path_s, signatures=0, rows=0)

    loaded_rows = 0
    loaded_signatures = Set{String}()
    lock(state.lock) do
        replace && empty!(state.history)
        for row in rows_in
            row isa AbstractDict || continue
            signature = strip(String(get(row, "signature", "")))
            isempty(signature) && continue
            route_name = strip(String(get(row, "route", "")))
            isempty(route_name) && continue
            route = Symbol(route_name)
            stats = _route_payload_stats(get(row, "stats", nothing))
            stats === nothing && continue

            bucket = get!(state.history, signature) do
                Dict{Symbol, OuterRouteStats}()
            end
            existing = get!(bucket, route) do
                OuterRouteStats()
            end
            existing.samples += stats.samples
            existing.successes += stats.successes
            existing.failures += stats.failures
            existing.elapsed_sum_s += stats.elapsed_sum_s
            existing.elapsed_sq_sum_s += stats.elapsed_sq_sum_s
            # Additive like every other field: `replace=false` merges a session's
            # own observations with the persisted ones, and an observation count
            # that did not merge would under-report the merged history it now
            # describes.
            existing.observations += stats.observations
            existing.campaigns += stats.campaigns
            loaded_rows += 1
            push!(loaded_signatures, signature)
        end
    end
    return (path=path_s, signatures=length(loaded_signatures), rows=loaded_rows)
end
