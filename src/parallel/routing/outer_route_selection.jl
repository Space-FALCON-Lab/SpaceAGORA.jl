@inline function _threads_or_none(threads_available::Bool)::Symbol
    return threads_available ? :threads : :none
end

@inline function _route_sat_bucket(n_sats::Int)::String
    if n_sats <= 1
        return "1"
    elseif n_sats <= 2
        return "2"
    elseif n_sats <= 4
        return "3_4"
    end
    return "5p"
end

@inline function _route_link_bucket(n_links::Int)::String
    if n_links <= 1
        return "1"
    elseif n_links <= 4
        return "2_4"
    elseif n_links <= 8
        return "5_8"
    end
    return "9p"
end

@inline function _route_max_link_bucket(n_links::Int)::String
    if n_links <= 1
        return "1"
    elseif n_links <= 2
        return "2"
    elseif n_links <= 4
        return "3_4"
    elseif n_links <= 8
        return "5_8"
    end
    return "9p"
end

@inline function _route_mission_bucket(mission_time_s::Float64)::String
    if mission_time_s <= 1800.0
        return "short"
    elseif mission_time_s <= 7200.0
        return "medium"
    end
    return "long"
end

@inline function _route_harmonics_bucket(L::Int)::String
    if L <= 0
        return "0"
    elseif L <= 10
        return "1_10"
    elseif L <= 20
        return "11_20"
    end
    return "21p"
end

@inline function _route_density_bucket(family::AbstractString)::String
    token = lowercase(strip(String(family)))
    if token in ("none", "vacuum", "noatmosphere")
        return "none"
    elseif token in ("gram_point", "gram")
        return "gram_pt"
    elseif token in ("gram_surrogate", "gram_offline_surrogate")
        return "gram_srg"
    elseif token in ("exponential", "exp")
        return "exp"
    elseif token in ("polynomialfit", "polyfit", "poly")
        return "poly"
    elseif token in ("nrlmsise00", "nrl")
        return "nrl"
    elseif isempty(token)
        return "unknown"
    end
    return replace(token, "|" => "_")
end

@inline function _route_solver_bucket(mode::AbstractString)::String
    token = lowercase(strip(String(mode)))
    if token in ("auto", "auto_stiff", "autostiff")
        return "auto"
    elseif token in ("rodas5p", "rodas")
        return "rodas"
    elseif token in ("tsit5",)
        return "tsit5"
    elseif startswith(token, "split_imex")
        return "split"
    elseif token in ("multirate",)
        return "mrate"
    elseif isempty(token)
        return "auto"
    end
    return replace(token, "|" => "_")
end

@inline function _route_interval_bucket(interval_s::Float64)::String
    if !isfinite(interval_s) || interval_s <= 0.0
        return "na"
    elseif interval_s <= 0.2
        return "ultra"
    elseif interval_s <= 1.0
        return "fast"
    elseif interval_s <= 10.0
        return "med"
    end
    return "slow"
end

@inline function _route_count_bucket(v::Int)::String
    if v <= 0
        return "0"
    elseif v == 1
        return "1"
    elseif v <= 3
        return "2_3"
    elseif v <= 6
        return "4_6"
    end
    return "7p"
end

@inline function _route_effector_cost_bucket(raw::AbstractString)::String
    token = lowercase(strip(String(raw)))
    if token in ("light", "medium", "heavy")
        return token
    elseif isempty(token)
        return "unknown"
    end
    return replace(token, "|" => "_")
end

@inline function _compat_outer_route_signature(f::OuterRouteFeatures)::String
    return join((
        "cat=$(f.category)",
        "sat=$(_route_sat_bucket(f.n_sats))",
        "links=$(_route_link_bucket(f.n_links))",
        "mission=$(_route_mission_bucket(f.mission_time_s))",
        "nbody=$(f.has_nbody ? "1" : "0")",
        "srp=$(f.has_srp ? "1" : "0")",
        "harm=$(_route_harmonics_bucket(f.harmonics_degree))",
        "ctrl=$(f.has_control ? "1" : "0")",
        "orient=$(f.orientation_on ? "1" : "0")"
    ), "|")
end

"""
    outer_route_signature(features) -> String

Return the canonical routing signature used to bucket simulation workloads for
outer-route policy and feedback.
"""
@inline function outer_route_signature(f::OuterRouteFeatures)::String
    return join((
        "cat=$(f.category)",
        "sat=$(_route_sat_bucket(f.n_sats))",
        "links=$(_route_link_bucket(f.n_links))",
        "maxlinks=$(_route_max_link_bucket(f.max_links_per_sat))",
        "mission=$(_route_mission_bucket(f.mission_time_s))",
        "nbody=$(f.has_nbody ? "1" : "0")",
        "srp=$(f.has_srp ? "1" : "0")",
        "harm=$(_route_harmonics_bucket(f.harmonics_degree))",
        "ctrl=$(f.has_control ? "1" : "0")",
        "orient=$(f.orientation_on ? "1" : "0")",
        "dens=$(_route_density_bucket(f.density_family))",
        "solver=$(_route_solver_bucket(f.solver_mode))",
        "dt=$(_route_interval_bucket(f.dt_max_orbit_s))",
        "ctrl_rate=$(_route_interval_bucket(f.control_rate_s))",
        "guid_rate=$(_route_interval_bucket(f.guidance_rate_s))",
        "nav_rate=$(_route_interval_bucket(f.navigation_rate_s))",
        "gram_srg=$(f.gram_surrogate_enabled ? "1" : "0")",
        "gram_grid=$(f.gram_static_grid_enabled ? "1" : "0")",
        "ctrl_eff=$(_route_count_bucket(f.control_effector_count))",
        "thermal=$(f.thermal_enabled ? "1" : "0")",
        "eff_cnt=$(_route_count_bucket(f.dynamic_effector_count))",
        "eff_cost=$(_route_effector_cost_bucket(f.effector_cost_class))"
    ), "|")
end

@inline function _outer_route_signature_hierarchy(f::OuterRouteFeatures)::Vector{String}
    full = outer_route_signature(f)
    mid = join((
        "cat=$(f.category)",
        "sat=$(_route_sat_bucket(f.n_sats))",
        "links=$(_route_link_bucket(f.n_links))",
        "maxlinks=$(_route_max_link_bucket(f.max_links_per_sat))",
        "mission=$(_route_mission_bucket(f.mission_time_s))",
        "nbody=$(f.has_nbody ? "1" : "0")",
        "srp=$(f.has_srp ? "1" : "0")",
        "harm=$(_route_harmonics_bucket(f.harmonics_degree))",
        "ctrl=$(f.has_control ? "1" : "0")",
        "orient=$(f.orientation_on ? "1" : "0")",
        "dens=$(_route_density_bucket(f.density_family))",
        "solver=$(_route_solver_bucket(f.solver_mode))",
        "thermal=$(f.thermal_enabled ? "1" : "0")",
        "eff_cost=$(_route_effector_cost_bucket(f.effector_cost_class))"
    ), "|")
    legacy = _compat_outer_route_signature(f)
    return unique(String[full, mid, legacy])
end

@inline function _route_elapsed_stats(
    stats::OuterRouteStats
)::NamedTuple{(:mean_s, :std_s), Tuple{Float64, Float64}}
    samples = max(1, stats.samples)
    mean_s = stats.elapsed_sum_s / samples
    mean_sq_s = stats.elapsed_sq_sum_s / samples
    variance_s = max(0.0, mean_sq_s - mean_s^2)
    return (mean_s=mean_s, std_s=sqrt(variance_s))
end

function _outer_route_stats_snapshot_internal(
    state::OuterRouteState,
    signature::String
)::Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate, :std_s, :campaigns), Tuple{Int, Float64, Float64, Float64, Int}}}
    lock(state.lock) do
        entry = get(state.history, signature, nothing)
        snap = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate, :std_s, :campaigns), Tuple{Int, Float64, Float64, Float64, Int}}}()
        if entry === nothing
            return snap
        end
        for (route, stats) in entry
            if stats.samples <= 0
                continue
            end
            success_rate = stats.successes / max(1, stats.samples)
            elapsed = _route_elapsed_stats(stats)
            snap[route] = (
                samples=stats.samples,
                mean_s=elapsed.mean_s,
                success_rate=success_rate,
                std_s=elapsed.std_s,
                campaigns=stats.campaigns
            )
        end
        return snap
    end
end

"""
    outer_route_stats_snapshot(state, signature) -> Dict

Return the aggregated routing statistics currently stored for a workload
signature.
"""
function outer_route_stats_snapshot(
    state::OuterRouteState,
    signature::String
)::Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}}
    base = _outer_route_stats_snapshot_internal(state, signature)
    snap = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), Tuple{Int, Float64, Float64}}}()
    for (route, info) in base
        snap[route] = (samples=info.samples, mean_s=info.mean_s, success_rate=info.success_rate)
    end
    return snap
end

@inline function _feature_is_lightweight(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool
    return f.n_sats <= t.outer_light_sat_threshold &&
        f.n_links <= t.outer_light_link_threshold &&
        f.mission_time_s <= t.outer_light_mission_threshold_s &&
        !f.has_nbody &&
        !f.has_control &&
        !f.orientation_on &&
        f.harmonics_degree == 0
end

@inline function _is_native_gram_point_density(f::OuterRouteFeatures)::Bool
    return _route_density_bucket(f.density_family) == "gram_pt" &&
        !f.gram_surrogate_enabled &&
        !f.gram_static_grid_enabled
end

@inline function _feature_heavy_for_process(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool
    if _is_native_gram_point_density(f)
        # Native GRAM point calls are lock-limited; prefer outer process isolation.
        return true
    end
    if _feature_is_lightweight(f, t)
        return false
    end
    return f.has_nbody || f.harmonics_degree >= 20 || f.mission_time_s > t.outer_light_mission_threshold_s
end

# Whether a Monte Carlo campaign is big enough that the process route is worth
# EXPLORING. Deliberately not the same question as whether it should be the
# default -- see _priority_outer_route_montecarlo.
@inline function _mc_process_worth_exploring(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool
    return f.montecarlo_samples >= t.mc_process_min_samples ||
        f.mission_time_s >= t.mc_process_min_mission_s
end

@inline function _priority_outer_route_montecarlo(
    f::OuterRouteFeatures,
    t::OuterRouteTuning;
    machine_class::Symbol,
    threads_available::Bool,
    parallel_enabled::Bool
)::Symbol
    if !parallel_enabled
        return :none
    end
    if f.montecarlo_samples <= 1
        return :none
    end
    # Threads, not process, and this is a measured reversal.
    #
    # This returned :process whenever the machine was medium/large AND either the
    # sample count cleared mc_process_min_samples (16) or the mission cleared
    # mc_process_min_mission_s (3600 s). Neither term is a proxy for PER-SAMPLE
    # COMPUTE, which is the only thing that decides whether a sample is worth
    # shipping to another process:
    #
    #   - Sample count measures total work AND total dispatch overhead. It rises
    #     on both sides of the comparison, so it cannot discriminate.
    #   - mission_time_s is SIMULATED seconds. A one-satellite, one-hour arc is
    #     3600 by this measure and 72 ms of actual compute.
    #
    # The result was that essentially every Monte Carlo campaign of 16+ samples
    # took the process route, and the process route lost every time. Measured at
    # 64 samples, 12 threads, median of 3 post-warm-up repeats, process against
    # threads:
    #
    #     montecarlo_multi_sat          0.038 s/sample   2.45x SLOWER
    #     independent_1sat_1hr          0.072 s/sample   2.90x SLOWER
    #     montecarlo_high_accuracy      0.072 s/sample   1.46x SLOWER
    #     montecarlo_mars_aerobraking   0.259 s/sample   1.41x SLOWER
    #     montecarlo_heavy_aerobraking  3.179 s/sample   1.69x SLOWER
    #
    # There is no crossover in that range. Even at 3.18 s of compute per sample
    # -- two orders of magnitude above the cheapest case -- Distributed's
    # per-campaign worker startup and per-sample serialisation still cost more
    # than they save. So there is no threshold on per-sample cost that would have
    # rescued the old rule either; the rule's premise was wrong, not its constants.
    #
    # Against the best static route this was worth +33% to +201% on every Monte
    # Carlo case in the paper's catalog, and the inner-axis ablations
    # (full_smart_nocalib / _noinner / _nopolicy / _innermodes_off / _nohints)
    # recovered none of it, which is what localised the cost here.
    #
    # Native GRAM is unaffected: _is_native_gram_point_density short-circuits to
    # the process route in both default_outer_route and outer_route_candidates
    # BEFORE this function is consulted, because there the process route is a
    # thread-safety requirement rather than a performance choice.
    #
    # Process stays in the CANDIDATE set (outer_route_candidates), so the bandit
    # can still discover it on a machine or workload where it wins -- it simply
    # has to earn the choice from measured feedback instead of being assumed.
    # This mirrors the no-regret floor on the inner axis: the heuristic's answer
    # is the default, and calibration may displace it only on evidence.
    #
    # WITH NO WORKER THREADS the comparison above does not apply and the answer
    # inverts. Everything measured against the process route measured it against
    # THREADS, and threads won every time. With one Julia thread there are no
    # threads to lose to: the comparison is process against serial, and a Monte
    # Carlo campaign is a set of independent samples, so the process pool is the
    # only parallelism available at all.
    #
    # Falling through to _threads_or_none returned :none here, which is fully
    # serial, and it was the single worst result in the benchmark set. At one
    # thread, 64 samples, the fastest fixed route is outer processes (R1_b) on
    # every Monte Carlo case measured, and both adaptive profiles trailed it by
    # most of a factor of two:
    #
    #     montecarlo_multi_sat        R1_b 0.316 s   R4 +89.0%   R5 +67.0%
    #     montecarlo_high_accuracy    R1_b 0.572 s   R4 +73.4%   R5 +77.5%
    #     montecarlo_mars_aerobraking R1_b 2.115 s   R4 +88.8%   R5 +86.6%
    #
    # The affordability and machine-class conditions are the same ones
    # outer_route_candidates uses to enumerate the process arm, so the default
    # can never name a route the selector did not offer.
    process_affordable = machine_class in (:large, :medium) && _mc_process_worth_exploring(f, t)
    if !threads_available
        return process_affordable ? :process : :none
    end
    # V2: give the outer axis to whichever route has more cores.
    #
    # The measurement behind the threads default above was taken with the
    # Distributed pool provisioned inside the timed window. The pool is now
    # persistent and warmed before the clock, and the 2026-09-02 paper run at
    # commit 3be5b0ec shows the opposite sign wherever the process route has at
    # least as many workers as the thread route has threads -- and only there:
    #
    #     montecarlo_heavy_aerobraking, 64 samples, 12 physical cores
    #     (threads, workers)   pinned threads   pinned process
    #     (12, 1)                  10.2 s          59.2 s
    #     ( 6, 2)                  13.9 s          31.5 s
    #     ( 4, 3)                  19.3 s          20.8 s
    #     ( 3, 4)                  24.1 s          16.2 s
    #     ( 2, 6)                  34.7 s          11.9 s
    #     ( 1,12)                  57.5 s           7.4 s
    #     independent_1sat_1hr, 256 samples
    #     (12,12)                   1.65 s          1.01 s
    #     ( 8,12)                   1.88 s          1.03 s
    #
    # Comparing worker count against thread count reproduces every row. The
    # shipped default chose threads on all of them and lost by up to 2.9x on
    # the wide-pool points; exploration could not rescue it because the guard
    # in _route_is_proven stops once threads beats serial, before process is
    # ever tried (see must_measure).
    if t.mc_route_by_core_budget && process_affordable &&
       t.process_max_workers >= max(1, t.outer_thread_budget)
        return :process
    end
    return :threads
end

"""
    default_outer_route(features; kwargs...) -> Symbol

Choose the default outer parallel route for a workload before adaptive
historical feedback is considered.
"""
function default_outer_route(
    f::OuterRouteFeatures;
    tuning::OuterRouteTuning=OuterRouteTuning(),
    machine_class::Symbol=:small,
    threads_available::Bool=true,
    parallel_enabled::Bool=true
)::Symbol
    if !parallel_enabled
        return :none
    end

    if _is_native_gram_point_density(f)
        return :process
    end

    if lowercase(strip(f.category)) == "montecarlo"
        return _priority_outer_route_montecarlo(
            f,
            tuning;
            machine_class=machine_class,
            threads_available=threads_available,
            parallel_enabled=parallel_enabled
        )
    end

    if _feature_is_lightweight(f, tuning)
        return :none
    end

    if f.n_sats >= tuning.inner_sat_threshold || f.n_links >= tuning.inner_link_threshold
        return :none
    end

    if tuning.spice_constellation_process_enabled &&
       f.n_sats >= tuning.spice_constellation_min_sats &&
       (f.has_nbody || f.has_srp)
        return machine_class in (:large, :medium) ? :process : _threads_or_none(threads_available)
    end

    if lowercase(strip(f.category)) == "satellite_scaling" && f.n_sats >= 4
        return _threads_or_none(threads_available)
    end

    if f.has_nbody || f.harmonics_degree >= 20
        return machine_class in (:large, :medium) ? :process : _threads_or_none(threads_available)
    end

    if machine_class in (:large, :medium)
        return :process
    end
    return _threads_or_none(threads_available)
end

"""
    outer_route_candidates(features; kwargs...) -> Vector{Symbol}

Return the feasible outer parallel routes for a workload under the current
tuning and machine assumptions.
"""
function outer_route_candidates(
    f::OuterRouteFeatures;
    tuning::OuterRouteTuning=OuterRouteTuning(),
    machine_class::Symbol=:small,
    threads_available::Bool=true,
    parallel_enabled::Bool=true
)::Vector{Symbol}
    if !parallel_enabled
        return Symbol[:none]
    end
    if _is_native_gram_point_density(f)
        return Symbol[:none, :process]
    end

    candidates = Symbol[:none]
    if threads_available
        push!(candidates, :threads)
    end
    # Candidacy is a weaker test than default-ness, and for Monte Carlo the two
    # are now different questions. _priority_outer_route_montecarlo no longer
    # returns :process (see its comment), so deriving candidacy from it would
    # drop the process route out of the bandit's arm set entirely and make the
    # reversal permanent and undiscoverable. Enumerate it here whenever the
    # campaign is large enough for the exploration to be affordable, and let
    # measured feedback decide.
    allow_process = if lowercase(strip(f.category)) == "montecarlo"
        f.montecarlo_samples > 1 &&
            machine_class in (:large, :medium) &&
            _mc_process_worth_exploring(f, tuning)
    else
        _feature_heavy_for_process(f, tuning)
    end
    if allow_process
        push!(candidates, :process)
    end
    return unique(candidates)
end

@inline function _route_ranked_candidates(candidates::Vector{Symbol}, default_route::Symbol)::Vector{Symbol}
    ranked = Symbol[]
    if default_route in candidates
        push!(ranked, default_route)
    end
    for route in (:threads, :none, :process)
        if route in candidates && !(route in ranked)
            push!(ranked, route)
        end
    end
    for route in candidates
        if !(route in ranked)
            push!(ranked, route)
        end
    end
    return ranked
end

@inline function _under_sampled_candidate(
    candidates::Vector{Symbol},
    snapshot,
    default_route::Symbol,
    min_samples::Int
)::Union{Nothing, Symbol}
    ranked = _route_ranked_candidates(candidates, default_route)
    for route in ranked
        info = get(snapshot, route, (samples=0, mean_s=Inf, success_rate=0.0, campaigns=0))
        # Campaigns, not samples: see OuterRouteStats.campaigns.
        if get(info, :campaigns, info.samples) < max(1, min_samples)
            return route
        end
    end
    return nothing
end

# Whether `route` has been tested AGAINST SOMETHING and won.
#
# Three conditions, and the third is the one that matters. The route needs
# enough observations to be trusted, no sampled alternative may beat it, and at
# least one alternative must actually have been sampled.
#
# Without that last clause the predicate is vacuous exactly when it is most
# dangerous: after a single campaign only the default has data, "no other route
# beats it" is trivially true because no other route has been tried, and the
# selector would stop exploring before it had learned anything. Comparing a
# measured route against nothing is not evidence that it is the right one -- it
# only means nothing has contradicted it yet.
#
# What this does buy, which is the point, is that once alternatives have been
# tried and lost the router stops paying to re-try them. That is the recurring
# cost persisted history is meant to eliminate: a restored signature carrying
# hundreds of observations across every route should exploit immediately, not
# re-run the round-robin its predecessor already paid for.
@inline function _route_is_proven(
    candidates::Vector{Symbol},
    snapshot,
    route::Symbol,
    min_samples::Int;
    # Arms that must have been measured before ANY arm counts as proven,
    # whatever proven_requires_all_candidates says. The route selector under
    # V2 names the parallel alternatives here: threads beating serial is no
    # evidence about process, and process was the arm never reached.
    must_measure::Vector{Symbol}=Symbol[],
)::Bool
    info = get(snapshot, route, nothing)
    info === nothing && return false
    # Campaigns, not samples: one campaign of N samples is one observation of
    # the route, and gating on N declared a route proven after a single, cold
    # look. See OuterRouteStats.campaigns.
    (get(info, :campaigns, info.samples) >= max(1, min_samples) && isfinite(info.mean_s)) || return false
    # An arm cannot be proven while an enumerated alternative has never been
    # measured at all.
    #
    # This used to require only that SOME alternative had been tried, and on a
    # width ladder that is satisfied by the worst possible one. The split
    # selector enumerates [w1, w2, w4, w8, w12]; once w12 and w1 both had a
    # timing and w12 was faster, the default was declared proven and w2, w4 and
    # w8 were never tried. Beating width 1 says nothing about width 8, and width
    # 8 is where the optimum actually sits: measured directly at 12 threads with
    # 64 samples, the full budget is 47% slower than width 8 on
    # montecarlo_multi_sat and 55% slower on montecarlo_high_accuracy.
    #
    # Requiring every enumerated candidate to have been measured costs a few
    # exploratory campaigns once per signature, against a standing error of that
    # size on every campaign afterwards.
    unmeasured_alternative = false
    tested_against_an_alternative = false
    # Only arms this selector enumerated. Scanning the whole snapshot was wrong
    # once the bucket held more than one selector's arms: after a single
    # campaign the split arm `split_threads_w4` sits beside the `:threads` route
    # arm carrying the same observation, so it registered as "an alternative
    # that was measured and did not win" and declared the route proven -- ending
    # exploration after one campaign, on evidence that was the route's own
    # timing recorded twice under two names.
    for other in candidates
        other === route && continue
        other_info = get(snapshot, other, nothing)
        if other_info === nothing || other_info.samples <= 0 || !isfinite(other_info.mean_s)
            other in must_measure && return false
            unmeasured_alternative = true
            continue
        end
        tested_against_an_alternative = true
        other_info.mean_s < info.mean_s && return false
    end
    if unmeasured_alternative && proven_requires_all_candidates()
        return false
    end
    return tested_against_an_alternative
end

# Whether ANY enumerated candidate is proven best -- the guard the selectors
# should be asking, where `_route_is_proven(default)` asks a narrower one.
#
# `_route_is_proven` returns false for a default as soon as a sampled
# alternative beats it (its third clause, deliberately). Gating forced
# exploration on that alone means a BEATEN default re-enables exploration, and
# `_under_sampled_candidate` then diverts the campaign to whichever ranked
# candidate has too few observations -- while the better-measured winner sits
# unconsulted. Measured shape: process at 0.90 s, threads at 0.10 s, `:none`
# never tried; the shipped selector runs `:none`.
#
# Proven-best for some candidate is the same predicate applied to that
# candidate: enough campaigns, nothing sampled beats it, at least one
# alternative tried. The exploration budget is spent only while no arm has
# earned the right to be exploited.
@inline function _any_candidate_proven(
    candidates::Vector{Symbol},
    snapshot,
    min_samples::Int;
    must_measure::Vector{Symbol}=Symbol[],
)::Bool
    for route in candidates
        _route_is_proven(candidates, snapshot, route, min_samples;
                         must_measure=must_measure) && return true
    end
    return false
end

@inline function _best_candidate(
    candidates::Vector{Symbol},
    snapshot
)::Union{Nothing, Symbol}
    best_route = nothing
    best_mean = Inf
    best_success_rate = -Inf
    for route in candidates
        info = get(snapshot, route, (samples=0, mean_s=Inf, success_rate=0.0))
        if info.samples <= 0 || !isfinite(info.mean_s)
            continue
        end
        if info.mean_s < best_mean - 1e-12
            best_route = route
            best_mean = info.mean_s
            best_success_rate = info.success_rate
        elseif isapprox(info.mean_s, best_mean; atol=1e-12, rtol=0.0) && info.success_rate > best_success_rate
            best_route = route
            best_success_rate = info.success_rate
        end
    end
    return best_route
end

@inline function _candidate_confidence_width(
    std_s::Float64,
    samples::Int,
    total_samples::Int,
    exploration_c::Float64
)::Float64
    samples <= 0 && return 0.0
    exploration = max(0.0, exploration_c)
    scaled = std_s * sqrt(log(max(2.0, Float64(total_samples) + 1.0)) / max(1.0, Float64(samples)))
    width = exploration * scaled
    return isfinite(width) ? max(0.0, width) : 0.0
end

@inline function _best_candidate_confidence(
    candidates::Vector{Symbol},
    snapshot,
    default_route::Symbol,
    exploration_c::Float64
)::NamedTuple{(:route, :confidence_s, :regret_s), Tuple{Union{Nothing, Symbol}, Float64, Float64}}
    total_samples = 0
    for route in candidates
        info = get(snapshot, route, nothing)
        info === nothing && continue
        total_samples += max(0, info.samples)
    end
    total_samples = max(1, total_samples)

    best_route = nothing
    best_score = Inf
    best_mean = Inf
    best_width = 0.0
    best_success_rate = -Inf
    best_observed_mean = Inf

    ranked = _route_ranked_candidates(candidates, default_route)
    for route in ranked
        info = get(snapshot, route, nothing)
        info === nothing && continue
        if info.samples <= 0 || !isfinite(info.mean_s)
            continue
        end
        best_observed_mean = min(best_observed_mean, info.mean_s)
        width = _candidate_confidence_width(info.std_s, info.samples, total_samples, exploration_c)
        score = info.mean_s - width
        if score < best_score - 1e-12
            best_route = route
            best_score = score
            best_mean = info.mean_s
            best_width = width
            best_success_rate = info.success_rate
        elseif isapprox(score, best_score; atol=1e-12, rtol=0.0) && info.success_rate > best_success_rate + 1e-12
            best_route = route
            best_mean = info.mean_s
            best_width = width
            best_success_rate = info.success_rate
        end
    end

    if best_route === nothing
        return (route=nothing, confidence_s=0.0, regret_s=0.0)
    end
    regret_s = isfinite(best_observed_mean) ? max(0.0, best_mean - best_observed_mean) : 0.0
    return (route=best_route, confidence_s=best_width, regret_s=regret_s)
end

"""
    select_outer_route!(state, features; kwargs...) -> Symbol

Select the outer parallel route for a workload, using adaptive history from
`state` when enabled.
"""
function select_outer_route!(
    state::OuterRouteState,
    f::OuterRouteFeatures;
    tuning::OuterRouteTuning=OuterRouteTuning(),
    machine_class::Symbol=:small,
    threads_available::Bool=true,
    parallel_enabled::Bool=true
)::Symbol
    default_route = default_outer_route(
        f;
        tuning=tuning,
        machine_class=machine_class,
        threads_available=threads_available,
        parallel_enabled=parallel_enabled
    )
    if !tuning.adaptive_enabled
        return default_route
    end

    candidates = outer_route_candidates(
        f;
        tuning=tuning,
        machine_class=machine_class,
        threads_available=threads_available,
        parallel_enabled=parallel_enabled
    )
    isempty(candidates) && return :none
    if !(default_route in candidates)
        default_route = first(candidates)
    end

    signature_chain = _outer_route_signature_hierarchy(f)
    signature = first(signature_chain)
    snapshot = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate, :std_s), Tuple{Int, Float64, Float64, Float64}}}()
    signature_used = signature
    for candidate_sig in signature_chain
        snap = _outer_route_stats_snapshot_internal(state, candidate_sig)
        if !isempty(snap)
            snapshot = snap
            signature_used = candidate_sig
            break
        end
    end
    chosen = default_route
    reason = "default"
    confidence_s = 0.0
    regret_s = 0.0
    if !isempty(snapshot)
        # Forced exploration is skipped once the default route has proven itself.
        #
        # `_under_sampled_candidate` guarantees every candidate `adaptive_min_samples`
        # observations before the selector will exploit anything. That is the right
        # behaviour cold, but it keeps costing after the answer is known: a
        # candidate set whose membership shifts between runs -- a different
        # machine_class, a changed threshold, or history restored for only some
        # routes -- re-triggers a guaranteed trial of a route the default already
        # beats. With persisted history this became a regression rather than a
        # cost: a signature carrying hundreds of observations of `:process` at
        # 0.1 s and `:threads` at 0.9 s would still divert to an unsampled
        # `:none`, which cold-start would never have chosen because
        # `default_outer_route` answers `:process` outright.
        #
        # So exploration stays mandatory only while the default is itself
        # unproven. Once it has the minimum samples and no sampled alternative
        # beats it, ranking falls through to UCB, which still explores -- through
        # the confidence width on candidates it has data for -- but is no longer
        # obliged to spend a full trial on every candidate it does not.
        # Under explore_until_any_proven the question is whether ANY arm has
        # earned exploitation, not whether the default has; see
        # _any_candidate_proven for the failure the narrower guard produces.
        default_proven = tuning.explore_until_any_proven ?
            _any_candidate_proven(candidates, snapshot, tuning.adaptive_min_samples;
                                  must_measure=Symbol[c for c in candidates if c !== :none]) :
            _route_is_proven(candidates, snapshot, default_route, tuning.adaptive_min_samples)
        explore = default_proven ?
            nothing :
            _under_sampled_candidate(candidates, snapshot, default_route, tuning.adaptive_min_samples)
        if !(explore === nothing)
            chosen = explore
            reason = "explore_hier"
        else
            best = _best_candidate_confidence(
                candidates,
                snapshot,
                default_route,
                tuning.adaptive_exploration_c
            )
            if !(best.route === nothing)
                chosen = best.route
                confidence_s = best.confidence_s
                regret_s = best.regret_s
                reason = chosen == default_route ? "default_ucb_hier" : "exploit_ucb_hier"
            end
        end
    end
    if tuning.trace
        println(
            "[outer-route] signature=$(signature) default=$(default_route) chosen=$(chosen) " *
            "reason=$(reason) signature_used=$(signature_used) candidates=$(join(string.(candidates), ',')) " *
            "confidence_s=$(round(confidence_s; digits=6)) regret_s=$(round(regret_s; digits=6))"
        )
    end
    return chosen
end

# ── Outer/inner split adaptation ──────────────────────────────────────────────
#
# How many outer workers to use, and therefore how many threads each one gets
# for inner parallelism. Until now this was arithmetic rather than a decision:
# the campaign runner took the widest outer split the route allowed and gave the
# inner budget whatever floor division left over, which for a threads route
# using the whole pool is exactly one thread -- inner parallelism disabled, not
# because that was measured to be right but because nothing chose otherwise.
#
# The two levels trade against each other. Widening the outer split shortens the
# queue but starves each sample of inner threads; narrowing it does the reverse.
# Which side wins depends on whether a single sample can use several threads
# productively, which is a property of the workload and not of the machine, so
# it cannot be settled by a fixed rule.
#
# The selector reuses the route bandit wholesale -- same per-signature bucket,
# same statistics, same confidence handling, same persistence -- by naming its
# arms `split_<route>_w<N>`. Route arms and split arms live side by side without
# colliding because each selector only ever scores the arms it enumerated.

@inline function _split_arm(route::Symbol, workers::Int)::Symbol
    return Symbol("split_", String(route), "_w", max(1, workers))
end

# Split arms live under their own signature namespace rather than beside the
# route arms.
#
# Sharing a bucket looked economical and was not. Anything that scans a bucket
# rather than an explicit candidate list sees both selectors' arms: _route_is_proven
# read one campaign's split arm as an independent alternative and declared the
# route proven after a single observation, and summing a snapshot's samples
# double-counts every campaign, because the split arm carries the same
# observation under a second name. Separating the namespace removes the whole
# class of interaction instead of patching each scan site as it is found.
const _SPLIT_SIGNATURE_PREFIX = "split|"

@inline function _split_signature_chain(f::OuterRouteFeatures)::Vector{String}
    return String[_SPLIT_SIGNATURE_PREFIX * sig for sig in _outer_route_signature_hierarchy(f)]
end

"""
    outer_split_candidates(route; budget, n_units, tuning) -> Vector{Int}

Worker counts worth trying for `route`, as a geometric ladder up to the widest
split the route and workload allow.

Geometric rather than exhaustive because the interesting differences are
multiplicative: 1, 2, 4, 8 covers the shape of the trade-off, while every
integer in between mostly re-measures the same regime at the cost of more
exploration. The widest split is always included, since it is what the previous
fixed rule would have chosen and the selector must be able to reproduce it.
"""
function outer_split_candidates(
    route::Symbol;
    budget::Int,
    n_units::Int,
    tuning::OuterRouteTuning=OuterRouteTuning()
)::Vector{Int}
    units = max(1, n_units)
    max_workers = if route === :process
        max(1, min(units, tuning.process_max_workers))
    elseif route === :threads
        max(1, min(units, max(1, budget)))
    else
        1
    end
    max_workers <= 1 && return Int[1]

    # The ladder starts at a quarter of the available width, not at 1.
    #
    # Every arm is an arm the selector must pay to measure, and on a Monte Carlo
    # campaign the very narrow ones cannot win: they run the same per-sample work
    # with less of the machine. Worse, they are actively harmful as candidates.
    # Ranking five arms from short, noisy campaigns is unreliable, and a narrow
    # arm that draws a lucky timing can be selected and held; measured on
    # montecarlo_mars_aerobraking, whose optimum IS the full width, opening the
    # whole ladder left the selector in a steady state of 0.057 s against 0.043 s
    # for the previous two-arm behaviour -- worse in the steady state, not merely
    # slower to converge.
    #
    # Declining to split at all is not lost by this: it is the `:none` route,
    # which the route selector enumerates separately and which this ladder would
    # only duplicate at w=1.
    # Filter the doubling grid rather than re-basing it: re-basing moves every
    # rung (a 12-wide budget became 3, 6, 12) and can drop the optimum, which
    # here is 8.
    floor_w = max(2, cld(max_workers, 4))
    out = Int[]
    w = 1
    while w < max_workers
        w >= floor_w && push!(out, w)
        w *= 2
    end
    push!(out, max_workers)
    return sort!(unique!(out))
end

"""
    select_outer_split!(state, features; route, budget, n_units, tuning) -> Int

Choose how many outer workers to run, given an already-selected `route`.

Cold, this returns the widest split, which is what the previous arithmetic rule
produced -- so an uncalibrated machine behaves exactly as before and the
adaptation can only improve on it. With history, the same UCB rule the route
selector uses picks the empirically fastest width, and the same
`adaptive_min_samples` guarantee explores each width before exploiting.

A process route always reports an inner budget of one thread per worker
regardless of width, since its workers are separate processes launched with a
single thread each; the width still matters there, because more processes is not
always faster once memory bandwidth or native-library instances dominate.
"""
function select_outer_split!(
    state::OuterRouteState,
    f::OuterRouteFeatures;
    route::Symbol,
    budget::Int,
    n_units::Int,
    tuning::OuterRouteTuning=OuterRouteTuning()
)::Int
    candidates = outer_split_candidates(route; budget=budget, n_units=n_units, tuning=tuning)
    length(candidates) <= 1 && return isempty(candidates) ? 1 : candidates[1]
    widest = candidates[end]
    tuning.adaptive_enabled || return widest

    arms = Symbol[_split_arm(route, w) for w in candidates]
    signature_chain = _split_signature_chain(f)
    snapshot = Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate, :std_s), Tuple{Int, Float64, Float64, Float64}}}()
    for candidate_sig in signature_chain
        snap = _outer_route_stats_snapshot_internal(state, candidate_sig)
        if any(a -> haskey(snap, a), arms)
            snapshot = snap
            break
        end
    end
    isempty(snapshot) && return widest

    widest_arm = _split_arm(route, widest)
    widest_proven = tuning.explore_until_any_proven ?
        _any_candidate_proven(arms, snapshot, tuning.adaptive_min_samples) :
        _route_is_proven(arms, snapshot, widest_arm, tuning.adaptive_min_samples)
    explore = widest_proven ?
        nothing :
        _under_sampled_candidate(arms, snapshot, widest_arm, tuning.adaptive_min_samples)
    if !(explore === nothing)
        for (w, arm) in zip(candidates, arms)
            arm === explore && return w
        end
        return widest
    end

    best = _best_candidate_confidence(arms, snapshot, widest_arm, tuning.adaptive_exploration_c)
    best.route === nothing && return widest
    for (w, arm) in zip(candidates, arms)
        arm === best.route && return w
    end
    return widest
end

"""
    record_outer_split_feedback!(state, features; route, workers, kwargs...)

Record the outcome of running `route` at `workers` wide, so later campaigns can
select the width empirically. Delegates to
[`record_outer_route_feedback!`](@ref) under the split arm, so the split shares
the route bandit's statistics and persistence rather than duplicating them.
"""
function record_outer_split_feedback!(
    state::OuterRouteState,
    f::OuterRouteFeatures;
    route::Symbol,
    workers::Int,
    successes::Int,
    failures::Int,
    elapsed_success_s::Float64=0.0,
    elapsed_success_sq_sum_s::Float64=NaN,
    tuning::OuterRouteTuning=OuterRouteTuning(),
    weight::Int=1,
)::Nothing
    return record_outer_route_feedback!(
        state, f;
        route=_split_arm(route, workers),
        successes=successes,
        failures=failures,
        elapsed_success_s=elapsed_success_s,
        elapsed_success_sq_sum_s=elapsed_success_sq_sum_s,
        tuning=tuning,
        signature_prefix=_SPLIT_SIGNATURE_PREFIX,
        weight=weight,
    )
end

"""
    outer_split_history_present(state, features, route) -> Bool

Whether any split arm for `route` has been observed under this workload's
split signature chain. False means the split selector would answer cold, which
is what the in-campaign race exists to replace.
"""
function outer_split_history_present(state::OuterRouteState, f::OuterRouteFeatures, route::Symbol)::Bool
    prefix = "split_" * String(route) * "_w"
    for sig in _split_signature_chain(f)
        snap = _outer_route_stats_snapshot_internal(state, sig)
        for arm in keys(snap)
            startswith(String(arm), prefix) && return true
        end
    end
    return false
end
