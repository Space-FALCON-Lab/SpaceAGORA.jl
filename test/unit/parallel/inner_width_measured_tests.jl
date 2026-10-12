using Test
using SpaceAGORA

# The predictive planner's inner width (the +bN of threads@w16+l0+b2), priced
# from campaign walls timed in this process. Shapes are the TRX50 X5 points:
# X5b mcgrid_8sat_16mc (n = 16, T = 32; b1 0.272 s, b2 0.246 s) and X5a
# mcgrid_16sat_8mc (n = 8, T = 32; b1 0.383 s, b4 0.345 s). No solves.

const IW = SpaceAGORA.SimulationCampaigns

_iw_cfg(T) = IW.PredictivePlannerConfig(
    margin = 0.15, guard_factor = 1.5, local_slots_max = T - 1, remote_overhead = 0.0,
    route_switch = false, heap_model = :locals, local_thrash = 3.0, measured_margin = 0.07)

function _iw_time!(sig, n, T, W, b, secs...)
    plan = IW._predictive_plan(:threads, W, 0, n, b == 1, nothing, 0.0; inner_budget = b)
    for s in secs
        IW._predictive_record_inner_wall!(sig, n, T, plan, s)
    end
end

function _iw_chosen(sig, n, T; pool = 0)
    curve = IW._predictive_measured_inner_curve(nothing, sig, n, T)
    p = IW.predictive_plan(n_samples = n, threads = T, process_workers = pool,
                           threads_candidate = true, local_slots_cap = T - 1,
                           constants = nothing, config = _iw_cfg(T), inner_curve = curve)
    return p, curve
end

@testset "measured inner width: X5b picks b2, X5a picks b4, flat keeps b1" begin
    # X5b: 0.272 -> 0.246 is a 9.6% edge, under the 15% margin but over 7%.
    _iw_time!("x5b", 16, 32, 16, 1, 0.9, 0.272, 0.275)   # first run cold; min kept
    _iw_time!("x5b", 16, 32, 16, 2, 0.33, 0.246)
    for pool in (0, 16)
        p, curve = _iw_chosen("x5b", 16, 32; pool = pool)
        @test curve.measured == [2]
        @test IW.predictive_plan_key(p.chosen) == "threads@w16+l0+b2"
        @test p.reason === :predicted_gain
        @test 0.07 <= p.gain < 0.15
        @test p.chosen.workers * p.chosen.inner_thread_budget <= 32
        @test all(q -> q.route !== :threads || q.workers * max(1, q.inner_thread_budget) <= 32, p.plans)
    end
    # The same measurement under the plain margin would not have deviated.
    curve = IW._predictive_measured_inner_curve(nothing, "x5b", 16, 32)
    strict = IW.PredictivePlannerConfig(margin = 0.15, guard_factor = 1.5, local_slots_max = 31,
                                        measured_margin = 0.15)
    q = IW.predictive_plan(n_samples = 16, threads = 32, process_workers = 0,
                           threads_candidate = true, local_slots_cap = 31, constants = nothing,
                           config = strict, inner_curve = curve)
    @test q.chosen.static_equivalent && q.reason === :margin_not_met

    # X5a: b4 at 8 consumers.
    _iw_time!("x5a", 8, 32, 8, 1, 0.383, 0.384)
    _iw_time!("x5a", 8, 32, 8, 4, 0.5, 0.345)
    p, _ = _iw_chosen("x5a", 8, 32)
    @test IW.predictive_plan_key(p.chosen) == "threads@w8+l0+b4"
    @test p.chosen.workers * p.chosen.inner_thread_budget <= 32

    # Flat: wider budget measured no faster than b1 -> b1.
    _iw_time!("flat", 16, 32, 16, 1, 0.272, 0.271)
    _iw_time!("flat", 16, 32, 16, 2, 0.4, 0.270)
    p, _ = _iw_chosen("flat", 16, 32)
    @test p.chosen.static_equivalent
    @test max(1, p.chosen.inner_thread_budget) == 1
end

@testset "no timing, no curve: the store-priced plan space is unchanged" begin
    @test IW._predictive_measured_inner_curve(nothing, "never-ran", 16, 32) === nothing
    base = IW.InnerSpeedupCurve(ones(4))
    # b = 1 not timed: the base curve is returned as is.
    @test IW._predictive_measured_inner_curve(base, "never-ran", 16, 32) === base
end

@testset "the probe times the enumerated wider width twice" begin
    sig = "probe-x5a"
    cfg = _iw_cfg(32)
    plan() = IW.predictive_plan(n_samples = 8, threads = 32, process_workers = 0,
                                threads_candidate = true, local_slots_cap = 31,
                                constants = nothing, config = cfg)
    probe() = IW._predictive_inner_probe(plan(), sig, 8, 32, cfg, nothing, IW.PredictiveCostTerms(), nothing)
    @test probe() === nothing                      # nothing timed yet
    _iw_time!(sig, 8, 32, 8, 1, 0.40)
    @test probe() === nothing                      # reference seen once: cold
    _iw_time!(sig, 8, 32, 8, 1, 0.383)
    p = probe()
    @test IW.predictive_plan_key(p) == "threads@w8+l0+b4"
    @test p.workers * p.inner_thread_budget <= 32
    _iw_time!(sig, 8, 32, 8, 4, 0.6)               # the cold run of the width
    @test IW.predictive_plan_key(probe()) == "threads@w8+l0+b4"
    _iw_time!(sig, 8, 32, 8, 4, 0.39)              # timed twice: never again
    @test probe() === nothing
    # One observation at either end is not a curve.
    _iw_time!("once", 16, 32, 16, 1, 0.272, 0.275)
    _iw_time!("once", 16, 32, 16, 2, 0.2)
    @test IW._predictive_measured_inner_curve(nothing, "once", 16, 32) === nothing
end
