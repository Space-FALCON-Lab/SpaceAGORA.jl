using Test
using SpaceAGORA

# A plan whose samples run on b > 1 threads is priced from the store rows
# measured the way those samples run: threads@W and the mixed local slots run
# under an outer split (outer=1 rows), the serial plan on the whole pool does
# not (outer=0 rows). Pooling the two priced a threads@W+b plan at a
# single-simulation speedup that the split's serialized RHS routes never reach.

const OC_SE = SpaceAGORA.SimulationEngine
const OC_SC = SpaceAGORA.SimulationCampaigns
const OC_TOKEN = OC_SE._RHS_CALIB_CODE_TOKEN

oc_sig(; budget, outer) =
    "v6|machine=testbox|budget=$(budget)|sats=5_8|effs=2|harm=0|eff=oc|dens=ExponentialAtmosphereModel|outer=$(outer)|code=$(OC_TOKEN)"
const OC_STEM = join(sort!(["machine=testbox", "sats=5_8", "effs=2", "harm=0", "eff=oc",
                            "dens=ExponentialAtmosphereModel", "code=$(OC_TOKEN)"]), "|")
oc_flat(a) = OC_SE._make_calib_flat_plan(a, :static)
oc_record(plan, ns) = OC_SE._rhs_sweep_timing_record(plan, Float64(ns), 4, 2, false)

function oc_reset_store!()
    lock(OC_SE._rhs_calib_lock) do
        empty!(OC_SE._rhs_calib_cache)
        OC_SE._rhs_calib_loaded[] = false
        OC_SE._rhs_calib_loaded_path[] = ""
    end
end

function oc_with_store(f)
    path = joinpath(mktempdir(), "rhs_calibration_outer.toml")
    withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => path) do
        oc_reset_store!()
        try
            # Single simulation: the flat route reaches 8x at eight threads.
            OC_SE._rhs_calib_store!(oc_sig(budget = 8, outer = 0), oc_flat(8), 1.0e6;
                                    timings = [oc_record(oc_flat(1), 8.0e6), oc_record(oc_flat(2), 4.0e6),
                                               oc_record(oc_flat(4), 2.0e6), oc_record(oc_flat(8), 1.0e6)])
            # Under an outer split the same shape barely gains.
            OC_SE._rhs_calib_store!(oc_sig(budget = 8, outer = 1), oc_flat(1), 8.0e6;
                                    timings = [oc_record(oc_flat(1), 8.0e6), oc_record(oc_flat(8), 7.6e6)])
            f()
        finally
            oc_reset_store!()
        end
    end
end

oc_cfg(T) = OC_SC.PredictivePlannerConfig(
    margin = 0.15, guard_factor = 1.5, local_slots_max = T - 1, remote_overhead = 0.0,
    route_switch = false, heap_model = :locals, local_thrash = 3.0)

oc_planning(; n, T, pool = 0, split, unsplit) = OC_SC.predictive_plan(
    n_samples = n, threads = T, process_workers = pool, threads_candidate = true,
    local_slots_cap = T - 1, constants = nothing, config = oc_cfg(T),
    inner_curve = split, inner_curve_unsplit = unsplit)

@testset "the curve reads only the rows measured the way a plan runs" begin
    oc_with_store() do
        single = OC_SE.rhs_inner_speedup_curve(OC_STEM; outer = false)
        split = OC_SE.rhs_inner_speedup_curve(OC_STEM; outer = true)
        pooled = OC_SE.rhs_inner_speedup_curve(OC_STEM)
        @test single[8] ≈ 8.0
        @test split[8] ≈ 8.0 / 7.6
        @test pooled[8] ≈ 8.0

        # The campaign asks for each by the stem its features registered.
        features = OC_SC.campaign_route_features(samples = 4)
        key = SpaceAGORA.ParallelProfiles.outer_route_signature(features)
        lock(OC_SC._CAMPAIGN_RHS_STEMS_LOCK) do
            OC_SC._CAMPAIGN_RHS_STEMS[key] = OC_STEM
        end
        try
            @test OC_SC._predictive_inner_curve(features; outer = true).speedup[8] ≈ 8.0 / 7.6
            @test OC_SC._predictive_inner_curve(features; outer = false).speedup[8] ≈ 8.0
        finally
            lock(OC_SC._CAMPAIGN_RHS_STEMS_LOCK) do
                delete!(OC_SC._CAMPAIGN_RHS_STEMS, key)
            end
        end
    end
end

@testset "a split plan is never priced from single-simulation rows" begin
    strong = OC_SC.InnerSpeedupCurve(Float64.(1:8))
    weak = OC_SC.InnerSpeedupCurve([1.0, 1.02, 1.03, 1.04, 1.04, 1.05, 1.05, 1.05])

    # Pooled pricing (the old behavior) chose threads@4 with eight threads each.
    old = oc_planning(n = 4, T = 32, split = strong, unsplit = strong)
    @test old.chosen.route === :threads && old.chosen.inner_thread_budget == 8

    # Priced from what a sample under the split measured, the b > 1 threads
    # plans lose; the serial plan on the whole pool keeps its own curve.
    new = oc_planning(n = 4, T = 32, split = weak, unsplit = strong)
    @test !(new.chosen.route === :threads && new.chosen.inner_thread_budget > 1)
    @test new.chosen.route === :none && new.chosen.inner_thread_budget == 32

    # No split measurements: no b > 1 split plan is offered at all.
    only_unsplit = oc_planning(n = 32, T = 16, pool = 8, split = nothing, unsplit = strong)
    @test !any(p -> p.route !== :none && p.inner_thread_budget > 1, only_unsplit.plans)
    # No unsplit measurements: no b > 1 serial plan.
    only_split = oc_planning(n = 4, T = 32, split = strong, unsplit = nothing)
    @test !any(p -> p.route === :none && p.inner_thread_budget > 1, only_split.plans)
    # Neither: the v1 plan space.
    v1 = oc_planning(n = 32, T = 16, pool = 8, split = nothing, unsplit = nothing)
    @test all(p -> p.inner_thread_budget <= 1, v1.plans)
end

@testset "one-thread plans are priced exactly as before" begin
    strong = OC_SC.InnerSpeedupCurve(Float64.(1:8))
    weak = OC_SC.InnerSpeedupCurve([1.0, 1.1, 1.2, 1.2])
    for (n, T, pool) in ((4, 32, 0), (32, 16, 8), (64, 24, 12), (11, 8, 4), (3, 4, 0))
        ref = Dict(OC_SC.predictive_plan_key(p) => p.makespan
                   for p in oc_planning(n = n, T = T, pool = pool, split = nothing, unsplit = nothing).plans)
        for (split, unsplit) in ((strong, strong), (weak, strong), (strong, nothing), (nothing, weak))
            planning = oc_planning(n = n, T = T, pool = pool, split = split, unsplit = unsplit)
            for p in planning.plans
                max(1, p.inner_thread_budget) == 1 || continue
                @test p.makespan == ref[OC_SC.predictive_plan_key(p)]
            end
        end
    end
end
