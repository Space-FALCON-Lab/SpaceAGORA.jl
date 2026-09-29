# The native-GRAM per-spacecraft memory model and the route decisions it drives.
#
# Everything here runs on fabricated machine facts -- the memory budget, the
# headroom and the per-worker package image are all supplied through the
# documented overrides -- so no GRAM build, no SPICE kernels and no knowledge of
# the host are needed. The measured constants themselves are sourced in
# docs/architecture/gram_memory_footprint.md; these tests assert the relations
# between them and the arithmetic they feed, not their literal values, so a
# re-measurement on another machine does not have to rewrite the suite.

using Test
using SpaceAGORA
const PPm = SpaceAGORA.ParallelProfiles

const MB = 1 << 20
const GB = 1 << 30

"""Machine facts for one profile: the budget is never the binding term here, so
the headroom alone decides, and the per-worker package image is pinned."""
function profile_env(; headroom_gb::Real, worker_image_gb::Real = 1.5, cores::Int = 12)
    return (
        "SPACEAGORA_PARALLEL_POLICY_V2" => "1",
        "SPACEAGORA_CORE_BUDGET" => string(cores),
        "SPACEAGORA_MEMORY_BUDGET_GB" => "100000",
        "SPACEAGORA_MEMORY_AVAILABLE_GB" => string(headroom_gb),
        "SPACEAGORA_PERF_WORKER_MEMORY_GB" => string(worker_image_gb),
    )
end

gram_features(n::Int) = PPm.OuterRouteFeatures(
    category = "montecarlo",
    n_sats = n,
    density_family = "gram_point",
    mission_time_s = 3600.0,
    montecarlo_samples = 64,
)

@testset "per-path constants" begin
    # Every native calling method carries a positive per-spacecraft charge, and
    # the unqualified default is the largest of them: a caller that cannot say
    # which path it is on must be charged the worst one.
    native = (:point, :freeze_per_step, :lookahead)
    for p in native
        @test PPm.gram_sat_memory_bytes(p) > 0
        @test PPm.gram_sat_memory_bytes(p) <= PPm._GRAM_SAT_MEMORY_BYTES
    end
    @test PPm._GRAM_SAT_MEMORY_BYTES == maximum(PPm.gram_sat_memory_bytes(p) for p in native)
    # The offline surrogate makes no native call and holds no per-spacecraft
    # native state, so it is charged nothing.
    @test PPm.gram_sat_memory_bytes(:surrogate) == 0
    # An unrecognized name is the default, never zero.
    @test PPm.gram_sat_memory_bytes(:something_new) == PPm._GRAM_SAT_MEMORY_BYTES

    # Whatever a re-measurement says, the model must stay far from the 90 MB the
    # constant carried before it was measured; a regression back to that order
    # of magnitude is what this bound catches.
    @test PPm._GRAM_SAT_MEMORY_BYTES < 30 * MB
end

@testset "native_gram_worker_extra_bytes" begin
    @test PPm.native_gram_worker_extra_bytes(0) == 0
    @test PPm.native_gram_worker_extra_bytes(-5) == 0
    @test PPm.native_gram_worker_extra_bytes(256) == 256 * PPm._GRAM_SAT_MEMORY_BYTES
    for p in (:point, :freeze_per_step, :lookahead, :surrogate)
        @test PPm.native_gram_worker_extra_bytes(1024; path = p) ==
              1024 * PPm.gram_sat_memory_bytes(p)
    end
    # The single env override applies to every path, as it did before the model
    # was split, so an operator who pins it still pins the whole model.
    withenv("SPACEAGORA_GRAM_SAT_MEMORY_MB" => "7") do
        for p in (:point, :freeze_per_step, :lookahead, :surrogate)
            @test PPm.native_gram_worker_extra_bytes(10; path = p) == 10 * 7 * MB
        end
    end
    withenv("SPACEAGORA_GRAM_SAT_MEMORY_MB" => "0") do
        @test PPm.native_gram_worker_extra_bytes(4096) == 0
    end
    for bad in ("-1", "abc")
        withenv("SPACEAGORA_GRAM_SAT_MEMORY_MB" => bad) do
            @test_throws ArgumentError PPm.native_gram_worker_extra_bytes(4)
        end
    end
end

@testset "headroom override supplies a machine other than the host" begin
    withenv(profile_env(headroom_gb = 223.5)...) do
        @test PPm.memory_headroom_bytes() == round(Int, 223.5 * GB)
    end
    # The budget still binds when it is the smaller term.
    withenv("SPACEAGORA_MEMORY_BUDGET_GB" => "2",
            "SPACEAGORA_MEMORY_AVAILABLE_GB" => "1000") do
        @test PPm.memory_headroom_bytes() <= 2 * GB
    end
    for bad in ("0", "-1", "abc")
        withenv("SPACEAGORA_MEMORY_AVAILABLE_GB" => bad) do
            @test_throws ArgumentError PPm.memory_headroom_bytes()
        end
    end
end

@testset "worker cap arithmetic on fabricated facts" begin
    # 52.5 GB of headroom, 1.5 GB per worker image: 35 workers with no workload
    # term at all.
    withenv(profile_env(headroom_gb = 52.5)...) do
        @test PPm.memory_worker_cap() == 35
        # A 1 GB workload term prices each worker at 2.5 GB: 21 workers.
        @test PPm.memory_worker_cap(extra_per_worker = GB) == 21
        # The old constant at 4096 spacecraft priced one worker at 360 GB, so
        # none fit; the measured model must leave the route affordable.
        old = 4096 * 90 * MB
        @test PPm.memory_worker_cap(extra_per_worker = old) == 0
        @test PPm.memory_worker_cap(
            extra_per_worker = PPm.native_gram_worker_extra_bytes(4096)) >= 2
    end
end

@testset "local slot cap on fabricated facts" begin
    withenv(profile_env(headroom_gb = 52.5)...) do
        extra = PPm.native_gram_worker_extra_bytes(1024)
        # With no workload term there is nothing measurable to reserve.
        @test PPm.memory_local_slot_cap(4) == typemax(Int) >> 1
        # Four workers at (1.5 GB + extra) each, the rest divided by one local
        # sample's working set.
        expected = fld(round(Int, 52.5 * GB) - 4 * (round(Int, 1.5 * GB) + extra), extra)
        @test PPm.memory_local_slot_cap(4; extra_per_worker = extra) == expected
        # Workers already alive are charged their workload term only, so more
        # local slots survive.
        @test PPm.memory_local_slot_cap(4; extra_per_worker = extra, resident = 4) >
              PPm.memory_local_slot_cap(4; extra_per_worker = extra)
    end
end

@testset "the route decision the constant drives" begin
    tuning() = PPm.OuterRouteTuning(memory_aware = true, mixed_dispatch = true)
    # This workstation's profile: 60 GB less the 10% OS reserve, less the
    # coordinator's own image, and a worker image pinned at 1.5 GB so the
    # arithmetic is the same on any host this runs on.
    withenv(profile_env(headroom_gb = 52.5, cores = 12)...) do
        PPm.refresh_machine_topology!()
        try
            for n in (256, 1024, 4096)
                f = gram_features(n)
                before = withenv("SPACEAGORA_GRAM_SAT_MEMORY_MB" => "90") do
                    PPm.effective_process_workers(f, tuning())
                end
                after = PPm.effective_process_workers(f, tuning())
                # The measured charge can only widen the pool, never narrow it.
                @test after >= before
                # And it leaves a usable pool at every size.
                @test after >= 2
            end
            # At 1024 and beyond the old constant priced one worker above the
            # whole budget, so the process route was withdrawn outright.
            for n in (1024, 4096)
                withenv("SPACEAGORA_GRAM_SAT_MEMORY_MB" => "90") do
                    @test PPm.effective_process_workers(gram_features(n), tuning()) < 2
                end
            end
        finally
            PPm.refresh_machine_topology!()
        end
    end
end
