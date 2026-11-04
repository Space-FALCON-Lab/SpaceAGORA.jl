using Test
using SpaceAGORA
const PPr = SpaceAGORA.ParallelProfiles
const SCamp = SpaceAGORA.SimulationCampaigns

@testset "memory_worker_cap does not re-charge workers that are already alive" begin
    withenv("SPACEAGORA_PERF_WORKER_MEMORY_GB" => "1.5", "SPACEAGORA_MEMORY_BUDGET_GB" => "") do
        cold = PPr.memory_worker_cap()
        @test PPr.memory_worker_cap(resident = 12) >= 12
        @test PPr.memory_worker_cap(resident = 12) >= cold
        # A workload term so large that not even the alive workers fit: as many as do, never more.
        huge = 1 << 40
        @test PPr.memory_worker_cap(resident = 12, extra_per_worker = huge) <= 12
        @test PPr.memory_worker_cap(resident = 0, extra_per_worker = huge) == 0
        # Local slots: alive workers are charged their workload term only.
        @test PPr.memory_local_slot_cap(12; extra_per_worker = 1 << 20, resident = 12) >=
              PPr.memory_local_slot_cap(12; extra_per_worker = 1 << 20)
    end
end

@testset "the tuning carries the resident count into the process cap" begin
    withenv("SPACEAGORA_PARALLEL_POLICY_V2" => "1", "SPACEAGORA_PERF_WORKER_MEMORY_GB" => "1.5") do
        t0 = PPr.OuterRouteTuning(memory_aware = true, process_workers_resident = 0)
        t12 = PPr.OuterRouteTuning(memory_aware = true, process_workers_resident = 12)
        @test t12.process_max_workers >= min(12, PPr.usable_core_budget())
        @test t12.process_max_workers >= t0.process_max_workers
        f = SCamp.campaign_route_features(samples = 64, n_sats = 16, density_family = "exponential", mission_time_s = 3600.0)
        @test PPr.effective_process_workers(f, t12) >= min(12, PPr.usable_core_budget())
    end
    # The campaign runner declares its pool's alive workers.
    t = SCamp._campaign_route_tuning()
    @test t.process_workers_resident == length(SpaceAGORA.campaign_process_pool().workers)
end
