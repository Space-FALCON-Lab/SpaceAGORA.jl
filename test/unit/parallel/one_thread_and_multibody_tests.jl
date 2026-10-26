using Test
using SpaceAGORA
const PP = SpaceAGORA.SimulationModel.ParallelPolicy
const AE = SpaceAGORA.SimulationModel.DynamicEffectors.AerodynamicEffectors

@testset "thread_policy_decision on one OS thread answers without reading the environment" begin
    d = PP.thread_policy_decision(8; mode = :auto, threshold = 1, source = :density_callback)
    @test d.use_threads == false || Base.Threads.nthreads() > 1
    @test d.allotment >= 1
    if Base.Threads.nthreads() == 1
        # The forced answer, and the same fields as the full path.
        @test d.use_threads == false && d.allotment == 1 && d.budget == 1 && d.desire == 1
        @test Set(keys(d)) == Set((:use_threads, :allotment, :budget, :mode, :threshold, :num_items, :adaptive_enabled, :desire))
        @test d.adaptive_enabled == false
    end
end

@testset "the multibody parallel mode is cached per solve and refreshed on demand" begin
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "off") do
        @test AE.refresh_multibody_parallel_mode!() == :off
        @test AE._multibody_parallel_mode() == :off
        withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "auto") do
            @test AE._multibody_parallel_mode() == :off            # cached: the solve has not restarted
            @test AE.refresh_multibody_parallel_mode!() == :auto  # the engine's per-solve refresh
            @test AE._multibody_parallel_mode() == :auto
        end
    end
    AE.refresh_multibody_parallel_mode!()
end

@testset "multibody callers share the neutral policy owner and per-solve cache" begin
    @test AE.refresh_multibody_parallel_mode! === PP.refresh_multibody_parallel_mode!
    @test AE._multibody_thread_decision === PP._multibody_thread_decision
    @test AE._MULTIBODY_MODE_CACHE === PP._MULTIBODY_MODE_CACHE
    @test SpaceAGORA.SimulationModel.DynamicEffectors.PerturbationEffectors._multibody_thread_decision === PP._multibody_thread_decision
    previous = PP._MULTIBODY_MODE_CACHE[]
    try
        withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "off") do
            @test PP.refresh_multibody_parallel_mode!() == :off
            @test AE._multibody_thread_decision(8) == (use_threads=false, allotment=1, mode=:off)
            withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
                    "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "4") do
                # The environment change must wait for the next solve refresh,
                # regardless of which historical namespace the caller uses.
                @test PP._multibody_thread_decision(8).mode == :off
                @test AE.refresh_multibody_parallel_mode!() == :auto
                @test PP._multibody_thread_decision(2) == (use_threads=false, allotment=1, mode=:auto)
                @test AE._multibody_thread_decision(2; heavy_work=false) == (use_threads=false, allotment=1, mode=:auto)
            end
        end
        withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on",
                "SPACEAGORA_MULTIBODY_MAX_THREADS" => "1") do
            PP.refresh_multibody_parallel_mode!()
            @test AE._multibody_thread_decision(1) == (use_threads=false, allotment=1, mode=:on)
            # The explicit cap must still win over forced-on mode.
            @test AE._multibody_thread_decision(8) == (use_threads=false, allotment=1, mode=:on)
        end
    finally
        PP._MULTIBODY_MODE_CACHE[] = previous
    end
end
