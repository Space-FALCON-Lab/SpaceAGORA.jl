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
