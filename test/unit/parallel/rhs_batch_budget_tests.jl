using Test
using SpaceAGORA
using Polyester

# The satellite_batch RHS fans a Polyester loop over the spacecraft. Its width
# must follow the solve's inner thread budget, not the machine's core count:
# under a threaded outer campaign every sample would otherwise use every core.

const SE = SpaceAGORA.SimulationEngine
const PPol = SpaceAGORA.SimulationModel.ParallelPolicy

@testset "_rhs_batch_workers follows the inner budget, bounded by physical cores" begin
    cores = Polyester.num_cores()
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => nothing) do
        # no budget advertised: the pool -- today's behaviour
        @test SE._rhs_batch_workers(nothing) == max(1, min(cores, PPol.effective_inner_thread_budget()))
    end
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "1") do
        @test SE._rhs_batch_workers(nothing) == 1
        @test SE._rhs_batch_minbatch(nothing, 32) == 32          # one batch: the whole range
        @test !SE._rhs_batch_parallel_enabled(nothing, 64) ||    # gate refuses a 1-wide batch
              SE._rhs_batch_workers(nothing) > 1
    end
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "3") do
        # effective_inner_thread_budget() caps the requested "3" at
        # Threads.nthreads() (see env_config.jl) before _rhs_batch_workers
        # caps it again at physical cores -- on a host with fewer than 3
        # threads there is no way to reach a 3-wide batch, so the expected
        # value must include that same nthreads cap rather than assuming
        # it never binds.
        expected = max(1, min(cores, PPol.effective_inner_thread_budget()))
        @test SE._rhs_batch_workers(nothing) == expected
        @test SE._rhs_batch_minbatch(nothing, 32) == cld(32, expected)
    end
end

# Polyester's closure conversion cannot capture testset-scope locals, so the
# probe loop lives in a function and writes into an array it is handed.
function _batch_thread_ids!(ids::Vector{Int})
    Polyester.@batch minbatch=1000 for i in eachindex(ids)
        ids[i] = Threads.threadid()
    end
    return ids
end

@testset "a Polyester @batch with one batch runs on the calling thread" begin
    # What makes minbatch >= n a serial loop rather than a 1-task dispatch.
    ids = _batch_thread_ids!(zeros(Int, 10))
    @test all(==(Threads.threadid()), ids)
end
