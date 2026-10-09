using Test
using SpaceAGORA
using Polyester
const PP_SBY = SpaceAGORA.SimulationModel.ParallelPolicy

# Regression (2026-10-09). The spin-barrier pool's workers keep spinning between
# dispatches. They used to poll with GC.safepoint() only, never yielding, so a
# task pinned to one of their threads could never run. Polyester's `@batch`
# pins its workers to threads: the RHS calibration sweep (predictive / R7
# profiles) timed a flat harmonics candidate on the spin pool and then a
# satellite_batch candidate through `@batch`, and the solve hung for good with
# two threads at 100% (TRX50 job 20261008-115111-1480699, 25 h, and reproduced
# on the workstation at 2 threads).
#
# Without the fix this test fails by timeout, not by hanging: it stops the
# spinners afterwards, which releases the pinned task.
@testset "Spin-barrier workers do not starve thread-pinned tasks" begin
    n = Threads.nthreads()
    if n < 2
        @test_skip n >= 2   # needs two or more default threads
    else
        # Which thread a spinner lands on is up to the scheduler, so at 2 threads one
        # trial can miss the thread Polyester needs. Several trials, 20 s each.
        for trial in 1:4
            hits = zeros(Int, 8n)
            # Creates the spin pool (nthreads-1 workers) and leaves its workers spinning.
            PP_SBY.threaded_foreach_worker_spin(:spin_yield_regression, n, n) do _w, _i
                nothing
            end
            done = Threads.Atomic{Bool}(false)
            pool = Threads.nthreads(:interactive) > 0 ? :interactive : :default
            batch = Threads.@spawn pool begin
                Polyester.@batch per=thread for i in 1:(8n)
                    hits[i] += 1
                end
                done[] = true
            end
            status = timedwait(() -> done[], 20.0; pollint=0.05)
            # Stop the spinners either way, so a failure cannot leave the process hung.
            PP_SBY._destroy_persistent_foreach_scope!(PP_SBY._active_policy_scope_id())
            wait(batch)
            @test status === :ok
            @test all(==(1), hits)
            status === :ok || break
        end
    end
end
