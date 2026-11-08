using Test
using SpaceAGORA
const PP_SBA = SpaceAGORA.SimulationModel.ParallelPolicy

# Regression (2026-11-08). The spin-barrier pool hands each dispatch to its
# workers through `pool.request::Base.RefValue{Any}`. The worker loop read the
# request's fields straight from that untyped value, so every worker boxed a copy
# of each field it read, on every dispatch. With the fused harmonics RHS the
# request carries a large closure, and this came to about 4 KB per worker per RHS
# call: a 64-satellite degree-50 constellation allocated 9.1 GB with the barrier
# on against 1.8 GB with it off (12 threads), and 16.9 GB at 24 threads. The
# worker now reads the request in a function specialized on its concrete type.
#
# The test compares a dispatch on two workers with one on every thread, so the
# fixed per-dispatch cost (the boxed request itself, the shared chunk counter)
# cancels and only per-worker allocation remains.
@testset "Spin-barrier dispatch does not allocate per worker" begin
    n = Threads.nthreads()
    if n < 3
        @test_skip n >= 3   # needs three or more default threads
    else
        acc = zeros(Int, n)
        captured = ntuple(k -> zeros(k), 24)   # a closure with many captures, like the RHS one
        f = (w, i) -> (@inbounds acc[w] += i + length(captured[1 + (i % 24)]); nothing)
        dispatch(workers) = PP_SBA.threaded_foreach_worker_spin(f, :spin_alloc_regression, n, workers)
        function bytes_per_dispatch(workers; reps=500)
            for _ in 1:50
                dispatch(workers)
            end
            # Least of several samples, so a stray allocation elsewhere cannot fail it.
            return minimum(1:5) do _
                (@allocated for _ in 1:reps
                    dispatch(workers)
                end) / reps
            end
        end
        try
            per_worker = (bytes_per_dispatch(n) - bytes_per_dispatch(2)) / (n - 2)
            # Before the fix this was 1536 B per worker for this closure at 4 threads.
            @test per_worker < 64
        finally
            PP_SBA._destroy_persistent_foreach_scope!(PP_SBA._active_policy_scope_id())
        end
    end
end
