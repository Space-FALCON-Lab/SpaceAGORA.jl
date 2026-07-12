# Canonical aggregator: no behavior ownership.
module ParallelProcess

using Distributed

include(joinpath(@__DIR__, "worker_pool.jl"))

export ProcessPool, campaign_process_pool, ensure_process_workers!, shutdown_process_pool!

end # module ParallelProcess
