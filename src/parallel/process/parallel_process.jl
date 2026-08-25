# Canonical aggregator: no behavior ownership.
module ParallelProcess

using Distributed

include(joinpath(@__DIR__, "worker_pool.jl"))
include(joinpath(@__DIR__, "density_service.jl"))

export ProcessPool, campaign_process_pool, ensure_process_workers!, shutdown_process_pool!
export density_process_pool, ensure_density_workers!, density_batch_remote
export density_service_partition, density_service_dispatch, install_worker_density_model!
export shutdown_density_workers!

end # module ParallelProcess
