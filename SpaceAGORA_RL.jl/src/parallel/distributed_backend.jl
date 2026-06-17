using Distributed

function setup_distributed_workers(n_workers::Integer)
    n = max(1, Int(n_workers))
    if n == 1
        return workers()
    end
    needed = n - length(workers())
    needed > 0 && addprocs(needed)
    return workers()[1:n]
end
