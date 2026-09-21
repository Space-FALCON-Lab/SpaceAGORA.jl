module ConstellationNetwork

export NetworkTopology, NetworkState, NetworkInputs, step_network!

Base.@kwdef struct NetworkTopology
    node_ids::Vector{Int} = Int[]
    adjacency::BitMatrix = falses(0, 0)
end

Base.@kwdef mutable struct NetworkState
    active_links::Int = 0
    queued_bits::Float64 = 0.0
    dropped_bits::Float64 = 0.0
end

Base.@kwdef struct NetworkInputs
    epoch_s::Float64 = 0.0
    dt_s::Float64 = 0.0
    requested_routes::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[]
end

function step_network!(::NetworkTopology, ::NetworkState, ::NetworkInputs)
    throw(ErrorException("Not implemented: step_network!(::NetworkTopology, ::NetworkState, ::NetworkInputs)"))
end

end # module ConstellationNetwork
