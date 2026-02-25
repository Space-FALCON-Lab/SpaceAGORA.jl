@kwdef mutable struct BaseThrusterModel <: AbstractThrusterModel
    thrust::Vector{Float64}
    direction::Vector{Float64}
    Δv::Vector{Float64}
    start_burn_time::Vector{Float64}
    stop_burn_time::Vector{Float64}
    Isp::Vector{Float64}
end
