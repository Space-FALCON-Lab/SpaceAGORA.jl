module Estimation

export EstimatorState, MeasurementPacket, estimate_state!

Base.@kwdef mutable struct EstimatorState
    x_hat::Vector{Float64} = Float64[]
    covariance::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    last_update_s::Float64 = 0.0
end

Base.@kwdef struct MeasurementPacket
    epoch_s::Float64 = 0.0
    values::Vector{Float64} = Float64[]
    source::String = ""
end

function estimate_state!(::EstimatorState, ::MeasurementPacket)
    throw(ErrorException("Not implemented: estimate_state!(::EstimatorState, ::MeasurementPacket)"))
end

end # module Estimation
