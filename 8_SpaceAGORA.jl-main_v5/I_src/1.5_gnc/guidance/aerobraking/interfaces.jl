abstract type AbstractAerobrakingStrategy end

struct EEdgStrategy <: AbstractAerobrakingStrategy end
struct TEdgStrategy <: AbstractAerobrakingStrategy end

Base.@kwdef struct AerobrakingGuidanceInput
    ip
    mission
    args
    index_ratio::Vector{Int} = Int[]
    state = nothing
    t::Float64 = 0.0
    position = 0
    current_position = 0
    gram_atmosphere = nothing
    heat_rate_control::Bool = false
    reevaluation_mode::Int = 1
    cnf = nothing
end

Base.@kwdef struct AerobrakingGuidanceOutput
    time_switch_1::Float64 = 0.0
    time_switch_2::Float64 = 0.0
    security_mode::Bool = false
end

function compute_aerobraking_guidance(::AbstractAerobrakingStrategy, ::AerobrakingGuidanceInput)
    throw(MethodError(compute_aerobraking_guidance, (AbstractAerobrakingStrategy, AerobrakingGuidanceInput)))
end
