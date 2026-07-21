const PAPER_ACTIONS_MPS = Float64[-1.0, -0.5, -0.3, -0.2, -0.1, -0.05, 0.0,
                                  0.05, 0.1, 0.2, 0.3, 0.5, 1.0]

struct AerobrakingAction
    index::Int
    delta_v_mps::Float64
    magnitude_mps::Float64
    phi_deg::Float64
    lowers_periapsis::Bool
    raises_periapsis::Bool
end

action_count() = length(PAPER_ACTIONS_MPS)
zero_action_index() = findfirst(==(0.0), PAPER_ACTIONS_MPS)

function action_from_index(index::Integer)
    1 <= index <= action_count() || throw(BoundsError(PAPER_ACTIONS_MPS, index))
    delta_v = PAPER_ACTIONS_MPS[index]
    magnitude = abs(delta_v)
    phi = delta_v <= 0 ? 0.0 : 180.0
    return AerobrakingAction(index, delta_v, magnitude, phi, delta_v < 0, delta_v > 0)
end

function action_from_delta_v(delta_v_mps::Real)
    delta_v = Float64(delta_v_mps)
    magnitude = abs(delta_v)
    phi = delta_v <= 0 ? 0.0 : 180.0
    return AerobrakingAction(0, delta_v, magnitude, phi, delta_v < 0, delta_v > 0)
end

function nearest_action_index(delta_v_mps::Real)
    _, idx = findmin(abs.(PAPER_ACTIONS_MPS .- Float64(delta_v_mps)))
    return idx
end
