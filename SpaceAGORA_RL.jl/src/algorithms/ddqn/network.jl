mutable struct QNetwork{M<:AbstractMatrix{Float32},V<:AbstractVector{Float32}}
    W1::M
    b1::V
    W2::M
    b2::V
    W3::M
    b3::V
end

mutable struct QNetworkGradients{M<:AbstractMatrix{Float32},V<:AbstractVector{Float32}}
    W1::M
    b1::V
    W2::M
    b2::V
    W3::M
    b3::V
end

function orthogonal_weight(rng::AbstractRNG, rows::Integer, cols::Integer; gain::Float32=Float32(sqrt(2)))
    rows_i = Int(rows)
    cols_i = Int(cols)
    if rows_i >= cols_i
        q = Matrix(qr(randn(rng, Float32, rows_i, cols_i)).Q)
        return gain .* q[:, 1:cols_i]
    else
        q = Matrix(qr(randn(rng, Float32, cols_i, rows_i)).Q)
        return gain .* transpose(q[:, 1:rows_i])
    end
end

function init_q_network(rng::AbstractRNG; input_dim::Int=9, hidden_dim::Int=1024, output_dim::Int=13)
    return QNetwork(
        orthogonal_weight(rng, hidden_dim, input_dim),
        zeros(Float32, hidden_dim),
        orthogonal_weight(rng, hidden_dim, hidden_dim),
        zeros(Float32, hidden_dim),
        orthogonal_weight(rng, output_dim, hidden_dim),
        zeros(Float32, output_dim),
    )
end

function Base.copy(net::QNetwork)
    return QNetwork(copy(net.W1), copy(net.b1), copy(net.W2), copy(net.b2),
                    copy(net.W3), copy(net.b3))
end

function copy_network!(dest::QNetwork, src::QNetwork)
    dest.W1 .= src.W1
    dest.b1 .= src.b1
    dest.W2 .= src.W2
    dest.b2 .= src.b2
    dest.W3 .= src.W3
    dest.b3 .= src.b3
    return dest
end

function zero_gradients_like(net::QNetwork)
    return QNetworkGradients(
        zero_array_like(net.W1), zero_array_like(net.b1),
        zero_array_like(net.W2), zero_array_like(net.b2),
        zero_array_like(net.W3), zero_array_like(net.b3),
    )
end

zero_network_like(net::QNetwork) =
    QNetwork(zero_array_like(net.W1), zero_array_like(net.b1),
             zero_array_like(net.W2), zero_array_like(net.b2),
             zero_array_like(net.W3), zero_array_like(net.b3))

relu_derivative(x) = x > 0 ? Float32(1) : Float32(0)

@inline _as_float32_array(values::AbstractArray{Float32}) = values
@inline _as_float32_array(values::AbstractArray{<:Real}) = Float32.(values)

function forward_cache(net::QNetwork, observations::AbstractMatrix{<:Real})
    x = _as_float32_array(observations)
    z1 = net.W1 * x .+ net.b1
    a1 = max.(z1, Float32(0))
    z2 = net.W2 * a1 .+ net.b2
    a2 = max.(z2, Float32(0))
    q = net.W3 * a2 .+ net.b3
    return x, z1, a1, z2, a2, q
end

predict_q(net::QNetwork, observations::AbstractMatrix{<:Real}) = forward_cache(net, observations)[end]
predict_q(net::QNetwork, observation::AbstractVector{<:Real}) =
    vec(predict_q(net, reshape(_as_float32_array(observation), :, 1)))

function _network_gradients_from_output_delta(
    net::QNetwork,
    cache::Tuple,
    dY::AbstractMatrix{Float32},
)
    x, z1, a1, z2, a2, q = cache
    size(dY) == size(q) || throw(DimensionMismatch("output delta must match network output size"))
    dW3 = dY * transpose(a2)
    db3 = vec(sum(dY; dims=2))
    da2 = transpose(net.W3) * dY
    dz2 = da2 .* Float32.(z2 .> 0)
    dW2 = dz2 * transpose(a1)
    db2 = vec(sum(dz2; dims=2))
    da1 = transpose(net.W2) * dz2
    dz1 = da1 .* Float32.(z1 .> 0)
    dW1 = dz1 * transpose(x)
    db1 = vec(sum(dz1; dims=2))

    return QNetworkGradients(dW1, db1, dW2, db2, dW3, db3)
end

function network_gradients_from_output_delta(net::QNetwork, observations::AbstractMatrix{Float32},
                                             dY::AbstractMatrix{Float32})
    return _network_gradients_from_output_delta(net, forward_cache(net, observations), dY)
end

function onehot_actions(actions::Vector{Int}, action_dim::Integer)
    encoded = zeros(Float32, Int(action_dim), length(actions))
    for (col, action) in pairs(actions)
        1 <= action <= action_dim || throw(BoundsError(1:action_dim, action))
        encoded[action, col] = 1f0
    end
    return encoded
end

function network_loss_and_gradients(net::QNetwork, observations::AbstractMatrix{Float32},
                                    actions::Vector{Int}, targets::AbstractVector{Float32};
                                    device::AbstractTrainingDevice=CPUTrainingDevice())
    cache = forward_cache(net, observations)
    q = cache[end]
    batch_size = size(observations, 2)
    length(actions) == batch_size || throw(DimensionMismatch("actions length must match batch size"))
    length(targets) == batch_size || throw(DimensionMismatch("targets length must match batch size"))
    encoded = to_device_array(device, onehot_actions(actions, size(q, 1)))
    target_row = reshape(to_device_array(device, _as_float32_array(targets)), 1, :)
    selected_q = sum(q .* encoded; dims=1)
    err = selected_q .- target_row
    loss = mean(abs2, err)
    dQ = (2f0 / Float32(batch_size)) .* err .* encoded
    return Float64(cpu_scalar(loss)), _network_gradients_from_output_delta(net, cache, dQ)
end

function gradient_norm(grads::QNetworkGradients)
    total = sum(abs2, grads.W1) + sum(abs2, grads.b1) +
            sum(abs2, grads.W2) + sum(abs2, grads.b2) +
            sum(abs2, grads.W3) + sum(abs2, grads.b3)
    return sqrt(Float64(cpu_scalar(total)))
end

function scale_gradients!(grads::QNetworkGradients, scale::Real)
    s = Float32(scale)
    grads.W1 .*= s
    grads.b1 .*= s
    grads.W2 .*= s
    grads.b2 .*= s
    grads.W3 .*= s
    grads.b3 .*= s
    return grads
end

function network_to_device(device::AbstractTrainingDevice, net::QNetwork)
    return QNetwork(
        to_device_array(device, net.W1),
        to_device_array(device, net.b1),
        to_device_array(device, net.W2),
        to_device_array(device, net.b2),
        to_device_array(device, net.W3),
        to_device_array(device, net.b3),
    )
end

cpu_network(net::QNetwork) = network_to_device(CPUTrainingDevice(), net)
