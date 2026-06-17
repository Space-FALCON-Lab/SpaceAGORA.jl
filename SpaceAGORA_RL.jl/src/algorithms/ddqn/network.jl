mutable struct QNetwork
    W1::Matrix{Float32}
    b1::Vector{Float32}
    W2::Matrix{Float32}
    b2::Vector{Float32}
    W3::Matrix{Float32}
    b3::Vector{Float32}
end

mutable struct QNetworkGradients
    W1::Matrix{Float32}
    b1::Vector{Float32}
    W2::Matrix{Float32}
    b2::Vector{Float32}
    W3::Matrix{Float32}
    b3::Vector{Float32}
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
        zeros(Float32, size(net.W1)), zeros(Float32, size(net.b1)),
        zeros(Float32, size(net.W2)), zeros(Float32, size(net.b2)),
        zeros(Float32, size(net.W3)), zeros(Float32, size(net.b3)),
    )
end

zero_network_like(net::QNetwork) =
    QNetwork(zeros(Float32, size(net.W1)), zeros(Float32, size(net.b1)),
             zeros(Float32, size(net.W2)), zeros(Float32, size(net.b2)),
             zeros(Float32, size(net.W3)), zeros(Float32, size(net.b3)))

relu_derivative(x) = x > 0 ? Float32(1) : Float32(0)

function forward_cache(net::QNetwork, observations::AbstractMatrix{<:Real})
    x = Float32.(observations)
    z1 = net.W1 * x .+ net.b1
    a1 = max.(z1, Float32(0))
    z2 = net.W2 * a1 .+ net.b2
    a2 = max.(z2, Float32(0))
    q = net.W3 * a2 .+ net.b3
    return x, z1, a1, z2, a2, q
end

predict_q(net::QNetwork, observations::AbstractMatrix{<:Real}) = forward_cache(net, observations)[end]
predict_q(net::QNetwork, observation::AbstractVector{<:Real}) =
    vec(predict_q(net, reshape(Float32.(observation), :, 1)))

function network_loss_and_gradients(net::QNetwork, observations::Matrix{Float32},
                                    actions::Vector{Int}, targets::Vector{Float32})
    x, z1, a1, z2, a2, q = forward_cache(net, observations)
    batch_size = size(x, 2)
    dQ = zeros(Float32, size(q))
    loss = Float32(0)
    for i in 1:batch_size
        err = q[actions[i], i] - targets[i]
        loss += err^2
        dQ[actions[i], i] = 2f0 * err / Float32(batch_size)
    end
    loss /= Float32(batch_size)

    dW3 = dQ * transpose(a2)
    db3 = vec(sum(dQ; dims=2))
    da2 = transpose(net.W3) * dQ
    dz2 = da2 .* Float32.(z2 .> 0)
    dW2 = dz2 * transpose(a1)
    db2 = vec(sum(dz2; dims=2))
    da1 = transpose(net.W2) * dz2
    dz1 = da1 .* Float32.(z1 .> 0)
    dW1 = dz1 * transpose(x)
    db1 = vec(sum(dz1; dims=2))

    return Float64(loss), QNetworkGradients(dW1, db1, dW2, db2, dW3, db3)
end

function gradient_norm(grads::QNetworkGradients)
    total = sum(abs2, grads.W1) + sum(abs2, grads.b1) +
            sum(abs2, grads.W2) + sum(abs2, grads.b2) +
            sum(abs2, grads.W3) + sum(abs2, grads.b3)
    return sqrt(total)
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
