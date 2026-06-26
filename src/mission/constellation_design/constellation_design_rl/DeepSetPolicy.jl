module DeepSetPolicy

using Flux
using LinearAlgebra
using Random

"""
    DualDeepSetPolicy

Dual DeepSet policy network for satellite seeding.
Separate encoders for satellites and clients with cross-attention.
"""
struct DualDeepSetPolicy
    # Satellite encoder (learns constellation patterns)
    satellite_φ::Chain      # Per-satellite: (6,) → (hidden,)
    satellite_ψ::Chain      # Global: (hidden,) → (embed,)
    
    # Client encoder (learns debris distribution patterns)  
    client_φ::Chain         # Per-client: (6,) → (hidden,)
    client_ψ::Chain         # Global: (hidden,) → (embed,)
    
    # Cross-attention (learns satellite-client relationships)
    cross_attention::Chain  # Combines sat+client embeddings
    
    # Policy head
    policy_head::Chain      # Output: 6 orbital elements (normalized [-1,1])
    
    # Configuration
    hidden_size::Int
    embed_size::Int
    device::Any
end

"""
    DualDeepSetPolicy(; hidden_size=32, embed_size=16, device=cpu)

Construct a DualDeepSetPolicy with configurable layer sizes.
"""
function DualDeepSetPolicy(; hidden_size::Int=32, embed_size::Int=16, device=cpu)
    # Satellite encoder
    satellite_φ = Chain(
        Dense(6, hidden_size, relu),
        Dense(hidden_size, hidden_size, relu)
    ) |> device
    
    satellite_ψ = Chain(
        Dense(hidden_size, embed_size, relu)
    ) |> device
    
    # Client encoder
    client_φ = Chain(
        Dense(6, hidden_size, relu),
        Dense(hidden_size, hidden_size, relu)
    ) |> device
    
    client_ψ = Chain(
        Dense(hidden_size, embed_size, relu)
    ) |> device
    
    # Cross-attention
    cross_attention = Chain(
        Dense(2 * embed_size, hidden_size, relu),
        Dense(hidden_size, hidden_size, relu)
    ) |> device
    
    # Policy head
    policy_head = Chain(
        Dense(hidden_size, hidden_size, relu),
        Dense(hidden_size, 6, tanh)  # Output normalized to [-1, 1]
    ) |> device
    
    return DualDeepSetPolicy(
        satellite_φ, satellite_ψ,
        client_φ, client_ψ,
        cross_attention,
        policy_head,
        hidden_size, embed_size, device
    )
end

"""
    (policy::DualDeepSetPolicy)(satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64}) -> Vector{Float64}

Forward pass through the DualDeepSetPolicy.
Returns action vector (6 orbital elements normalized to [-1, 1]).
"""
function (policy::DualDeepSetPolicy)(satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64})
    # Encode satellites
    n_sats = size(satellite_orbitals, 1)
    if n_sats > 0
        sat_features = hcat([policy.satellite_φ(satellite_orbitals[i:i, :]') for i in 1:n_sats]...)
        sat_embedding = policy.satellite_ψ(mean(sat_features, dims=2))
    else
        sat_embedding = zeros(policy.embed_size, 1) |> policy.device
    end
    
    # Encode clients
    n_clients = size(client_orbitals, 1)
    if n_clients > 0
        client_features = hcat([policy.client_φ(client_orbitals[i:i, :]') for i in 1:n_clients]...)
        client_embedding = policy.client_ψ(mean(client_features, dims=2))
    else
        client_embedding = zeros(policy.embed_size, 1) |> policy.device
    end
    
    # Cross-attention
    combined = vcat(sat_embedding, client_embedding)
    cross_features = policy.cross_attention(combined)
    
    # Policy head
    action = policy.policy_head(cross_features)
    
    return vec(action)
end

"""
    get_embeddings(policy::DualDeepSetPolicy, satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64}) -> NamedTuple

Extract intermediate embeddings for visualization and logging.
"""
function get_embeddings(policy::DualDeepSetPolicy, satellite_orbitals::Matrix{Float64}, client_orbitals::Matrix{Float64})
    # Encode satellites
    n_sats = size(satellite_orbitals, 1)
    if n_sats > 0
        sat_φ_output = hcat([policy.satellite_φ(satellite_orbitals[i:i, :]') for i in 1:n_sats]...)
        sat_embedding = policy.satellite_ψ(mean(sat_φ_output, dims=2))
    else
        sat_φ_output = zeros(policy.hidden_size, 1) |> policy.device
        sat_embedding = zeros(policy.embed_size, 1) |> policy.device
    end
    
    # Encode clients
    n_clients = size(client_orbitals, 1)
    if n_clients > 0
        client_φ_output = hcat([policy.client_φ(client_orbitals[i:i, :]') for i in 1:n_clients]...)
        client_embedding = policy.client_ψ(mean(client_φ_output, dims=2))
    else
        client_φ_output = zeros(policy.hidden_size, 1) |> policy.device
        client_embedding = zeros(policy.embed_size, 1) |> policy.device
    end
    
    # Cross-attention
    combined = vcat(sat_embedding, client_embedding)
    cross_attn_output = policy.cross_attention(combined)
    
    return (
        sat_φ_output = sat_φ_output,
        client_φ_output = client_φ_output,
        cross_attn_output = cross_attn_output,
        sat_embedding = sat_embedding,
        client_embedding = client_embedding
    )
end

"""
    check_gpu_availability() -> Bool, Any

Check if GPU (CUDA) is available and return device.
Returns (is_available, device).
"""
function check_gpu_availability()
    try
        if @isdefined(CUDA) && CUDA.has_cuda()
            @info "GPU detected, using CUDA for training"
            return true, gpu
        else
            @warn "GPU not available, falling back to CPU training"
            return false, cpu
        end
    catch
        @warn "CUDA not available, falling back to CPU training"
        return false, cpu
    end
end

export DualDeepSetPolicy, check_gpu_availability, get_embeddings

end # module DeepSetPolicy
