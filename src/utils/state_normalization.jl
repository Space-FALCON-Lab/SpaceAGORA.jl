"""
    state_normalization.jl

Provides utilities for normalizing integration states in SpaceAGORA simulations.
This module handles normalization of different state vector components including:
- Quaternion normalization (enforces unit norm constraint)
- Position vector normalization (optional scaling)
- Velocity vector normalization (optional scaling)
- Full state normalization workflows

## Key Functions:
- `normalize_quaternion!` - Normalize quaternion component to unit norm
- `normalize_state_component!` - Normalize a vector state component
- `StateNormalizationConfig` - Configuration for normalization strategy
- `create_state_normalization_callback` - Generate normalization event callbacks
"""

module StateNormalization
    using LinearAlgebra
    using StaticArrays
    using DifferentialEquations

    export StateNormalizationConfig, NormalizationStrategy
    export normalize_quaternion!, normalize_state_component!, normalize_full_state!
    export create_state_normalization_callback, create_quaternion_normalization_event
    export QuaternionNormalizationStrategy

    """
    NormalizationStrategy

    Abstract type for different state normalization strategies.
    """
    abstract type NormalizationStrategy end

    """
    QuaternionNormalizationStrategy <: NormalizationStrategy

    Strategy for normalizing quaternion states to maintain unit-norm constraint.

    ## Fields:
    - `enabled::Bool` - Whether quaternion normalization is active
    - `tolerance::Float64` - Deviation from unit norm that triggers normalization
    - `normalize_every_step::Bool` - If true, normalize at every step; if false, normalize on-demand
    """
    struct QuaternionNormalizationStrategy <: NormalizationStrategy
        enabled::Bool
        tolerance::Float64
        normalize_every_step::Bool
    end

    """
    StateNormalizationConfig

    Configuration for state normalization behavior during integration.

    ## Fields:
    - `quaternion_strategy::QuaternionNormalizationStrategy` - Configuration for quaternion normalization
    - `verbose::Bool` - Print statistics about normalization events
    - `track_statistics::Bool` - Track normalization statistics during integration
    """
    struct StateNormalizationConfig
        quaternion_strategy::QuaternionNormalizationStrategy
        verbose::Bool
        track_statistics::Bool

        function StateNormalizationConfig(;
            quaternion_enabled::Bool=true,
            quaternion_tolerance::Float64=1e-10,
            quaternion_every_step::Bool=false,
            verbose::Bool=false,
            track_statistics::Bool=false
        )
            quat_strat = QuaternionNormalizationStrategy(
                quaternion_enabled,
                quaternion_tolerance,
                quaternion_every_step
            )
            new(quat_strat, verbose, track_statistics)
        end
    end

    """
    normalize_quaternion!(q::AbstractVector{<:Real})

    Normalize a quaternion vector to unit norm in-place.

    A quaternion q = [q_x, q_y, q_z, q_w] must satisfy ||q|| = 1.
    This function rescales the vector to enforce this constraint.

    ## Arguments:
    - `q::AbstractVector{<:Real}` - Quaternion vector of length 4, modified in-place

    ## Returns:
    - `norm::Float64` - Original norm before normalization (useful for diagnostics)

    ## Example:
    ```julia
    q = [0.1, 0.2, 0.3, 0.9]
    original_norm = normalize_quaternion!(q)
    @assert abs(norm(q) - 1.0) < 1e-14
    ```
    """
    function normalize_quaternion!(q::AbstractVector{<:Real})
        q_norm = norm(q)
        if q_norm > 0.0
            q ./= q_norm
        else
            # Degenerate case - set to identity quaternion
            q[1:3] .= 0.0
            q[4] = 1.0
        end
        return q_norm
    end

    """
    normalize_quaternion!(q::SVector{4, <:Real})

    Normalize a quaternion SVector to unit norm (returns new SVector).

    ## Arguments:
    - `q::SVector{4, <:Real}` - Quaternion vector of length 4

    ## Returns:
    - `q_normalized::SVector{4, Float64}` - Normalized quaternion
    - `original_norm::Float64` - Original norm (as tuple if unpacked as q_norm, original_norm = ...)
    """
    function normalize_quaternion!(q::SVector{4, T}) where T <: Real
        q_norm = norm(q)
        if q_norm > 0.0
            return q ./ q_norm, q_norm
        else
            return SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), q_norm
        end
    end

    """
    normalize_state_component!(component::AbstractVector{<:Real})

    Generic function to normalize any vector state component to unit norm.

    Useful for normalizing angular momentum vectors, unit direction vectors, etc.

    ## Arguments:
    - `component::AbstractVector{<:Real}` - Vector to normalize (modified in-place)

    ## Returns:
    - `original_norm::Float64` - Norm before normalization
    """
    function normalize_state_component!(component::AbstractVector{<:Real})
        comp_norm = norm(component)
        if comp_norm > 0.0
            component ./= comp_norm
        end
        return comp_norm
    end

    """
    normalize_full_state!(state, config::StateNormalizationConfig)

    Normalize all applicable components of the integration state vector.

    ## Arguments:
    - `state` - Integration state (ComponentVector or SVector, depending on implementation)
    - `config::StateNormalizationConfig` - Configuration specifying what to normalize

    ## Behavior:
    - If `config.quaternion_strategy.enabled` is true, normalizes state.q to unit norm
    - Stores statistics if `config.track_statistics` is true
    - Prints info if `config.verbose` is true
    """
    function normalize_full_state!(state, config::StateNormalizationConfig)
        if config.quaternion_strategy.enabled && haskey(state, :q)
            orig_norm = normalize_quaternion!(state.q)
            if config.verbose
                @info "Normalized quaternion from norm=$orig_norm to norm=$(norm(state.q))"
            end
        end
    end

    """
    create_quaternion_normalization_event(config::QuaternionNormalizationStrategy)

    Create a DiscreteCallback for quaternion normalization during integration.

    ## Arguments:
    - `config::QuaternionNormalizationStrategy` - Configuration for quaternion normalization

    ## Returns:
    - `callback::DiscreteCallback` - Event callback for use with DifferentialEquations solvers

    ## Usage:
    ```julia
    quat_config = QuaternionNormalizationStrategy(enabled=true, tolerance=1e-10)
    callback = create_quaternion_normalization_event(quat_config)
    solve(prob, <solver>, callback=callback)
    ```
    """
    function create_quaternion_normalization_event(config::QuaternionNormalizationStrategy)
        if !config.enabled
            return nothing
        end

        # Condition: check if quaternion deviates from unit norm
        function quaternion_condition(y, t, integrator)
            # Handle both ComponentVector and SVector cases
            q = if haskey(y, :q)
                y.q  # ComponentVector case
            else
                # Try accessing via indices - for legacy SVector case
                y[7:10]  # Assumes indices 7-10 contain quaternion
            end

            if config.normalize_every_step
                # Always trigger normalization at every step
                return true
            else
                # Only trigger if tolerance exceeded
                return abs(norm(q) - 1.0) > config.tolerance
            end
        end

        # Affect: normalize the quaternion
        function quaternion_affect!(integrator)
            if haskey(integrator.u, :q)
                normalize_quaternion!(integrator.u.q)
            else
                # Legacy case: normalize indices 7-10
                q_slice = @view integrator.u[7:10]
                normalize_quaternion!(q_slice)
            end
        end

        return DiscreteCallback(quaternion_condition, quaternion_affect!)
    end

    """
    create_state_normalization_callback(config::StateNormalizationConfig)

    Create a complete normalization callback for the integration state.

    ## Arguments:
    - `config::StateNormalizationConfig` - Full configuration for state normalization

    ## Returns:
    - `callback::Union{DiscreteCallback, Nothing}` - Event callback or nothing if disabled

    ## Example:
    ```julia
    norm_config = StateNormalizationConfig(
        quaternion_enabled=true,
        quaternion_tolerance=1e-11,
        quaternion_every_step=false,
        verbose=true,
        track_statistics=true
    )
    callback = create_state_normalization_callback(norm_config)
    ```
    """
    function create_state_normalization_callback(config::StateNormalizationConfig)
        quat_callback = create_quaternion_normalization_event(config.quaternion_strategy)
        return quat_callback
    end

    """
    get_state_deviation_metrics(state, reference_state=nothing)::Dict{String, Float64}

    Compute deviation metrics for the state to diagnose normalization issues.

    ## Arguments:
    - `state` - Current integration state
    - `reference_state` - Optional reference state for comparison

    ## Returns:
    - `metrics::Dict{String, Float64}` - Dictionary with deviation metrics
      - `quaternion_norm_deviation` - |norm(q) - 1.0|
      - `quaternion_norm` - norm(q)
      - `angular_velocity_norm` - norm(ω)
    """
    function get_state_deviation_metrics(state, reference_state=nothing)::Dict{String, Float64}
        metrics = Dict{String, Float64}()

        if haskey(state, :q)
            q_norm = norm(state.q)
            metrics["quaternion_norm"] = q_norm
            metrics["quaternion_norm_deviation"] = abs(q_norm - 1.0)
        end

        if haskey(state, :ω)
            metrics["angular_velocity_norm"] = norm(state.ω)
        end

        return metrics
    end

end  # module StateNormalization
