module OrbitalElementTrackingScheduling

using ForwardDiff
using JuMP
using LinearAlgebra
using Random
using SCIP

include("classical_element_sensitivity.jl")
using .ClassicalElementSensitivity: classical_element_jacobian

export wrap_to_pi,
       element_tracking_error,
       scaled_target,
       build_element_reachability,
       build_path_output_reachability,
       solve_miqp_tracking,
       classical_elements_from_rv,
       element_jacobian_inertial,
       element_sensitivity_rtn,
       element_jacobian_inertial_ad,
       element_sensitivity_rtn_ad,
       embed_target_element_matrix,
       build_toy_tracking_case

wrap_to_pi(angle::Real) = atan(sin(angle), cos(angle))

function element_tracking_error(
    nominal_elements::AbstractVector{<:Real},
    target_elements::AbstractVector{<:Real},
)
    length(nominal_elements) == 6 ||
        throw(ArgumentError("Nominal elements must contain [a, e, i, Omega, omega, nu]."))
    length(target_elements) == 6 ||
        throw(ArgumentError("Target elements must contain [a, e, i, Omega, omega, nu]."))
    all(isfinite, nominal_elements) && all(isfinite, target_elements) ||
        throw(ArgumentError("Orbital elements must be finite."))

    difference = zeros(Float64, 5)
    difference[1] = target_elements[1] - nominal_elements[1]
    difference[2] = target_elements[2] - nominal_elements[2]
    for index in 3:5
        difference[index] = wrap_to_pi(target_elements[index] - nominal_elements[index])
    end
    return difference
end

function scaled_target(
    element_difference::AbstractVector{<:Real},
    element_scales::AbstractVector{<:Real},
)
    !isempty(element_difference) || throw(ArgumentError("The tracking target cannot be empty."))
    length(element_scales) == length(element_difference) ||
        throw(ArgumentError("The tracking target and scales must have equal length."))
    all(scale -> isfinite(scale) && scale > 0.0, element_scales) ||
        throw(ArgumentError("Element scales must be finite and positive."))
    return Float64.(element_difference) ./ Float64.(element_scales)
end

function build_element_reachability(ad_list, bd_edge_list, scaled_terminal_output)
    interval_count = length(ad_list)
    interval_count > 0 || throw(ArgumentError("The dynamics history cannot be empty."))
    length(bd_edge_list) == interval_count ||
        throw(ArgumentError("A and B histories must have equal length."))

    output = Matrix{Float64}(scaled_terminal_output)
    output_count, state_count = size(output)
    output_count > 0 || throw(ArgumentError("The scaled output must have at least one row."))
    edge_count = size(bd_edge_list[1], 2)
    edge_count > 0 || throw(ArgumentError("At least one candidate edge is required."))
    reachability = zeros(Float64, output_count, edge_count, interval_count)
    sensitivity = copy(output)

    for interval in interval_count:-1:1
        ad = Matrix{Float64}(ad_list[interval])
        bd = Matrix{Float64}(bd_edge_list[interval])
        size(ad) == (state_count, state_count) ||
            throw(ArgumentError("Every A matrix must be square with the output state dimension."))
        size(bd) == (state_count, edge_count) ||
            throw(ArgumentError("Every B matrix must have a consistent state and edge dimension."))
        reachability[:, :, interval] .= sensitivity * bd
        sensitivity = sensitivity * ad
    end
    return reachability
end

function build_path_output_reachability(ad_list, bd_edge_list, scaled_path_outputs)
    interval_count = length(ad_list)
    interval_count > 0 || throw(ArgumentError("The dynamics history cannot be empty."))
    length(bd_edge_list) == interval_count ||
        throw(ArgumentError("A and B histories must have equal length."))
    length(scaled_path_outputs) == interval_count + 1 ||
        throw(ArgumentError("Path outputs are required at every state sample."))

    first_output = Matrix{Float64}(scaled_path_outputs[1])
    output_count, state_count = size(first_output)
    output_count > 0 || throw(ArgumentError("Path outputs must have at least one row."))
    edge_count = size(bd_edge_list[1], 2)
    edge_count > 0 || throw(ArgumentError("At least one candidate edge is required."))
    reachability = zeros(
        Float64,
        output_count,
        edge_count,
        interval_count,
        interval_count + 1,
    )

    for stage in 2:(interval_count + 1)
        sensitivity = Matrix{Float64}(scaled_path_outputs[stage])
        size(sensitivity) == (output_count, state_count) ||
            throw(ArgumentError("Every path output matrix must have consistent dimensions."))
        for interval in (stage - 1):-1:1
            ad = Matrix{Float64}(ad_list[interval])
            bd = Matrix{Float64}(bd_edge_list[interval])
            size(ad) == (state_count, state_count) ||
                throw(ArgumentError("Every A matrix must match the path-output state dimension."))
            size(bd) == (state_count, edge_count) ||
                throw(ArgumentError("Every B matrix must have a consistent state and edge dimension."))
            reachability[:, :, interval, stage] .= sensitivity * bd
            sensitivity = sensitivity * ad
        end
    end
    return reachability
end

function incident_edge_indices(edges, satellite)
    return [index for index in eachindex(edges) if satellite in edges[index]]
end

function project_terminal_limits!(values, feasible, edges, terminal_count)
    edge_count, interval_count = size(values)
    values .*= feasible
    for _ in 1:4, interval in 1:interval_count, satellite in eachindex(terminal_count)
        incident = incident_edge_indices(edges, satellite)
        limit = terminal_count[satellite]
        isempty(incident) && continue
        limit == 0 && (values[incident, interval] .= 0.0; continue)
        vector = values[incident, interval]
        sum(vector) <= limit && continue
        lower = minimum(vector) - 1.0
        upper = maximum(vector)
        for _ in 1:50
            threshold = (lower + upper) / 2.0
            if sum(clamp.(vector .- threshold, 0.0, 1.0)) > limit
                lower = threshold
            else
                upper = threshold
            end
        end
        values[incident, interval] .= clamp.(vector .- upper, 0.0, 1.0)
    end
    return values
end

function binary_tracking_warm_start(gamma, target, weights, feasible, edges, terminal_count)
    output_count, edge_count, interval_count = size(gamma)
    weighted_gamma = Diagonal(sqrt.(weights)) * reshape(gamma, output_count, :)
    weighted_target = sqrt.(weights) .* target
    lipschitz = opnorm(weighted_gamma)^2
    iszero(lipschitz) && return zeros(Float64, edge_count, interval_count)

    relaxed = zeros(Float64, edge_count, interval_count)
    accelerated = copy(relaxed)
    momentum = 1.0
    for _ in 1:200_000
        gradient = reshape(
            transpose(weighted_gamma) *
            (weighted_gamma * vec(accelerated) - weighted_target),
            edge_count,
            interval_count,
        )
        next_relaxed = project_terminal_limits!(
            clamp.(accelerated .- gradient ./ lipschitz, 0.0, 1.0),
            feasible,
            edges,
            terminal_count,
        )
        next_momentum = (1.0 + sqrt(1.0 + 4.0 * momentum^2)) / 2.0
        accelerated = next_relaxed .+
            ((momentum - 1.0) / next_momentum) .* (next_relaxed .- relaxed)
        relaxed = next_relaxed
        momentum = next_momentum
    end

    objective(command) = sum(
        weights .* (reshape(gamma, output_count, :) * vec(command) .- target) .^ 2,
    )
    best = zeros(Float64, edge_count, interval_count)
    best_objective = objective(best)
    rng = MersenneTwister(42)
    for _ in 1:500_000
        candidate = Float64.(rand(rng, edge_count, interval_count) .< relaxed) .* feasible
        for interval in 1:interval_count, satellite in eachindex(terminal_count)
            incident = incident_edge_indices(edges, satellite)
            while sum(candidate[incident, interval]) > terminal_count[satellite]
                active = [edge for edge in incident if candidate[edge, interval] == 1.0]
                scores = relaxed[active, interval] .* rand(rng, length(active))
                candidate[active[argmin(scores)], interval] = 0.0
            end
        end
        candidate_objective = objective(candidate)
        if candidate_objective < best_objective
            best .= candidate
            best_objective = candidate_objective
        end
    end
    return best
end

function schedule_is_feasible(command, feasible, edges, terminal_count)
    all(command .<= feasible) || return false
    for interval in axes(command, 2), satellite in eachindex(terminal_count)
        incident = incident_edge_indices(edges, satellite)
        sum(command[incident, interval]) <= terminal_count[satellite] || return false
    end
    return true
end

function local_binary_warm_start(
    gamma,
    target,
    weights,
    feasible,
    edges,
    terminal_count,
    reference,
    max_hamming_distance,
)
    output_count, edge_count, interval_count = size(gamma)
    current = reference .* feasible
    objective(command) = sum(
        weights .* (reshape(gamma, output_count, :) * vec(command) .- target) .^ 2,
    )
    hamming(command) = sum(command .!= reference)
    current_objective = objective(current)
    variable_count = edge_count * interval_count
    while hamming(current) < max_hamming_distance
        best = current
        best_objective = current_objective
        for first in 1:variable_count
            for second in first:variable_count
                candidate = copy(current)
                candidate[first] = 1.0 - candidate[first]
                if second != first
                    candidate[second] = 1.0 - candidate[second]
                end
                hamming(candidate) <= max_hamming_distance || continue
                schedule_is_feasible(candidate, feasible, edges, terminal_count) || continue
                candidate_objective = objective(candidate)
                if candidate_objective < best_objective
                    best = candidate
                    best_objective = candidate_objective
                end
            end
        end
        best_objective < current_objective || break
        current = best
        current_objective = best_objective
    end
    return current
end

function solve_miqp_tracking(
    reachability::AbstractArray{<:Real, 3},
    desired_output::AbstractVector{<:Real},
    feasible,
    edge_satellites;
    nterm = nothing,
    Qdiag = nothing,
    rho_active::Real = 1.0e-6,
    reference_u = nothing,
    max_hamming_distance::Union{Nothing, Int} = nothing,
    relative_gap::Real = 0.0,
    time_limit_s::Union{Nothing, Real} = nothing,
    simultaneous_edge_pairs = Tuple{Int, Int}[],
    balanced_edge_pairs = Tuple{Int, Int}[],
    rho_pair_imbalance::Real = 0.0,
    path_reachability = nothing,
    path_offset = nothing,
    Qdiag_path = nothing,
    path_bounds_scaled = nothing,
    path_integration_weights = nothing,
    Qdiag_integrated_path = nothing,
    integrated_path_bounds_scaled = nothing,
    silent::Bool = false,
)
    output_count, edge_count, interval_count = size(reachability)
    output_count > 0 || throw(ArgumentError("Reachability must have at least one output row."))
    length(desired_output) == output_count ||
        throw(ArgumentError("The desired output dimension does not match reachability."))
    size(feasible) == (edge_count, interval_count) ||
        throw(ArgumentError("Feasibility must be edge_count by interval_count."))
    length(edge_satellites) == edge_count ||
        throw(ArgumentError("One endpoint pair is required for every edge."))
    all(edge -> length(edge) == 2 && edge[1] != edge[2] && minimum(edge) >= 1, edge_satellites) ||
        throw(ArgumentError("Every edge must connect two distinct positive satellite indices."))
    all(value -> value == 0 || value == 1 || value == false || value == true, feasible) ||
        throw(ArgumentError("Feasibility entries must be binary."))
    all(pair -> 1 <= pair[1] <= edge_count && 1 <= pair[2] <= edge_count && pair[1] != pair[2], simultaneous_edge_pairs) ||
        throw(ArgumentError("Simultaneous edge-pair indices must identify two distinct candidate edges."))
    all(pair -> 1 <= pair[1] <= edge_count && 1 <= pair[2] <= edge_count && pair[1] != pair[2], balanced_edge_pairs) ||
        throw(ArgumentError("Balanced edge-pair indices must identify two distinct candidate edges."))
    all(isfinite, desired_output) && all(isfinite, reachability) ||
        throw(ArgumentError("Reachability and the desired output must be finite."))

    satellite_count = maximum(maximum(edge) for edge in edge_satellites)
    terminal_count = isnothing(nterm) ? ones(Int, satellite_count) : Int.(nterm)
    length(terminal_count) >= satellite_count ||
        throw(ArgumentError("Terminal counts do not cover every satellite in the edge list."))
    all(>=(0), terminal_count) || throw(ArgumentError("Terminal counts must be nonnegative."))

    weights = isnothing(Qdiag) ? ones(Float64, output_count) : Float64.(Qdiag)
    length(weights) == output_count ||
        throw(ArgumentError("One tracking weight is required per terminal output."))
    all(value -> isfinite(value) && value >= 0.0, weights) ||
        throw(ArgumentError("Tracking weights must be finite and nonnegative."))
    isfinite(rho_active) && rho_active >= 0.0 ||
        throw(ArgumentError("Active-link regularization must be finite and nonnegative."))
    isfinite(rho_pair_imbalance) && rho_pair_imbalance >= 0.0 ||
        throw(ArgumentError("Pair-imbalance regularization must be finite and nonnegative."))
    isfinite(relative_gap) && relative_gap >= 0.0 ||
        throw(ArgumentError("Relative MIP gap must be finite and nonnegative."))
    isnothing(time_limit_s) || (isfinite(time_limit_s) && time_limit_s > 0.0) ||
        throw(ArgumentError("Time limit must be finite and positive."))
    isnothing(reference_u) == isnothing(max_hamming_distance) || throw(ArgumentError(
        "reference_u and max_hamming_distance must be provided together.",
    ))
    reference = if isnothing(reference_u)
        nothing
    else
        size(reference_u) == (edge_count, interval_count) || throw(ArgumentError(
            "The reference schedule must be edge_count by interval_count.",
        ))
        all(value -> value == 0 || value == 1, reference_u) || throw(ArgumentError(
            "The reference schedule must be binary.",
        ))
        max_hamming_distance >= 0 || throw(ArgumentError(
            "The maximum Hamming distance must be nonnegative.",
        ))
        Float64.(reference_u)
    end

    path_gamma = isnothing(path_reachability) ? nothing : Float64.(path_reachability)
    path_output_count = 0
    path_stage_count = 0
    path_baseline = zeros(Float64, 0, 0)
    path_weights = Float64[]
    path_bounds = nothing
    integration_weights = nothing
    integrated_path_weights = Float64[]
    integrated_path_bounds = nothing
    if !isnothing(path_gamma)
        ndims(path_gamma) == 4 ||
            throw(ArgumentError("Path reachability must have four dimensions."))
        path_output_count, path_edge_count, path_interval_count, path_stage_count =
            size(path_gamma)
        path_output_count > 0 && path_stage_count > 0 ||
            throw(ArgumentError("Path reachability must contain outputs and stages."))
        path_edge_count == edge_count && path_interval_count == interval_count ||
            throw(ArgumentError("Path and terminal reachability must share edge and interval dimensions."))
        all(isfinite, path_gamma) ||
            throw(ArgumentError("Path reachability must be finite."))
        path_baseline = isnothing(path_offset) ?
            zeros(Float64, path_output_count, path_stage_count) : Float64.(path_offset)
        size(path_baseline) == (path_output_count, path_stage_count) ||
            throw(ArgumentError("Path offset dimensions must match path outputs and stages."))
        all(isfinite, path_baseline) || throw(ArgumentError("Path offsets must be finite."))
        path_weights = isnothing(Qdiag_path) ?
            ones(Float64, path_output_count) : Float64.(Qdiag_path)
        length(path_weights) == path_output_count ||
            throw(ArgumentError("One path weight is required per path output."))
        all(value -> isfinite(value) && value >= 0.0, path_weights) ||
            throw(ArgumentError("Path weights must be finite and nonnegative."))
        if !isnothing(path_bounds_scaled)
            path_bounds = Float64.(path_bounds_scaled)
            length(path_bounds) == path_output_count ||
                throw(ArgumentError("One path bound is required per path output."))
            all(value -> isfinite(value) && value > 0.0, path_bounds) ||
                throw(ArgumentError("Path bounds must be finite and positive."))
        end
        if !isnothing(path_integration_weights)
            integration_weights = Float64.(path_integration_weights)
            length(integration_weights) == path_stage_count || throw(ArgumentError(
                "One integration weight is required per path stage.",
            ))
            all(value -> isfinite(value) && value >= 0.0, integration_weights) ||
                throw(ArgumentError("Path integration weights must be finite and nonnegative."))
            isapprox(sum(integration_weights), 1.0; atol = 1.0e-12, rtol = 1.0e-12) ||
                throw(ArgumentError("Path integration weights must sum to one."))
            integrated_path_weights = isnothing(Qdiag_integrated_path) ?
                zeros(Float64, path_output_count) : Float64.(Qdiag_integrated_path)
            length(integrated_path_weights) == path_output_count || throw(ArgumentError(
                "One integrated-path weight is required per path output.",
            ))
            all(value -> isfinite(value) && value >= 0.0, integrated_path_weights) ||
                throw(ArgumentError("Integrated-path weights must be finite and nonnegative."))
            if !isnothing(integrated_path_bounds_scaled)
                integrated_path_bounds = Float64.(integrated_path_bounds_scaled)
                length(integrated_path_bounds) == path_output_count || throw(ArgumentError(
                    "One integrated-path bound is required per path output.",
                ))
                all(value -> isfinite(value) && value > 0.0, integrated_path_bounds) ||
                    throw(ArgumentError("Integrated-path bounds must be finite and positive."))
            end
        elseif !isnothing(Qdiag_integrated_path) || !isnothing(integrated_path_bounds_scaled)
            throw(ArgumentError("Integrated-path options require path integration weights."))
        end
    elseif !isnothing(path_offset) || !isnothing(Qdiag_path) || !isnothing(path_bounds_scaled) ||
           !isnothing(path_integration_weights) || !isnothing(Qdiag_integrated_path) ||
           !isnothing(integrated_path_bounds_scaled)
        throw(ArgumentError("Path options require path reachability."))
    end

    gamma = Float64.(reachability)
    target = Float64.(desired_output)
    model = Model(SCIP.Optimizer)
    silent && set_silent(model)
    set_optimizer_attribute(model, "numerics/feastol", 1.0e-10)
    set_optimizer_attribute(model, "limits/gap", float(relative_gap))
    set_optimizer_attribute(model, "limits/absgap", 0.0)
    isnothing(time_limit_s) ||
        set_optimizer_attribute(model, "limits/time", float(time_limit_s))
    @variable(model, u[1:edge_count, 1:interval_count], Bin)
    warm_start = if isnothing(reference)
        binary_tracking_warm_start(
            gamma,
            target,
            weights,
            feasible,
            edge_satellites,
            terminal_count,
        )
    else
        local_binary_warm_start(
            gamma,
            target,
            weights,
            feasible,
            edge_satellites,
            terminal_count,
            reference,
            max_hamming_distance,
        )
    end
    warm_start_output = reshape(gamma, output_count, :) * vec(warm_start)
    warm_start_tracking_cost = sum(weights .* (warm_start_output .- target) .^ 2)
    set_start_value.(u, warm_start)

    for edge in 1:edge_count, interval in 1:interval_count
        @constraint(model, u[edge, interval] <= feasible[edge, interval])
    end
    for interval in 1:interval_count, satellite in 1:satellite_count
        incident = incident_edge_indices(edge_satellites, satellite)
        @constraint(
            model,
            sum(u[edge, interval] for edge in incident) <= terminal_count[satellite],
        )
    end
    for (first_edge, second_edge) in simultaneous_edge_pairs, interval in 1:interval_count
        @constraint(model, u[first_edge, interval] == u[second_edge, interval])
    end
    if !isnothing(reference)
        @constraint(
            model,
            sum(
                reference[edge, interval] == 1.0 ?
                    1.0 - u[edge, interval] : u[edge, interval]
                for edge in 1:edge_count, interval in 1:interval_count
            ) <= max_hamming_distance,
        )
    end

    @expression(
        model,
        terminal_output[output_index = 1:output_count],
        sum(
            gamma[output_index, edge, interval] * u[edge, interval]
            for edge in 1:edge_count, interval in 1:interval_count
        ),
    )
    @expression(
        model,
        tracking_cost,
        sum(
            weights[output_index] * (terminal_output[output_index] - target[output_index])^2
            for output_index in 1:output_count
        ),
    )
    @expression(
        model,
        active_cost,
        float(rho_active) * sum(u[edge, interval] for edge in 1:edge_count, interval in 1:interval_count),
    )
    @expression(
        model,
        pair_imbalance_cost,
        float(rho_pair_imbalance) * sum(
            (u[first_edge, interval] - u[second_edge, interval])^2
            for (first_edge, second_edge) in balanced_edge_pairs, interval in 1:interval_count
        ),
    )
    eccentricity_path_cost = @expression(model, 0.0)
    integrated_path_cost = @expression(model, 0.0)
    if !isnothing(path_gamma)
        @expression(
            model,
            path_output[component = 1:path_output_count, stage = 1:path_stage_count],
            path_baseline[component, stage] + sum(
                path_gamma[component, edge, interval, stage] * u[edge, interval]
                for edge in 1:edge_count, interval in 1:interval_count
            ),
        )
        eccentricity_path_cost = @expression(
            model,
            sum(
                path_weights[component] * path_output[component, stage]^2
                for component in 1:path_output_count, stage in 1:path_stage_count
            ),
        )
        if !isnothing(path_bounds)
            for component in 1:path_output_count, stage in 1:path_stage_count
                @constraint(model, -path_bounds[component] <= path_output[component, stage])
                @constraint(model, path_output[component, stage] <= path_bounds[component])
            end
        end
        if !isnothing(integration_weights)
            @expression(
                model,
                integrated_path_output[component = 1:path_output_count],
                sum(
                    integration_weights[stage] * path_output[component, stage]
                    for stage in 1:path_stage_count
                ),
            )
            integrated_path_cost = @expression(
                model,
                sum(
                    integrated_path_weights[component] * integrated_path_output[component]^2
                    for component in 1:path_output_count
                ),
            )
            if !isnothing(integrated_path_bounds)
                for component in 1:path_output_count
                    @constraint(
                        model,
                        -integrated_path_bounds[component] <= integrated_path_output[component],
                    )
                    @constraint(
                        model,
                        integrated_path_output[component] <= integrated_path_bounds[component],
                    )
                end
            end
        end
    end
    objective_multiplier = 1.0
    @objective(
        model,
        Min,
        objective_multiplier * (
            tracking_cost + eccentricity_path_cost + integrated_path_cost + active_cost +
            pair_imbalance_cost
        ),
    )
    optimize!(model)
    has_values(model) || error(
        "Orbital-element tracking MIQP returned no feasible incumbent " *
        "(termination=$(termination_status(model)), primal=$(primal_status(model))).",
    )

    u_value = round.(Float64, value.(u))
    output_value = reshape(gamma, output_count, :) * vec(u_value)
    tracking_cost_value = sum(weights .* (output_value .- target) .^ 2)
    active_cost_value = float(rho_active) * sum(u_value)
    pair_imbalance_cost_value = float(rho_pair_imbalance) * sum(
        (
            (u_value[first_edge, interval] - u_value[second_edge, interval])^2
            for (first_edge, second_edge) in balanced_edge_pairs, interval in 1:interval_count
        );
        init = 0.0,
    )
    path_output_value = zeros(Float64, path_output_count, path_stage_count)
    eccentricity_path_cost_value = 0.0
    integrated_path_output_value = Float64[]
    integrated_path_cost_value = 0.0
    if !isnothing(path_gamma)
        for stage in 1:path_stage_count
            path_output_value[:, stage] .= path_baseline[:, stage] .+
                reshape(path_gamma[:, :, :, stage], path_output_count, :) * vec(u_value)
        end
        eccentricity_path_cost_value = sum(
            path_weights .* vec(sum(abs2, path_output_value; dims = 2)),
        )
        if !isnothing(integration_weights)
            integrated_path_output_value = path_output_value * integration_weights
            integrated_path_cost_value = sum(
                integrated_path_weights .* integrated_path_output_value .^ 2,
            )
        end
    end
    solver_objective_value = try
        objective_value(model)
    catch
        tracking_cost_value + eccentricity_path_cost_value + integrated_path_cost_value +
            active_cost_value + pair_imbalance_cost_value
    end
    solver_objective_bound = try
        objective_bound(model)
    catch
        NaN
    end
    solver_relative_gap = try
        JuMP.relative_gap(model)
    catch
        Inf
    end
    return (
        model = model,
        status = termination_status(model),
        primal_status = primal_status(model),
        u = u_value,
        warm_start = warm_start,
        warm_start_output = warm_start_output,
        warm_start_tracking_cost = warm_start_tracking_cost,
        output = output_value,
        desired_output = target,
        tracking_error = output_value .- target,
        tracking_cost = tracking_cost_value,
        path_output = path_output_value,
        eccentricity_path_cost = eccentricity_path_cost_value,
        integrated_path_output = integrated_path_output_value,
        integrated_path_cost = integrated_path_cost_value,
        active_cost = active_cost_value,
        pair_imbalance_cost = pair_imbalance_cost_value,
        objective = tracking_cost_value + eccentricity_path_cost_value +
            integrated_path_cost_value + active_cost_value + pair_imbalance_cost_value,
        solver_objective = solver_objective_value,
        objective_bound = solver_objective_bound,
        relative_gap = solver_relative_gap,
        objective_multiplier = objective_multiplier,
    )
end

function classical_elements_from_rv(x::AbstractVector{T}; mu::Real = 1.0) where {T<:Real}
    length(x) == 6 || throw(ArgumentError("A Cartesian state must contain three positions and three velocities."))
    mu > 0.0 || throw(ArgumentError("The gravitational parameter must be positive."))
    mu_t = one(T) * mu
    position = x[1:3]
    velocity = x[4:6]
    radius = norm(position)
    speed = norm(velocity)
    angular_momentum_vector = cross(position, velocity)
    angular_momentum = norm(angular_momentum_vector)
    inertial_normal = [zero(T), zero(T), one(T)]
    node_vector = cross(inertial_normal, angular_momentum_vector)
    eccentricity_vector =
        ((speed^2 - mu_t / radius) .* position .- dot(position, velocity) .* velocity) ./ mu_t
    eccentricity = norm(eccentricity_vector)
    energy = 0.5 * speed^2 - mu_t / radius
    semimajor_axis = -mu_t / (2.0 * energy)
    inclination = acos(clamp(angular_momentum_vector[3] / angular_momentum, -one(T), one(T)))
    raan = atan(node_vector[2], node_vector[1])
    angular_momentum_direction = angular_momentum_vector / angular_momentum
    argument_of_periapsis = atan(
        dot(cross(node_vector, eccentricity_vector), angular_momentum_direction),
        dot(node_vector, eccentricity_vector),
    )
    true_anomaly = atan(
        dot(cross(eccentricity_vector, position), angular_momentum_direction),
        dot(eccentricity_vector, position),
    )
    return [semimajor_axis, eccentricity, inclination, raan, argument_of_periapsis, true_anomaly]
end

element_jacobian_inertial(position, velocity; mu=1.0) =
    classical_element_jacobian(position, velocity; mu=mu, include_true_anomaly=true)

function element_sensitivity_rtn(position, velocity, rtn_frame; mu=1.0)
    rotation = [rtn_frame zeros(3, 3); zeros(3, 3) rtn_frame]
    return element_jacobian_inertial(position, velocity; mu=mu) * rotation
end

# Automatic differentiation is retained only as an independent validation reference.
function element_jacobian_inertial_ad(
    position::AbstractVector{<:Real},
    velocity::AbstractVector{<:Real};
    mu::Real = 1.0,
)
    length(position) == 3 && length(velocity) == 3 ||
        throw(ArgumentError("Position and velocity must each have three components."))
    state = vcat(Float64.(position), Float64.(velocity))
    return ForwardDiff.jacobian(local_state -> classical_elements_from_rv(local_state; mu = mu), state)
end

function element_sensitivity_rtn_ad(
    position::AbstractVector{<:Real},
    velocity::AbstractVector{<:Real},
    rtn_frame::AbstractMatrix{<:Real};
    mu::Real = 1.0,
)
    size(rtn_frame) == (3, 3) || throw(ArgumentError("The RTN frame must be 3 by 3."))
    inertial_jacobian = element_jacobian_inertial_ad(position, velocity; mu = mu)
    rtn_to_inertial = zeros(Float64, 6, 6)
    rtn_to_inertial[1:3, 1:3] .= rtn_frame
    rtn_to_inertial[4:6, 4:6] .= rtn_frame
    return inertial_jacobian * rtn_to_inertial
end

function embed_target_element_matrix(
    satellite_output::AbstractMatrix{<:Real},
    target_satellite::Int,
    satellite_count::Int,
)
    satellite_count >= 1 || throw(ArgumentError("At least one satellite is required."))
    1 <= target_satellite <= satellite_count ||
        throw(ArgumentError("The target satellite index is invalid."))
    size(satellite_output, 2) == 6 ||
        throw(ArgumentError("A satellite element output matrix must have six columns."))
    output = zeros(Float64, size(satellite_output, 1), 6 * satellite_count)
    columns = (6 * (target_satellite - 1) + 1):(6 * target_satellite)
    output[:, columns] .= satellite_output
    return output
end

function build_toy_tracking_case()
    nominal_elements = [
        7378.1363,
        0.0010,
        deg2rad(30.0),
        deg2rad(40.0),
        deg2rad(10.0),
        deg2rad(0.0),
    ]
    target_elements = [
        nominal_elements[1] + 100.0,
        nominal_elements[2] + 2.0e-5,
        nominal_elements[3] + deg2rad(0.5),
        nominal_elements[4] + deg2rad(0.6),
        nominal_elements[5] - deg2rad(0.5),
        nominal_elements[6] + deg2rad(120.0),
    ]
    element_scales = [
        100.0,
        2.0e-5,
        deg2rad(0.5),
        deg2rad(0.6),
        deg2rad(0.5),
    ]
    desired_output = scaled_target(
        element_tracking_error(nominal_elements, target_elements),
        element_scales,
    )
    edges = [(2, 1), (3, 1), (1, 2), (1, 3), (2, 3)]
    interval_count = 5
    reachability = zeros(Float64, 5, length(edges), interval_count)
    reachability[1, 1, 1] = 1.0
    reachability[2, 2, 2] = 1.0
    reachability[3, 3, 3] = 1.0
    reachability[4, 4, 4] = 1.0
    reachability[5, 5, 5] = -1.0
    known_schedule = Matrix{Float64}(I, length(edges), interval_count)
    return (
        nominal_elements = nominal_elements,
        target_elements = target_elements,
        element_scales = element_scales,
        desired_output = desired_output,
        edges = edges,
        reachability = reachability,
        feasible = ones(Int, length(edges), interval_count),
        nterm = ones(Int, 3),
        known_schedule = known_schedule,
    )
end

end
