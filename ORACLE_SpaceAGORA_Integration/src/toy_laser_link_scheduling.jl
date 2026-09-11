module ToyLaserLinkScheduling

using LinearAlgebra
using JuMP
using HiGHS

export sat_base,
       rR,
       rT,
       rN,
       vR,
       vT,
       vN,
       build_toy_problem,
       transition_product,
       compute_edge_time_gains,
       compute_gramian_scores,
       known_solution,
       optimize_signed_terminal_change,
       optimize_signed_semimajor_axis_increase,
       optimize_gramian_schedule,
       simulate_schedule,
       directed_pair_geometry_data,
       edge_input_matrix_rtn,
       normalized_semimajor_axis,
       semimajor_axis_sensitivity_rtn,
       embed_satellite_output_row

sat_base(i::Int) = 6 * (i - 1)
rR(i::Int) = sat_base(i) + 1
rT(i::Int) = sat_base(i) + 2
rN(i::Int) = sat_base(i) + 3
vR(i::Int) = sat_base(i) + 4
vT(i::Int) = sat_base(i) + 5
vN(i::Int) = sat_base(i) + 6

function unitvec(n::Int, index::Int)
    vector = zeros(Float64, n)
    vector[index] = 1.0
    return vector
end

function build_toy_problem(; H::Int = 5)
    H >= 1 || throw(ArgumentError("Horizon must contain at least one step."))

    nsat = 3
    nx = 6 * nsat
    edges = [(1, 2), (1, 3), (2, 3)]
    edge_count = length(edges)
    ad_list = [Matrix{Float64}(I, nx, nx) for _ in 1:H]

    b12 = unitvec(nx, vT(1)) - unitvec(nx, vT(2))
    b13 = unitvec(nx, vR(1)) - unitvec(nx, vR(3))
    b23 = unitvec(nx, vT(2)) - unitvec(nx, vT(3))
    bd = hcat(b12, b13, b23)
    bd_list = [copy(bd) for _ in 1:H]

    ca = zeros(Float64, 1, nx)
    ca[1, rR(1)] = 2.0
    ca[1, vT(1)] = 2.0

    return (
        Ns = nsat,
        nx = nx,
        Ne = edge_count,
        H = H,
        edges = edges,
        Ad_list = ad_list,
        Bd_list = bd_list,
        Ca = ca,
        feasible = ones(Int, edge_count, H),
        nterm = ones(Int, nsat),
    )
end

function transition_product(ad_list, kf::Int, jplus::Int)
    isempty(ad_list) && throw(ArgumentError("Dynamics history cannot be empty."))
    1 <= jplus <= kf || throw(ArgumentError("Require 1 <= jplus <= kf."))
    kf - 1 <= length(ad_list) || throw(ArgumentError("Transition exceeds dynamics history."))
    nx = size(ad_list[1], 1)
    transition = Matrix{Float64}(I, nx, nx)
    for index in jplus:(kf - 1)
        transition = Matrix{Float64}(ad_list[index]) * transition
    end
    return transition
end

function compute_edge_time_gains(ad_list, bd_list, terminal_output, horizon::Int)
    length(ad_list) >= horizon || throw(ArgumentError("A history is shorter than the horizon."))
    length(bd_list) >= horizon || throw(ArgumentError("B history is shorter than the horizon."))
    terminal_state_index = horizon + 1
    edge_count = size(bd_list[1], 2)
    gains = zeros(Float64, edge_count, horizon)
    for step in 1:horizon
        transition = transition_product(ad_list, terminal_state_index, step + 1)
        for edge in 1:edge_count
            gains[edge, step] = dot(
                vec(terminal_output),
                transition * bd_list[step][:, edge],
            )
        end
    end
    return gains
end

compute_gramian_scores(gains::AbstractMatrix{<:Real}) = Float64.(gains) .^ 2

function known_solution(problem)
    command = zeros(Float64, problem.Ne, problem.H)
    command[1, :] .= 1.0
    objective = 2.0 * problem.H
    return command, objective
end

function incident_edges(edges, satellite)
    return [edge for edge in eachindex(edges) if satellite in edges[edge]]
end

function optimize_signed_terminal_change(
    gains,
    feasible,
    edges;
    nterm,
    desired_direction::Real,
)
    edge_count, horizon = size(gains)
    isfinite(desired_direction) && !iszero(desired_direction) ||
        throw(ArgumentError("Desired direction must be a finite nonzero scalar."))
    direction = sign(float(desired_direction))
    directed_gains = direction .* Float64.(gains)
    objective_scale = maximum(abs, directed_gains)
    objective_scale = iszero(objective_scale) ? 1.0 : objective_scale
    scaled_gains = directed_gains ./ objective_scale
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, command[1:edge_count, 1:horizon], Bin)
    @objective(
        model,
        Max,
        sum(scaled_gains[edge, step] * command[edge, step] for edge in 1:edge_count, step in 1:horizon),
    )
    for edge in 1:edge_count, step in 1:horizon
        @constraint(model, command[edge, step] <= feasible[edge, step])
    end
    for step in 1:horizon, satellite in eachindex(nterm)
        incident = incident_edges(edges, satellite)
        @constraint(
            model,
            sum(command[edge, step] for edge in incident) <= nterm[satellite],
        )
    end
    optimize!(model)
    is_solved_and_feasible(model) || error("Signed scheduling MILP did not solve to feasibility.")
    command_value = value.(command)
    signed_change = sum(gains .* command_value)
    signed_change = iszero(signed_change) ? 0.0 : signed_change
    return (
        u = command_value,
        objective = signed_change,
        directed_objective = direction * signed_change,
        desired_direction = direction,
        status = termination_status(model),
    )
end

function optimize_signed_semimajor_axis_increase(
    gains,
    feasible,
    edges;
    nterm,
)
    return optimize_signed_terminal_change(
        gains,
        feasible,
        edges;
        nterm = nterm,
        desired_direction = 1.0,
    )
end

function optimize_gramian_schedule(scores, feasible, edges; nterm)
    edge_count, horizon = size(scores)
    objective_scale = maximum(abs, scores)
    objective_scale = iszero(objective_scale) ? 1.0 : objective_scale
    scaled_scores = scores ./ objective_scale
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, command[1:edge_count, 1:horizon], Bin)
    @objective(
        model,
        Max,
        sum(scaled_scores[edge, step] * command[edge, step] for edge in 1:edge_count, step in 1:horizon),
    )
    for edge in 1:edge_count, step in 1:horizon
        @constraint(model, command[edge, step] <= feasible[edge, step])
    end
    for step in 1:horizon, satellite in eachindex(nterm)
        incident = incident_edges(edges, satellite)
        @constraint(
            model,
            sum(command[edge, step] for edge in incident) <= nterm[satellite],
        )
    end
    optimize!(model)
    is_solved_and_feasible(model) || error("Gramian scheduling MILP did not solve to feasibility.")
    command_value = value.(command)
    return (
        u = command_value,
        objective = sum(scores .* command_value),
        status = termination_status(model),
    )
end

function simulate_schedule(problem, command::AbstractMatrix{<:Real})
    size(command) == (problem.Ne, problem.H) ||
        throw(ArgumentError("Command matrix must be Ne by H."))
    states = zeros(Float64, problem.nx, problem.H + 1)
    for step in 1:problem.H
        states[:, step + 1] .=
            problem.Ad_list[step] * states[:, step] +
            problem.Bd_list[step] * command[:, step]
    end
    return states
end

function directed_pair_geometry_data(pair_data, directed_edges)
    return [
        begin
            i, j = directed_edge
            i != j || throw(ArgumentError("A directed edge must connect distinct satellites."))
            source_index = findfirst(
                data -> data.pair == (i, j) || data.pair == (j, i),
                pair_data,
            )
            isnothing(source_index) && throw(ArgumentError("No geometry exists for directed edge $(directed_edge)."))
            source = pair_data[source_index]
            orientation = source.pair == (i, j) ? 1.0 : -1.0
            merge(
                source,
                (
                    pair = (i, j),
                    khat = orientation .* Float64.(source.khat),
                ),
            )
        end for directed_edge in directed_edges
    ]
end

function edge_input_matrix_rtn(
    frames,
    pair_data;
    acceleration_gain::Real,
    acceleration_scale::Real,
)
    nsat = length(frames)
    acceleration_scale > 0.0 || throw(ArgumentError("Acceleration scale must be positive."))
    all(size(frame) == (3, 3) for frame in frames) ||
        throw(ArgumentError("Every RTN frame must be 3 by 3."))
    input_matrix = zeros(Float64, 6 * nsat, length(pair_data))
    for (edge, data) in enumerate(pair_data)
        i, j = data.pair
        direction_bar =
            (float(data.zeta) * float(acceleration_gain) / float(acceleration_scale)) .*
            Float64.(data.khat)
        velocity_rows_i = (6 * (i - 1) + 4):(6 * i)
        velocity_rows_j = (6 * (j - 1) + 4):(6 * j)
        # Treat (i,j) as a directed beam: i recoils and receiver j is pushed along k_ij.
        input_matrix[velocity_rows_i, edge] .= -transpose(frames[i]) * direction_bar
        input_matrix[velocity_rows_j, edge] .= transpose(frames[j]) * direction_bar
    end
    return input_matrix
end

function normalized_semimajor_axis(rbar::AbstractVector, vbar::AbstractVector)
    energy = 0.5 * dot(vbar, vbar) - 1.0 / norm(rbar)
    return -1.0 / (2.0 * energy)
end

function semimajor_axis_sensitivity_rtn(
    rbar::AbstractVector,
    vbar::AbstractVector,
    frame::AbstractMatrix;
    step::Real = 1.0e-6,
)
    step > 0.0 || throw(ArgumentError("Finite-difference step must be positive."))
    sensitivity = zeros(Float64, 1, 6)
    for component in 1:6
        perturbation = zeros(Float64, 6)
        perturbation[component] = step
        dr = frame * perturbation[1:3]
        dv = frame * perturbation[4:6]
        plus = normalized_semimajor_axis(rbar + dr, vbar + dv)
        minus = normalized_semimajor_axis(rbar - dr, vbar - dv)
        sensitivity[1, component] = (plus - minus) / (2.0 * step)
    end
    return sensitivity
end

function embed_satellite_output_row(satellite_row, target_satellite::Int, nsat::Int)
    size(satellite_row) == (1, 6) || throw(ArgumentError("Satellite output row must be 1 by 6."))
    1 <= target_satellite <= nsat || throw(ArgumentError("Target satellite index is invalid."))
    output = zeros(Float64, 1, 6 * nsat)
    columns = (6 * (target_satellite - 1) + 1):(6 * target_satellite)
    output[1, columns] .= vec(satellite_row)
    return output
end

end
