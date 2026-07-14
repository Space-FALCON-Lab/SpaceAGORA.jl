function seed_solution_for_save_csv!(
    solution::Solution;
    n_bodies::Int,
    n_reaction_wheels::Int,
    n_thrusters::Int,
    base::Float64=1.0,
    closed_form::Bool=false
)
    solution.physical_properties.α = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.β = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_rate = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_load = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.rw_h = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.rw_τ = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.thruster_forces = [Float64[] for _ in 1:n_thrusters]

    function _push_sample!(field)
        if field isa Vector{Float64}
            push!(field, base)
        elseif field isa Vector{Int64}
            push!(field, 1)
        elseif field isa Vector{Vector{Float64}}
            for subfield in field
                push!(subfield, base)
            end
        end
        return nothing
    end

    for group in (solution.orientation, solution.physical_properties, solution.performance, solution.forces)
        for fname in fieldnames(typeof(group))
            _push_sample!(getfield(group, fname))
        end
    end

    if closed_form
        solution.closed_form.t_cf = [base + 10.0]
        solution.closed_form.h_cf = [base + 20.0]
        solution.closed_form.γ_cf = [base + 30.0]
        solution.closed_form.v_cf = [base + 40.0]
    end

    return solution
end
