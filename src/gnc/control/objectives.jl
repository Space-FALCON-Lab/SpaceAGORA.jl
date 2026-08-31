#=
"""
    Objective metadata and solution outputs.

    ET and MED use the same condensed matrices. This file keeps the mode labels
    and the solved trajectory outputs in one place.
"""
=#
objective_kind(::TargetEnergyMode) = :target_energy
objective_kind(::MaxEnergyDepletionMode) = :max_energy_depletion

objective_label(::TargetEnergyMode) = "target_energy"
objective_label(::MaxEnergyDepletionMode) = "max_deplete"

function terminal_energy_affine_terms(problem::AerobrakingMPCProblem)
    S_E = build_output_selector(problem.N, problem.ny; yidx=4)
    r_E = zeros(problem.N)
    r_E[end] = 1.0
    H_Erow = vec(r_E' * S_E * problem.H)
    E_off = (r_E' * S_E * (problem.Mx * problem.δX0 .+ vec(problem.Ybar')))[1]
    return H_Erow, E_off
end

function commanded_area(problem::AerobrakingMPCProblem, ΔU)
    return fill(problem.Abar_m2, problem.N) .+ ΔU
end

function solution_from_area_delta(
    mode::AerobrakingMPCMode,
    problem::AerobrakingMPCProblem,
    ΔU::AbstractVector{<:Real};
    ok::Bool,
    objective::Real=NaN,
    slacks::AbstractVector{<:Real}=Float64[],
    solver_status::Symbol=ok ? :solved : :not_solved,
)
    H_Erow, E_off = terminal_energy_affine_terms(problem)
    ΔU_vec = Float64.(collect(ΔU))
    Y_pred = predict_outputs(problem.H, problem.Mx, problem.δX0, problem.Ybar, problem.ny, ΔU_vec)
    return AerobrakingMPCSolution(
        mode=objective_kind(mode),
        ok=Bool(ok),
        ΔU=ΔU_vec,
        slacks=Float64.(collect(slacks)),
        objective=Float64(objective),
        terminal_energy=dot(H_Erow, ΔU_vec) + E_off,
        commanded_area_m2=commanded_area(problem, ΔU_vec),
        predicted_outputs=Y_pred,
        solver_status=solver_status,
    )
end
