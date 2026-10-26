# Objective definitions shared by energy targeting and maximum depletion.
objective_kind(::TargetEnergyMode) = :target_energy
objective_kind(::MaxEnergyDepletionMode) = :max_energy_depletion

# The predicted end-of-pass energy is affine in the area increments:
# E_N = H_Erow ⋅ ΔU + E_off, from the last node's energy output (index 4).
function terminal_energy_affine_terms(problem::AerobrakingMPCProblem)
    row = (problem.N - 1) * problem.ny + 4
    H_Erow = vec(problem.H[row, :])
    E_off = (problem.Mx[row, :]' * problem.δX0) + problem.Ybar[problem.N, 4]
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
