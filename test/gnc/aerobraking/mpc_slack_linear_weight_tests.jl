using Test
using LinearAlgebra
using SpaceAGORA
using SpaceAGORA.SimulationModel.DynamicsKS

# Linear slack term of the QP: zero by default, and with a weight above the
# multiplier of the softened drag limit the slack vanishes where the limit can be met.

function _slw_config(; slack_linear_weight=0.0)
    AerobrakingMPCConfig(
        mode=MaxEnergyDepletionMode(), bus_reference_area_m2=5.0, controllable_area_m2=5.0,
        mass_kg=500.0, drag_coefficient=2.2, qdot_max_w_cm2=1.0, heat_load_max_j_cm2=1.0,
        drag_max_n=10.0, area_slew_max_m2_s=0.2, use_constraints=true,
        use_slew_constraint=false, use_qdot_constraint=false, use_heat_load_constraint=false,
        use_drag_constraint=true, target_energy_mj_kg=0.0, area_weight=1.0e-6,
        area_slew_weight=0.0, slack_weight=1.0e3, target_energy_weight=0.0,
        max_depletion_energy_weight=1.0, osqp_eps_abs=1.0e-9, osqp_eps_rel=1.0e-9,
        osqp_max_iter=200_000, slack_linear_weight=slack_linear_weight)
end

# Condensed problem: the drag of every node is 13 N at the nominal area (10 m^2)
# and changes by 1 N per m^2, against a 10 N limit, so the limit is met at 7 m^2.
# The end-of-pass energy falls by 1 J/kg per m^2 of area at every node, so maximum
# depletion pushes the area up against the drag limit. The multiplier of each
# normalized drag row is 2.
function _slw_problem(N)
    ny = 4
    H = zeros(N * ny, N)
    for k in 1:N
        H[(k - 1) * ny + 2, k] = 1.0
        H[(N - 1) * ny + 4, k] = -1.0
    end
    Ybar = zeros(N, ny)
    Ybar[:, 2] .= 13.0
    Ybar[:, 4] .= 1.0
    return AerobrakingMPCProblem(
        params=KSPropagationParams(Re=6_378_137.0, μ=3.986004418e14, J2=0.0, Ω=0.0),
        H=H, Mx=zeros(N * ny, 9), δX0=zeros(9), N=N, ny=ny,
        t=collect(range(0.0, step=1.0, length=N)),
        Ybar=Ybar, Xbar=zeros(N, 9), Abar_m2=10.0)
end

@testset "linear slack term: default keeps the quadratic-only cost" begin
    @test _slw_config().slack_linear_weight == 0.0
    problem = _slw_problem(5)
    solution = solve_mpc_qp(problem, _slw_config())
    @test solution.ok
    # Quadratic penalty alone: slack = multiplier / (2 slack_weight) = 1e-3 on every node.
    @test all(isapprox.(solution.slacks, 1.0e-3; atol=1.0e-4))
    @test all(solution.predicted_outputs[:, 2] .> 10.0)
end

@testset "linear slack term: the limit is met where it can be" begin
    problem = _slw_problem(5)
    solution = solve_mpc_qp(problem, _slw_config(slack_linear_weight=20.0))
    @test solution.ok
    @test maximum(abs.(solution.slacks)) < 1.0e-5
    @test all(isapprox.(solution.commanded_area_m2, 7.0; atol=1.0e-4))
    @test maximum(solution.predicted_outputs[:, 2]) <= 10.0 + 1.0e-4
    @test_throws ArgumentError solve_mpc_qp(problem, _slw_config(slack_linear_weight=-1.0))
end
