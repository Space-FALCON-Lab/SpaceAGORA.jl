using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

@testset "RPO nonlinear two-body comparison plant" begin
    guidance = SpaceAGORA.SimulationModel.GuidanceHooks
    n = 0.0011
    orbit = guidance.rpo_two_body_reference_orbit(n)
    @test sqrt(orbit.μ / orbit.radius_m^3) ≈ n
    @test_throws ArgumentError guidance.rpo_init_two_body_plant(zeros(6), 0.0)

    # Two different circular orbits have an analytic relative trajectory.
    # A large separation makes a mistaken HCW plant distinguishable.
    a = orbit.radius_m
    b = a + 100_000.0
    nc = sqrt(orbit.μ / b^3)
    initial = [b - a, 0.0, 0.0, 0.0, (nc - n) * b, 0.0]
    plant = guidance.rpo_init_two_body_plant(initial, n)
    duration = 600.0
    actual = guidance.rpo_step_two_body!(plant, zeros(3), duration)
    phase = (nc - n) * duration
    expected = [b * cos(phase) - a, b * sin(phase), 0.0,
                -b * (nc - n) * sin(phase), b * (nc - n) * cos(phase), 0.0]
    @test actual[1:3] ≈ expected[1:3] atol=1e-5 rtol=0
    @test actual[4:6] ≈ expected[4:6] atol=1e-8 rtol=0
    @test plant.state[1:3] ≈ a .* [cos(n * duration), sin(n * duration), 0.0] atol=1e-5 rtol=0
    radial_hcw = (4 - 3cos(n * duration)) * initial[1] +
        2 / n * (1 - cos(n * duration)) * initial[5]
    @test abs(actual[1] - radial_hcw) > 1.0

    # Thrust affects only the chaser. Halving the integration step preserves
    # sub-millimetre accuracy, including the RTN-to-inertial thrust rotation.
    initial = [10.0, -4.0, 2.0, 0.01, 0.02, -0.01]
    forced = guidance.rpo_init_two_body_plant(initial, n)
    refined = guidance.rpo_init_two_body_plant(initial, n)
    coast = guidance.rpo_init_two_body_plant(initial, n)
    command = [0.01, -0.005, 0.003]
    actual = guidance.rpo_step_two_body!(forced, command, 20.0)
    fine = guidance.rpo_step_two_body!(refined, command, 20.0; max_step_s=0.025)
    unforced = guidance.rpo_step_two_body!(coast, zeros(3), 20.0)
    @test actual ≈ fine atol=1e-6 rtol=0
    @test forced.state[1:6] == coast.state[1:6]
    @test norm(actual[1:3] - unforced[1:3]) > 1.0

    # Exercise the comparison evaluator itself, not only the plant helper.
    geometry = SpaceAGORA.SimulationModel.RPOReferenceGeometry(
        SpaceAGORA.SimulationModel.RPOStationGeometry(reshape([0.0, 1e6, 0.0], 3, 1)),
    )
    config = guidance.rpo_740_mpc_final_pso_config(
        retime_dt_s=0.5, retime_a_max_mps2=0.01,
    )
    tracking = guidance.RPOLQMPCTrackingSettings(dt_s=0.5, horizon=4, settle_time_s=1.0)
    path = [10.0 11.0; 0.0 0.0; 0.0 0.0]
    result = guidance.rpo_track_retimed_path_lqmpc(path, path[:, end], geometry, config, tracking)
    replay = guidance.rpo_init_two_body_plant(result.x_hist[:, 1], tracking.mean_motion_radps)
    for k in axes(result.u_hist, 2)
        state = guidance.rpo_step_two_body!(replay, result.u_hist[:, k], tracking.dt_s)
        @test state ≈ result.x_hist[:, k + 1] atol=1e-10 rtol=0
    end
    @test all(isfinite, result.x_hist)
end
