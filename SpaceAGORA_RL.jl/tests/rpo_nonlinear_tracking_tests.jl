using Test
using LinearAlgebra
using SpaceAGORA_RL

# Load the optional backend before compiling the testset body.
SpaceAGORA_RL._spaceagora_rpo_modules()

@testset "Actuator-aware RPO uses nonlinear realized-thrust feedback" begin
    modules = SpaceAGORA_RL._spaceagora_rpo_modules()
    guidance = modules.guidance
    navigation = modules.navigation
    geometry = navigation.RPOReferenceGeometry(
        navigation.RPOStationGeometry(reshape([0.0, 1e6, 0.0], 3, 1)),
    )
    tracking = guidance.RPOLQMPCTrackingSettings(dt_s=0.1, horizon=4, settle_time_s=0.0)
    config = RPOHyPRRLConfig(safe_distance_m=0.0)
    scenario = RPOHyPRRLScenario(
        start_rtn=[10.0, 0.0, 0.0], goal_rtn=[11.0, 0.0, 0.0],
        geometry=geometry, tracking_settings=tracking,
        pso_config=guidance.rpo_740_mpc_final_pso_config(),
    )
    t_ref = collect(0.0:tracking.dt_s:2.0)
    r_ref = repeat(scenario.start_rtn, 1, length(t_ref))
    r_ref[1, :] .+= t_ref ./ 2
    v_ref = zeros(size(r_ref))
    identity = [0.0, 0.0, 0.0, 1.0]
    result = SpaceAGORA_RL._coupled_rpo_tracking(
        modules, t_ref, r_ref, v_ref, scenario.goal_rtn,
        [0.0, 1.0], hcat(identity, identity), tracking, scenario, config,
    )
    plant = guidance.rpo_init_two_body_plant(result.state_history[:, 1], tracking.mean_motion_radps)
    for k in axes(result.realized_acceleration_history, 2)
        state = guidance.rpo_step_two_body!(plant, result.realized_acceleration_history[:, k], tracking.dt_s)
        @test state ≈ result.state_history[:, k + 1] atol=1e-10 rtol=0
    end
    @test result.propellant_used_kg > 0.0
    @test norm(result.command_history - result.realized_acceleration_history) > 1e-8
    @test all(isfinite, result.state_history)
end
