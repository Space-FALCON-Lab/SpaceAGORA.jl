@testset "T-EDG Strategy Contract" begin
    @test isdefined(SimulationModel.GuidanceHooks, :target_planning)
    @test isdefined(SimulationModel.GuidanceHooks, :control_solarpanels_targeting_num_int)
    @test isdefined(SimulationModel.GuidanceHooks, :control_solarpanels_targeting_closed_form)

    cnf_state = (time_switch_1 = 12.0, time_switch_2 = 34.0, security_mode = false)
    input = SimulationModel.GuidanceHooks.AerobrakingGuidanceInput(
        ip=nothing,
        mission=nothing,
        args=(; cnf = cnf_state),
        cnf=cnf_state,
    )
    out = SimulationModel.GuidanceHooks.compute_aerobraking_guidance(
        SimulationModel.GuidanceHooks.TEdgStrategy(),
        input,
    )
    @test out.time_switch_1 == 12.0
    @test out.time_switch_2 == 34.0
    @test out.security_mode == false
end
