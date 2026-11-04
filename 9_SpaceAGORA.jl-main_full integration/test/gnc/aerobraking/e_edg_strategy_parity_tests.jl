@testset "E-EDG Strategy Contract" begin
    @test isdefined(SimulationModel.GuidanceHooks, :compute_e_edg_guidance_window!)
    @test hasmethod(
        SimulationModel.GuidanceHooks.compute_aerobraking_guidance,
        (SimulationModel.GuidanceHooks.EEdgStrategy, SimulationModel.GuidanceHooks.AerobrakingGuidanceInput),
    )

    for fn in (
        :switch_calculation_with_integration,
        :switch_calculation,
        :second_time_switch_recalc_with_integration,
        :second_time_switch_recalc,
        :security_mode,
    )
        @test isdefined(SimulationModel.GuidanceHooks, fn)
    end
end
