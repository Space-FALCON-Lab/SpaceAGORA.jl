@testset "Aerobraking Policy Selector Stub" begin
    cfg = SimulationModel.AerobrakingPolicyConfig()
    selector = SimulationModel.DefaultAerobrakingPolicySelector()

    input_default = SimulationModel.GuidanceHooks.AerobrakingGuidanceInput(
        ip=nothing,
        mission=nothing,
        args=Dict{Symbol, Any}(),
    )
    @test SimulationModel.select_strategy(selector, cfg, input_default) == SimulationModel.E_EDG

    input_t = SimulationModel.GuidanceHooks.AerobrakingGuidanceInput(
        ip=nothing,
        mission=nothing,
        args=Dict{Symbol, Any}(:aerobraking_strategy => "t_edg"),
    )
    @test SimulationModel.select_strategy(selector, cfg, input_t) == SimulationModel.T_EDG

    drl = SimulationModel.DRLPolicyAdapterStub()
    @test_throws ErrorException SimulationModel.select_strategy(drl, cfg, input_default)
end
