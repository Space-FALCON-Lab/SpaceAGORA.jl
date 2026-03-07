@testset "Aerobraking Policy Selector Stub" begin
    cfg = SimulationModel.AerobrakingPolicyConfig()
    cfg_t = SimulationModel.AerobrakingPolicyConfig(default_strategy=SimulationModel.T_EDG)
    selector = SimulationModel.DefaultAerobrakingPolicySelector()

    input_default = SimulationModel.GuidanceHooks.AerobrakingGuidanceInput(
        ip=nothing,
        mission=nothing,
        args=nothing,
    )
    @test SimulationModel.select_strategy(selector, cfg, input_default) == SimulationModel.E_EDG
    @test SimulationModel.select_strategy(selector, cfg_t, input_default) == SimulationModel.T_EDG

    input_t = SimulationModel.GuidanceHooks.AerobrakingGuidanceInput(
        ip=nothing,
        mission=nothing,
        args=Dict{Symbol, Any}(:aerobraking_strategy => "t_edg"),
    )
    @test SimulationModel.select_strategy(selector, cfg, input_t) == SimulationModel.E_EDG
    @test SimulationModel.select_strategy(selector, cfg_t, input_t) == SimulationModel.T_EDG

    drl = SimulationModel.DRLPolicyAdapterStub()
    @test_throws ErrorException SimulationModel.select_strategy(drl, cfg, input_default)
end
