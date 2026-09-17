using Test
using SpaceAGORA

@testset "shared example imports preserve caller module aliases" begin
    common_path = normpath(joinpath(@__DIR__, "..", "..", "examples", "common.jl"))
    aliases = quote
        const SimulationEngine = SpaceAGORA.SimulationEngine
        const SimulationModel = SpaceAGORA.SimulationModel
        const RuntimeServices = SpaceAGORA.RuntimeServices
        const SM = SimulationModel
    end

    for predeclared in (false, true)
        @testset "aliases declared before include: $predeclared" begin
            # Separate caller modules match examples and analytical studies.
            # Including common.jl must not load GRAM or start a simulation.
            probe_module = Module(gensym(:ExampleImports))
            Core.eval(probe_module, :(using SpaceAGORA))
            predeclared && Core.eval(probe_module, aliases)
            Base.include(probe_module, common_path)

            # AnalyticalPerturbationModels declares SM after its common include.
            # An imported binding rejects this otherwise valid local constant.
            Core.eval(probe_module, aliases)
            @test probe_module.SimulationEngine === SpaceAGORA.SimulationEngine
            @test probe_module.SimulationModel === SpaceAGORA.SimulationModel
            @test probe_module.RuntimeServices === SpaceAGORA.RuntimeServices
            @test probe_module.SM === SpaceAGORA.SimulationModel
            @test probe_module.run_simulation === SpaceAGORA.run_simulation
            @test probe_module.quat_mult === SpaceAGORA.SimulationModel.quat_mult
            @test probe_module.make_example_config === SpaceAGORA.TelemetryVerification.make_example_config
            @test probe_module.make_three_body_spacecraft === SpaceAGORA.TelemetryVerification.make_three_body_spacecraft
            @test probe_module.run_and_report === SpaceAGORA.TelemetryVerification.run_and_report
        end
    end
end
