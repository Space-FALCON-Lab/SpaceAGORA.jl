module MissionHistoryContractTests
using Test, SpaceAGORA
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO, "test", "support", "mission_history_contract.jl"))
include(joinpath(REPO, "scripts", "sync_companion_support.jl"))
@testset "Core mission-history contract and generated-source parity" begin
    @test isempty(CompanionSupportSync.sync_support(REPO))
    MissionHistoryContract.check(SpaceAGORA.SimulationModel.SceneVisualization._validated_mission_history)
end
end
