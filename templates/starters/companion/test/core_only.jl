# Run in a fresh process: SpaceAGORA alone loads none of this companion, and a
# compatibility request for another contract version is refused.
using Test, SpaceAGORA

@testset "core-only session" begin
    loaded() = any(package -> package.name == "StarterName", keys(Base.loaded_modules))
    @test !loaded()
    services = SpaceAGORA.ExternalPropagationServices
    @test_throws services.CompatibilityError services.require_version(v"0.0.1")
    @test !loaded()
end
