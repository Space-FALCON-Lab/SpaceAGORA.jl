using Test, Serialization, SpaceAGORA
const ES = SpaceAGORA.EDGServices
@testset "Core loads without the EDG companion" begin
    @test !ES.available()
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAEDGExt) === nothing
    @test_throws ES.EDGUnavailableError ES.require_provider()
    @test all(id.name != "SpaceAGORAEDG" for id in keys(Base.loaded_modules))
    @test ES.CONTRACT_VERSION == v"3.0.0"
    @test isdefined(SpaceAGORA, :run_simulation)
end


include("fixtures.jl")
function serialized_snapshot(value)
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end
@testset "An EDG run without its provider fails before mutating the mission" begin
    mktempdir() do directory
        output = joinpath(directory, "not-created")
        configuration = missing_provider_configuration(output)
        state = only(configuration.guidance_model.guidance_effectors).state
        state.last_alpha_rad[1] = 0.123
        state.selected_mode[1] = :targeting
        before_state = serialized_snapshot(state)
        before_vehicle = serialized_snapshot(configuration.dynamics_model.spacecraft)
        @test_throws ES.EDGUnavailableError run_simulation(configuration; isolate_state=false)
        @test !ispath(output)
        @test serialized_snapshot(state) == before_state
        @test serialized_snapshot(configuration.dynamics_model.spacecraft) == before_vehicle
        @test only(configuration.control_model.control_effectors).state === state
        @test !ES.available()
    end
end
