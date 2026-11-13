using Test

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PROVIDER = Symbol(only(ARGS))
const ORIGINAL_LOAD_PATH = copy(LOAD_PATH)
const EXAMPLE, STOP = PROVIDER == :SpaceAGORAEDG ?
    ("AGORA_Odyssey_Energy_Depletion_Control_Test.jl", "\nsetup_gram_example!()") :
    ("Earth_RPO_CubeSat_MPC.jl", "\nconst DEFAULT_STATION_DIMS_M")
const PATH = joinpath(ROOT, "examples", EXAMPLE)
# Exercise the real imports without running a mission or querying native GRAM.
Base.include_string(@__MODULE__, first(split(read(PATH, String), STOP)), PATH)

using JSON
@testset "fresh provider process selects plotting dependencies" begin
    @test v"0.21" <= pkgversion(JSON) < v"0.22"
    @test realpath(Base.active_project()) == realpath(joinpath(ROOT, "examples", "Project.toml"))
    @test LOAD_PATH == ORIGINAL_LOAD_PATH
    @test isdefined(@__MODULE__, PROVIDER)
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAPlotlyJSExt) !== nothing
    if PROVIDER == :SpaceAGORAEDG
        @test SpaceAGORA.EDGServices.available()
    end
end
