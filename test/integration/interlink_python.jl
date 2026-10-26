# Both fresh-process load orders must activate the optional matching backend.
using Test
length(ARGS) == 1 || error("Choose core-first or python-first")
if only(ARGS) == "core-first"
    using SpaceAGORA
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt) === nothing
    using PythonCall
elseif only(ARGS) == "python-first"
    using PythonCall
    using SpaceAGORA
else
    error("Choose core-first or python-first")
end
@test pkgversion(PythonCall) >= v"0.9.36"
@test Base.get_extension(SpaceAGORA, :SpaceAGORAPythonCallExt) !== nothing
include(joinpath(@__DIR__, "..", "unit", "simulation", "interlink_tests.jl"))
include(joinpath(@__DIR__, "..", "unit", "simulation", "interlink_python_tests.jl"))
