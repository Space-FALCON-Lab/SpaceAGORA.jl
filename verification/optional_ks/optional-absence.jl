# Run with JULIA_LOAD_PATH=@:@stdlib and the separately resolved absence environment.
include("baseline-fixtures.jl")
using Test
@testset "ordinary route with KS package absent" begin
    @test LOAD_PATH==["@","@stdlib"]
    @test Base.find_package("KSPropagator")===nothing
    @test Base.get_extension(SpaceAGORA,:SpaceAGORAKSPropagatorExt)===nothing
    result=run_simulation(baseline_config(:two_body,:dp8;duration=2.);return_solution=true)
    @test result.t[end]==2.
    @test all(isfinite,result.u[end])
    @test Base.get_extension(SpaceAGORA,:SpaceAGORAKSPropagatorExt)===nothing
    caught=try
        run_simulation(baseline_config(:two_body,:ks;duration=2.);return_solution=true)
        nothing
    catch err
        err
    end
    @test caught isa ArgumentError
    @test occursin("requires the optional KSPropagator package",sprint(showerror,caught))
    @test Base.get_extension(SpaceAGORA,:SpaceAGORAKSPropagatorExt)===nothing
end
