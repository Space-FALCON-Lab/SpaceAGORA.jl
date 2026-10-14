using Test
const smoke=joinpath(@__DIR__,"..","smoke","ci_examples_suite_smoke.jl")
const routing=Module(:ExampleRoutingWitness)
Base.include(ex -> ex == :(examples = list_examples()) ? :(examples = String[]) : ex,
    routing,smoke)
@testset "Example commands isolate the optional companion" begin
    for name in ("Earth_RPO_CubeSat_MPC.jl", "Earth_RPO_CubeSat_MPC_Batch.jl",
            "Earth_RPO_CubeSat_MPC_PlannerComparison.jl", "Earth_RPO_CubeSat_MPC_Replanning.jl",
            "Robot_Arm_Planner_Cloth_Demo.jl")
        cmd=routing.example_command(joinpath(routing.EXAMPLES_DIR,name))
        @test "--project=$(joinpath(routing.EXAMPLES_DIR,"rpo_planner_env"))" in cmd.exec
        @test all(!startswith(arg,"--load") && arg!="-L" for arg in cmd.exec)
    end
    @test routing.example_project("Earth_Torque_Free_Test.jl")==routing.REPO_ROOT
    mktempdir() do dir
        # Use the actual smoke child process to prove an ordinary example cannot
        # see the companion, even after the parent has loaded it elsewhere.
        probe=joinpath(dir,"baseline_core_witness.jl")
        write(probe,"using SpaceAGORA; @assert !hypr_available(); @assert Base.find_package(\"SpaceAGORAHYPR\") === nothing; println(\"baseline_core_only_ok\")")
        ok,has_nan,output=routing.run_example(probe)
        @test ok
        @test !has_nan
        @test occursin("baseline_core_only_ok",output)
    end
end
