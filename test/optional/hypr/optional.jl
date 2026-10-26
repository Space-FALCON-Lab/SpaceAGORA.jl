using Test, Pkg
order=only(ARGS)
if order=="core-first"
    using SpaceAGORA
    @test !SpaceAGORA.hypr_available()
    const before_type=SpaceAGORA.RPOPSOConfig
    const before_function=SpaceAGORA.SimulationModel.rpo_pso_plan_path
    @test_throws SpaceAGORA.HYPRUnavailableError SpaceAGORA.make_rpo_configuration(
        planner=SpaceAGORA.HYPRRPOPlanner(SpaceAGORA.RPOPSOConfig()))
    using HYPR
    @test before_type===SpaceAGORA.RPOPSOConfig
    @test before_function===SpaceAGORA.SimulationModel.rpo_pso_plan_path
elseif order=="hypr-first"
    using HYPR
    @test Base.get_extension(HYPR, :HYPRSpaceAGORAExt) === nothing
    using SpaceAGORA
else
    error("Expected core-first or hypr-first")
end
@testset "Direct companion imports preserve planner identities" begin
    @test SpaceAGORA.hypr_available()
    @test_throws MethodError SpaceAGORA.SimulationModel.rpo_pso_plan_path()
    @test any(p.name=="HYPR" for p in values(Pkg.dependencies()))
    @test all(p.name!="SpaceAGORAHYPR" for p in values(Pkg.dependencies()))
    adapter = Base.get_extension(HYPR, :HYPRSpaceAGORAExt)
    @test adapter isa Module
    @test adapter.initialized()
    @test all(isdefined(adapter, name) for name in (:SwarmPolicy, :RPO, :RobotArm, :PlannerAdapter))
    S=SpaceAGORA.SimulationModel
    @test SpaceAGORA.RPOPSOConfig===S.RPOPSOConfig===S.GuidanceHooks.RPOPSOConfig
    @test SpaceAGORA.RobotArmHYPRConfig===S.RobotArmPlanning.RobotArmHYPRConfig
    @test SpaceAGORA.RobotArmHYPRResult===S.RobotArmPlanning.RobotArmHYPRResult
    @test parentmodule(S.GuidanceHooks.rpo_pso_plan_path)===S.GuidanceHooks
    @test parentmodule(S.RobotArmPlanning.plan_robot_arm_motion_hypr)===S.RobotArmPlanning
    @test any(m.module===adapter.RPO for m in methods(S.rpo_pso_plan_path))
    @test parentmodule(S.HYPRUtils.hypr_path_length)===S.HYPRUtils
    @test parentmodule(S.GuidanceHooks.rpo_retime_samples)===S.GuidanceHooks
end
include(joinpath(@__DIR__,"..","..","unit","gnc","rpo_planner_adapter_tests.jl"))
# One order runs the wider integration campaign; both run the seeded adapters.
if order=="core-first"
    for path in ("rpo_port_tests.jl", "gnc/rpo_planner_compatibility_tests.jl", "gnc/rpo_planner_lifecycle_tests.jl",
            "gnc/rpo_ownership_reproducibility_tests.jl", "gnc/rpo_hypr_manuscript_tests.jl",
            "gnc/rpo_final_isparo_reconstruction_tests.jl",
            "gnc/shared_sampling_tests.jl", "gnc/shared_metrics_tests.jl", "robotics/runtests.jl")
        Base.include(Module(gensym(:HYPRIntegration)),joinpath(@__DIR__,"..","..","unit",path))
    end
end
