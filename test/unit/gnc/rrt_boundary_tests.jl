using Test, Random, LinearAlgebra, StaticArrays

# A real standalone load: no HYPR configuration, objective or optimizer module.
module StandaloneRRT
    using Random, LinearAlgebra, StaticArrays
    const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
    include(joinpath(ROOT, "src/gnc/shared/path_geometry.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/sampling_settings.jl"))
    include(joinpath(ROOT, "src/gnc/shared/rpo/path_sampling.jl"))
    include(joinpath(ROOT, "src/gnc/rrt/tree_operations.jl"))
    include(joinpath(ROOT, "src/gnc/rrt/rrt_connect.jl"))
end

@testset "RRT executes with explicit policies and no HYPR" begin
    R = StandaloneRRT
    for name in (:RPOPSOConfig, :rpo_pso_config, :rpo_pso_bounds,
                 :rpo_post_refine_path, :rpo_normalized_path_cost_components)
        @test !isdefined(R, name)
    end
    start = SVector(-2.0, 0.0, 0.0)
    goal = SVector(2.0, 0.0, 0.0)
    bounds = (SVector(-3.0, -3.0, -3.0), SVector(3.0, 3.0, 3.0))
    score(p) = 3.0 * R.rpo_path_length(p)
    components(p) = (total=score(p), policy=:standalone_length)
    function outside_sphere(a, b)
        d = b - a
        t = dot(d, d) == 0.0 ? 0.0 : clamp(-dot(a, d) / dot(d, d), 0.0, 1.0)
        return norm(a + t * d) >= 0.5
    end
    for (planner, settings) in (
        (R.rpo_rrt_connect_plan_path, R.RPORRTConnectSettings(n_iters=300, shortcut_iters=8)),
        (R.rpo_rrt_star_plan_path, R.RPORRTStarSettings(n_iters=300, shortcut_iters=8)))
        refinements = Ref(0)
        refine(p) = (refinements[] += 1; (copy(p), score(p), true))
        common = (bounds=bounds, settings=settings, evaluate_components=components,
            evaluate_cost=score)
        for seed in (741, 742)
            rng = MersenneTwister(seed)
            expected_rng = copy(rng)
            direct = planner(start, goal, nothing; common..., rng=rng,
                edge_is_safe=(a,b)->true, refine_path=refine)
            @test direct.path == hcat(start, goal)
            @test direct.cost === 12.0
            @test direct.components.policy === :standalone_length
            @test direct.iterations == 0 && direct.path_found
            @test refinements[] == 0
            @test rand(rng) == rand(expected_rng)
            @test !hasproperty(direct, :config)
        end
        for seed in (741, 742)
            detour = planner(start, goal, nothing; common..., rng=MersenneTwister(seed),
                edge_is_safe=outside_sphere)
            @test detour.path_found
            @test detour.iterations > 0
            @test detour.path[:,1] == start && detour.path[:,end] == goal
            @test all(outside_sphere(detour.path[:,j],detour.path[:,j+1]) for j in 1:size(detour.path,2)-1)
            @test detour.cost == score(detour.path)
            @test !detour.refinement_improved
            @test refinements[] == 0
        end
        refined = planner(start, goal, nothing; common..., rng=MersenneTwister(741),
            edge_is_safe=outside_sphere, refine_path=refine)
        @test refined.path_found && refined.refinement_improved
        @test refinements[] == 1
        @test refined.cost == score(refined.path)
        failed = planner(start, goal, nothing; common..., rng=MersenneTwister(741),
            edge_is_safe=(a,b)->false)
        @test !failed.path_found
        @test failed.path == hcat(start, goal) # Legacy diagnostic only.
        @test failed.iterations == settings.n_iters
        @test_throws ErrorException planner(start, goal, nothing; common...,
            edge_is_safe=(a,b)->error("constraint callback failed"))
    end
end

@testset "RRT boundary rejects definition-time HYPR coupling" begin
    root = StandaloneRRT.ROOT
    path = joinpath(root, "src/gnc/rrt/rrt_connect.jl")
    source = read(path, String)
    function load_rrt(text)
        target = Module(gensym(:IsolatedRRT))
        Core.eval(target, :(using Random, LinearAlgebra, StaticArrays))
        Base.include_string(target, text, path)
        target
    end
    @test !isdefined(load_rrt(source), :RPOPSOConfig)
    rejected = try
        load_rrt(source * "\nrrt_forbidden(cfg::RPOPSOConfig) = cfg\n")
        nothing
    catch e
        e
    end
    @test rejected isa LoadError && rejected.error isa UndefVarError
    @test rejected isa LoadError && rejected.error isa UndefVarError &&
        rejected.error.var === :RPOPSOConfig
end
