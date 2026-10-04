using Test, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const EXAMPLE = joinpath(ROOT, "examples", "rpo_planner_env")

@testset "Documented example setup preserves tracked files" begin
    tracked = filter(!isempty, split(read(`git -C $ROOT ls-files -z -- examples/rpo_planner_env`, String), '\0'))
    @test "examples/rpo_planner_env/Project.template.toml" in tracked
    @test "examples/rpo_planner_env/Project.toml" ∉ tracked
    @test "examples/rpo_planner_env/Manifest.toml" ∉ tracked
    before = Dict(path => read(joinpath(ROOT, path)) for path in tracked)
    root_before = Dict(name => read(joinpath(ROOT, name)) for name in ("Project.toml", "Manifest.toml"))
    for attempt in 1:2
        withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do
            run(`$(Base.julia_cmd()) --startup-file=no --project=$EXAMPLE $(joinpath(EXAMPLE, "setup.jl"))`)
        end
        @test Dict(path => read(joinpath(ROOT, path)) for path in tracked) == before
        @test Dict(name => read(joinpath(ROOT, name)) for name in keys(root_before)) == root_before
        for name in ("Project.toml", "Manifest.toml")
            path = joinpath(EXAMPLE, name)
            @test isfile(path)
            @test success(`git -C $ROOT check-ignore -q -- $path`)
        end
        project = TOML.parsefile(joinpath(EXAMPLE, "Project.toml"))
        @test all(haskey(project["deps"], name) for name in ("SpaceAGORA", "SpaceAGORAHYPR", "HYPR"))
    end
end
