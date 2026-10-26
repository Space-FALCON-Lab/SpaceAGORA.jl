using Test, TOML
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "setup_visualizer.jl"))
const VI = VisualizerInstallation
@testset "Explicit immutable visualizer installation" begin
    mktempdir() do parent
        root = joinpath(parent, "core"); mkpath(joinpath(root, "packages"))
        project = joinpath(root, "Project.toml")
        write(project, "name = \"SpaceAGORA\"\n[deps]\nJSON = \"682c06a0-de6a-54ab-a142-c8b1cf79cde6\"\n")
        write(joinpath(root, "Manifest.toml"), "# unchanged sentinel\n")
        pin = joinpath(root, "packages", "VisualizerSource.toml")
        values = Dict("url"=>"https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl.git",
            "rev"=>repeat("a",40), "uuid"=>string(VI.VISUALIZER_UUID), "version"=>"0.2.0")
        writepin(d) = open(io -> TOML.print(io, d), pin, "w")
        writepin(values)
        spec = VI.source_spec(root; local_path="")
        @test spec.rev == repeat("a",40)
        @test spec.url == values["url"]
        for (key, bad) in (("rev","main"),("rev","v0.2.0"),("rev",repeat("g",40)),
                ("uuid","wrong"),("version","0.9.0"),("url","https://example.com/repo.git"))
            writepin(merge(values, Dict(key=>bad)))
            @test_throws ErrorException VI.source_spec(root; local_path="")
        end
        writepin(values)
        env = joinpath(parent,"separate")
        plan = VI.installation_plan(env; root, local_path="")
        @test !ispath(env)
        @test plan.project["deps"]["SpaceAGORAVisualizer"] == string(VI.VISUALIZER_UUID)
        @test plan.project["sources"]["SpaceAGORAVisualizer"]["rev"] == values["rev"]
        @test plan.project["sources"]["SpaceAGORA"]["path"] == root
        core = VI.installation_plan(env; root, with_visualizer=false)
        @test !haskey(core.project["deps"], "SpaceAGORAVisualizer")
        @test !haskey(core.project["sources"], "SpaceAGORAVisualizer")
        @test core.spec === nothing
        dev=joinpath(parent,"renderer");mkpath(dev)
        devproject=joinpath(dev,"Project.toml")
        write(devproject,"name = \"SpaceAGORAVisualizer\"\nuuid = \"$(VI.VISUALIZER_UUID)\"\nversion = \"0.2.0\"\n")
        @test VI.source_spec(root; local_path=dev).path == dev
        @test VI.installation_plan(env; root,local_path=dev).project["sources"]["SpaceAGORAVisualizer"] == Dict("path"=>dev)
        @test_throws ErrorException VI.setup(dev;root,local_path=dev)
        write(devproject,replace(read(devproject,String),"0.2.0"=>"0.9.0"))
        @test_throws ErrorException VI.source_spec(root;local_path=dev)
        before=(read(project),read(joinpath(root,"Manifest.toml")),readdir(root))
        alias=joinpath(parent,"alias");symlink(root,alias;dir_target=true)
        active=Base.active_project()
        for candidate in (root,root*"/",joinpath(root,"."),alias)
            @test_throws ErrorException VI.setup(candidate;root,with_visualizer=false)
            @test (read(project),read(joinpath(root,"Manifest.toml")),readdir(root)) == before
            @test Base.active_project() == active
        end
        mkpath(env)
        open(io->TOML.print(io,plan.project),joinpath(env,"Project.toml"),"w")
        @test_throws ErrorException VI.installation_plan(env;root,with_visualizer=false)
    end
end
