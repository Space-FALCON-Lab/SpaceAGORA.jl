using Test, TOML
include(joinpath(@__DIR__, "..", "..", "scripts", "setup_hypr.jl"))
@testset "Immutable HYPR installation contract" begin
    mktempdir() do root
        folder=joinpath(root,"packages");mkpath(folder)
        pin=joinpath(folder,"HYPRSource.toml")
        write(pin,"url = \"https://github.com/Space-FALCON-Lab/HYPR.jl.git\"\nrev = \"$(repeat("a",40))\"\nuuid = \"$(HYPRInstallation.HYPR_UUID)\"\nversion = \"0.5.0\"\n")
        spec=HYPRInstallation.source_spec(root;local_path="")
        @test spec.rev == repeat("a",40)
        @test spec.url == "https://github.com/Space-FALCON-Lab/HYPR.jl.git"
        for bad in ("main", "v0.1.0", "abc123", repeat("g",40))
            text=read(pin,String);write(pin,replace(text,repeat("a",40)=>bad))
            @test_throws ErrorException HYPRInstallation.source_spec(root;local_path="")
            write(pin,text)
        end
        # HYPRServices 3.0.0 refuses the 0.4 and earlier providers and any version it does not list
        for unsupported in ("0.1.0", "0.1.1", "0.2.0", "0.3.0", "0.4.0", "0.9.0")
            old=read(pin,String);write(pin,replace(old,"0.5.0"=>unsupported))
            @test_throws ErrorException HYPRInstallation.source_spec(root;local_path="")
            write(pin,old)
        end
        dev=joinpath(root,"dev");mkpath(dev)
        write(joinpath(dev,"Project.toml"),"uuid = \"$(HYPRInstallation.HYPR_UUID)\"\nversion = \"0.5.0\"\n")
        @test HYPRInstallation.source_spec(root;local_path=dev).path == dev
        dev_project=joinpath(dev,"Project.toml")
        write(dev_project,replace(read(dev_project,String),"0.5.0"=>"0.4.0"))
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path=dev)
        write(dev_project,replace(read(dev_project,String),"0.4.0"=>"0.3.0"))
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path=dev)
        write(dev_project,replace(read(dev_project,String),"0.3.0"=>"0.2.0"))
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path=dev)
        write(dev_project,replace(read(dev_project,String),"0.2.0"=>"0.1.1"))
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path=dev)
        write(joinpath(dev,"Project.toml"),"uuid = \"$(HYPRInstallation.HYPR_UUID)\"\nversion = \"0.9.0\"\n")
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path=dev)
    end
    @test_throws ErrorException HYPRInstallation.setup(HYPRInstallation.ROOT)
end

@testset "Root aliases refuse before any installation write" begin
    mktempdir() do parent
        root = joinpath(parent, "checkout")
        mkpath(root)
        project = joinpath(root, "Project.toml")
        manifest = joinpath(root, "Manifest.toml")
        write(project, "name = \"RootGuardFixture\"\nuuid = \"d271120b-0ade-4aeb-b5fc-5e636b5dd1d3\"\n[deps]\n")
        write(manifest, "# Sentinel manifest: refusal must not reach Pkg.\n")
        before = (read(project), read(manifest), readdir(root))
        link = joinpath(parent, "checkout-link")
        symlink(root, link; dir_target=true)
        active = Base.active_project()
        cd(parent) do
            for root_spelling in (root, root * "/"), environment in (
                    root, root * "/", joinpath(root, "."),
                    "checkout", joinpath("checkout", "..", "checkout"), link, link * "/")
                err = try
                    HYPRInstallation.setup(environment; root=root_spelling, with_hypr=false)
                    nothing
                catch caught
                    caught
                end
                @test err isa ErrorException
                @test occursin("Use a separate HYPR project", sprint(showerror, err))
                @test (read(project), read(manifest), readdir(root)) == before
                @test Base.active_project() == active
            end
        end
    end
end

@testset "Retired loader manifest migration preserves weak identities" begin
    mktempdir() do directory
        manifest = joinpath(directory, "Manifest.toml")
        core_uuid = "afbfb69f-5c0b-4832-b760-43725dff8540"
        loader_uuid = "4096c904-79d0-4d9a-9f5d-f1d6c7a2f1b4"
        original = Dict("manifest_format" => "2.0", "julia_version" => string(VERSION),
            "deps" => Dict(
                "SpaceAGORA" => [Dict("uuid" => core_uuid, "version" => "0.5.1")],
                "HYPR" => [Dict("uuid" => string(HYPRInstallation.HYPR_UUID), "version" => "0.4.0",
                    "weakdeps" => ["SpaceAGORA", "SpaceAGORAHYPR"],
                    "extensions" => Dict("HYPRSpaceAGORAExt" => "SpaceAGORA"))],
                "SpaceAGORAHYPR" => [Dict("uuid" => loader_uuid, "version" => "0.2.0",
                    "path" => joinpath(directory, "removed-loader"), "deps" => ["SpaceAGORA", "HYPR"])]))
        open(manifest, "w") do io; TOML.print(io, original; sorted=true); end
        @test !ispath(original["deps"]["SpaceAGORAHYPR"][1]["path"])
        @test HYPRInstallation.retire_loader_manifest!(manifest) === nothing
        migrated = TOML.parsefile(manifest)
        @test !haskey(migrated["deps"], "SpaceAGORAHYPR")
        @test migrated["deps"]["SpaceAGORA"] == original["deps"]["SpaceAGORA"]
        @test migrated["deps"]["HYPR"][1]["weakdeps"] == Dict(
            "SpaceAGORA" => core_uuid, "SpaceAGORAHYPR" => loader_uuid)
        @test migrated["deps"]["HYPR"][1]["extensions"] == original["deps"]["HYPR"][1]["extensions"]
        @test length(HYPRInstallation.Pkg.Types.read_manifest(manifest)) == 2
        once = read(manifest)
        @test HYPRInstallation.retire_loader_manifest!(manifest) === nothing
        @test read(manifest) == once
        # A package that really imports the old loader needs its own migration.
        # Do not silently erase that dependency or rewrite its implementation.
        for dependencies in (["SpaceAGORAHYPR"], Dict("SpaceAGORAHYPR" => loader_uuid))
            blocked = deepcopy(original)
            blocked["deps"]["Consumer"] = [Dict("uuid" => "7296cf7c-9f96-44d3-93c0-ef879936dd03",
                "version" => "1.0.0", "deps" => dependencies)]
            open(manifest, "w") do io; TOML.print(io, blocked; sorted=true); end
            before = read(manifest)
            err = try HYPRInstallation.retire_loader_manifest!(manifest); nothing catch caught; caught end
            @test err isa ErrorException
            @test occursin("Consumer still depends on the retired", sprint(showerror, err))
            @test read(manifest) == before
            project = joinpath(directory, "Project.toml")
            write(project, "[deps]\nSpaceAGORAHYPR = \"$loader_uuid\"\n[preferences]\nkeep = true\n")
            project_before = read(project)
            @test_throws ErrorException HYPRInstallation.setup(directory; with_hypr=false)
            @test read(project) == project_before
            @test read(manifest) == before
        end
    end
end
