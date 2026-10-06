using Test, TOML
include(joinpath(@__DIR__, "..", "scripts", "setup_hypr.jl"))
@testset "Immutable HYPR installation contract" begin
    mktempdir() do root
        folder=joinpath(root,"packages","SpaceAGORAHYPR");mkpath(folder)
        pin=joinpath(folder,"HYPRSource.toml")
        write(pin,"url = \"https://github.com/Space-FALCON-Lab/HYPR.jl.git\"\nrev = \"$(repeat("a",40))\"\nuuid = \"$(HYPRInstallation.HYPR_UUID)\"\nversion = \"0.1.0\"\n")
        spec=HYPRInstallation.source_spec(root;local_path="")
        @test spec.rev == repeat("a",40)
        @test spec.url == "https://github.com/Space-FALCON-Lab/HYPR.jl.git"
        for bad in ("main", "v0.1.0", "abc123", repeat("g",40))
            text=read(pin,String);write(pin,replace(text,repeat("a",40)=>bad))
            @test_throws ErrorException HYPRInstallation.source_spec(root;local_path="")
            write(pin,text)
        end
        old=read(pin,String)
        write(pin,replace(old,"0.1.0"=>"0.1.1"))
        @test HYPRInstallation.source_spec(root;local_path="").rev == repeat("a",40)
        write(pin,old)
        old=read(pin,String);write(pin,replace(old,"0.1.0"=>"0.9.0"))
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path="")
        write(pin,old)
        dev=joinpath(root,"dev");mkpath(dev)
        write(joinpath(dev,"Project.toml"),"uuid = \"$(HYPRInstallation.HYPR_UUID)\"\nversion = \"0.1.0\"\n")
        @test HYPRInstallation.source_spec(root;local_path=dev).path == dev
        dev_project=joinpath(dev,"Project.toml")
        write(dev_project,replace(read(dev_project,String),"0.1.0"=>"0.1.1"))
        @test HYPRInstallation.source_spec(root;local_path=dev).path == dev
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
