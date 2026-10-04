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
        old=read(pin,String);write(pin,replace(old,"0.1.0"=>"0.9.0"))
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path="")
        write(pin,old)
        dev=joinpath(root,"dev");mkpath(dev)
        write(joinpath(dev,"Project.toml"),"uuid = \"$(HYPRInstallation.HYPR_UUID)\"\nversion = \"0.1.0\"\n")
        @test HYPRInstallation.source_spec(root;local_path=dev).path == dev
        write(joinpath(dev,"Project.toml"),"uuid = \"$(HYPRInstallation.HYPR_UUID)\"\nversion = \"0.9.0\"\n")
        @test_throws ErrorException HYPRInstallation.source_spec(root;local_path=dev)
    end
    @test_throws ErrorException HYPRInstallation.setup(HYPRInstallation.ROOT)
end
