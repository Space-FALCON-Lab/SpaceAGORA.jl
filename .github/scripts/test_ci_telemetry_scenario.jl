using Test, TOML
include(joinpath(@__DIR__, "ci_telemetry_scenario.jl"))
@testset "CI telemetry partitions preserve every original scenario" begin
    path=normpath(joinpath(@__DIR__,"../../test/telemetry_benchmark_manifest.toml"))
    original=TOML.parsefile(path)
    names=[x["name"] for x in original["scenarios"]]
    @test Set(names)==Set(["odyssey","vex","earth_gmat"])
    @test length(names)==3
    combined=[]
    mktempdir() do dir
        for name in names
            dest=joinpath(dir,name*".toml")
            telemetry_scenario_main([path,name,dest])
            selected=TOML.parsefile(dest)
            @test length(selected["scenarios"])==1
            @test selected["scenarios"][1]==only(filter(x->x["name"]==name,original["scenarios"]))
            @test Dict(k=>v for (k,v) in selected if k!="scenarios")==Dict(k=>v for (k,v) in original if k!="scenarios")
            append!(combined,selected["scenarios"])
        end
        @test combined==original["scenarios"]
        @test_throws ErrorException telemetry_scenario_main([path,"unknown",joinpath(dir,"bad.toml")])
        @test !isfile(joinpath(dir,"bad.toml"))
    end
    @test_throws ErrorException telemetry_scenario_main([path,names[1],path])
    @test TOML.parsefile(path)==original
    duplicate=deepcopy(original);push!(duplicate["scenarios"],duplicate["scenarios"][1])
    @test_throws ErrorException select_telemetry_scenario(duplicate,names[1])
end
