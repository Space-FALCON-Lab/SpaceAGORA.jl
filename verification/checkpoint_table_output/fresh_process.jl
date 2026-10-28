using Test,TOML
length(ARGS)==1 || error("usage: checkpoint-fresh-process.jl NEW_OUTPUT_DIRECTORY")
out=abspath(only(ARGS));isdir(out) && !isempty(readdir(out)) && error("Refusing previous evidence overwrite")
mkpath(out);project=dirname(Base.active_project())
worker=get(ENV,"CHECKPOINT_TABLE_WORKER",joinpath(@__DIR__,"worker.jl"))
function child(phase,directory,duration,rate,interval)
    mkpath(directory)
    open(joinpath(directory,phase*".log"),"w") do io
        run(pipeline(`$(Base.julia_cmd()) --startup-file=no --project=$project $worker $phase $directory $duration $rate $interval`;stdout=io,stderr=io))
    end
end
child("direct",joinpath(out,"direct"),1.04,.1,.37)
child("interrupt",joinpath(out,"cold"),1.04,.1,.37)
child("resume",joinpath(out,"cold"),1.04,.1,.37)
short=nextfloat(.3)
child("checkpoint",joinpath(out,"extended"),short,.1,short)
child("resume",joinpath(out,"extended"),1.04,.1,short)
direct=TOML.parsefile(joinpath(out,"direct/direct.toml"))
cold=TOML.parsefile(joinpath(out,"cold/resume.toml"))
short_result=TOML.parsefile(joinpath(out,"extended/checkpoint.toml"))
extended=TOML.parsefile(joinpath(out,"extended/resume.toml"))
# Literal expected nominal epochs make the intended decimal cadence explicit.
required=[0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.04]
@testset "Cold process and mission extension preserve table cadence" begin
    @test direct["time"]==required
    for result in (cold,extended)
        @test isempty(setdiff(required,result["time"]))
        @test length(unique(result["time"]))==length(result["time"])
        @test first(result["time"])==0.
        @test last(result["time"])==1.04
        @test result["columns"]==direct["columns"]
        @test result["max_affine_position_error_m"]<=3e-5
    end
    @test short_result["time"]==[0.,.1,.2,.3,short]
    @test count(==(.3),extended["time"])==1
    @test count(==(short),extended["time"])==1
    @test .37 in cold["time"] && .74 in cold["time"]
end
println("FRESH_CHECKPOINT_AND_EXTENSION_SUCCESS")
