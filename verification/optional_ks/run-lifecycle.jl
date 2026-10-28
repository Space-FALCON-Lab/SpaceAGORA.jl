# Every worker and checker is a fresh process; no inherited model objects/caches.
using TOML, SHA
length(ARGS)==1 || error("usage: run-lifecycle.jl NEW_OUTPUT_DIRECTORY")
output=abspath(only(ARGS))
isdir(output) && !isempty(readdir(output)) && error("Refusing to overwrite prior evidence: $output")
mkpath(output)
project=dirname(Base.active_project())
files=Dict{String,String}()
for name in ("lifecycle-inputs.toml","baseline-fixtures.jl","lifecycle-fixtures.jl",
             "lifecycle-worker.jl","lifecycle-reference.jl","lifecycle-direct-check.jl",
             "lifecycle-checkpoint-check.jl","run-lifecycle.jl")
    files[name]=bytes2hex(sha256(read(joinpath(@__DIR__,name))))
end
for name in ("Project.toml","Manifest.toml")
    files["environment/"*name]=bytes2hex(sha256(read(joinpath(project,name))))
end
for (pkg,names) in (("SpaceAGORA",("src/orchestration/simulation/callbacks/event_callbacks.jl",
    "ext/SpaceAGORAKSPropagatorExt.jl","src/orchestration/simulation/engine/ks_solver.jl")),
    ("KSPropagator",("src/integrators/cartesian_bridge.jl",)))
    root=dirname(dirname(Base.find_package(pkg)))
    for name in names
        files[pkg*"/"*name]=bytes2hex(sha256(read(joinpath(root,name))))
    end
end
open(joinpath(output,"source.toml"),"w") do io
    TOML.print(io,Dict("julia"=>string(VERSION),"project"=>project,"files"=>files,
        "scope"=>"Frozen direct,segmented,cold-resume lifecycle, independent Cartesian refinement; no performance claim"))
end
function child(script,args...;log)
    cmd=`$(Base.julia_cmd()) --startup-file=no --project=$project $(joinpath(@__DIR__,script)) $args`
    open(joinpath(output,log),"w") do io
        run(pipeline(cmd;stdout=io,stderr=io))
    end
    println("COMPLETED ",log);flush(stdout)
end
function testcase(prefix,coincident)
    extra=coincident ? ["coincident"] : String[]
    for mode in ("dp8","ks")
        child("lifecycle-worker.jl",mode,"direct",joinpath(output,prefix*"-"*mode),extra...;
            log=prefix*"-"*mode*"-direct.log")
    end
    if !coincident
        for mode in ("dp8","ks")
            child("lifecycle-worker.jl",mode,"checkpoint",joinpath(output,"main-checkpoint-"*mode);
                log="main-"*mode*"-checkpoint.log")
        end
    end
    dir=joinpath(output,prefix*"-resume-ks")
    child("lifecycle-worker.jl","ks","interrupt",dir,extra...;log=prefix*"-interrupt.log")
    child("lifecycle-worker.jl","ks","resume",dir,extra...;log=prefix*"-resume.log")
    child("lifecycle-direct-check.jl",prefix,output;log=prefix*"-direct-check.log")
    child("lifecycle-checkpoint-check.jl",prefix,output,output;log=prefix*"-checkpoint-check.log")
end
@sync begin
    @async testcase("main",false)
    @async testcase("coincident",true)
end
println("LIFECYCLE_SUITE_SUCCESS ",output)
