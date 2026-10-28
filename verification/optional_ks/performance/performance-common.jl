ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="false"
ENV["SPACEAGORA_GRAM_TRACK_CACHE_MODE"]="off"
using SpaceAGORA,TOML,LinearAlgebra,SHA
include(joinpath(pkgdir(SpaceAGORA),"verification/optional_ks/baseline-fixtures.jl"))
BLAS.set_num_threads(1)
Threads.nthreads()==1 && BLAS.get_num_threads()==1 || error("Use one Julia and BLAS thread")
const BENCH_TIMES=collect(0.:60.:86400.)
const BENCH_ROOT=get(ENV,"KS_QUALIFICATION_ROOT",@__DIR__)
const BENCH_KEYS=("position_m","velocity_m_s")
filehash(path)=bytes2hex(sha256(read(path)))
function physical_sources()
    sources=Dict{String,String}()
    for (name,root) in (("SpaceAGORA",pkgdir(SpaceAGORA)),("KSPropagator",dirname(dirname(Base.find_package("KSPropagator")))))
        for base in (name=="SpaceAGORA" ? ("src","ext") : ("src",))
            for (dir,_,files) in walkdir(joinpath(root,base)), file in files
                path=joinpath(dir,file)
                sources[name*"/"*relpath(path,root)]=filehash(path)
            end
        end
        sources[name*"/Project.toml"]=filehash(joinpath(root,"Project.toml"))
    end
    return sources
end
function experiment_identity()
    gate=TOML.parsefile(joinpath(@__DIR__,"performance-correctness-gate.toml"))
    source=physical_sources()
    source==gate["physical_source_sha256"] || error("Current source differs from qualified physical tree")
    evidence=Dict(name=>filehash(joinpath(BENCH_ROOT,name)) for name in keys(gate["evidence_sha256"]))
    evidence==gate["evidence_sha256"] || error("Correctness evidence changed")
    environment=dirname(Base.active_project())
    manifest=TOML.parsefile(joinpath(environment,"Manifest.toml"))
    only(manifest["deps"]["KSPropagator"])["repo-rev"]==gate["ks_installed_commit"] || error("Wrong KS source revision")
    versions=Dict(name=>join([string(get(entry,"version","stdlib")) for entry in entries],",") for (name,entries) in manifest["deps"])
    scripts=Dict(name=>filehash(joinpath(@__DIR__,name)) for name in
        ("performance-common.jl","performance-calibration.jl","performance-worker.jl","performance-summary.jl","run-performance-process.py","performance-plan.md","performance-correctness-gate.toml"))
    fixture=Dict(name=>filehash(joinpath(pkgdir(SpaceAGORA),"verification/optional_ks",name)) for name in ("baseline-fixtures.jl","baseline-inputs.toml"))
    return Dict("physical_sources"=>source,"correctness_evidence"=>evidence,"scripts"=>scripts,"fixtures"=>fixture,
        "environment"=>Dict(n=>filehash(joinpath(environment,n)) for n in ("Project.toml","Manifest.toml")),
        "resolved_versions"=>versions,"julia"=>string(VERSION),"blas"=>string(BLAS.get_config()),
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "physical_settings"=>Dict(k=>get(ENV,k,"") for k in ("SPACEAGORA_DENSITY_FREEZE_PER_STEP","SPACEAGORA_GRAM_TRACK_CACHE_MODE")),
        "sample_grid"=>"0:60:86400 seconds; 1441 samples; max1m/.001m/s; mass1e-9kg; heat frozen budget")
end
function solution_contract(sol,mode)
    complete=sol.t[end]==86400. && SE.SciMLBase.successful_retcode(sol)
    result=Dict{String,Any}("complete"=>complete,"endpoint_s"=>Float64(sol.t[end]),"retcode"=>string(sol.retcode),"clock_applicable"=>mode===:ks)
    if mode===:ks
        clock=sol.alg.diagnostics.max_clock_residual_s
        result["clock_residual_s"]=clock
        result["complete"] &= isfinite(clock) && clock<=1e-7
    end
    return result
end
function bench_metrics(sol,reference)
    pe=0.;ve=0.;me=0.;qr=0.
    for t in BENCH_TIMES
        u=Vector(sol(t));r=Vector(reference(t))
        all(isfinite,u) && all(isfinite,r) || error("Nonfinite sample")
        pe=max(pe,norm(u[1:3]-r[1:3]));ve=max(ve,norm(u[4:6]-r[4:6]))
        me=max(me,abs(u[7]-r[7]));qr=max(qr,abs(u[8]-r[8])/max(1e-6,1e-5*abs(r[8])))
    end
    return Dict("position_m"=>pe,"velocity_m_s"=>ve,"mass_kg"=>me,"heat_budget_ratio"=>qr)
end
bench_pass(m)=all(isfinite,values(m)) && m["position_m"]<=1. && m["velocity_m_s"]<=.001 && m["mass_kg"]<=1e-9 && m["heat_budget_ratio"]<=1.
positive_errors(m)=all(isfinite(m[key]) && m[key]>0 for key in BENCH_KEYS)
bench_match(m,k)=bench_pass(m) && positive_errors(m) && positive_errors(k) && all(.5<=m[key]/k[key]<=2. for key in BENCH_KEYS)
bench_score(m,k)=positive_errors(m) && positive_errors(k) ? sum(log(m[key]/k[key])^2 for key in BENCH_KEYS) : Inf
