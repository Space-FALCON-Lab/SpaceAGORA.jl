ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="false"
ENV["SPACEAGORA_GRAM_TRACK_CACHE_MODE"]="off"
length(ARGS)==4 || error("usage: performance-worker.jl CALIBRATION_TOML CASE METHOD NEW_OUTPUT_DIRECTORY")
const METHOD=Symbol(ARGS[3]);METHOD in (:ks,:dp8) || error("unknown method")
package_load=@timed begin
    @eval using SpaceAGORA
    METHOD===:ks && (@eval using KSPropagator)
end
include("performance-common.jl")
const CASE=Symbol(ARGS[2]);CASE in (:two_body,:J2_plus_drag) || error("unknown case")
cal=TOML.parsefile(ARGS[1]);cal["status"]=="matched accuracy selected; no speed measurements" || error("Incomplete calibration")
const SELECTED=cal["selected"][string(CASE)]
const TARGET=cal["KS_current"][string(CASE)]["accuracy"]
const IDENTITY=experiment_identity();IDENTITY==cal["identity"] || error("Calibration identity mismatch")
const OUT=abspath(ARGS[4]);isdir(OUT) && !isempty(readdir(OUT)) && error("Refusing to overwrite timing")
mkpath(OUT)
const RECORDS=Dict{String,Any}[]
rtol=METHOD===:ks ? 1e-12 : SELECTED["rtol"]
atol=METHOD===:ks ? 1e-13 : SELECTED["atol"]
dtmax=METHOD===:ks ? 5. : SELECTED["dtmax_s"]
config=baseline_config(CASE,METHOD;rtol,atol,dtmax)
timing_record(t)=Dict("wall_s"=>t.time,"compile_s"=>t.compile_time,"recompile_s"=>t.recompile_time,"allocated_bytes"=>t.bytes,"gc_s"=>t.gctime)
function persist(status)
    IDENTITY==experiment_identity() || error("Experiment identity drift")
    open(joinpath(OUT,"timing.toml"),"w") do io
        TOML.print(io,Dict("status"=>status,"case"=>string(CASE),"method"=>string(METHOD),"identity"=>IDENTITY,
            "hardware"=>Dict("cpu"=>Sys.cpu_info()[1].model,"kernel"=>string(Sys.KERNEL),"machine"=>Sys.MACHINE,"hardware_threads"=>Sys.CPU_THREADS,"total_memory_bytes"=>Sys.total_memory()),
            "calibration_sha256"=>filehash(ARGS[1]),"package_load"=>timing_record(package_load),"rtol"=>rtol,"atol"=>atol,"dtmax_s"=>dtmax,"runs"=>RECORDS))
    end
end
persist("before first use")
GC.gc()
# Dynamic evaluation keeps first-call compilation inside the measured boundary.
# This first-use record includes eval overhead; warm calls do not use eval.
first_run=@timed @eval run_simulation($config;return_solution=true)
ref=full_reference_solution(CASE;rtol=2.5e-13,atol=2.5e-13,dtmax=10.)
function record_run(timing,kind,index)
    c=solution_contract(timing.value,METHOD)
    row=merge(timing_record(timing),Dict{String,Any}("kind"=>kind,"index"=>index,"contract"=>c,"matched_accuracy"=>false))
    push!(RECORDS,row)
    if c["complete"]
        m=bench_metrics(timing.value,ref);row["accuracy"]=m;row["matched_accuracy"]=bench_match(m,TARGET)
    end
    persist("in progress")
    println("TIMING ",CASE," ",METHOD," ",row);flush(stdout)
    row["matched_accuracy"] || error("Run fails complete physical endpoint/calibrated accuracy gates")
    return row
end
record_run(first_run,"first use",0)
first_run=nothing
stable=Ref(0)
for index in 1:5
    GC.gc()
    timed=@timed run_simulation(config;return_solution=true)
    row=record_run(timed,"warm candidate",index)
    row["stable_execution"]=row["compile_s"]<=.01*row["wall_s"]
    row["stable_execution"] && (stable[]+=1)
    persist("in progress")
    timed=nothing
    stable[]>=3 && break
end
persist(stable[]>=3 ? "complete worker; paired accuracy summary still required" : "residual compilation; no execution ratio")
stable[]>=3 || error("Fewer than three stable warm calls")
println("MATCHED_ACCURACY_WORKER_SUCCESS ",CASE," ",METHOD)
