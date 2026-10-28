using TOML
length(ARGS)==2 || error("usage: performance-summary.jl PERFORMANCE_OUTPUT_DIRECTORY NEW_SUMMARY_TOML")
root=abspath(ARGS[1]);out=abspath(ARGS[2]);isfile(out) && error("Refusing to overwrite summary")
function summarize(root,out)
records=Dict{String,Any}();identity=nothing;calibration_id=nothing
median3(xs)=sort(xs)[2]
for case in ("two_body","J2_plus_drag")
    methods=Dict{String,Any}();runs=Dict{String,Any}()
    for mode in ("ks","dp8")
        r=TOML.parsefile(joinpath(root,case*"-"*mode,"timing.toml"))
        monitor=TOML.parsefile(joinpath(root,case*"-"*mode*"-monitor.toml"))
        monitor["status"]=="terminal clean observations" && monitor["exit_code"]==0 && monitor["background_activity_valid"] || error("Worker background/process gate invalid")
        monitor["command"][end-2]==case && monitor["command"][end-1]==mode && abspath(monitor["command"][end])==joinpath(root,case*"-"*mode) || error("Misassociated process monitor")
        if isnothing(calibration_id);calibration_id=r["calibration_sha256"];else;calibration_id==r["calibration_sha256"] || error("Mixed calibrations");end
        r["case"]==case && r["method"]==mode || error("Misplaced case/method evidence")
        r["status"]=="complete worker; paired accuracy summary still required" || error("Worker incomplete")
        if isnothing(identity)
            identity=r["identity"]
        else
            identity==r["identity"] || error("Mixed experiment identities")
        end
        all(x->x["matched_accuracy"] && x["contract"]["complete"],r["runs"]) || error("Unqualified timed run")
        stable=filter(x->get(x,"stable_execution",false),r["runs"])
        length(stable)==3 || error("Exactly three stable repeats required")
        runs[mode]=stable
        wall=[x["wall_s"] for x in stable]
        methods[mode]=Dict("median_wall_s"=>median3(wall),"min_wall_s"=>minimum(wall),"max_wall_s"=>maximum(wall),"stable_runs"=>stable,"first_use_and_warmups"=>filter(x->!get(x,"stable_execution",false),r["runs"]),"package_load"=>r["package_load"],"hardware"=>r["hardware"])
    end
    # Every stable cross-pair must match BOTH actual achieved error metrics.
    ratios=Dict{String,Float64}()
    for key in ("position_m","velocity_m_s")
        a=[r["accuracy"][key] for r in runs["ks"]];b=[r["accuracy"][key] for r in runs["dp8"]]
        worst=max(maximum(a)/minimum(b),maximum(b)/minimum(a))
        isfinite(worst) && worst<=2. || error("Actual paired $key differs by more than2x; no matched-accuracy speed ratio")
        ratios[key]=worst
    end
    records[case]=Dict("methods"=>methods,"worst_cross_pair_error_ratio"=>ratios,"KS_over_DP8_median_wall_ratio"=>methods["ks"]["median_wall_s"]/methods["dp8"]["median_wall_s"])
end
open(out,"w") do io;TOML.print(io,Dict("status"=>"complete bounded matched-accuracy comparison","identity"=>identity,"calibration_sha256"=>calibration_id,"cases"=>records));end
println("ACTUAL_PAIRED_ACCURACY_SUMMARY_SUCCESS")

end
summarize(root,out)
