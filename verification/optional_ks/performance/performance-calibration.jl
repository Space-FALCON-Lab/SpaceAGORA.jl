# Accuracy-only calibration; no speed measurements.
include("performance-common.jl")
length(ARGS)==2 || error("usage: performance-calibration.jl BASELINE_RESULTS_TOML NEW_OUTPUT_DIRECTORY")
function calibrate(baseline_path,out)
    baseline=TOML.parsefile(baseline_path)
    baseline_path==joinpath(BENCH_ROOT,"baseline-attempt1/results.toml") || error("Use the source-bound baseline record")
    count(r->r["preset"]=="fine" && r["passed"],baseline["rows"])==4 || error("Four fine prerequisites required")
    isdir(out) && !isempty(readdir(out)) && error("Refusing to overwrite calibration")
    mkpath(out)
    identity=experiment_identity()
    rows=Dict{String,Any}[];selected=Dict{String,Any}();ks=Dict{String,Any}()
    function persist(status)
        identity==experiment_identity() || error("Experiment drift; frozen identity will not be relabeled")
        open(joinpath(out,"calibration.toml"),"w") do io
            TOML.print(io,Dict("status"=>status,"identity"=>identity,"rows"=>rows,"selected"=>selected,"KS_current"=>ks))
        end
    end
    for case in (:two_body,:J2_plus_drag)
        ref=full_reference_solution(case;rtol=2.5e-13,atol=2.5e-13,dtmax=10.)
        # Recompute the fine target on the CURRENT installed source, retaining the historical verdict separately.
        current=run_simulation(baseline_config(case,:ks;rtol=1e-12,atol=1e-13,dtmax=5.);return_solution=true)
        contract=solution_contract(current,:ks);contract["complete"] || error("Current KS calibration incomplete")
        k=bench_metrics(current,ref);bench_pass(k) && positive_errors(k) || error("Current fine KS fails benchmark prerequisites")
        prior=only(filter(r->r["case"]==string(case) && r["preset"]=="fine",baseline["rows"]))
        ks[string(case)]=Dict("accuracy"=>k,"contract"=>contract,"historical_accuracy"=>prior["independent_error"])
        persist("calibration in progress")
        println("CURRENT_KS_ACCURACY ",case," ",k);flush(stdout)
        tried=Float64[];matches=Dict{String,Any}[]
        function attempt(rtol)
            push!(tried,rtol)
            row=Dict{String,Any}("case"=>string(case),"rtol"=>rtol,"atol"=>rtol/10,"dtmax_s"=>600.,"status"=>"started","matched_accuracy"=>false)
            push!(rows,row);persist("calibration in progress")
            try
                normal=run_simulation(baseline_config(case,:dp8;rtol,atol=rtol/10,dtmax=600.);return_solution=true)
                c=solution_contract(normal,:dp8);row["contract"]=c
                if c["complete"]
                    m=bench_metrics(normal,ref);row["accuracy"]=m;row["passed_common_target"]=bench_pass(m)
                    row["matched_accuracy"]=bench_match(m,k);row["log_score"]=bench_score(m,k)
                    row["status"]="complete"
                    row["matched_accuracy"] && push!(matches,row)
                else
                    row["status"]="incomplete solution; rejected before interpolation"
                end
            catch err
                row["status"]="execution error";row["error"]=sprint(showerror,err)
                persist("calibration aborted on recorded exception")
                rethrow() # Preserve unexpected errors; never hide them as useful accuracy points.
            end
            persist("calibration in progress")
            println("CALIBRATION ",case," ",row);flush(stdout)
        end
        for rtol in (1e-10,3e-10,1e-9,3e-9,1e-8,3e-8,1e-7,3e-7,1e-6)
            attempt(rtol)
        end
        for _ in 1:3
            !isempty(matches) && break
            prior_grid=sort(copy(tried))
            for (a,b) in zip(prior_grid[1:end-1],prior_grid[2:end]);attempt(sqrt(a*b));end
        end
        if isempty(matches)
            persist("no matched candidate; no speed ratio justified")
            error("No accuracy match for $case")
        end
        selected[string(case)]=first(sort(matches;by=r->r["log_score"]))
        persist("calibration in progress")
    end
    persist("matched accuracy selected; no speed measurements")
    println("MATCHED_ACCURACY_CALIBRATION_SUCCESS")
end
calibrate(abspath(ARGS[1]),abspath(ARGS[2]))
