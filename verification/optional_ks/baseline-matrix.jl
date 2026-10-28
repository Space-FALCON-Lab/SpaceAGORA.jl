# Run with --project=integration/environment. Frozen inputs are baseline-inputs.toml.
ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="false"
ENV["SPACEAGORA_GRAM_TRACK_CACHE_MODE"]="off"
include("baseline-fixtures.jl")
using KSPropagator, TOML
const OUT=get(ENV,"KS_BASELINE_OUTPUT",joinpath(@__DIR__,"baseline-attempt1"))
isdir(OUT) && !isempty(readdir(OUT)) && error("Refusing to overwrite previous baseline evidence: $OUT")
mkpath(OUT)
const TIMES=collect(0.:60.:86400.)
function metrics(sol,reference)
    pe=0.;ve=0.;me=0.;qe=0.;qratio=0.
    for t in TIMES
        u=Vector(sol(t));r=Vector(reference(t))
        pe=max(pe,norm(u[1:3]-r[1:3]));ve=max(ve,norm(u[4:6]-r[4:6]))
        if length(u)>=8 && length(r)>=8
            me=max(me,abs(u[7]-r[7]));qe=max(qe,abs(u[8]-r[8]))
            qratio=max(qratio,abs(u[8]-r[8])/max(1e-6,1e-5*abs(r[8])))
        end
    end
    return Dict("position_m"=>pe,"velocity_m_s"=>ve,"mass_kg"=>me,
        "heat_J_cm2"=>qe,"heat_budget_ratio"=>qratio)
end
function raw_save(path,sol)
    open(path,"w") do io
        println(io,"physical_time_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,mass_kg,heat_J_cm2")
        for t in TIMES
            println(io,join(vcat(t,Vector(sol(t))),','))
        end
    end
end
function run_matrix()
    records=Dict{String,Any}[]
    for case in (:two_body,:J2_only,:drag_only,:J2_plus_drag)
        println("BASELINE_BEGIN ",case);flush(stdout)
        reference_coarse=full_reference_solution(case;rtol=1e-12,atol=1e-12,dtmax=20.)
        reference=full_reference_solution(case;rtol=2.5e-13,atol=2.5e-13,dtmax=10.)
        refinement=metrics(reference_coarse,reference)
        refgood=refinement["position_m"]<=.01 && refinement["velocity_m_s"]<=1e-5 &&
            refinement["heat_budget_ratio"]<=.01
        analytical=case===:two_body ? metrics(reference,analytic_kepler) : Dict{String,Any}()
        if case===:two_body
            refgood &= analytical["position_m"]<=.01 && analytical["velocity_m_s"]<=1e-5
        end
        raw_save(joinpath(OUT,"$(case)-reference.csv"),reference)
        println("REFERENCE ",case," refined=",refgood," ",refinement," analytic=",analytical);flush(stdout)
        # Normal SpaceAGORA DP8 route, same configured force models/inputs.
        normal=run_simulation(baseline_config(case,:dp8;rtol=2.5e-13,atol=2.5e-13,dtmax=10.);return_solution=true)
        normal_metrics=metrics(normal,reference)
        raw_save(joinpath(OUT,"$(case)-normal-DP8.csv"),normal)
        println("NORMAL_DP8 ",case," ",normal_metrics);flush(stdout)
        for (preset,rtol,atol,dtmax) in (("coarse",1e-8,1e-9,20.),("medium",1e-10,1e-11,10.),("fine",1e-12,1e-13,5.))
            println("KS_BEGIN ",case," ",preset);flush(stdout)
            sol=nothing
            try
                sol=run_simulation(baseline_config(case,:ks;rtol,atol,dtmax);return_solution=true)
            catch err
                record=Dict{String,Any}("case"=>string(case),"preset"=>preset,"passed"=>false,
                    "status"=>"execution error", "error"=>sprint(showerror,err),
                    "reference_qualified"=>refgood,"reference_refinement"=>refinement,
                    "normal_route_error"=>normal_metrics,"rtol"=>rtol,"atol"=>atol,"dtmax_s"=>dtmax)
                push!(records,record)
                open(joinpath(OUT,"results.toml"),"w") do io
                    TOML.print(io,Dict("status"=>"in progress","julia"=>string(VERSION),"rows"=>records))
                end
                println("KS_ERROR ",case," ",preset," ",record["error"]);flush(stdout)
                continue
            end
            m=metrics(sol,reference);parity=metrics(sol,normal)
            clock=sol.alg.diagnostics.max_clock_residual_s
            passed=refgood && m["position_m"]<=1. && m["velocity_m_s"]<=.001 &&
                m["mass_kg"]<=1e-9 && m["heat_budget_ratio"]<=1. && clock<=1e-7 &&
                sol.t[end]==86400.
            raw_save(joinpath(OUT,"$(case)-$(preset)-KS.csv"),sol)
            record=Dict{String,Any}("case"=>string(case),"preset"=>preset,"passed"=>passed,
                "reference_qualified"=>refgood,"independent_error"=>m,"normal_route_error"=>normal_metrics,
                "KS_vs_normal"=>parity,"reference_refinement"=>refinement,"reference_analytic"=>analytical,
                "clock_residual_s"=>clock,"attempts"=>sol.alg.diagnostics.attempts,
                "nonlinear_failures"=>sol.alg.diagnostics.nonlinear_failures,"rhs_calls_steps_only"=>sol.alg.diagnostics.rhs_calls,
                "rtol"=>rtol,"atol"=>atol,"dtmax_s"=>dtmax)
            push!(records,record)
            open(joinpath(OUT,"results.toml"),"w") do io
                TOML.print(io,Dict("status"=>"in progress","julia"=>string(VERSION),"rows"=>records))
            end
            println("KS_RESULT ",case," ",preset," passed=",passed," ",m," clock=",clock);flush(stdout)
        end
    end
    open(joinpath(OUT,"results.toml"),"w") do io
        TOML.print(io,Dict("status"=>"complete bounded matrix; review required","julia"=>string(VERSION),"rows"=>records))
    end
    println("BASELINE_MATRIX_TERMINAL rows=",length(records));flush(stdout)
end
run_matrix()
