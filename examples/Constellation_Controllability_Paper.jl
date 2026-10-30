ENV["GKSwstype"]="100"
using SpaceAGORA, LinearAlgebra, CSV, DataFrames
const PAPER_PLOTS=!("--no-plots" in ARGS || "--check" in ARGS)
if PAPER_PLOTS
    using Plots
end
BLAS.set_num_threads(1)

function run_paper_constellation(;check=false)
    case=paper_walker_constellation()
    orbits=check ? .05 : 5.;n=check ? 12 : 1200
    times=collect(range(0.,orbits*case.period_s;length=n+1))
    out=joinpath(@__DIR__,"..","outputs",check ? "paper_constellation_interlinks_smoke" : "paper_constellation_interlinks")
    mkpath(out)
    reference=replay_laser_schedule(case,times,zeros(Int,length(case.edges),n))
    nominal=laser_element_history(case,times,reference)
    CSV.write(joinpath(out,"initial_conditions.csv"),DataFrame(case.initial_elements,
        ["satellite","shell","plane","slot","a_m","e","i_rad","RAAN_rad","omega_rad","M_rad"]))
    scale=[1.,1.,180/pi,180/pi,180/pi];names=["a_m","e","i_deg","RAAN_deg","omega_deg"]
    rows=NamedTuple[]
    for objective in (:a,:i,:omega,:combined)
        result=plan_laser_schedule(case,times;objective,reference)
        replay=replay_laser_schedule(case,times,result.command)
        actual=laser_element_history(case,times,replay);delta=actual-nominal
        delta[3:5,:].=atan.(sin.(delta[3:5,:]),cos.(delta[3:5,:]))
        folder=joinpath(out,"maximize_$objective");mkpath(folder)
        header=vcat(["time_s","time_orbits"],"actual_".*names,"predicted_".*names,
            "controlled_".*names,"reference_".*names)
        data=hcat(times,times/case.period_s,(scale.*delta)',(scale.*result.predicted)',
            (scale.*actual)',(scale.*nominal)')
        CSV.write(joinpath(folder,"history.csv"),DataFrame(data,header))
        CSV.write(joinpath(folder,"schedule.csv"),DataFrame(hcat(times[1:end-1],times[2:end],result.command'),
            vcat(["start_time_s","end_time_s"],["$(s)_to_$t" for (s,t) in case.edges])))
        CSV.write(joinpath(folder,"terminal_gains.csv"),DataFrame(reduce(vcat,[permutedims(
            vcat(k,e,result.available[e,k],result.scores[e,k],scale.*result.gains[:,e,k]))
            for k in 1:n for e in eachindex(case.edges)]),vcat(["interval","edge","available","score"],names)))
        violations=0
        for k in 1:n,t in range(times[k],times[k+1];length=21)[1:end-1]
            e=findfirst(!iszero,result.command[:,k]);e===nothing && continue
            violations+=!laser_geometry(replay(t),case;t)[e].available
        end
        for j in 1:5
            push!(rows,(objective=string(objective),element=names[j],predicted=scale[j]*result.predicted[j,end],
                nonlinear=scale[j]*delta[j,end],active_intervals=sum(result.command),sampled_geometry_violations=violations))
        end
        @assert isapprox(sum(result.scores.*result.command),dot(result.weights,result.predicted[:,end]);rtol=1e-8,atol=1e-8)
        if !check && PAPER_PLOTS
            labels=["Delta a (m)","Delta e","Delta i (deg)","Delta RAAN (deg)","Delta omega (deg)"]
            panels=[plot(times/case.period_s,scale[j]*delta[j,:];label="Nonlinear KS",
                ylabel=labels[j],left_margin=10Plots.mm) for j in 1:5]
            for j in 1:5
                plot!(panels[j],times/case.period_s,scale[j]*result.predicted[j,:];label="Linear",linestyle=:dash)
            end
            plot!(panels[end];xlabel="Time (orbits)")
            savefig(plot(panels...;layout=(5,1),size=(750,1000)),joinpath(folder,"element_comparison.pdf"))
            p=plot(;xlabel="Time (orbits)",ylabel="Emitter ID",legend=false,yticks=2:27,
                size=(750,450),left_margin=8Plots.mm)
            for k in 1:n
                e=findfirst(!iszero,result.command[:,k]);e===nothing && continue
                id=case.edges[e][1]
                plot!(p,times[k:k+1]/case.period_s,[id,id];linewidth=3,color=:black)
            end
            savefig(p,joinpath(folder,"schedule.pdf"))
        end
        println("$objective: linear da=$(result.predicted[1,end]) m; nonlinear da=$(delta[1,end]) m")
        flush(stdout)
    end
    CSV.write(joinpath(out,"terminal_comparison.csv"),DataFrame(rows))
    println("Results: $(abspath(out))")
end

run_paper_constellation(;check="--check" in ARGS)
