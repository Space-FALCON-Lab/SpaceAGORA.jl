using SpaceAGORA, JuMP, HiGHS, LinearAlgebra, CSV, DataFrames, Printf
BLAS.set_num_threads(1)

# Battery, thermal and data-constrained terminal-gain scheduling on the Walker-27
# constellation: the eight cases of Constellation-Controllability's
# scripts/walker27_resource_scheduling.jl. Each schedule is checked against the
# closed-form resource replay and replayed nonlinearly in KS; resources are also
# replayed on the actual (laser-on) orbit's environment. Every hardware coefficient
# comes from the parts walker27_resource_constellation mounts (resource_model).
#   julia --project=examples examples/Constellation_Resource_Scheduling.jl [--check]

const DISPLAY=[1.,1.,180/pi,180/pi,180/pi]  # a in m, e, angles in deg
const NAMES=["a_m","e","i_deg","RAAN_deg","omega_deg"]

function run_constellation_resource_scheduling(;check=false,time_limit=check ? 10. : 120.)
    case=walker27_resource_constellation();orbits=check ? .25 : 5.
    fine=collect(range(0.,orbits*case.period_s;length=round(Int,240orbits)+1))
    out=joinpath(@__DIR__,"..","outputs",check ? "constellation_resource_scheduling_smoke" : "constellation_resource_scheduling")
    mkpath(out)
    # The five-orbit delivery requirement cannot fit a quarter orbit.
    p=check ? resource_model(case;terminal_bits=vcat(15e9,fill(.61e9,26))) : resource_model(case)
    specs=[("unconstrained_a",false,false,false,:a),("battery_a",true,false,false,:a),
        ("thermal_a",false,true,false,:a),("data_a",false,false,true,:a),
        ("all_a",true,true,true,:a),("all_i",true,true,true,:i),
        ("all_omega",true,true,true,:omega),("all_joint",true,true,true,:combined)]
    check && (specs=specs[[1,5]])
    reference=replay_laser_schedule(case,fine,zeros(Int,length(case.edges),length(fine)-1))
    summary=NamedTuple[];terminal=NamedTuple[]
    for (name,battery,thermal,data,objective) in specs
        r=plan_resource_schedule(case,fine;objective,battery,thermal,data,parameters=p,reference,time_limit)
        if !r.feasible
            push!(summary,(;case=name,status=r.status));println("$name: $(r.status)");continue
        end
        states=resource_replay(r.u,r.v,r.flow,r.downlink,case.edges,r.times,r.env,p;solar=r.solar)
        battery && @assert maximum(abs,states.energy-r.energy)<1e-5
        thermal && @assert maximum(abs,states.temperature-r.temperature)<1e-5
        data && @assert maximum(abs,states.queue-r.queue)<1e-5
        delivered=dot(r.downlink,diff(r.times))
        @assert abs(sum(states.queue[:,end])+delivered-sum(p.initial_Gbit)-last(r.times)*sum(p.generation_Gbps))<1e-5
        actual=replay_laser_schedule(case,fine,r.command)
        d=laser_element_history(case,[last(fine)],actual)-laser_element_history(case,[last(fine)],reference)
        d[3:5].=atan.(sin.(d[3:5]),cos.(d[3:5]))
        # Resources on the environment the commanded (not the reference) orbit sees.
        actual_states=resource_replay(r.u,r.v,r.flow,r.downlink,case.edges,r.times,
            resource_environment(case,r.times,actual),p;solar=states.solar)
        violation(x)=max(0.,x...)
        push!(summary,(;case=name,status=r.status,relative_gap=r.gap,solve_s=r.seconds,
            active_blocks=sum(r.u),objective=r.objective,upper_bound=r.bound,
            min_energy_kWh=minimum(states.energy),max_temperature_K=maximum(states.temperature),
            max_queue_Gbit=maximum(states.queue),delivered_Gbit=delivered,
            actual_energy_violation_kWh=violation((maximum(p.reserve_kWh.-actual_states.energy),maximum(actual_states.energy.-p.capacity_kWh))),
            actual_temperature_violation_K=violation((maximum(actual_states.temperature.-p.maximum_K),maximum(p.minimum_K.-actual_states.temperature)))))
        for j in 1:5
            push!(terminal,(;case=name,element=NAMES[j],linear=DISPLAY[j]*r.predicted[j],ks_nonlinear=DISPLAY[j]*d[j]))
        end
        CSV.write(joinpath(out,"schedule_$name.csv"),DataFrame(hcat(r.times[1:end-1],r.times[2:end],r.u'),
            vcat(["start_s","end_s"],["$(s)_to_$t" for (s,t) in case.edges])))
        @printf("%-16s %-10s gap=%.3g active=%3d  linear Δa=%9.3f m  KS Δa=%9.3f m  min E=%.3f kWh  max T=%.2f K\n",
            name,r.status,r.gap,sum(r.u),r.predicted[1],d[1],minimum(states.energy),maximum(states.temperature))
        flush(stdout)
        CSV.write(joinpath(out,"summary.csv"),vcat((DataFrame([s]) for s in summary)...;cols=:union))
        CSV.write(joinpath(out,"terminal_comparison.csv"),DataFrame(terminal))
    end
    CSV.write(joinpath(out,"summary.csv"),vcat((DataFrame([s]) for s in summary)...;cols=:union))
    println("Results: $(abspath(out))")
end

run_constellation_resource_scheduling(;check="--check" in ARGS)
