using SpaceAGORA, LinearAlgebra, CSV, DataFrames, Printf
BLAS.set_num_threads(1)

# Long-horizon behaviour of the fixed-reference terminal-gain scheduler and the
# receding-horizon MPC around it. Every accuracy figure compares the linear
# prediction with a nonlinear KS replay of the same schedule.
#  1. horizon_sweep.csv:   one plan per horizon, linear vs nonlinear terminal change.
#  2. geometry_drift.csv:  for the longest one-shot plan, how far each commanded
#                          link's replayed force direction is from the planned one.
#  3. mpc_summary.csv and mpc_segments.csv: MPC executing one orbit or the whole
#     horizon per plan, its per-segment prediction check, and one continuous
#     replay of the executed schedule.
# Gains are against the no-laser trajectory, averaged over the last orbit: J2
# short-period terms make terminal osculating a unreliable once targets drift.

const SCALE=[1.,1.,180/pi,180/pi,180/pi]
element_index(o)=o==:a ? 1 : o==:i ? 3 : 5
wrap!(d)=(d[3:5,:].=atan.(sin.(d[3:5,:]),cos.(d[3:5,:]));d)

function run_constellation_laser_mpc(;check=false)
    case=walker27_constellation()
    per_orbit=check ? 24 : 240;orbits=check ? 3. : 100.;horizon=check ? 48 : 1200
    sweep=check ? [.5,1.,2.] : [1.,5.,10.,20.,40.,80.,100.]
    objectives=check ? (:a,) : (:a,:i,:omega)
    out=joinpath(@__DIR__,"..","outputs",check ? "constellation_laser_mpc_smoke" : "constellation_laser_mpc")
    mkpath(out)
    grid(orb)=collect(range(0.,orb*case.period_s;length=round(Int,orb*per_orbit)+1))

    sweep_rows=NamedTuple[];drift_rows=NamedTuple[]
    for orb in sweep
        times=grid(orb);reference=replay_laser_schedule(case,times,zeros(Int,length(case.edges),length(times)-1))
        for objective in objectives
            j=element_index(objective)
            plan=plan_laser_schedule(case,times;objective,reference)
            replay=replay_laser_schedule(case,times,plan.command)
            d=wrap!(laser_element_history(case,[last(times)],replay)-laser_element_history(case,[last(times)],reference))
            linear=SCALE[j]*plan.predicted[j,end];nonlinear=SCALE[j]*d[j,end]
            push!(sweep_rows,(;orbits=orb,objective=string(objective),linear,nonlinear,
                relative_error=abs(linear-nonlinear)/abs(nonlinear)))
            @printf("one-shot %6.1f orbits %-6s linear=%.6g nonlinear=%.6g\n",orb,objective,linear,nonlinear)
            orb==last(sweep) || continue
            for o in 0:ceil(Int,orb)-1
                samples=NTuple{4,Float64}[]
                for k in o*per_orbit+1:min((o+1)*per_orbit,length(times)-1)
                    e=findfirst(!iszero,plan.command[:,k]);e===nothing && continue
                    t=(times[k]+times[k+1])/2
                    g0=laser_geometry(reference(t),case;t)[e];g1=laser_geometry(replay(t),case;t)[e]
                    f0,f1=g0.target_force,g1.target_force
                    push!(samples,(Float64(!g1.available),rad2deg(acos(clamp(dot(f0,f1)/(norm(f0)*norm(f1)),-1,1))),
                        norm(f1)/norm(f0),abs(g1.distance-g0.distance)))
                end
                isempty(samples) && continue
                mean(i)=sum(s[i] for s in samples)/length(samples)
                push!(drift_rows,(;objective=string(objective),orbit=o,unavailable_fraction=mean(1),
                    mean_force_angle_deg=mean(2),max_force_angle_deg=maximum(s[2] for s in samples),
                    mean_force_ratio=mean(3),mean_range_change_m=mean(4)))
            end
        end
    end
    CSV.write(joinpath(out,"horizon_sweep.csv"),DataFrame(sweep_rows))
    CSV.write(joinpath(out,"geometry_drift.csv"),DataFrame(drift_rows))

    times=grid(orbits);last_orbit=times[end-per_orbit:end]
    nominal=laser_element_history(case,last_orbit,replay_laser_schedule(case,times,zeros(Int,length(case.edges),length(times)-1)))
    summary_rows=NamedTuple[];segment_rows=NamedTuple[]
    for objective in objectives,execute in (per_orbit,horizon)
        j=element_index(objective)
        r=run_laser_mpc(case,times;horizon,execute,objective)
        continuous=replay_laser_schedule(case,times,r.command)
        d=wrap!(laser_element_history(case,last_orbit,continuous)-nominal)
        y=continuous(last(times))
        restart_gap=maximum(norm(r.final_cartesian[1:3,s]-ks_state_to_cartesian(view(y,10s-9:10s)).position_ii_m)
            for s in axes(r.final_cartesian,2))
        predicted=[SCALE[j]*s.predicted[j] for s in r.segments];actual=[SCALE[j]*s.actual[j] for s in r.segments]
        # Relative error only where the segment moved the element appreciably.
        moved=abs.(actual).>=.01maximum(abs.(actual));err=sort(abs.(predicted-actual)[moved]./abs.(actual[moved]))
        push!(summary_rows,(;objective=string(objective),horizon,execute,replans=length(r.segments),
            active_intervals=sum(r.command),gain_last_orbit_mean=SCALE[j]*sum(d[j,:])/size(d,2),
            gain_terminal=SCALE[j]*d[j,end],segment_rel_error_median=err[cld(end,2)],
            segment_rel_error_max=maximum(err),chained_vs_continuous_position_m=restart_gap))
        for (s,p,a) in zip(r.segments,predicted,actual)
            push!(segment_rows,(;objective=string(objective),execute,first_interval=s.first_interval,predicted=p,actual=a))
        end
        @printf("MPC %-6s execute=%4d gain(last-orbit mean)=%.6g  segment error median=%.2e max=%.2e  restart gap=%.2e m\n",
            objective,execute,summary_rows[end].gain_last_orbit_mean,summary_rows[end].segment_rel_error_median,
            summary_rows[end].segment_rel_error_max,restart_gap)
        flush(stdout)
    end
    CSV.write(joinpath(out,"mpc_summary.csv"),DataFrame(summary_rows))
    CSV.write(joinpath(out,"mpc_segments.csv"),DataFrame(segment_rows))
    println("Results: $(abspath(out))")
end

run_constellation_laser_mpc(;check="--check" in ARGS)
