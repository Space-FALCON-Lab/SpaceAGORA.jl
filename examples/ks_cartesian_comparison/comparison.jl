module NativeCoordinateComparison
using SpaceAGORA, OrdinaryDiffEq, DiffEqCallbacks, LinearAlgebra, CSV, DataFrames
using OrdinaryDiffEqLowOrderRK: DP5
const CS=SpaceAGORA.ConstellationScheduling
const KS=SpaceAGORA.SimulationModel.DynamicsKS
const FT=SpaceAGORA.SimulationModel.FrameTransforms
const OUTPUT=joinpath(@__DIR__,"..","..","outputs","ks_cartesian_comparison")

function ten_initial(e)
    case=walker27_constellation();ids=[1,4,7,10,13,15,18,21,24,27]
    elements=copy(case.initial_elements[ids,:]);elements[:,6].=e
    # Preserve the original comparison's perigee distances as eccentricity changes.
    elements[:,5].*=.99/(1-e)
    x=zeros(6,10)
    for j in 1:10
        a,_,i,node,arg,M=elements[j,5:10];E=M
        for _ in 1:15;E-=(E-e*sin(E)-M)/(1-e*cos(E));end
        nu=2atan(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2))
        r,v=FT.orbitalelemtorv(CS.SVector(a,e,i,node,arg,nu,0.),case.config.planet)
        x[:,j].=vcat(r,v)
    end
    (;x,elements,planet=case.config.planet,period=2pi*sqrt(elements[1,5]^3/case.config.planet.μ))
end

to_ks(x,p;t=0.)=vcat([Vector(KS.cartesian_to_ks_state(x[1:3,j],x[4:6,j],p;elapsed_time_s=t)) for j in axes(x,2)]...)
function to_cart(y)
    reduce(hcat,[begin
        q=KS.ks_state_to_cartesian(view(y,10j-9:10j))
        vcat(q.position_ii_m,q.velocity_ii_m)
    end for j in 1:length(y)÷10])
end

function cart_rhs!(du,u,p,t)
    for j in 1:length(u)÷6
        k=6j-6;r=view(u,k+1:k+3)
        du[k+1:k+3].=view(u,k+4:k+6)
        du[k+4:k+6].=-p.μ*r/norm(r)^3+KS.ks_j2_acceleration_si(r,p)
    end
end

function ks_rhs!(du,u,p,s)
    rho_target=sum(abs2,view(u,1:4))
    for j in 1:length(u)÷10
        rows=10j-9:10j;y=view(u,rows)
        du[rows].=KS.nonlinear_ks_rhs(y,p)*(rho_target/sum(abs2,view(y,1:4)))
    end
end

function solve_nonlinear(c,route,tf,alg;tol=1e-10,step=0.)
    if route==:cartesian
        prob=ODEProblem(cart_rhs!,vec(c.x),(0.,tf),c.planet)
        sol=solve(prob,alg;reltol=tol,abstol=tol,adaptive=step==0.,
            dt=step==0. ? 1. : step,maxiters=100000000)
    else
        # The target clock is the common independent-variable map for every satellite.
        y=to_ks(c.x,c.planet);a=c.elements[1,5]
        callback=ContinuousCallback((u,s,integrator)->u[10]-tf,
            integrator->terminate!(integrator);abstol=1e-9,reltol=0.,save_positions=(true,true))
        prob=ODEProblem(ks_rhs!,y,(0.,2tf/(a*(1-c.elements[1,6]))),c.planet)
        sol=solve(prob,alg;reltol=tol,abstol=tol,adaptive=step==0.,
            dt=(step==0. ? 1. : step)/a,callback,maxiters=100000000)
    end
    string(sol.retcode) in ("Success","Terminated") || error("Propagation failed: $(sol.retcode)")
    sol
end

function physical_sample(sol,route,t)
    route==:cartesian && return reshape(sol(t),6,:)
    lo=first(sol.t);hi=last(sol.t)
    for _ in 1:60
        mid=(lo+hi)/2
        sol(mid;idxs=10)<t ? (lo=mid) : (hi=mid)
    end
    to_cart(sol((lo+hi)/2))
end

function sample_trajectory(sol,route,times,nsat)
    X=Array{Float64}(undef,6,nsat,length(times))
    for (k,t) in enumerate(times)
        X[:,:,k].=physical_sample(sol,route,t)
    end
    X
end

function errors(X,Y)
    dp=[norm(X[1:3,j,k]-Y[1:3,j,k]) for j in axes(X,2),k in axes(X,3)]
    dv=[norm(X[4:6,j,k]-Y[4:6,j,k]) for j in axes(X,2),k in axes(X,3)]
    (;max_position_m=maximum(dp),rms_position_m=sqrt(sum(abs2,dp)/length(dp)),
        max_velocity_m_s=maximum(dv),rms_velocity_m_s=sqrt(sum(abs2,dv)/length(dv)))
end

function measure_nonlinear(c,route,tf,alg,tol,step,times,Y)
    @timed begin
        solve_stats=@timed solve_nonlinear(c,route,tf,alg;tol,step)
        sol=solve_stats.value
        X=sample_trajectory(sol,route,times,size(c.x,2))
        err=errors(X,Y)
        clock=route==:cartesian ? 0. : maximum(abs.(sol.u[end][10:10:end].-tf))
        (;solve_stats,err,clock)
    end
end

function nonlinear_runs(;horizons=(1,5),steps=(5.,20.,100.),tols=(1e-8,1e-10,1e-12),out=OUTPUT)
    mkpath(out);rows=NamedTuple[];truth_rows=NamedTuple[]
    for e in (.01,.3),orbits in horizons
        c=ten_initial(e);tf=orbits*c.period;times=collect(range(0.,tf;length=240orbits+1))
        CSV.write(joinpath(out,"initial_e$(e).csv"),DataFrame(c.elements,
            ["original_id","shell","plane","slot","a_m","e","i_rad","RAAN_rad","omega_rad","M_rad"]))
        truth=solve_nonlinear(c,:cartesian,tf,Vern9();tol=2e-14)
        refined=solve_nonlinear(c,:cartesian,tf,Vern9();tol=1e-14)
        Y=sample_trajectory(refined,:cartesian,times,size(c.x,2))
        TY=sample_trajectory(truth,:cartesian,times,size(c.x,2))
        push!(truth_rows,(;e,orbits,errors(TY,Y)...))
        for route in (:cartesian,:ks_target_s), (name,alg) in (("Tsit5",Tsit5()),("DP5",DP5())),
            (mode,value) in vcat([("adaptive",tol) for tol in tols],[("fixed",s) for s in steps])
            tol=mode=="adaptive" ? value : 1e-12;step=mode=="fixed" ? value : 0.
            warm_times=collect(range(0.,min(tf,100.);length=3))
            WY=sample_trajectory(refined,:cartesian,warm_times,size(c.x,2))
            measure_nonlinear(c,route,min(tf,100.),alg,tol,step,warm_times,WY)
            GC.gc()
            total=measure_nonlinear(c,route,tf,alg,tol,step,times,Y)
            s=total.value.solve_stats
            push!(rows,(;e,orbits,route=string(route),integrator=name,mode,setting=value,
                solve_s=s.time,total_s=total.time,solve_gc_s=s.gctime,total_gc_s=total.gctime,
                solve_bytes=s.bytes,total_bytes=total.bytes,steps=length(s.value.t)-1,
                clock_error_s=total.value.clock,total.value.err...))
            CSV.write(joinpath(out,"nonlinear.csv"),DataFrame(rows))
        end
        CSV.write(joinpath(out,"truth_refinement.csv"),DataFrame(truth_rows))
        println("Nonlinear comparison finished e=$e, $orbits orbits");flush(stdout)
    end
end

"""Direct Cartesian replay, using exactly the planner's native interlink force API."""
function cart_replay(case,times,command)
    model=deepcopy(case.interlinks);c=case.config
    function rhs!(du,u,p,t)
        cart_rhs!(du,u,c.planet,t)
        k=clamp(searchsortedlast(times,t),1,size(command,2))
        if any(!iszero,view(command,:,k))
            x=reshape(u,6,:)
            native=(;sc=[(;pos=x[1:3,j],vel=x[4:6,j],mass=c.mass_kg) for j in axes(x,2)])
            ctx=CS.link_context(case,t)
            CS.update_availability!(model,case.spacecraft,native;context=ctx)
            model.graph.active.=(command[:,k].==1).&model.graph.available
            for j in axes(x,2)
                du[6j-2:6j].+=CS.interlink_force_on_spacecraft(model,native,j;
                    spacecraft=case.spacecraft,context=ctx)/c.mass_kg
            end
        end
    end
    sol=solve(ODEProblem(rhs!,vec(case.initial_cartesian),(first(times),last(times))),Tsit5();
        reltol=c.tolerance,abstol=c.tolerance,dtmax=minimum(diff(times))/4,tstops=times[2:end-1])
    string(sol.retcode)=="Success" || error("Cartesian replay failed")
    sol
end

function cart_elements(sol,t,case)
    x=sol(t);collect(FT.rvtoorbitalelement(CS.SVector{3}(x[1:3]),CS.SVector{3}(x[4:6]),case.config.planet)[1:5])
end
function element_delta(a,b)
    d=a-b;d[3:5].=atan.(sin.(d[3:5]),cos.(d[3:5]));d
end

"""Target-only inertial ZOH predictions on the selected binary schedule.

KS retains its native fictitious-time tangent and input lift. Physical-time
switch boundaries require removing the clock perturbation after each interval.
Helpers enter the input geometry on the no-laser reference, as in the original
first-order scheduler. This is not a coupled active-reference model.
"""
function inertial_predictions(case,times,plan;substeps=1)
    p=case.config.planet;reference=plan.reference
    a=case.initial_elements[1,5]
    D=Diagonal(vcat(fill(sqrt(a),4),fill(sqrt(p.μ),4),p.μ/a,sqrt(a^3/p.μ)))
    dy=zeros(10);dx=zeros(6);rows=NamedTuple[]
    for k in 1:length(times)-1
      for sub in 1:substeps
        dt=(times[k+1]-times[k])/substeps
        t0=times[k]+(sub-1)*dt;t1=times[k]+sub*dt;tm=(t0+t1)/2
        full=reference(tm);ym=full[1:10];y0=reference(t0)[1:10];y1=reference(t1)[1:10]
        x=KS.ks_state_to_cartesian(ym);cm=CS.cartesian_model(x.position_ii_m,x.velocity_ii_m,p)
        geom=CS.laser_geometry(full,case;t=tm)
        accelerations=reduce(hcat,[g.target_force/case.config.mass_kg for g in geom])
        selected=accelerations*plan.command[:,k]
        F,V=CS.zoh(cm.A_cartesian,reshape(cm.B_cartesian*selected,6,1),dt)
        dx=F*dx+vec(V)
        As=Matrix(KS.evaluate_ks_dynamics_and_jacobians(ym,p;compute_input_jacobian=Val(false)).state)
        Bs=reshape(KS.ks_acceleration_input_matrix(ym)*selected,10,1)
        ds=dt/6*(1/sum(abs2,y0[1:4])+4/sum(abs2,ym[1:4])+1/sum(abs2,y1[1:4]))
        Fs,Vs=CS.zoh(D\As*D,D\Bs,ds)
        dy=D*(Fs*(D\dy)+vec(Vs))
        # A fixed-s sensitivity includes an arrival-time perturbation. Shift it
        # along the nominal flow to the requested physical activation boundary.
        dy-=KS.nonlinear_ks_rhs(y1,p)/sum(abs2,y1[1:4])*dy[10]
        endmap=KS.ks_physical_tangent(y1,p);xe=KS.ks_state_to_cartesian(y1)
        C=CS.orbital_element_jacobian(xe.position_ii_m,xe.velocity_ii_m,p.μ)
        kz=C*endmap.H*dy;cz=C*dx;scale=[1.,1.,180/pi,180/pi,180/pi]
        if sub==substeps
          for (j,element) in enumerate(("a_m","e","i_deg","RAAN_deg","omega_deg"))
            push!(rows,(;time_orbits=t1/case.period_s,element,ks_native_s=scale[j]*kz[j],
                cartesian_inertial=scale[j]*cz[j],ks_clock_perturbation_s=dy[10]))
          end
        end
      end
    end
    DataFrame(rows)
end

function scheduling_runs(;orbits=5,intervals=240orbits,out=OUTPUT)
    mkpath(out);case=walker27_constellation();times=collect(range(0.,orbits*case.period_s;length=intervals+1))
    off=zeros(Int,length(case.edges),intervals)
    reference=replay_laser_schedule(case,times,off)
    cart_reference=cart_replay(case,times,off)
    # Warm compilation before timing either implementation.
    for route in (:ks,:cartesian)
        plan_laser_schedule(case,times[1:3];reference,linearization=route)
    end
    rows=NamedTuple[];timings=NamedTuple[];scale=[1.,1.,180/pi,180/pi,180/pi]
    names=["a_m","e","i_deg","RAAN_deg","omega_deg"]
    for objective in (:a,:i,:omega,:combined)
        plans=Dict()
        for route in (:ks,:cartesian)
            GC.gc();stats=@timed plan_laser_schedule(case,times;reference,objective,linearization=route)
            plans[route]=stats.value
            push!(timings,(;objective=string(objective),route=string(route),stage="planning",
                wall_s=stats.time,gc_s=stats.gctime,bytes=stats.bytes))
        end
        kp=plans[:ks];cp=plans[:cartesian]
        differences=count(k->kp.command[:,k]!=cp.command[:,k],1:intervals)
        folder=joinpath(out,"maximize_$objective");mkpath(folder)
        for (route,plan) in ((:ks,kp),(:cartesian,cp))
            CSV.write(joinpath(folder,"schedule_$route.csv"),DataFrame(
                hcat(times[1:end-1],times[2:end],plan.command'),
                vcat(["start_s","end_s"],["$(s)_to_$t" for (s,t) in case.edges])))
        end
        # Replay the KS-selected schedule in both coordinates, even if a tied choice differs.
        ks_stats=@timed replay_laser_schedule(case,times,kp.command)
        cart_stats=@timed cart_replay(case,times,kp.command)
        for (route,st) in ((:ks,ks_stats),(:cartesian,cart_stats))
            push!(timings,(;objective=string(objective),route=string(route),stage="replay",
                wall_s=st.time,gc_s=st.gctime,bytes=st.bytes))
        end
        history=NamedTuple[]
        for (k,t) in enumerate(times)
            kd=element_delta(CS.elements(view(ks_stats.value(t),1:10),case.config),
                CS.elements(view(reference(t),1:10),case.config))
            cd=element_delta(cart_elements(cart_stats.value,t,case),cart_elements(cart_reference,t,case))
            for j in 1:5
                record=(;objective=string(objective),element=names[j],time_orbits=t/case.period_s,
                    ks_linear=scale[j]*kp.predicted[j,k],cartesian_linear=scale[j]*cp.predicted[j,k],
                    ks_nonlinear=scale[j]*kd[j],cartesian_nonlinear=scale[j]*cd[j])
                push!(history,record)
                k==length(times) && push!(rows,(;record...,different_intervals=differences,
                    max_gain_difference=maximum(abs.(kp.gains[j,:,:]-cp.gains[j,:,:]))*scale[j]))
            end
        end
        CSV.write(joinpath(folder,"history.csv"),DataFrame(history))
        refinement=DataFrame()
        for substeps in (1,4,16)
            inertial=inertial_predictions(case,times,kp;substeps)
            substeps==1 && CSV.write(joinpath(folder,"inertial_zoh.csv"),inertial)
            lastrows=copy(inertial[end-4:end,:]);lastrows.substeps=fill(substeps,5)
            append!(refinement,lastrows)
        end
        CSV.write(joinpath(folder,"inertial_refinement.csv"),refinement)
        CSV.write(joinpath(out,"terminal.csv"),DataFrame(rows))
        CSV.write(joinpath(out,"scheduling_timing.csv"),DataFrame(timings))
        println("Scheduling comparison finished $objective ($differences different intervals)");flush(stdout)
    end
end

"""Receding-horizon MPC (`run_laser_mpc`) with KS versus Cartesian linearization.

Defaults follow examples/Constellation_Laser_MPC.jl: five-orbit planning horizon,
one orbit executed per replan. Segment errors compare each plan's linear delta a
with the nonlinear KS replay of the executed commands. Each executed schedule is
also replayed continuously in KS and Cartesian against that route's no-laser run.
"""
function mpc_runs(;orbits=5,horizon=1200,execute=240,out=OUTPUT)
    mkpath(out);case=walker27_constellation();n=round(Int,240orbits)
    times=collect(range(0.,orbits*case.period_s;length=n+1));off=zeros(Int,length(case.edges),n)
    reference=replay_laser_schedule(case,times,off);cart_reference=cart_replay(case,times,off)
    for route in (:ks,:cartesian)  # warm compilation
        run_laser_mpc(case,times[1:5];horizon=4,execute=2,linearization=route)
    end
    rows=NamedTuple[]
    for objective in (:a,:i,:omega,:combined)
        commands=Dict()
        for route in (:ks,:cartesian)
            GC.gc();stats=@timed run_laser_mpc(case,times;horizon,execute,objective,linearization=route)
            r=stats.value;commands[route]=r.command
            predicted=[s.predicted[1] for s in r.segments];actual=[s.actual[1] for s in r.segments]
            ks=replay_laser_schedule(case,times,r.command);cart=cart_replay(case,times,r.command)
            tf=last(times)
            ks_da=element_delta(CS.elements(view(ks(tf),1:10),case.config),CS.elements(view(reference(tf),1:10),case.config))[1]
            cart_da=element_delta(cart_elements(cart,tf,case),cart_elements(cart_reference,tf,case))[1]
            push!(rows,(;objective=string(objective),route=string(route),orbits,horizon,execute,
                replans=length(r.segments),wall_s=stats.time,gc_s=stats.gctime,bytes=stats.bytes,
                segment_max_abs_error_a_m=maximum(abs.(predicted-actual)),
                segment_sum_predicted_a_m=sum(predicted),segment_sum_actual_a_m=sum(actual),
                ks_nonlinear_a_m=ks_da,cartesian_nonlinear_a_m=cart_da,
                active_intervals=sum(r.command),different_intervals=0))
        end
        d=count(k->commands[:ks][:,k]!=commands[:cartesian][:,k],1:n)
        rows[end]=merge(rows[end],(;different_intervals=d));rows[end-1]=merge(rows[end-1],(;different_intervals=d))
        CSV.write(joinpath(out,"mpc.csv"),DataFrame(rows))
        println("MPC comparison finished $objective, $orbits orbits ($d different intervals)");flush(stdout)
    end
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    using .NativeCoordinateComparison
    NativeCoordinateComparison.BLAS.set_num_threads(1)
    mode=isempty(ARGS) ? "short" : ARGS[1]
    mode in ("short","full","nonlinear","schedule","extended","mpc") || error("Choose short, full, nonlinear, schedule, extended or mpc")
    if mode=="extended"
        NativeCoordinateComparison.nonlinear_runs(;horizons=(1,5,50,100,200),
            out=joinpath(NativeCoordinateComparison.OUTPUT,"extended_horizons"))
    end
    if mode in ("short","full","nonlinear")
        NativeCoordinateComparison.nonlinear_runs(;horizons=mode=="full" ? (1,5,50,100,200) : (1,5))
    end
    if mode in ("short","full","schedule")
        NativeCoordinateComparison.scheduling_runs()
    end
    mode=="mpc" && NativeCoordinateComparison.mpc_runs()
end
