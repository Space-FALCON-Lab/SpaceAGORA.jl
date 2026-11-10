function _cartesian_states(state)
    x=zeros(6,length(state)÷10)
    for s in axes(x,2)
        c=ks_state_to_cartesian(view(state,slice(s)))
        x[:,s].=vcat(c.position_ii_m,c.velocity_ii_m)
    end
    x
end

"""Receding-horizon MPC over the fixed-reference terminal-gain scheduler.

At each step, `plan_laser_schedule` plans over the next `horizon` intervals of
`times` (fewer near the end) from the current nonlinear state, with its all-off
reference from that state. The first `execute` intervals are applied in the
nonlinear KS replay, and planning restarts from that replay's final state.

Each entry of `segments` holds the segment's first interval, the planner's
linear change in the five elements at the end of the executed intervals and the
nonlinear change produced by the same commands. Both changes are relative to
that segment's all-off reference. `final_cartesian` is the chained replay's
state at `last(times)`, one 6-vector column per spacecraft.
"""
function run_laser_mpc(case,times;horizon,execute,objective=:a,linearization=:ks)
    n=length(times)-1
    1<=execute<=horizon || throw(ArgumentError("Need 1 <= execute <= horizon."))
    command=zeros(Int,length(case.edges),n);segments=NamedTuple[]
    segment_case=case;k=0
    while k<n
        m=min(execute,n-k)
        window=times[k+1:min(k+horizon,n)+1];executed=times[k+1:k+m+1]
        plan=plan_laser_schedule(segment_case,window;objective,linearization)
        u=plan.command[:,1:m];command[:,k+1:k+m].=u
        y=replay_laser_schedule(segment_case,executed,u)(last(executed))
        actual=elements(view(y,slice(1)),case.config)-
            elements(view(plan.reference(last(executed)),slice(1)),case.config)
        actual[3:5].=atan.(sin.(actual[3:5]),cos.(actual[3:5]))
        push!(segments,(;first_interval=k+1,predicted=plan.predicted[:,m+1],actual))
        segment_case=merge(segment_case,(;initial_cartesian=_cartesian_states(y)))
        k+=m
    end
    return (;command,segments,final_cartesian=segment_case.initial_cartesian,times=collect(times))
end
