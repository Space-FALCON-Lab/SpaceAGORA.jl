isdefined(@__MODULE__,:orbit_initial) || include("baseline-fixtures.jl")
# Independent eight-state Cartesian model: no SpaceAGORA force/control hooks.
function lifecycle_reference_rhs!(dy,y,phase,t)
    reference_rhs!(dy,y,:J2_only,t)
    r=@view y[1:3];v=@view y[4:6];mass=y[7]
    rho=RHO_REF*exp(-(norm(r)-RE-H_REF)/SCALE_HEIGHT)
    relative=v-cross([0.,0.,SPIN],r)
    dy[4:6].-=.5*rho*(2.2*phase.area)/mass*norm(relative).*relative
    dy[4:6].+=1e-7*sin(t/700.)*(r/7e6)+[phase.held_force_N/mass,2e-8*cos(t/900.),0.]
    if phase.burning
        dy[4:6].+=(.2/mass)*(v/norm(v))
        dy[7]=-.2/(300.0*9.80665)
    else
        dy[7]=0.
    end
    dy[8]=reference_heat(y)
    return nothing
end
function lifecycle_reference(;coincident=false,rtol=2.5e-13,atol=2.5e-13,dtmax=5.)
    duration=coincident ? 442. : 6000.
    actions=coincident ? [221.] : [611.25,1879.125,3599.75]
    boundaries=sort!(unique!(vcat(0.,collect(13.:13.:prevfloat(duration)),actions,duration)))
    r,v=orbit_initial();state=vcat(r,v,500.,0.)
    solutions=Any[];roots=Float64[]
    for i in 1:(length(boundaries)-1)
        left,right=boundaries[i:i+1]
        area=coincident ? 10. : (left>=3599.75 ? 8. : (left>=611.25 ? 15. : 10.))
        # Boundaries include exact integer13s ticks. This held value is fixed
        # throughout all stages, including the incoming segment's right edge.
        tick=13.0*floor(left/13.)
        phase=(area=area,held_force_N=.02*sin(tick/1000.),burning=!coincident && 611.25<=left<1879.125)
        problem=SE.ODEProblem(lifecycle_reference_rhs!,state,(left,right),phase)
        sol=SE.solve(problem,SE.DP8();reltol=rtol,abstol=atol,dtmax=dtmax,
            save_everystep=true,dense=true)
        push!(solutions,sol)
        if sol.u[1][3]>0. && sol.u[end][3]<=0.
            lo,hi=left,right
            for _ in 1:80
                mid=lo+(hi-lo)/2
                (mid==lo || mid==hi) && break
                if sol(mid)[3]>0.
                    lo=mid
                else
                    hi=mid
                end
            end
            push!(roots,lo+(hi-lo)/2)
        end
        state=copy(sol.u[end])
        if right==(coincident ? 221. : 1879.125)
            state[5]+=.05
        end
    end
    function at(t;continuity=:right)
        t==duration && return copy(state)
        i=searchsortedlast(boundaries,Float64(t))
        continuity===:left && i>1 && t==boundaries[i] && (i-=1)
        1<=i<=length(solutions) || throw(BoundsError(boundaries,t))
        return solutions[i](t)
    end
    return (;at,boundaries,solutions,roots,final=state)
end
lifecycle_mass_oracle(t)=500.0 -(.2/(300.0*9.80665))*clamp(t-611.25,0.,1267.875)
