module SpaceAGORAKSPropagatorExt
using SpaceAGORA, KSPropagator
import OrdinaryDiffEqCore as OC
const SE = SpaceAGORA.SimulationEngine

function require_capability()
    isdefined(KSPropagator,:coupled_cartesian_step_doubling) &&
        isdefined(KSPropagator,:CoupledMidpointFailure) || throw(ArgumentError(
            "This KSPropagator source does not provide the qualified coupled stepping API. Install the documented integration source; legacy 0.1 sources are insufficient."))
    return nothing
end

mutable struct Diagnostics
    attempts::Int
    nonlinear_failures::Int
    rhs_calls::Int
    max_clock_residual_s::Float64
end
struct CoupledKS <: OC.OrdinaryDiffEqAdaptiveAlgorithm
    mu::Float64
    position_indices::NTuple{3,Int}
    velocity_indices::NTuple{3,Int}
    diagnostics::Diagnostics
end
OC.alg_order(::CoupledKS) = 2
OC.alg_adaptive_order(::CoupledKS) = 2

function algorithm(prob,args)
    require_capability()
    z=zero(prob.u0)
    z.sc[1].pos .= 1.
    positions=Tuple(findall(!iszero,Vector(z)))
    fill!(z,0.);z.sc[1].vel .= 1.
    velocities=Tuple(findall(!iszero,Vector(z)))
    length(positions)==3 && length(velocities)==3 || throw(ArgumentError("KS requires a three-dimensional Cartesian state."))
    return CoupledKS(Float64(args.environment_model.planet.μ),positions,velocities,Diagnostics(0,0,0,0.))
end

mutable struct Cache{U,R} <: OC.OrdinaryDiffEqMutableCache
    u::U
    uprev::U
    tmp::U
    k::R
    fsalfirst::R
    error::U
end
OC.full_cache(c::Cache) = (c.u,c.uprev,c.tmp,c.k,c.fsalfirst,c.error)
OC.get_fsalfirstlast(c::Cache,u) = (c.fsalfirst,c.k)
function OC.alg_cache(alg::CoupledKS,u,rate_prototype,::Type{UT},::Type{BT},::Type{TT},
        uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true},verbose) where {UT,BT,TT}
    return Cache(u,uprev,zero(u),zero(rate_prototype),zero(rate_prototype),zero(u))
end

function OC.initialize!(integrator,c::Cache)
    integrator.kshortsize=2
    integrator.fsalfirst=c.fsalfirst
    integrator.fsallast=c.k
    resize!(integrator.k,2)
    integrator.k[1]=c.fsalfirst;integrator.k[2]=c.k
    SE._trial_cartesian_rhs!(c.fsalfirst,integrator.uprev,integrator.p,Float64(integrator.t))
    OC.increment_nf!(integrator.stats,1)
    return nothing
end

function OC.perform_step!(integrator,c::Cache,repeat_step=false)
    (;uprev,u,p,t,dt,alg)=integrator
    alg.diagnostics.attempts+=1
    calls=Ref(0)
    function rhs!(dy,y,time)
        calls[]+=1
        SE._trial_cartesian_rhs!(dy,y,p,Float64(time))
    end
    # Even after a callback or repeated candidate at the same time, the start
    # derivative must correspond to this physical state, not a time-only cache.
    rhs!(integrator.fsalfirst,uprev,t)
    result=nothing
    try
        if !p.is_active[1]
            all(iszero,integrator.fsalfirst) || throw(ArgumentError("Inactive KS spacecraft has nonzero dynamics outside the initial support contract."))
            u .= uprev
            fill!(c.error,0.)
            OC.set_EEst!(integrator,0.)
        else
            result=KSPropagator.coupled_cartesian_step_doubling(rhs!,uprev,t,dt,alg.mu;
                position_indices=alg.position_indices,velocity_indices=alg.velocity_indices)
            u .= result.state
            alg.diagnostics.max_clock_residual_s=max(alg.diagnostics.max_clock_residual_s,result.clock_error)
            for i in eachindex(c.error)
                at=integrator.opts.abstol isa Number ? integrator.opts.abstol : integrator.opts.abstol[i]
                rt=integrator.opts.reltol isa Number ? integrator.opts.reltol : integrator.opts.reltol[i]
                denominator=at+rt*max(abs(uprev[i]),abs(u[i]))
                e=result.error_estimate[i]
                c.error[i]=denominator>0 ? e/denominator : (iszero(e) ? 0. : Inf)
            end
            OC.set_EEst!(integrator,integrator.opts.internalnorm(c.error,t))
        end
        # Refresh accepted-candidate state, force and heat/output buffers. The
        # outer SciML controller alone decides acceptance and runs callbacks.
        rhs!(integrator.fsallast,u,t+dt)
    catch err
        if err isa KSPropagator.CoupledMidpointFailure
            alg.diagnostics.nonlinear_failures+=1
            alg.diagnostics.nonlinear_failures<=3 && @debug "KS candidate nonlinear failure" t dt exception=err
            integrator.force_stepfail=true
            OC.set_EEst!(integrator,Inf)
        else
            rethrow()
        end
    finally
        OC.increment_nf!(integrator.stats,calls[])
        alg.diagnostics.rhs_calls+=calls[]
    end
    return nothing
end
end
