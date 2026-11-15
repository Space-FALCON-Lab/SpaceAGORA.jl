module ContactHybridTests
using Test, LinearAlgebra, StaticArrays, SpaceAGORA
include("../helpers/contact_fixture.jl")
const C=SpaceAGORA.ContactLoadServices
const Z=SVector(0.,0.,0.)
struct HybridContact <: C.AbstractHybridContactModel
    kind::Symbol
end
mutable struct HybridRuntime <: C.AbstractHybridContactRuntime
    kind::Symbol
    projected::Bool
    released::Bool
    accepted::Vector{Float64}
    release_time::Float64
    closed::Bool
end
const CREATED=HybridRuntime[]
C.contact_preflight(::HybridContact,::NTuple{2,C.ContactBody})=nothing
function C.contact_prepare(m::HybridContact,::NTuple{2,C.ContactBody})
    rt=HybridRuntime(m.kind,false,false,Float64[],NaN,false)
    push!(CREATED,rt);rt
end
C.contact_close!(rt::HybridRuntime)=(rt.closed=true;nothing)
function C.contact_evaluate!(rt::HybridRuntime,states::NTuple{2,C.ContactState},t::Float64)
    f=rt.released ? SVector(.1,0.,0.) : Z
    (C.ContactWrench(-f,Z),C.ContactWrench(f,Z))
end
C.contact_event_value(rt::HybridRuntime,states,t::Float64)=rt.released ? 1. : .003-t
function C.contact_accept!(rt::HybridRuntime,states::NTuple{2,C.ContactState},t::Float64)
    push!(rt.accepted,t)
    if !rt.released && t>=.003-1e-12
        rt.released=true;rt.release_time=t
    end
    if rt.kind==:late && (t<.003-1e-12 || length(rt.accepted)<2 || rt.accepted[end-1]!=t)
        return (C.ContactVelocityCorrection(Z,Z),C.ContactVelocityCorrection(Z,Z))
    end
    if !rt.projected
        rt.projected=true
        relative=states[2].velocity_m_s[1]-states[1].velocity_m_s[1]
        delta=relative/21
        rt.kind==:energy && (delta=-delta)
        a=C.ContactVelocityCorrection(SVector(delta,0.,0.),Z)
        b=C.ContactVelocityCorrection(SVector(-20delta,0.,0.),Z)
        rt.kind==:momentum && (a=C.ContactVelocityCorrection(Z,Z))
        return (a,b)
    end
    (C.ContactVelocityCorrection(Z,Z),C.ContactVelocityCorrection(Z,Z))
end
@testset "accepted hybrid transitions and conservative projection" begin
    base=fixture()
    sol=run_simulation(C.with_contact_model(base,HybridContact(:valid));return_solution=true,visualization=false)
    rt=last(CREATED)
    @test rt.closed && rt.projected && rt.released
    @test first(rt.accepted)==0.
    @test issorted(rt.accepted)
    @test rt.release_time≈.003 atol=1e-10
    @test last(sol.t)==.01
    run=sol.prob.p.shared_buffers.contact_runtime[]
    @test run.projection_count[]==1
    @test only(run.projections).time_s==0.
    @test run.projection_work_J[]≈-20/21*.25^2/2 atol=1e-15
    @test only(run.projections).work_J==run.projection_work_J[]
    for kind in (:momentum,:energy)
        @test_throws ArgumentError run_simulation(C.with_contact_model(base,HybridContact(kind));return_solution=true,visualization=false)
        @test last(CREATED).closed
    end
end
@testset "capture after a continuous event saves the corrected state" begin
    sol=run_simulation(C.with_contact_model(fixture(),HybridContact(:late));return_solution=true,visualization=false)
    rt=last(CREATED);run=sol.prob.p.shared_buffers.contact_runtime[]
    @test rt.projected && rt.closed
    @test only(run.projections).time_s≈.003 atol=1e-10
    @test abs(sol(rt.release_time;continuity=:right).sc[1].att_v[1,1])<1e-14
end
end
