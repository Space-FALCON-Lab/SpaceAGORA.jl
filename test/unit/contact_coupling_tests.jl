module ContactCouplingTests
using Test, LinearAlgebra, StaticArrays, SpaceAGORA
const C = SpaceAGORA.ContactLoadServices
const S = SpaceAGORA.SimulationEngine
const SM = SpaceAGORA.SimulationModel
const Z = SVector(0.,0.,0.)
const Q = SVector(0.,0.,0.,1.)

struct FakeContact <: C.AbstractContactModel
    mode::Symbol
end
mutable struct FakeRuntime <: C.AbstractContactRuntime
    mode::Symbol
    closed::Bool
    observations::Vector{Tuple{NTuple{2,C.ContactState},Float64}}
end
const CREATED = FakeRuntime[]
C.contact_preflight(::FakeContact, bodies::NTuple{2,C.ContactBody}) = nothing
function C.contact_prepare(model::FakeContact, bodies::NTuple{2,C.ContactBody})
    runtime = FakeRuntime(model.mode, false, Tuple{NTuple{2,C.ContactState},Float64}[])
    push!(CREATED,runtime)
    return runtime
end
C.contact_close!(rt::FakeRuntime) = (rt.closed=true; nothing)
function C.contact_evaluate!(rt::FakeRuntime, states::NTuple{2,C.ContactState}, t::Float64)
    rt.closed && error("closed runtime")
    push!(rt.observations,(states,t))
    rt.mode === :throw && error("intentional provider failure")
    rt.mode === :badtype && return nothing
    rt.mode === :nan && return (C.ContactWrench(SVector(NaN,0.,0.),Z),C.ContactWrench(Z,Z))
    rt.mode === :unbalanced && return (C.ContactWrench(Z,Z),C.ContactWrench(SVector(10.,0.,0.),Z))
    rt.mode === :zero && return (C.ContactWrench(Z,Z),C.ContactWrench(Z,Z))
    f = SVector(10.,0.,0.)
    point = states[1].position_m + SVector(.5,.2,0.)
    return (C.ContactWrench(-f,cross(point-states[1].position_m,-f)),
        C.ContactWrench(f,cross(point-states[2].position_m,f)))
end

function fixture(; h=2.0^-12, duration=.01, position=SVector(7e6,0.,0.))
    bus = Link(root=true,m=20.,dims=MVector(1.,.8,.6))
    node = CompliantTopologyNode(:bumper; mass_kg=1.,inertia_body_kg_m2=.004,
        position=SVector(.12,0.,0.),velocity=SVector(-.25,0.,0.))
    edge = CompliantTopologyEdge(:mount,0,1; parent_point_body=Z, child_point_body=SVector(-.11,0.,0.),
        k_translation_n_m=200.,c_translation_n_s_m=.2,k_rotation_n_m_rad=2.,c_rotation_n_m_s_rad=.02)
    attachment = CompliantAttachment(model=build_compliant_topology([node],[edge]),link=bus,mount_point=SVector(.5,0.,0.))
    sc = SpacecraftModel(links=[bus],root=bus,prop_mass=0.,inertia_tensor=bus.inertia,
        initial_condition=CartesianInitialCondition(position,SVector(.25/21,0.,0.)),attachments=[attachment])
    base = make_example_config(planet=SM.make_no_gram_planet(:earth),spacecraft=sc,
        mission_time=duration,initial_time=InitialTime(year=2020,month=1,day=1),
        dynamic_effectors=(),ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=true,keplerian=false,verbose=false,results=false,
        solver_config=SolverConfig(solver_mode=:tsit5))
    return SM.SimConfig._with_configuration(base;
        mission_configuration=SM.SimConfig._with_fields(base.mission_configuration; data_rate=duration/2),
        integration_tolerances=IntegrationTolerances(dt_max_orbit=h,dt_max_atmosphere=h,
            reltol_orbit=1e-10,abstol_orbit=1e-11,reltol_quaternion=1e-10,abstol_quaternion=1e-11))
end

@testset "experimental contact service and refusals" begin
    SpaceAGORA.Conformance.test_service_facade(C; version=v"0.1.0")
    @test_throws ArgumentError C.ContactBody(-1.,Matrix{Float64}(I,3,3))
    @test_throws ArgumentError C.ContactBody(1.,zeros(3,3))
    args=C.with_contact_model(fixture(),FakeContact(:zero))
    @test S._validate_run_models!(args,:tsit5) === nothing
    copywith=SM.SimConfig._with_configuration
    # Existing interlink examples pass a bare scheduling policy by keyword.
    # The contact field must preserve that convenience and every provider.
    bare=copywith(args;scheduling_policy_model=SemimajorAxisRatePolicy())
    @test bare.scheduling_policy_model isa SchedulingPolicyModel
    @test bare.contact_models==args.contact_models
    @test isempty(copywith(fixture();scheduling_policy_model=SemimajorAxisRatePolicy()).contact_models)
    @test_throws ArgumentError S._validate_run_models!(args,:rodas5p)
    @test_throws ArgumentError S._validate_run_models!(copywith(args;contact_models=(:invalid,)),:tsit5)
    @test_throws ArgumentError S._validate_run_models!(copywith(args;contact_models=(FakeContact(:zero),FakeContact(:zero))),:tsit5)
    @test_throws ArgumentError S._validate_run_models!(copywith(args;simulation_settings=SimulationSettings(checkpoint_enabled=true)),:tsit5)
    @test_throws ArgumentError S._validate_run_models!(C.with_contact_model(fixture(h=.1),FakeContact(:zero)),:tsit5)
    @test_throws ArgumentError run_simulation(args;extra_callbacks=(:unsupported,))
end

@testset "smooth loads precede contact; each COM wrench enters once" begin
    args=C.with_contact_model(fixture(),FakeContact(:force))
    u=S.build_initial_conditions(args)
    p=S._prepare_run_runtime!(args,u,:tsit5)
    owner=Ref{Union{Nothing,C.ContactRunState}}(nothing)
    S._initialize_contact_runtime!(p,u,owner)
    rt=owner[].runtime
    try
        du=zero(u); S.spacecraft_dynamics!(du,u,p,0.)
        states,t=last(rt.observations)
        @test t==0.
        @test states[1].position_m==Z
        @test states[2].position_m≈SVector(.62,0.,0.)
        @test states[2].velocity_m_s≈SVector(-.25,0.,0.)
        @test states[1].smooth_force_n≈SVector(1.95,0.,0.)
        @test states[1].smooth_force_n + states[2].smooth_force_n≈Z
        @test du.sc[1].vel≈SVector((1.95-10)/20,0.,0.)
        @test du.sc[1].att_v[:,1]+du.sc[1].vel≈SVector(10-1.95,0.,0.)
        @test args.dynamics_model.spacecraft[1].inertia_tensor*du.sc[1].ω≈SVector(0.,0.,2.)
        @test .004du.sc[1].att_ω[:,1]≈SVector(0.,0.,-2.)
        before=copy(u)
        other=copy(u);other.sc[1].att_r[1,1]+=.02
        S.spacecraft_dynamics!(zero(u),other,p,2.)
        again=zero(u);S.spacecraft_dynamics!(again,u,p,0.)
        @test again==du
        @test u==before
        @test first(rt.observations)==last(rt.observations)
        owner[].busy[]=true
        @test_throws ArgumentError S.spacecraft_dynamics!(zero(u),u,p,0.)
        owner[].busy[]=false
    finally
        C.contact_close!(rt)
    end
    @test rt.closed
end

@testset "zero contact identity, independent runtimes and cleanup" begin
    args=fixture()
    plain=run_simulation(args;return_solution=true,visualization=false)
    coupled=C.with_contact_model(args,FakeContact(:zero))
    a=run_simulation(coupled;return_solution=true,visualization=false)
    ra=last(CREATED)
    b=run_simulation(coupled;return_solution=true,visualization=false)
    rb=last(CREATED)
    @test a.u==plain.u && a.t==plain.t
    @test a.u==b.u && a.t==b.t
    @test ra!==rb && ra.closed && rb.closed
    @test !isempty(ra.observations) && ra.observations==rb.observations
    for mode in (:throw,:badtype,:nan,:unbalanced)
        @test_throws Exception run_simulation(C.with_contact_model(args,FakeContact(mode));return_solution=true,visualization=false)
        @test last(CREATED).closed
    end
end
end
