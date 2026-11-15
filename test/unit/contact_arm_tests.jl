module ContactArmTests
using Test, LinearAlgebra, StaticArrays, SpaceAGORA
const C=SpaceAGORA.ContactLoadServices
const S=SpaceAGORA.SimulationEngine
const SM=SpaceAGORA.SimulationModel
const CA=SM.ClothRobotArmDynamics
const CM=SM.ClothMultibody
const Z=SVector(0.,0.,0.)
const Q=SVector(0.,0.,0.,1.)

struct ArmContact <: C.AbstractContactModel
    force::SVector{3,Float64}
end
mutable struct ArmRuntime <: C.AbstractContactRuntime
    force::SVector{3,Float64}
    observations::Vector{Any}
end
C.contact_preflight(::ArmContact,::NTuple{N,C.ContactBody}) where {N}=nothing
C.contact_prepare(c::ArmContact,::NTuple{N,C.ContactBody}) where {N}=ArmRuntime(c.force,Any[])
C.contact_close!(::ArmRuntime)=nothing
function C.contact_evaluate!(rt::ArmRuntime,states::NTuple{N,C.ContactState},t::Float64) where {N}
    push!(rt.observations,(states,t))
    point=last(states).position_m+SVector(.04,.03,0.)
    ntuple(N) do i
        f=i==1 ? -rt.force : i==N ? rt.force : Z
        C.ContactWrench(f,cross(point-states[i].position_m,f))
    end
end

function fixture(;duration=.01)
    bus=Link(root=true,m=20.,dims=MVector(1.,.8,.6))
    arm=SM.default_robot_arm_model(link_lengths_m=(.3,.25,.2),link_radii_m=(.03,.025,.02),
        link_masses_kg=(.7,.5,.3),mount_offset_body=(.5,0.,0.))
    pose=SM.RobotArmBasePose(Z,Q)
    plan=SM.RobotArmPlanBuffer(arm,pose,[.2,-.3,.1]).plan
    effector=SM.RobotArmReactionEffector(plan=plan,k_translation_n_m=200.,
        c_translation_n_s_m=.2,k_rotation_n_m_rad=2.,c_rotation_n_m_s_rad=.02)
    sc=SpacecraftModel(links=[bus],root=bus,prop_mass=0.,inertia_tensor=bus.inertia,
        initial_condition=CartesianInitialCondition(SVector(7e6,0.,0.),Z))
    args=make_example_config(planet=SM.make_no_gram_planet(:earth),spacecraft=sc,
        mission_time=duration,initial_time=InitialTime(year=2020,month=1,day=1),
        dynamic_effectors=(effector,),ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=true,keplerian=false,verbose=false,results=false,
        solver_config=SolverConfig(solver_mode=:tsit5))
    SM.SimConfig._with_configuration(args;
        mission_configuration=SM.SimConfig._with_fields(args.mission_configuration;data_rate=duration/2),
        integration_tolerances=IntegrationTolerances(dt_max_orbit=2.0^-12,dt_max_atmosphere=2.0^-12,
            reltol_orbit=1e-10,abstol_orbit=1e-11,reltol_quaternion=1e-10,abstol_quaternion=1e-11))
end

@testset "all compliant arm links participate before derivatives" begin
    base=fixture(); args=C.with_contact_model(base,ArmContact(SVector(-.4,.3,.2)))
    @test S._validate_run_models!(args,:tsit5)===nothing
    bodies=S._contact_bodies(args)
    @test length(bodies)==4
    @test [b.mass_kg for b in bodies]==[20.,.7,.5,.3]
    u=S.build_initial_conditions(args)
    u.sc[1].arm_r[2,2]+=.001
    u.sc[1].arm_v[:,3].=[.01,-.02,.03]
    p=S._prepare_run_runtime!(args,u,:tsit5)
    plain=S._prepare_run_runtime!(base,u,:tsit5)
    owner=Ref{Union{Nothing,C.ContactRunState}}(nothing)
    S._initialize_contact_runtime!(p,u,owner)
    try
        du=zero(u); baseline=zero(u)
        S.spacecraft_dynamics!(baseline,u,plain,0.)
        S.spacecraft_dynamics!(du,u,p,0.)
        states,t=last(owner[].runtime.observations)
        @test length(states)==4 && t==0.
        @test norm(sum(s.smooth_force_n for s in states))<1e-12
        moment=sum(s.smooth_torque_nm+cross(s.position_m, s.smooth_force_n) for s in states)
        @test norm(moment)<1e-10
        w=C.contact_evaluate!(owner[].runtime,states,0.)
        @test du.sc[1].vel-baseline.sc[1].vel≈w[1].force_n/20 atol=1e-12
        for i in 1:3
            b=bodies[i+1]; q=states[i+1].quaternion
            @test du.sc[1].arm_v[:,i]-baseline.sc[1].arm_v[:,i]≈w[i+1].force_n/b.mass_kg atol=1e-12
            @test du.sc[1].arm_ω[:,i]-baseline.sc[1].arm_ω[:,i]≈b.inertia_body_kg_m2\(CM._rot(q)'*w[i+1].torque_nm) atol=1e-11
        end
        before=copy(u); altered=copy(u);altered.sc[1].arm_r[1,2]+=.02
        S.spacecraft_dynamics!(zero(u),altered,p,1.)
        replay=zero(u);S.spacecraft_dynamics!(replay,u,p,0.)
        @test replay==du && u==before
    finally
        C.contact_close!(owner[].runtime)
    end
end

@testset "zero contact preserves the existing robot-arm route" begin
    args=fixture()
    a=run_simulation(args;return_solution=true,visualization=false)
    b=run_simulation(C.with_contact_model(args,ArmContact(Z));return_solution=true,visualization=false)
    @test a.t==b.t
    @test a.u==b.u
end
@testset "multiple attachments subtract the final bus acceleration once" begin
    bus=Link(root=true,m=20.,dims=MVector(1.,.8,.6))
    attachments=map((.2,-.2)) do offset
        node=CompliantTopologyNode(:bumper;mass_kg=1.,inertia_body_kg_m2=.004,position=SVector(.12,0.,0.))
        edge=CompliantTopologyEdge(:mount,0,1;parent_point_body=Z,child_point_body=SVector(-.11,0.,0.),
            k_translation_n_m=200.,c_translation_n_s_m=.2,k_rotation_n_m_rad=2.,c_rotation_n_m_s_rad=.02)
        CompliantAttachment(model=build_compliant_topology([node],[edge]),link=bus,mount_point=SVector(.5,offset,0.))
    end
    sc=SpacecraftModel(links=[bus],root=bus,prop_mass=0.,inertia_tensor=bus.inertia,
        initial_condition=CartesianInitialCondition(SVector(7e6,0.,0.),Z),attachments=collect(attachments))
    base=SM.SimConfig._with_configuration(fixture();dynamics_model=DynamicsModel([sc],()))
    args=C.with_contact_model(base,ArmContact(SVector(-.4,.3,.2)))
    @test S._validate_run_models!(args,:tsit5)===nothing
    @test [b.mass_kg for b in S._contact_bodies(args)]==[20.,1.,1.]
    u=S.build_initial_conditions(args);p=S._prepare_run_runtime!(args,u,:tsit5)
    plain=S._prepare_run_runtime!(base,u,:tsit5);owner=Ref{Union{Nothing,C.ContactRunState}}(nothing)
    S._initialize_contact_runtime!(p,u,owner)
    try
        du=zero(u);baseline=zero(u)
        S.spacecraft_dynamics!(baseline,u,plain,0.);S.spacecraft_dynamics!(du,u,p,0.)
        states,_=last(owner[].runtime.observations);loads=C.contact_evaluate!(owner[].runtime,states,0.)
        bus_delta=du.sc[1].vel-baseline.sc[1].vel
        @test bus_delta≈loads[1].force_n/20 atol=1e-12
        for i in 1:2
            @test du.sc[1].att_v[:,i]-baseline.sc[1].att_v[:,i]+bus_delta≈loads[i+1].force_n atol=1e-12
            @test .004*(du.sc[1].att_ω[:,i]-baseline.sc[1].att_ω[:,i])≈loads[i+1].torque_nm atol=1e-12
        end
    finally
        C.contact_close!(owner[].runtime)
    end
end
end
