module LinkAttitudeCheckpointTests
using Test, SpaceAGORA, StaticArrays, Serialization
const SE=SpaceAGORA.SimulationEngine
const C=SpaceAGORA.Conformance

function vehicle()
    make_three_body_spacecraft(bus_dims=(1.,2.,3.),panel_dims=(0.01,2.,1.),
        bus_mass=50.,panel_mass_each=5.,panel_offset_y=2.,
        ic=InitialCondition(ra=6_500_000.,rp=6_500_000.,i=0.,Ω=0.,ω=0.,ν=0.))
end
args_for(sc)=(dynamics_model=(spacecraft=[sc],),)
checkpoint(entry)=(runtime_state=(link_attitudes=entry,),)
from_file(x)=(io=IOBuffer();serialize(io,x);seekstart(io);deserialize(io))

@testset "Non-root attitude checkpoint ownership" begin
    source=vehicle(); panel=SolarPanelAngleOfAttackControlModel((2,3))
    apply_panel_angle_of_attack!(panel,source,0.27,2.)
    snapshot=from_file(SE._link_attitude_checkpoint_state(args_for(source)))
    source.links[2].q[1]=0.123
    @test snapshot.spacecraft[1].links[1].q != Tuple(source.links[2].q)
    fresh=vehicle();rootq=copy(fresh.root.q)
    entry=SE._check_link_attitude_checkpoint_state(args_for(fresh),checkpoint(snapshot))
    @test SE._restore_link_attitude_checkpoint_state!(args_for(fresh),entry) === nothing
    for record in snapshot.spacecraft[1].links
        @test Tuple(fresh.links[record.index].q) == record.q
        @test fresh.links[record.index].α == record.alpha_rad
    end
    @test fresh.root.q == rootq
    before=[(Tuple(l.q),l.α) for l in fresh.links]
    invalids=(merge(snapshot,(version=2,)),
        merge(snapshot,(spacecraft=[],)),
        merge(snapshot,(spacecraft=[merge(snapshot.spacecraft[1],(id=99,))],)),
        merge(snapshot,(spacecraft=[merge(snapshot.spacecraft[1],(links=snapshot.spacecraft[1].links[1:1],))],)),
        merge(snapshot,(spacecraft=[merge(snapshot.spacecraft[1],(links=[snapshot.spacecraft[1].links[1],
            merge(snapshot.spacecraft[1].links[2],(q=(NaN,0.,0.,1.),))],))],)))
    for invalid in invalids
        @test_throws ArgumentError SE._check_link_attitude_checkpoint_state(args_for(fresh),checkpoint(invalid))
        @test [(Tuple(l.q),l.α) for l in fresh.links] == before
    end
    swapped=vehicle();swapped.links[2],swapped.links[3]=swapped.links[3],swapped.links[2]
    @test_throws ArgumentError SE._check_link_attitude_checkpoint_state(args_for(swapped),checkpoint(snapshot))
    @test_logs (:warn,r"no non-root link attitudes") SE._check_link_attitude_checkpoint_state(args_for(fresh),(;))
    @test [(Tuple(l.q),l.α) for l in fresh.links] == before
end

@testset "Hinged panel exact geometry and next rate-limited tick" begin
    function hinged()
        sc=vehicle()
        add_actuator!(sc,sc.links[2],PanelHinge(axis=SVector(0.,1.,0.),max_rate_rad_s=0.1,angle_rad=0.2))
        sc
    end
    sc=hinged();panel=SolarPanelAngleOfAttackControlModel((2,3))
    apply_panel_angle_of_attack!(panel,sc,0.,0.)
    apply_panel_angle_of_attack!(panel,sc,0.,2.)
    geometry=SE._link_attitude_checkpoint_state(args_for(sc))
    callbacks=SE._callback_checkpoint_state(args_for(sc))
    fresh=hinged();p=(shared_buffers=(callback_resume=Ref(false),),)
    entry=SE._check_link_attitude_checkpoint_state(args_for(fresh),checkpoint(from_file(geometry)))
    SE._restore_callback_checkpoint_state!(args_for(fresh),(runtime_state=(callbacks=from_file(callbacks),),),p)
    SE._restore_link_attitude_checkpoint_state!(args_for(fresh),entry)
    @test fresh.links[2].q == sc.links[2].q
    @test fresh.links[2].α == sc.links[2].α
    @test fresh.links[2].hinges[1].angle_rad == sc.links[2].hinges[1].angle_rad
    @test fresh.links[2].hinges[1].last_update_s == 2.
    apply_panel_angle_of_attack!(panel,fresh,0.,3.)
    @test fresh.links[2].hinges[1].angle_rad ≈ 0.5
end

# Its initializer deliberately writes the configured geometry. A resumed run
# must reinstate the checkpoint geometry after that reset and before the RHS.
struct PanelTicks <: AbstractControlEffectorModel end
function SpaceAGORA.initialize_control!(::PanelTicks,u,p,t)
    apply_panel_angle_of_attack!(SolarPanelAngleOfAttackControlModel((2,3)),p.args.dynamics_model.spacecraft[1],pi/2,t)
end
function SpaceAGORA.update_control!(::PanelTicks,u,p,t::Float64,i::Int64)
    apply_panel_angle_of_attack!(SolarPanelAngleOfAttackControlModel((2,3)),p.args.dynamics_model.spacecraft[i],0.2+0.01t,t)
end
SpaceAGORA.control_wrench(::PanelTicks,x::StateSample,p,t::Float64)=(SVector(0.,0.,0.),SVector(0.,0.,0.))
SpaceAGORA.spacecraft_scope(::PanelTicks)=BoundSpacecraft((1,))
@testset "Resume between control ticks after an initializer resets geometry" begin
    for hinged in (false,true)
    planet=make_no_gram_planet(:earth);control=PanelTicks();sc=vehicle()
    if hinged
        add_actuator!(sc,sc.links[2],PanelHinge(axis=SVector(0.,1.,0.),max_rate_rad_s=0.1,angle_rad=0.2))
    end
    base=make_example_config(planet=planet,spacecraft=sc,mission_time=6.,
        initial_time=InitialTime(year=2020,month=1,day=1),
        density_model=PolynomialFitAtmosphereModel(planet),ephemerides_model=SimpleEphemeridesModel(),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),AerodynamicCommandedAreaDragModel()),
        results=false,verbose=false,orientation_sim=false,solver_config=SolverConfig(solver_mode=:dp8))
    args=SimulationConfiguration(; (name=>getproperty(base,name) for name in propertynames(base) if name!=:control_model)...,
        control_model=ControlModel(control_effectors=(control,),control_rates=[2.]))
    check=C._gnc_checkpoint_round_trip(control,:control_model,args)
    println(check)
    @test check.status === :pass
    end
end
end
