module ControlForecastFixture
using SpaceAGORA, StaticArrays, Serialization, TOML, LinearAlgebra, Test, Distributed, SHA
import SpaceAGORA.TelemetryVerification: make_example_config
const SM=SpaceAGORA.SimulationModel
const L=SpaceAGORA.SimulationLifecycle
const SE=SpaceAGORA.SimulationEngine

mutable struct ProbeControl <: AbstractControlEffectorModel
    directory::String
    captured::Bool
    held::Float64
    ticks::Vector{Float64}
end
SpaceAGORA.spacecraft_scope(::ProbeControl)=BoundSpacecraft((1,))
SpaceAGORA.held_command_time(c::ProbeControl,i)=isempty(c.ticks) ? -Inf : last(c.ticks)
SpaceAGORA.control_wrench(c::ProbeControl,x,p,t)=(SVector(0.,c.held,0.),SVector(0.,0.,0.))
L.checkpoint_state(c::ProbeControl)=(captured=c.captured,held=c.held,ticks=copy(c.ticks))
L.restore_checkpoint_state!(c::ProbeControl,s)=(c.captured=s.captured;c.held=s.held;c.ticks=copy(s.ticks);nothing)
L.initialize_control!(c::ProbeControl,u,p,t)=(c.captured=false;c.held=0.;empty!(c.ticks);nothing)
function SpaceAGORA.update_control!(c::ProbeControl,u,p,t::Float64,i::Int64)
    if !c.captured && t>=0.6
        c.captured=true
        snapshot=capture_control_forecast(u,p,t,c)
        serialize(joinpath(c.directory,"snapshot.bin"),snapshot)
    end
    push!(c.ticks,t)
    # Non-idempotent command detects both a skipped and a repeated tick.
    c.held+=10t
    nothing
end

function parent(directory)
    ispath(directory) && error("Refusing overwrite")
    mkpath(directory)
    root=SM.Link(root=true,m=100.,dims=MVector(1.,1.,1.))
    ic=CartesianInitialCondition(SVector(7e6,0.,0.),SVector(0.,7500.,0.))
    sc=SpacecraftModel(links=[root],root=root,prop_mass=0.,inertia_tensor=root.inertia,initial_condition=ic,id=1)
    base=make_example_config(planet=SM.Earth(),spacecraft=sc,mission_time=3.,
        initial_time=InitialTime(year=2020,month=1,day=1),dynamic_effectors=(),
        ephemerides_model=SimpleEphemeridesModel(),orientation_sim=false,keplerian=false,EI_km=1.,verbose=false,results=false)
    args=SM.SimConfig._with_configuration(base;
        control_model=ControlModel(control_effectors=(ProbeControl(directory,false,0.,Float64[]),),control_rates=[0.2]),
        solver_config=SolverConfig(solver_mode=:split_imex),
        simulation_settings=SimulationSettings(results=false,verbose=false,generate_plots=false,results_directory=directory))
    out=run_simulation(args;return_results=true,return_solution=true)
    snapshot=deserialize(joinpath(directory,"snapshot.bin"))
    candidate=only(snapshot.configuration.control_model.control_effectors)
    @test_throws ArgumentError run_control_forecast(snapshot,candidate;duration_s=2.4,sample_period_s=0.1)
    @test !only(args.control_model.control_effectors).captured
    @test isempty(only(args.control_model.control_effectors).ticks)
    @test candidate.ticks==[0.2,0.4]
    @test snapshot.clock_origin_s==0.
    @test snapshot.checkpoint.t==3*0.2
    serialize(joinpath(directory,"parent.bin"),(state=collect(out.solution.u[end]),final_time=out.solution.t[end]))
    # Check owned-worker cleanup even when the user computation fails.
    before=Set(Distributed.workers())
    @test_throws ErrorException with_control_forecast_worker(snapshot) do evaluate
        error("deliberate forecast client failure")
    end
    @test Set(Distributed.workers())==before
    @test_throws ArgumentError SE._check_forecast_data((state=(channel=RemoteChannel(()->Channel{Int}(1)),),))
    @test_throws ArgumentError SE._check_forecast_data((state=Future(),))
    @test_throws ArgumentError SE._check_forecast_data(Channel{Int}(1))
    @test_throws ArgumentError SE._check_forecast_data((callback=identity,))
    @test_throws RemoteException with_control_forecast_worker(snapshot;
            packages=(Base.PkgId(Base.UUID("ffffffff-ffff-ffff-ffff-ffffffffffff"),"MissingForecastTestPackage"),)) do evaluate
        error("must not enter client")
    end
    @test Set(Distributed.workers())==before
    @test_throws InterruptException with_control_forecast_worker(snapshot) do evaluate
        throw(InterruptException())
    end
    @test Set(Distributed.workers())==before
    println("parent and failure cleanup passed")
end

function child(directory)
    snapshot=deserialize(joinpath(directory,"snapshot.bin"))
    candidate=only(snapshot.configuration.control_model.control_effectors)
    parent=deserialize(joinpath(directory,"parent.bin"))
    function digest(x)
        io=IOBuffer();serialize(io,x);bytes2hex(sha256(take!(io)))
    end
    original=digest(snapshot)
    result=run_control_forecast(snapshot,candidate;duration_s=3-snapshot.checkpoint.t,sample_period_s=0.1)
    repeat=run_control_forecast(snapshot,candidate;duration_s=3-snapshot.checkpoint.t,sample_period_s=0.1)
    @test digest(snapshot)==original
    @test result.retcode==repeat.retcode=="Success"
    @test result.final_time_s==repeat.final_time_s==parent.final_time==3.
    @test result.start_time_s==snapshot.checkpoint.t
    @test result.clock_origin_s==snapshot.clock_origin_s
    @test norm(result.final_state[1:3]-parent.state[1:3])<1e-5
    @test norm(result.final_state[4:6]-parent.state[4:6])<1e-8
    @test result.final_state==repeat.final_state
    @test isequal(result.table,repeat.table)
    @test_throws ArgumentError run_control_forecast(snapshot,candidate;duration_s=-1.,sample_period_s=0.1)
    # Offset the clock by 0.05 seconds and seed its two prior commands. The
    # Cartesian starting state is arbitrary; integrate the remaining piecewise
    # constant acceleration analytically instead of borrowing another solver.
    offset_candidate=deepcopy(candidate)
    offset_candidate.ticks=[0.25,0.45]
    offset_candidate.held=10sum(offset_candidate.ticks)
    checkpoint=merge(snapshot.checkpoint,(t=0.65,
        runtime_state=merge(snapshot.checkpoint.runtime_state,
            (periodic_ticks=(version=2,time_s=0.65,pending=[(clock=1,period_s=0.2)]),)),))
    offset=ControlForecastSnapshot(snapshot.configuration,checkpoint,0.05,snapshot.owner_pid,snapshot.kernels,snapshot.environment)
    result_offset=run_control_forecast(offset,offset_candidate;duration_s=2.4,sample_period_s=0.1)
    position=copy(checkpoint.u[1:3]);velocity=copy(checkpoint.u[4:6]);held=offset_candidate.held
    times=vcat([0.05+i*0.2 for i in 3:14],[3.05])
    for i in 1:length(times)-1
        held+=10times[i]
        dt=times[i+1]-times[i]
        acceleration=[0.,held/100.,0.]
        position .+= velocity.*dt .+ acceleration.*(0.5dt^2)
        velocity .+= acceleration.*dt
    end
    @test result_offset.clock_origin_s==0.05
    @test result_offset.final_time_s≈3.05
    @test norm(result_offset.final_state[1:3]-position)<1e-5
    @test norm(result_offset.final_state[4:6]-velocity)<1e-8
    open(joinpath(directory,"comparison.toml"),"w") do io
        TOML.print(io,Dict("position_difference_m"=>norm(result.final_state[1:3]-parent.state[1:3]),
            "velocity_difference_m_s"=>norm(result.final_state[4:6]-parent.state[4:6]),
            "repeat_state_exact"=>result.final_state==repeat.final_state,
            "snapshot_unchanged"=>digest(snapshot)==original,
            "start_time_s"=>result.start_time_s,"final_time_s"=>result.final_time_s,
            "nonzero_origin_position_difference_m"=>norm(result_offset.final_state[1:3]-position),
            "nonzero_origin_velocity_difference_m_s"=>norm(result_offset.final_state[4:6]-velocity)))
    end
    println("forecast repeatability and uninterrupted-state checks passed")
end

function kernel_replay()
    planets=SM.Planets
    spice=SE.SPICE
    mktempdir() do directory
        paths=[joinpath(directory,"first.tpc"),joinpath(directory,"second.tpc")]
        for (index,path) in enumerate(paths)
            write(path,"KPL/PCK\n\\begindata\nFORECAST_TEST_VALUE = ( $index )\n\\begintext\n")
        end
        kernels=[(path=path,kind="TEXT",sha256=bytes2hex(sha256(read(path)))) for path in paths[[1,2,1]]]
        planets._furnsh_once(paths[2])
        generation=planets._SPICE_KERNEL_POOL_GENERATION[]
        key=SE._srp_ephemeris_reuse_key("earth",0.,1.,1.)
        SE._SRP_EPHEMERIS_REUSE_CACHE[key]=SM.SRPSunEphemerisCache([0.,1.],[SVector(1.,0.,0.),SVector(1.,0.,0.)])
        planets._EARTH_CACHE[("test","test")]=SM.Earth()
        SE._initialize_control_forecast_kernels!(kernels)
        @test planets._SPICE_KERNEL_POOL_GENERATION[]>generation
        @test isempty(planets._EARTH_CACHE)
        @test isempty(SE._SRP_EPHEMERIS_REUSE_CACHE)
        @test spice.ktotal("ALL")==3
        @test [spice.kdata(i,"ALL",8192,8192)[1] for i in 1:3]==paths[[1,2,1]]
        @test spice.gdpool("FORECAST_TEST_VALUE";start=1,room=1)==[1.]
        write(paths[2],"changed contents")
        @test_throws ArgumentError SE._initialize_control_forecast_kernels!(kernels)
        @test spice.ktotal("ALL")==3
        @test spice.gdpool("FORECAST_TEST_VALUE";start=1,room=1)==[1.]
        rm(paths[2])
        @test_throws ArgumentError SE._initialize_control_forecast_kernels!(kernels)
        @test spice.ktotal("ALL")==3
    end
    nothing
end

end
