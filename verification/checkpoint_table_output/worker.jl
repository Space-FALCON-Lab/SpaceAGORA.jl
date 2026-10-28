# Standalone ordinary-route reproduction; no KS package or source is loaded.
using SpaceAGORA,Test,TOML
import SpaceAGORA.TelemetryVerification: make_example_config
const S=SpaceAGORA.SimulationModel
const SC=S.SimulationCallbacks
function config(directory;checkpoint=false,resume=false,duration=1.04,rate=.1,interval=.37)
    root=S.Link(root=true,m=100.,dims=S.MVector(1.,1.,1.))
    ic=S.CartesianInitialCondition(S.SVector(4e9,0.,0.),S.SVector(2.,-3.,1.))
    sc=S.SpacecraftModel(links=[root],root=root,prop_mass=1.,inertia_tensor=root.inertia,initial_condition=ic,id=1)
    cfg=make_example_config(planet=S.Earth(),spacecraft=sc,mission_time=duration,
        initial_time=S.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.),
        dynamic_effectors=(),density_model=S.NoAtmosphereModel(),ephemerides_model=S.SimpleEphemeridesModel(),
        orientation_sim=false,keplerian=false,EI_km=1.,verbose=false,results=false)
    return S.SimConfig._with_configuration(cfg;
        simulation_settings=S.SimulationSettings(results=false,verbose=false,generate_plots=false,normalize=false,
            save_csv=false,results_directory=directory,checkpoint_enabled=checkpoint,
            checkpoint_interval_s=interval,checkpoint_directory=joinpath(directory,"checkpoints"),resume_from_checkpoint=resume),
        mission_configuration=S.MissionConfiguration(mission_type=S.MissionTime,keplerian=false,
            mission_time=duration,orientation_sim=false,data_rate=rate))
end
length(ARGS)==5 || error("usage: worker.jl PHASE DIRECTORY DURATION RATE CHECKPOINT_INTERVAL")
phase=Symbol(ARGS[1]);directory=abspath(ARGS[2]);duration=parse(Float64,ARGS[3])
rate=parse(Float64,ARGS[4]);interval=parse(Float64,ARGS[5]);mkpath(directory)
cfg=config(directory;checkpoint=phase!==:direct,resume=phase===:resume,duration,rate,interval)
function interrupt_factory(args)
    SC.DiscreteCallback((u,t,int)->t==.5,int->error("deliberate fresh-process interruption");
        initialize=(cb,u,t,int)->begin
            t<.5<=int.p.shared_buffers.solve_segment_end_time[] && SC.add_tstop!(int,.5)
            nothing
        end,save_positions=(false,false))
end
try
    result=run_simulation(cfg;return_results=true,return_solution=false,
        priority_callback_factories=phase===:interrupt ? (interrupt_factory,) : ())
    phase===:interrupt && error("expected interruption was missing")
    times=result.table.time
    errors=[sqrt((result.table.sc1_pos_1[i]-(4e9+2t))^2+
                 (result.table.sc1_pos_2[i]+3t)^2+(result.table.sc1_pos_3[i]-t)^2)
            for (i,t) in enumerate(times)]
    record=Dict("phase"=>string(phase),"time"=>times,"rows"=>length(times),
        "max_affine_position_error_m"=>maximum(errors),"columns"=>string.(propertynames(result.table)),
        "duration_s"=>duration,"rate_s"=>rate,"interval_s"=>interval,"julia"=>string(VERSION))
    open(joinpath(directory,string(phase)*".toml"),"w") do io;TOML.print(io,record);end
    println("CHECKPOINT_WORKER_SUCCESS ",phase," rows=",length(times))
catch err
    if phase===:interrupt && occursin("deliberate fresh-process interruption",sprint(showerror,err))
        println("CHECKPOINT_WORKER_EXPECTED_INTERRUPTION")
    else
        rethrow()
    end
end
