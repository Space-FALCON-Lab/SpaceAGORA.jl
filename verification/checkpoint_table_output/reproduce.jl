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
length(ARGS)==1 || error("usage: checkpoint-table-probe.jl NEW_OUTPUT_DIRECTORY")
out=abspath(only(ARGS));isdir(out) && !isempty(readdir(out)) && error("Refusing prior evidence overwrite")
mkpath(out)
direct=run_simulation(config(joinpath(out,"direct"));return_results=true,return_solution=false)
segmented=run_simulation(config(joinpath(out,"segmented");checkpoint=true);return_results=true,return_solution=false)
function interrupt_factory(args)
    return SC.DiscreteCallback((u,t,int)->t==.5,int->error("deliberate checkpoint reproduction interruption");
        initialize=(cb,u,t,int)->begin
            t<.5<=int.p.shared_buffers.solve_segment_end_time[] && SC.add_tstop!(int,.5)
            nothing
        end,save_positions=(false,false))
end
cold_dir=joinpath(out,"cold-resume")
try
    run_simulation(config(cold_dir;checkpoint=true);return_results=true,return_solution=false,
        priority_callback_factories=(interrupt_factory,))
    error("expected interruption was missing")
catch err
    occursin("deliberate checkpoint reproduction interruption",sprint(showerror,err)) || rethrow()
end
# A separate fresh-process replay is added by the full qualification suite;
# this minimal ordinary reproduction intentionally tests the returned full prefix.
resumed=run_simulation(config(cold_dir;checkpoint=true,resume=true);return_results=true,return_solution=false)
grid=vcat(collect(range(0.,1.04;step=.1)),1.04)
records=Dict{String,Any}()
for (name,result) in (("direct",direct),("segmented",segmented),("resumed",resumed))
    times=result.table.time
    records[name]=Dict("times"=>times,"missing_required"=>setdiff(grid,times),
        "rows"=>length(times),"all_unique"=>length(unique(times))==length(times))
end
open(joinpath(out,"observed.toml"),"w") do io;TOML.print(io,records);end
@testset "Ordinary checkpoint table retains mission grid" begin
    @test Base.get_extension(SpaceAGORA,:SpaceAGORAKSPropagatorExt)===nothing
    for name in ("direct","segmented","resumed")
        @test isempty(records[name]["missing_required"])
        @test records[name]["all_unique"]
        @test first(records[name]["times"])==0.
        @test last(records[name]["times"])==1.04
    end
end
println("ORDINARY_TABLE_PROBE_SUCCESS")
