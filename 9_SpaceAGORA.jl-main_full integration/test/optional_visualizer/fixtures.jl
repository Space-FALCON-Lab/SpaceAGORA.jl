using SpaceAGORA, StaticArrays, DataFrames, Arrow, JSON, Base64
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft
const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization

function visualizer_args(output)
    craft = make_three_body_spacecraft(bus_dims=(1.0,1.0,1.0), panel_dims=(0.01,1.0,0.5),
        bus_mass=100.0, panel_mass_each=1.0, panel_offset_y=1.0,
        ic=SM.InitialCondition(ra=4.5e6,rp=3.8e6,i=0.0,ω=0.0,Ω=0.0,ν=0.0),
        prop_mass=0.0, id=1)
    return make_example_config(planet=make_no_gram_planet(:mars), spacecraft=craft,
        mission_time=0.2, initial_time=SM.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),
        dynamic_effectors=(), density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false,
        keplerian=true, verbose=false, results=true, results_directory=output)
end

function saved_visualizer_fixture(output)
    mkpath(output)
    args = visualizer_args(output)
    scene = build_visualization_scene(args; rotation_max_samples=2)
    prefix = joinpath(output, "simulation_results")
    write_visualization_scene(prefix * "_scene.json", scene)
    table = DataFrame(time=[0.0,1.0,2.0],sc1_pos_1=[4e6,4e6+10,4e6+20],
        sc1_pos_2=[0.0,2.0,4.0],sc1_pos_3=[0.0,3.0,6.0],
        sc1_q_1=zeros(3),sc1_q_2=zeros(3),sc1_q_3=zeros(3),sc1_q_4=ones(3),
        sc1_metric=[3.0,2.0,1.0])
    Arrow.write(prefix * ".feather", table)
    return (; args, scene, prefix, table)
end

function embedded_visualizer_payload(path)
    html=read(path,String)
    captured=match(r"window\.SPACEAGORA_VIEWER = (.*?);\n</script>"s,html)
    captured === nothing && error("HTML lacks its viewer payload")
    return JSON.parse(captured.captures[1])
end
