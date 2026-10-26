# Run from an environment containing SpaceAGORA:
# julia --project=. examples/AGORA_Mission_Workflow.jl output/coast-001
# Choose an absent or empty destination for each run.
using SpaceAGORA

function build_mission_coast(; duration_s=300.0)
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=45.0, ω=0.0, Ω=0.0, ν=0.0),
        prop_mass=0.0, id=1)
    return make_example_config(; planet, spacecraft, mission_time=Float64(duration_s),
        initial_time=InitialTime(year=2024, month=1, day=1),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(), ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, verbose=false,
        solver_config=SolverConfig(solver_mode=:tsit5, parallel=false))
end

plan = MissionPlan("Earth coast", build_mission_coast;
    phases=(MissionPhase(:coast; window_s=(0.0, 300.0)),))
destination = isempty(ARGS) ? joinpath("output", "mission-coast-001") : only(ARGS)
mission_run = run_mission(plan; output_dir=destination)
println("Saved run to: ", mission_run.directory)
println("Viewer payload: ", mission_run.payload_path)
println("Inspect the table through mission_run.results.table, or open the directory with SpaceAGORAVisualizer.")
