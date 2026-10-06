const SM = SimulationModel

@inline _example_smoke_enabled() = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
@inline _example_smoke_results_enabled() = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE_RESULTS", "0") == "1"

@inline function _example_smoke_mission_time(default_time::Float64)
    raw = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME", "120.0")
    parsed = tryparse(Float64, raw)
    if parsed === nothing || !(parsed > 0.0)
        return min(default_time, 120.0)
    end
    return min(default_time, parsed)
end

function _example_smoke_args(args::SM.SimulationConfiguration)
    if !_example_smoke_enabled()
        return args
    end

    mc = args.mission_configuration
    ss = args.simulation_settings
    keep_results = _example_smoke_results_enabled()
    mc_smoke = SM.MissionConfiguration(
        mission_type=mc.mission_type,
        keplerian=mc.keplerian,
        number_of_orbits=max(1, min(mc.number_of_orbits, 1)),
        mission_time=_example_smoke_mission_time(mc.mission_time),
        orientation_sim=mc.orientation_sim,
        num_steps_to_save=max(50, min(mc.num_steps_to_save, 200)),
        data_rate=mc.data_rate
    )
    ss_smoke = SM.SimulationSettings(
        results=keep_results,
        verbose=false,
        # Keep smoke outputs local to the current run directory to avoid cross-run collisions.
        results_directory=joinpath(pwd(), "output"),
        generate_plots=false,
        generate_filenames=ss.generate_filenames,
        normalize=false,
        save_csv=keep_results
    )

    return SM.SimConfig._with_configuration(args;
        simulation_settings=ss_smoke,
        mission_configuration=mc_smoke,
    )
end

using ..SimulationModel.ExampleConfiguration: make_three_body_spacecraft, make_example_config

function run_and_report(args::SM.SimulationConfiguration)
    args_eff = _example_smoke_args(args)
    t = @elapsed run_simulation(args_eff)
    csv_path = joinpath(args_eff.simulation_settings.results_directory, "simulation_results.csv")
    saved_csv_path = nothing
    if args_eff.simulation_settings.results && isfile(csv_path)
        df = CSV.read(csv_path, DataFrame)
        println("Saved $(nrow(df)) samples to $(abspath(csv_path))")
        saved_csv_path = csv_path
    end
    println("COMPUTATIONAL TIME = $(t) s")
    return saved_csv_path
end
