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
    mc_smoke = SM.SimConfig._with_fields(mc;
        number_of_orbits=max(1, min(mc.number_of_orbits, 1)),
        mission_time=_example_smoke_mission_time(mc.mission_time),
        num_steps_to_save=max(50, min(mc.num_steps_to_save, 200)),
    )
    output_override = strip(get(ENV, "SPACEAGORA_CLI_OUTPUT_DIR", ""))
    # A smoke run writes at most its results table: no plots, scene or
    # checkpoints. Every other setting, physics included, is the example's own.
    ss_smoke = SM.SimConfig._with_fields(ss;
        results=keep_results,
        verbose=false,
        # Honor an explicit CLI destination; otherwise retain per-run smoke isolation.
        results_directory=isempty(output_override) ? joinpath(pwd(), "output") : abspath(output_override),
        generate_plots=false,
        normalize=false,
        save_csv=keep_results,
        save_visualization_scene=false,
        checkpoint_enabled=false,
        resume_from_checkpoint=false,
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
    # The long layout writes its spacecraft table as simulation_results_spacecraft.csv.
    long_layout = args_eff.simulation_settings.results_layout === :long
    csv_path = joinpath(args_eff.simulation_settings.results_directory,
        long_layout ? "simulation_results_spacecraft.csv" : "simulation_results.csv")
    saved_csv_path = nothing
    if args_eff.simulation_settings.results && isfile(csv_path)
        df = CSV.read(csv_path, DataFrame)
        println("Saved $(nrow(df)) $(long_layout ? "spacecraft rows" : "samples") to $(abspath(csv_path))")
        saved_csv_path = csv_path
    end
    println("COMPUTATIONAL TIME = $(t) s")
    return saved_csv_path
end
