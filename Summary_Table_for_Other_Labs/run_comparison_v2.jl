# julia --startup-file=no --project=2_SpaceAGORA.jl run_comparison_v2.jl [--smoke] [--video]

using Pkg

#################
# Project paths #
#################
const WORKSPACE_ROOT = @__DIR__
const SPACEAGORA_ROOT = joinpath(WORKSPACE_ROOT, "2_SpaceAGORA.jl")
const PROTOTYPE_ROOT = joinpath(WORKSPACE_ROOT, "1_Kuang's Prototype Code")

Pkg.activate(SPACEAGORA_ROOT; io = devnull)
using CSV, DataFrames, TOML, Printf

###################
# Shared settings #
###################
# Defaults for the shared runners; batch overrides are in comparison_settings_v2 below.
const VIDEO_OPTIONS = (
    show_earth    = true,
    helper_trails = false,
    target_trail  = true,
)

const TEST_CASE = (
    helpers                = 10,
    helper_altitude_km      = 1000.0,
    target_altitude_km      = 1050.0,
    helper_inclination_deg = 0.0,
    target_inclination_deg = 0.0,
    target_ecc             = 0.0,
    target_nu_deg          = 0.0,
    orbits                 = 10.0,

    laser_range_km         = 200.0,
    laser_power_w          = 10_000.0,
    magnification          = 100.0,
    mass_kg                = 227.0,

    prototype_schedule     = "gve_sma",
    spaceagora_schedule    = "gve_sma",
    spaceagora_dt_max_s     = 10.0,
)

###################
# Case validation #
###################
function validate_case(settings)
    settings.helpers isa Integer && settings.helpers >= 1 || error("helpers must be a positive integer")
    for key in (:helper_altitude_km, :target_altitude_km, :orbits, :laser_range_km, :mass_kg, :spaceagora_dt_max_s)
        value = getproperty(settings, key)
        isfinite(value) && value > 0 || error("$key must be finite and positive")
    end
    for key in (:laser_power_w, :magnification)
        value = getproperty(settings, key)
        isfinite(value) && value >= 0 || error("$key must be finite and nonnegative")
    end
    0 <= settings.target_ecc < 1 || error("target_ecc must be in [0, 1)")
    for key in (:helper_inclination_deg, :target_inclination_deg, :target_nu_deg)
        isfinite(getproperty(settings, key)) || error("$key must be finite")
    end
    gve_schedules = ("gve_sma", "gve_ecc", "gve_inc", "gve_raan", "gve_argp")
    settings.prototype_schedule in ("none", gve_schedules...) || error("Unsupported prototype schedule")
    settings.spaceagora_schedule in ("naive_next_entering", "positive_along_track", gve_schedules...) ||
        error("Unsupported SpaceAGORA schedule")
end

#######################
# Simulation duration #
#######################
function resolve_duration(settings; smoke = false)
    reference_radius_m = 6_378_137.0
    reference_mu = 3.986004418e14
    target_period_s = 2pi * sqrt((reference_radius_m + settings.target_altitude_km * 1e3)^3 / reference_mu)
    duration_s = settings.orbits * target_period_s
    isfinite(duration_s) && duration_s > 0 || error("Derived duration must be finite and positive")
    return merge(settings, (; target_period_s, duration_s = smoke ? min(60.0, duration_s) : duration_s))
end

####################
# Prototype runner #
####################
function run_prototype(settings, smoke; video = false,
                       video_duration_seconds = smoke ? 2.0 : 100.0,
                       video_options = VIDEO_OPTIONS)
    ### 1. Define orbits ###
    helper_orbits = [
        (
            a_m   = R_EARTH + settings.helper_altitude_km * 1e3,
            e     = 0.0,
            i_deg = settings.helper_inclination_deg,
            Ω_deg = 0.0,
            ω_deg = 0.0,
            ν_deg = 360.0 * (helper - 1) / settings.helpers,
        ) for helper in 1:settings.helpers
    ]
    target_orbit = (
        a_m   = R_EARTH + settings.target_altitude_km * 1e3,
        e     = settings.target_ecc,
        i_deg = settings.target_inclination_deg,
        Ω_deg = 0.0,
        ω_deg = 0.0,
        ν_deg = settings.target_nu_deg,
    )
    orbits = vcat(helper_orbits, [target_orbit])

    ### 2. Laser / cavity parameters and output paths ###
    cavity = Dict{Tuple{Int,Int},Dict{Symbol,Any}}(
        (helper, settings.helpers + 1) => Dict(:B => settings.magnification, :Pin => settings.laser_power_w)
        for helper in 1:settings.helpers
    )
    paths = scenario_paths(joinpath(PROTOTYPE_ROOT, "output"), settings, settings.duration_s;
        source   = "prototype",
        smoke,
        schedule = Symbol(settings.prototype_schedule),
    )

    ### 3. Run and verify the simulation ###
    solution, params, _, _ = mktempdir() do image_dir
        run_open_cavity_multi(orbits;
            mass_kg       = settings.mass_kg,
            Pm            = zeros(length(orbits), length(orbits)),
            cavity,
            use_los       = true,
            min_range     = 0.0,
            max_range     = settings.laser_range_km * 1e3,
            T_seconds     = settings.duration_s,
            helper_num    = settings.helpers,
            use_J2        = true,
            useDrag       = false,
            gve_schedule  = Symbol(settings.prototype_schedule),
            record_events = true,
            output_interval_s = get(settings, :output_interval_s, 10.0),
            verbose       = false,
            result_plots  = false,
            IMG_DIR       = image_dir,
        )
    end
    string(solution.retcode) == "Success" || error("Prototype failed: $(solution.retcode)")
    solution.t == output_times(settings.duration_s, get(settings, :output_interval_s, 10.0)) || error("Prototype trajectory grid mismatch")

    ### 4. Save results and extract encounters ###
    feather = save_timeseries_feather(solution, params; feather_dir = paths.feather)
    extract_feather_encounters(feather; output_dir = paths.csv)
    result = Dict{String,Any}(
        "csv_dir"         => paths.csv,
        "feather_dir"     => paths.feather,
        "duration_s"      => solution.t[end],
        "retcode"         => string(solution.retcode),
        "solver"          => "Vern9",
        "mu"              => MU,
        "earth_radius_m"  => R_EARTH,
        "light_speed_mps" => C,
    )

    ### 5. Optional animation ###
    if video
        video_path = joinpath(dirname(dirname(paths.feather)), "videos", basename(paths.feather), "animation.mp4")
        animate_all_satellites_3d_smooth_helper_target(solution, params, settings.helpers;
            output_file      = video_path,
            duration_seconds = video_duration_seconds,
            video_options...,
        )
        result["video_path"] = video_path
    end
    return result
end

#####################
# SpaceAGORA runner #
#####################
function run_spaceagora(settings, smoke)
    ### 1. Build the shared case options ###
    planet = make_no_gram_planet(:earth)
    period = 2pi * sqrt((planet.Rp_e + settings.target_altitude_km * 1e3)^3 / planet.μ)
    common = (; (key => getproperty(settings, key) for key in
        (:helpers, :helper_altitude_km, :target_altitude_km, :helper_inclination_deg,
         :target_inclination_deg, :target_ecc, :target_nu_deg, :laser_range_km,
         :laser_power_w, :magnification, :mass_kg))...)
    opts = OracleCase2Options(;
        common...,
        orbits     = settings.duration_s / period,
        schedule   = Symbol(settings.spaceagora_schedule),
        dt_max_s   = settings.spaceagora_dt_max_s,
        output_interval_s = get(settings, :output_interval_s, 10.0),
        beta       = 1.0,
        eta        = 1.0,
        smoke,
        output_dir = joinpath(SPACEAGORA_ROOT, "output"),
    )

    ### 2. Run and verify the simulation ###
    result = run_open_cavity_case_native(opts)
    string(result.summary.retcode) == "Success" || error("SpaceAGORA failed: $(result.summary.retcode)")
    isapprox(result.sol.t[end], settings.duration_s; atol=1e-8, rtol=0) || error("SpaceAGORA duration mismatch")

    ### 3. Collect saved result paths and model metadata ###
    paths = scenario_paths(opts.output_dir, opts, period * opts.orbits; source = "spaceagora", smoke)
    return Dict(
        "csv_dir"         => paths.csv,
        "feather_dir"     => result.results_dir,
        "duration_s"      => result.sol.t[end],
        "retcode"         => string(result.summary.retcode),
        "solver"          => string(result.summary.solver),
        "mu"              => planet.μ,
        "earth_radius_m"  => planet.Rp_e,
        "light_speed_mps" => 299792458.0,
    )
end

####################
# Analysis helpers #
####################
function read_analysis(run_info)
    directory = joinpath(run_info["csv_dir"], "analysis")
    extrema = CSV.read(joinpath(directory, "extrema.csv"), DataFrame)
    encounters = CSV.read(joinpath(directory, "encounters.csv"), DataFrame)
    state_count = length(CSV.File(joinpath(directory, "encounter_states.csv"); select=[:time_s]))
    return (; extrema, encounters, state_count)
end

format_value(value) = ismissing(value) ? "unavailable" : value isa AbstractFloat ? @sprintf("%.9g", value) : string(value)

function markdown_table(io, headers, rows)
    println(io, "| ", join(headers, " | "), " |")
    println(io, "| ", join(fill("---", length(headers)), " | "), " |")
    for row in rows
        println(io, "| ", join((replace(format_value(value), "|"=>"\\|") for value in row), " | "), " |")
    end
    println(io)
end

function paired_encounters(table, prefix)
    ordered = sort(table, [:sc_a, :sc_b, :start_time_s])
    counts = Dict{Tuple{Int,Int},Int}()
    ordered.occurrence = [begin
        pair = (row.sc_a, row.sc_b)
        counts[pair] = get(counts, pair, 0) + 1
    end for row in eachrow(ordered)]
    select!(ordered, :sc_a, :sc_b, :occurrence, :start_time_s, :end_time_s,
        :geometry_duration_s, :laser_on_s, :start_clipped, :end_clipped)
    rename!(ordered, [name=>prefix*name for name in names(ordered) if name ∉ ("sc_a", "sc_b", "occurrence")])
    return ordered
end

#####################
# Comparison report #
#####################
function write_report(path, settings, runs; smoke = false, saved_results = false)
    ### 1. Match metrics and encounter intervals ###
    baseline, candidate = read_analysis(runs[1]), read_analysis(runs[2])
    metrics = outerjoin(select(baseline.extrema, :metric, :extremum, :value=>:prototype),
        select(candidate.extrema, :metric, :extremum, :value=>:spaceagora);
        on=[:metric, :extremum], validate=(true, true))
    sort!(metrics, [:metric, :extremum])
    intervals = outerjoin(paired_encounters(baseline.encounters, "prototype_"),
        paired_encounters(candidate.encounters, "spaceagora_");
        on=[:sc_a, :sc_b, :occurrence], validate=(true, true))
    sort!(intervals, [:sc_a, :sc_b, :occurrence])

    ### 2. Write scenario, saved results, and comparison tables ###
    open(path, "w") do io
        println(io, "# ORACLE Prototype and SpaceAGORA Comparison\n")
        println(io, saved_results ? "Generated from existing saved outputs; no simulations were rerun. Solver return codes and model wall times are unavailable unless recorded separately.\n" :
            smoke ? "Smoke verification run; not a full-duration equivalence study.\n" : "Fresh runs of both models using the shared case below.\n")
        println(io, "## 1. Scenario Conditions\n")
        markdown_table(io, ["Setting", "Value"], [(string(key), value) for (key, value) in pairs(settings)])
        println(io, "Requested orbits count initial target Keplerian periods using the prototype reference Earth radius (6378137 m) and mu (3.986004418e14 m^3/s^2). Both models use the same derived duration in seconds, not a counted-revolution stopping condition. Smoke runs cap this duration at 60 seconds. Folder names retain the derived duration in seconds.\n")
        println(io, "J2 enabled; drag disabled; output every $(get(settings, :output_interval_s, 10.0)) seconds plus the final endpoint. Helpers are circular and evenly spaced in anomaly; RAAN and argument of periapsis are zero. Altitudes specify initial semimajor axis minus each model's Earth radius. SpaceAGORA beta and eta are 1. ORACLE here means Kuang's prototype.\n")
        markdown_table(io, ["Model", "Status", "Solver", "Simulated duration (s)", "Model run wall time (s)", "Earth radius (m)", "mu (m^3/s^2)", "Light speed (m/s)"],
            [(name, run_info["retcode"], run_info["solver"], run_info["duration_s"], get(run_info, "run_wall_time_s", missing), run_info["earth_radius_m"], run_info["mu"], run_info["light_speed_mps"])
             for (name, run_info) in zip(("ORACLE prototype", "SpaceAGORA"), runs)])
        println(io, "Model run wall time includes setup, simulation, recording, analysis, and any requested video rendering. It excludes worker startup and imports, but includes first-call compilation; it is not a warmed solver-only benchmark.\n")
        println(io, "## 2. Saved Analysis Results\n")
        for (name, run_info, analysis) in zip(("ORACLE Prototype", "SpaceAGORA"), runs, (baseline, candidate))
            println(io, "### $name\n")
            if haskey(run_info, "video_path")
                relative = replace(relpath(run_info["video_path"], dirname(abspath(path))), '\\'=>'/')
                println(io, "Video: [animation.mp4]($(replace(relative, " "=>"%20")))\n")
            end
            for (label, directory) in (("Feather", run_info["feather_dir"]), ("Analysis CSVs", joinpath(run_info["csv_dir"], "analysis")))
                relative = replace(relpath(directory, dirname(abspath(path))), '\\'=>'/')
                println(io, "$label: [$(relative)]($(replace(relative, " "=>"%20")))/\n")
            end
            println(io, "Geometric encounters: $(nrow(analysis.encounters)); saved encounter-state rows: $(analysis.state_count).\n")
            markdown_table(io, ["Metric", "Extremum", "Value", "sc_a", "sc_b", "Time (s)"],
                [(row.metric, row.extremum, row.value, row.sc_a, row.sc_b, row.time_s) for row in eachrow(analysis.extrema)])
            markdown_table(io, ["Encounter", "sc_a", "sc_b", "Start (s)", "End (s)", "Geometry (s)", "Laser on (s)", "Start clipped", "End clipped"],
                [(row.encounter_id, row.sc_a, row.sc_b, row.start_time_s, row.end_time_s, row.geometry_duration_s,
                  row.laser_on_s, row.start_clipped, row.end_clipped) for row in eachrow(analysis.encounters)])
        end
        println(io, "## 3. Comparison\n")
        println(io, "Differences are SpaceAGORA minus prototype. Percent differences use the absolute prototype value; zero or unavailable references have no percentage.\n")
        markdown_table(io, ["Metric", "Extremum", "Prototype", "SpaceAGORA", "Difference", "Difference (%)"],
            [(row.metric, row.extremum, row.prototype, row.spaceagora, row.spaceagora-row.prototype,
              ismissing(row.prototype) || iszero(row.prototype) ? missing : 100*(row.spaceagora-row.prototype)/abs(row.prototype)) for row in eachrow(metrics)])
        println(io, "Encounters are matched by unordered spacecraft pair and chronological occurrence, not by run-local IDs. Unmatched occurrences show unavailable differences; a missed encounter can shift subsequent occurrence matching.\n")
        markdown_table(io, ["sc_a", "sc_b", "Occurrence", "Entry difference (s)", "Exit difference (s)", "Geometry difference (s)", "Laser-on difference (s)"],
            [(row.sc_a, row.sc_b, row.occurrence, row.spaceagora_start_time_s-row.prototype_start_time_s,
              row.spaceagora_end_time_s-row.prototype_end_time_s,
              row.spaceagora_geometry_duration_s-row.prototype_geometry_duration_s,
              row.spaceagora_laser_on_s-row.prototype_laser_on_s) for row in eachrow(intervals)])
        println(io, "No numerical acceptance tolerance was specified: these tables quantify agreement, not a physical-equivalence verdict. Motion extrema use only encounter states and entry/exit boundaries; unrelated pairs and out-of-encounter times are excluded. Clipped encounters cover only the recorded portion. Earth/light-speed constants, solvers, LOS checks, and scheduling differ. SpaceAGORA uses accepted-step scheduling and a velocity-kick callback; prototype forces are evaluated in the ODE. The report reads saved analysis CSVs after both worker processes exit; it does not compare full trajectories or equate the two delta-v diagnostics.\n")
    end
    return path
end

###############
# Batch cases #
###############
# t denotes target altitude (km); i denotes target inclination (degrees).
# Case 3 duplicates case 2 and is run only once.
const COMPARISON_CASES_V2 = (
    (id="case1_N20_t1050_i0_gve_sma", helpers=20, target_altitude_km=1050.0, inclination=0.0, target_nu_deg=0.0, schedule="gve_sma"),
    (id="case2_N20_t1000_i5_gve_sma", helpers=20, target_altitude_km=1000.0, inclination=5.0, target_nu_deg=0.0, schedule="gve_sma"),
    (id="case4_N100_t1050_i0_gve_sma", helpers=100, target_altitude_km=1050.0, inclination=0.0, target_nu_deg=0.0, schedule="gve_sma"),
    (id="case5_N20_t1000_i5_gve_inc", helpers=20, target_altitude_km=1000.0, inclination=5.0, target_nu_deg=0.0, schedule="gve_inc"),
    (id="case6_N20_t1150_i0_gve_sma", helpers=20, target_altitude_km=1150.0, inclination=0.0, target_nu_deg=0.0, schedule="gve_sma"),
    (id="case7_N20_t1000_i0_gve_sma", helpers=20, target_altitude_km=1000.0, inclination=0.0, target_nu_deg=1.0, schedule="gve_sma"),
)

##################
# Batch settings #
##################
function comparison_settings_v2(case; smoke=false)
    ### 1. Apply case-specific and batch-wide overrides ###
    settings = merge(TEST_CASE, (
        helpers=case.helpers,
        helper_altitude_km=1000.0,
        helper_inclination_deg=0.0,
        target_altitude_km=case.target_altitude_km,
        target_inclination_deg=case.inclination,
        target_nu_deg=case.target_nu_deg,
        orbits=50.0,
        output_interval_s=10.0,
        prototype_schedule=case.schedule,
        spaceagora_schedule=case.schedule,
    ))

    ### 2. Validate and derive the common duration ###
    validate_case(settings)
    return resolve_duration(settings; smoke)
end

##################
# Worker imports #
##################
# Each child process loads only its model; animation dependencies are optional.
if length(ARGS) == 4 && ARGS[1] == "--worker-v2"
    if ARGS[2] == "prototype"
        include("1_Kuang's Prototype Code/test16_feather.jl")
        if get(TOML.parsefile(ARGS[3]), "video", false)
            include("3_Kuang's Prototype Code_for_animation/comparison_animation.jl")
        end
    elseif ARGS[2] == "spaceagora"
        include("2_SpaceAGORA.jl/ORACLE/run_case2_laser_links.jl")
    else
        error("Unknown worker: $(ARGS[2])")
    end
end

##############
# Batch main #
##############
function main_v2(arguments=ARGS)
    ### 1. Execute and time a single worker request ###
    if length(arguments) == 4 && arguments[1] == "--worker-v2"
        config = TOML.parsefile(arguments[3])
        settings = (; (Symbol(key)=>value for (key, value) in config["case"])...)
        elapsed = @elapsed begin
            result = if arguments[2] == "prototype"
                get(config, "video", false) ? run_animated_comparison(settings, config["smoke"]) :
                    run_prototype(settings, config["smoke"])
            else
                run_spaceagora(settings, config["smoke"])
            end
        end
        result["run_wall_time_s"] = elapsed
        open(arguments[4], "w") do io
            TOML.print(io, result)
        end
        return nothing
    end

    ### 2. Parse batch options ###
    all(argument -> argument in ("--smoke", "--video"), arguments) ||
        error("Usage: julia run_comparison_v2.jl [--smoke] [--video]")
    smoke = "--smoke" in arguments
    video = "--video" in arguments
    suffix = smoke ? "_smoke" : ""

    ### 3. Initialize the batch index ###
    index_path = joinpath(WORKSPACE_ROOT, "comparison_summary_v2$(suffix).md")
    open(index_path, "w") do io
        println(io, "# $(length(COMPARISON_CASES_V2))-Case Comparison\n")
        println(io, "Requested case 3 duplicates case 2 and is run only once. Each case has one target and helpers at 1000 km with zero inclination. Both models use the same scheduler and 50 initial target orbital periods.\n")
        println(io, "Trajectory recording interval: 10 seconds plus the final endpoint. Animation: $(video ? "enabled" : "disabled"). Simulation settings otherwise inherit TEST_CASE from run_comparison_v2.jl.\n")
        smoke && println(io, "Smoke validation: 60 simulated seconds per model.\n")
        println(io, "Wall times include model setup, simulation, recording, analysis, and any requested rendering; worker startup/imports are excluded, first-call compilation is included.\n")
    end

    ### 4. Run each case sequentially in separate model processes ###
    for case in COMPARISON_CASES_V2
        settings = comparison_settings_v2(case; smoke)
        println("\n$(case.id): $(settings.helpers) helpers, one target, $(settings.duration_s) seconds")

        ### 4.1. Share configuration and collect worker results ###
        runs = mktempdir() do directory
            config_path = joinpath(directory, "case.toml")
            open(config_path, "w") do io
                TOML.print(io, Dict("case"=>Dict(string(key)=>value for (key, value) in pairs(settings)), "smoke"=>smoke, "video"=>video))
            end
            results = Dict{String,Any}[]
            for model in ("prototype", "spaceagora")
                result_path = joinpath(directory, model*".toml")
                println("Running $model...")
                run(`$(Base.julia_cmd()) --startup-file=no --project=$SPACEAGORA_ROOT $(@__FILE__) --worker-v2 $model $config_path $result_path`)
                push!(results, TOML.parsefile(result_path))
            end
            results
        end
        ### 4.2. Write the case report and append its index entry ###
        report_path = joinpath(WORKSPACE_ROOT, "comparison_v2_$(case.id)$(suffix).md")
        write_report(report_path, settings, runs; smoke)
        open(index_path, "a") do io
            println(io, "## $(case.id)\n")
            println(io, "Target: $(settings.target_altitude_km) km, $(case.inclination) deg; helpers: $(case.helpers); scheduler: $(case.schedule).\n")
            println(io, "[Comparison report]($(basename(report_path)))\n")
            markdown_table(io, ["Model", "Simulated duration (s)", "Model run wall time (s)"],
                [(name, result["duration_s"], result["run_wall_time_s"])
                 for (name, result) in zip(("ORACLE prototype", "SpaceAGORA"), runs)])
            if haskey(runs[1], "video_path")
                video_path = replace(relpath(runs[1]["video_path"], WORKSPACE_ROOT), '\\'=>'/')
                println(io, "[Prototype animation]($(replace(video_path, " "=>"%20")))\n")
            end
        end
        println("Report saved to $report_path")
    end
    println("Batch index saved to $index_path")
    return index_path
end

####################
# Script entrypoint #
####################
if abspath(PROGRAM_FILE) == @__FILE__
    main_v2()
end