include("run_comparison_v2.jl")
include("1_Kuang's Prototype Code/extract_feather_encounters.jl")
using Test

const SUMMARY_DIRECTORY = joinpath(WORKSPACE_ROOT, "4_Comparasion_Summary")

function saved_bundle(case, model, interval)
    root = model == "prototype" ? PROTOTYPE_ROOT : SPACEAGORA_ROOT
    pattern = "$(model)_N$(case.helpers)_h1000km_t$(Int(case.target_altitude_km))km_ih0_it$(Int(case.inclination))_"
    names = filter(name -> startswith(name, pattern) && endswith(name, "_$(case.schedule)_J2true_dt$(interval)s"),
        readdir(joinpath(root, "output", "feather")))
    name = only(names)
    return (feather=joinpath(root, "output", "feather", name), csv=joinpath(root, "output", "CSV", name))
end

function check_saved_grid(table, interval, duration)
    expected = collect(0.0:Float64(interval):duration)
    last(expected) < duration && push!(expected, duration)
    @test length(table.time) == length(expected)
    @test table.time[1:end-1] == expected[1:end-1]
    @test isapprox(last(table.time), duration; atol=1e-8, rtol=0)
    @test parse(Float64, Arrow.getmetadata(table)["output_interval_s"]) == interval
    @test Arrow.getmetadata(table)["target_id"] == "21"
    for spacecraft in 1:21, field in ("pos", "vel"), component in 1:3
        @test all(value -> !ismissing(value) && isfinite(value),
            getproperty(table, Symbol("sc$(spacecraft)_$(field)_$(component)")))
    end
end

function compare_saved_model(io, case, model, duration)
    fine_paths, coarse_paths = saved_bundle(case, model, 1), saved_bundle(case, model, 10)
    fine = Arrow.Table(joinpath(fine_paths.feather, "trajectory.feather"))
    coarse = Arrow.Table(joinpath(coarse_paths.feather, "trajectory.feather"))
    @testset "$(case.id) $model coverage" begin
        check_saved_grid(fine, 1, duration)
        check_saved_grid(coarse, 10, duration)
    end
    indices = [searchsortedfirst(fine.time, time) for time in coarse.time[1:end-1]]
    push!(indices, length(fine.time))
    @test fine.time[indices[1:end-1]] == coarse.time[1:end-1]
    @test isapprox(last(fine.time), last(coarse.time); atol=1e-8, rtol=0)
    println(io, "### $model\n")
    for (label, paths) in (("dt1s", fine_paths), ("dt10s", coarse_paths))
        relative = replace(relpath(paths.feather, SUMMARY_DIRECTORY), " "=>"%20")
        println(io, "[$label source Feather bundle]($relative/)\n")
    end
    markdown_table(io, ["Cadence", "Rows", "Start (s)", "End (s)", "Grid and finite-state checks"],
        [(label, length(table.time), first(table.time), last(table.time), "PASS")
         for (label, table) in (("dt1s", fine), ("dt10s", coarse))])
    differences = []
    for (label, pattern) in (("Position components (m)", r"^sc\d+_pos_\d$"),
                            ("Velocity components (m/s)", r"^sc\d+_vel_\d$"),
                            ("Accumulated delta-v components (m/s)", r"^dv_.*_accumulated$"))
        columns = filter(name -> occursin(pattern, String(name)), propertynames(fine))
        @test !isempty(columns)
        shared = maximum(maximum(abs.(getproperty(coarse, name) .- getproperty(fine, name)[indices])) for name in columns)
        endpoint = maximum(abs(last(getproperty(coarse, name))-last(getproperty(fine, name))) for name in columns)
        push!(differences, (label, shared, endpoint))
    end
    markdown_table(io, ["Quantity", "Maximum absolute shared-time difference", "Maximum absolute final difference"], differences)
    mismatches = count(!isequal(first_value, second_value) for (first_value, second_value) in
        zip(fine.laser_active_helper[indices], coarse.laser_active_helper))
    println(io, "Shared timestamps: $(length(indices)); active-helper mismatches: $mismatches.\n")
    event_rows = []
    for file in ("geometry_encounters.feather", "laser_on.feather")
        fine_events = DataFrame(Arrow.Table(joinpath(fine_paths.feather, file)))
        coarse_events = DataFrame(Arrow.Table(joinpath(coarse_paths.feather, file)))
        push!(event_rows, (file, nrow(fine_events), nrow(coarse_events), isequal(fine_events, coarse_events)))
    end
    markdown_table(io, ["Event table", "dt1s rows", "dt10s rows", "Exactly equal (all fields)"], event_rows)
    mktempdir() do directory
        fine_analysis = extract_feather_encounters(fine_paths.feather; output_dir=directory)
        coarse_analysis = read_analysis(Dict("csv_dir"=>coarse_paths.csv))
        println(io, "Encounter analysis tables exactly equal: $(isequal(fine_analysis.encounters, coarse_analysis.encounters)). " *
            "Encounter-state rows: $(nrow(fine_analysis.samples)) (dt1s), $(coarse_analysis.state_count) (dt10s).\n")
        metrics = innerjoin(select(fine_analysis.extrema, :metric, :extremum, :value=>:fine),
            select(coarse_analysis.extrema, :metric, :extremum, :value=>:coarse);
            on=[:metric, :extremum], validate=(true, true))
        @test nrow(metrics) == nrow(fine_analysis.extrema) == nrow(coarse_analysis.extrema)
        markdown_table(io, ["Metric", "Extremum", "dt1s", "dt10s", "dt10s minus dt1s"],
            [(row.metric, row.extremum, row.fine, row.coarse, row.coarse-row.fine) for row in eachrow(metrics)])
    end
    return Dict("csv_dir"=>coarse_paths.csv, "feather_dir"=>coarse_paths.feather,
        "duration_s"=>last(coarse.time), "retcode"=>"Endpoint verified; return code unavailable",
        "solver"=>missing, "earth_radius_m"=>model == "prototype" ? 6378137.0 : 6378136.6,
        "mu"=>model == "prototype" ? 3.986004418e14 : 3.98600436e14,
        "light_speed_mps"=>model == "prototype" ? 3.0e8 : 299792458.0)
end

function compare_saved_cadences()
    report = joinpath(SUMMARY_DIRECTORY, "dt1s_vs_dt10s.md")
    open(report, "w") do io
        println(io, "# Saved dt1s versus dt10s Results\n")
        println(io, "Existing outputs from both models are compared within each model, for cases 1 and 2. " *
            "No dynamics were rerun and no trajectory was downsampled. Target altitude is labelled `t`; helper altitude is labelled `h`.\n")
        println(io, "Checks cover complete output grids, all 21 spacecraft positions and velocities, shared timestamps including the final endpoint, " *
            "active helpers, full recorded event tables, and encounter extrema. Endpoint tolerance is 1e-8 s; interior timestamps match exactly. " *
            "Differences below are measurements, not acceptance tolerances.\n")
        println(io, "Missing dt1s analysis CSVs are reconstructed from their original Feather bundles in temporary directories using the existing extractor. " *
            "The dt10s analysis CSVs are read from the saved output folders. Original outputs are unchanged. " *
            "Runtime and solver return codes cannot be recovered from these Feather bundles. Constants in the generated case summaries use the existing model configuration.\n")
        for case in COMPARISON_CASES_V2[1:2]
            println(io, "## $(case.id)\n")
            settings = comparison_settings_v2(case)
            runs = [compare_saved_model(io, case, model, settings.duration_s) for model in ("prototype", "spaceagora")]
            name = "comparison_v2_$(case.id)_dt10s.md"
            write_report(joinpath(SUMMARY_DIRECTORY, name), settings, runs; saved_results=true)
            println(io, "[dt10s model comparison]($name)\n")
            GC.gc()
        end
        println(io, "## Interpretation\n")
        println(io, "Zero position/velocity differences mean identical saved dynamics at the shared times, not necessarily identical sampled extrema. " *
            "A one-second grid can capture motion extrema missed by a ten-second grid. Prototype accumulated delta-v uses trapezoidal integration " *
            "on saved samples, so its diagnostic can differ even when shared trajectory states and event intervals agree. " *
            "This comparison does not claim the two different models are physically equivalent.\n")
    end
    println("Saved cadence comparison: $report")
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_saved_cadences()
end