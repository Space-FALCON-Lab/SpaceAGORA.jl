using CSV, DataFrames, SHA

function summarize_deliverable_extrema(root=@__DIR__)
    deliverables = joinpath(root, "6_Deliverables")
    folders = filter(isdir, readdir(deliverables; join=true))
    tables = map(folders) do folder
        report = only(filter(name -> occursin(r"^comparison_v2_case\d+_.*\.md$", name),
            readdir(folder)))
        case_id = parse(Int, match(r"_case(\d+)_", report).captures[1])
        table = CSV.read(joinpath(folder, "extrema_analysis.csv"), DataFrame)
        insertcols!(table, 1, :test_case => fill(case_id, nrow(table)),
            :scenario => fill(basename(folder), nrow(table)))
        table
    end
    combined = sort!(vcat(tables...), :test_case)
    summary = combined[Int[], :]
    for group in groupby(combined, [:metric, :extremum])
        values = collect(skipmissing(group.value))
        isempty(values) && continue
        extremum = first(group.extremum)
        extremum in ("min", "max") || error("Unknown extremum: $extremum")
        selected = extremum == "min" ? minimum(values) : maximum(values)
        append!(summary, filter(row -> !ismissing(row.value) && row.value == selected, group))
    end
    path = joinpath(deliverables, "extrema_analysis_all_test_cases.csv")
    CSV.write(path, summary)
    @assert isequal(CSV.read(path, DataFrame), summary)
    println("Verified cross-case extrema from ", length(tables), " scenarios: ", path)
    return summary
end

function export_deliverables(root=@__DIR__)
    prototype_root = joinpath(root, "1_Kuang's Prototype Code", "output")
    spaceagora_root = joinpath(root, "2_SpaceAGORA.jl", "output", "CSV")
    scenarios = filter(name -> startswith(name, "prototype_") &&
        isdir(joinpath(prototype_root, "feather", name)),
        readdir(joinpath(prototype_root, "feather")))
    isempty(scenarios) && error("No prototype scenarios found")
    reports = filter(path -> endswith(path, ".md"),
        readdir(joinpath(root, "4_Comparasion_Summary"); join=true))
    bundles = map(scenarios) do prototype_scenario
        scenario = replace(prototype_scenario, r"^prototype_" => "")
        report = only(filter(path -> occursin("spaceagora_" * scenario, read(path, String)), reports))
        source = joinpath(spaceagora_root, "spaceagora_" * scenario)
        encounters_path = joinpath(source, "sim_output", "encounter_record.csv")
        extrema_path = joinpath(source, "analysis", "extrema.csv")
        encounters = CSV.read(encounters_path, DataFrame)
        extrema = CSV.read(extrema_path, DataFrame)
        longest = only(eachrow(filter(row -> row.metric == "encounter_duration_s" &&
            row.extremum == "max", extrema)))
        selected = only(eachrow(filter(row -> row.encounter_id == longest.encounter_id, encounters)))
        @assert selected.duration_s == longest.value == maximum(encounters.duration_s)
        trajectory_path = joinpath(source, "analysis", "encounter_states",
            "encounter_$(longest.encounter_id).csv")
        trajectory = CSV.read(trajectory_path, DataFrame)
        insertcols!(trajectory, 1, :encounter_id => fill(longest.encounter_id, nrow(trajectory)))
        video_root = joinpath(prototype_root, "videos", prototype_scenario)
        copies = [
            joinpath(source, "sim_output", "laser_link_record.csv") => "laser_link_record.csv",
            extrema_path => "extrema_analysis.csv",
            joinpath(video_root, prototype_scenario * ".mp4") => prototype_scenario * ".mp4",
            joinpath(video_root, prototype_scenario * ".png") => prototype_scenario * ".png",
        ]
        for (path, _) in copies
            isfile(path) || error("Missing source: $path")
        end
        encounter_table = select(encounters, :encounter_id, :sc_a, :sc_b,
            :start_time_s, :end_time_s, :duration_s)
        (; scenario, encounter_table, trajectory, copies, report, encounter_id=longest.encounter_id)
    end
    for bundle in bundles
        destination = joinpath(root, "6_Deliverables", bundle.scenario)
        mkpath(destination)
        for (filename, table) in (
            "encounter_record.csv" => bundle.encounter_table,
            "Maximum-Duration Encounter Trajectory.csv" => bundle.trajectory,
        )
            path = joinpath(destination, filename)
            CSV.write(path, table)
            @assert isequal(CSV.read(path, DataFrame), table)
        end
        for (source, filename) in bundle.copies
            destination_path = joinpath(destination, filename)
            cp(source, destination_path; force=true)
            @assert open(sha256, source) == open(sha256, destination_path)
        end
        report_content = replace(read(bundle.report, String), "](../" => "](../../")
        report_path = joinpath(destination, basename(bundle.report))
        write(report_path, report_content)
        @assert read(report_path, String) == report_content
        println("VERIFIED ", bundle.scenario, ": ", nrow(bundle.encounter_table),
            " encounters; maximum-duration encounter ", bundle.encounter_id,
            "; ", nrow(bundle.trajectory), " trajectory rows; 7 deliverables")
    end
    println("Exported and verified ", length(bundles), " scenario folders in 6_Deliverables")
    summarize_deliverable_extrema(root)
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_deliverables()
end