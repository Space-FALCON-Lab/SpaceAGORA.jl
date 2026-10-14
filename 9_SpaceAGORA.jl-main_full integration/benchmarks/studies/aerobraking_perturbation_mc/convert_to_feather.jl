# Converts trajectory_with_active_force.csv → .feather for all cases in a run directory.
# Usage: julia --project=. benchmarks/studies/aerobraking_perturbation_mc/convert_to_feather.jl <run_dir>
#        Omit <run_dir> to use the most recent timestamped run in output/aerobraking_perturbation_mc/.
# Options:
#   --overwrite   Replace existing .feather files

using Arrow
using CSV
using DataFrames
using Printf

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_BASE = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")

function _to_arrow_safe!(df::DataFrame)
    for name in names(df)
        col = df[!, name]
        if any(v -> v isa Symbol, col)
            df[!, name] = [v isa Symbol ? String(v) : v for v in col]
        end
    end
    return df
end

function convert_run(run_dir::String; overwrite::Bool=false)
    csv_files = String[]
    for (root, _, files) in walkdir(run_dir)
        for f in files
            f == "trajectory_with_active_force.csv" && push!(csv_files, joinpath(root, f))
        end
    end
    sort!(csv_files)
    isempty(csv_files) && (println("No trajectory_with_active_force.csv files found in $run_dir"); return 0)
    println("Found $(length(csv_files)) CSV files")

    converted = 0
    skipped = 0
    for (i, csv_path) in enumerate(csv_files)
        feather_path = replace(csv_path, ".csv" => ".feather")
        if isfile(feather_path) && !overwrite
            skipped += 1
            continue
        end
        df = CSV.read(csv_path, DataFrame)
        _to_arrow_safe!(df)
        Arrow.write(feather_path, df)
        converted += 1
        (converted % 25 == 0 || converted == 1) &&
            @printf("[%3d/%d] %s\n", i, length(csv_files), basename(dirname(csv_path)))
        flush(stdout)
    end
    println("Done: converted=$(converted), skipped=$(skipped)")
    return converted
end

function _most_recent_run(base::String)
    isdir(base) || error("Output directory not found: $base")
    entries = filter(
        d -> isdir(joinpath(base, d)) && occursin(r"^\d{8}_\d{6}$", d),
        readdir(base),
    )
    isempty(entries) && error("No timestamped run directories found in $base")
    return joinpath(base, last(sort(entries)))
end

function main(args=ARGS)
    args = collect(args)
    overwrite = "--overwrite" in args
    positional = filter(a -> !startswith(a, "-"), args)
    run_dir = isempty(positional) ? _most_recent_run(DEFAULT_OUTPUT_BASE) : abspath(positional[1])
    isdir(run_dir) || error("Not a directory: $run_dir")
    println("Run directory: $run_dir")
    convert_run(run_dir; overwrite=overwrite)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
