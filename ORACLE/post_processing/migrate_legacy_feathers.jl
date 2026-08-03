#!/usr/bin/env julia
# Moves legacy feather/toml output files that sit directly inside a T*s/ folder
# (i.e. not yet in a schedule subfolder) into a new naive_next_entering/ subfolder,
# so the directory layout matches what the updated 9_Runners.jl now produces.
#
# Safe to re-run: skips T*s/ folders that have no bare feather/toml files.
#
# Usage (from repo root):
#   julia --project=. ORACLE/post_processing/migrate_legacy_feathers.jl

const REPO_ROOT     = normpath(joinpath(@__DIR__, "..", ".."))
const PAPER_ROOT    = joinpath(REPO_ROOT, "output", "paper_plot_mode")
const SCHEDULE_DIRS = ("naive_next_entering", "positive_along_track")
const MOVE_EXTS     = (".feather", ".toml", ".csv")

function migrate()
    isdir(PAPER_ROOT) || error("Paper output directory not found: $PAPER_ROOT")

    moved_total = 0

    # Walk every T*s/ folder inside the expected hierarchy
    for (root, dirs, files) in walkdir(PAPER_ROOT)
        # Only act on directories whose name looks like T<digits>s
        basename(root) == "summary.csv" && continue
        b = basename(root)
        startswith(b, "T") && endswith(b, "s") || continue

        # Skip if this folder is itself already a schedule subfolder
        parent_b = basename(dirname(root))
        parent_b in SCHEDULE_DIRS && continue

        # Find files directly in this T*s/ folder that match the target extensions
        legacy_files = filter(f -> any(endswith(f, ext) for ext in MOVE_EXTS), files)
        isempty(legacy_files) && continue

        dest_dir = joinpath(root, "naive_next_entering")
        mkpath(dest_dir)

        for f in legacy_files
            src = joinpath(root, f)
            dst = joinpath(dest_dir, f)
            if isfile(dst)
                println("  SKIP (already exists): $dst")
                continue
            end
            mv(src, dst)
            println("  MOVED: $src\n      → $dst")
            moved_total += 1
        end
    end

    println("\nDone. $moved_total file(s) moved into naive_next_entering/ subfolders.")
end

migrate()
