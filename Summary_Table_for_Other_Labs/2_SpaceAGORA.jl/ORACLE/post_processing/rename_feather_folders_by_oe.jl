#!/usr/bin/env julia
# Walks every simulation_results.feather under output/paper_plot_mode, computes
# the starting eccentricity (e) and true anomaly (ν) of sc1, and renames the
# containing schedule folder to:
#
#   {schedule_base}_e{e:.4f}_nu{nu_deg:.4f}
#
# where {schedule_base} is the existing folder name with any trailing _e{…}
# suffix stripped so it is never doubled.
#
# Dry-run by default — prints what would change without touching the filesystem.
# Pass --execute to perform the renames and update each manifest's absolute path.
#
# Usage:
#   julia --project=. ORACLE/post_processing/rename_feather_folders_by_oe.jl
#   julia --project=. ORACLE/post_processing/rename_feather_folders_by_oe.jl --execute
#   julia --project=. ORACLE/post_processing/rename_feather_folders_by_oe.jl --base-dir output/paper_plot_mode/h1000km_t850km

using Arrow
using LinearAlgebra
using Printf
using TOML

const MU = 3.986004418e14   # Earth gravitational parameter [m³/s²]
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# ── CLI flags ─────────────────────────────────────────────────────────────────
const DRY_RUN  = !("--execute" in ARGS)
# By default only rename folders that already carry an _e{…} suffix
# (the new-style runs).  Pass --all to also process legacy bare-name folders
# like gve_inc, gve_sma, naive_next_entering, etc.
const ALL_FOLDERS = "--all" in ARGS

_bd_idx  = findfirst(==("--base-dir"), ARGS)
const BASE_DIR = if _bd_idx !== nothing
    p = ARGS[_bd_idx + 1]
    isabspath(p) ? p : joinpath(REPO_ROOT, p)
else
    joinpath(REPO_ROOT, "output", "paper_plot_mode")
end

# ── Orbital element helpers ───────────────────────────────────────────────────
function compute_e_and_nu(pos::Vector{Float64}, vel::Vector{Float64})
    r     = norm(pos)
    h_vec = cross(pos, vel)
    e_vec = cross(vel, h_vec) / MU - pos / r
    e     = norm(e_vec)
    if e < 1e-6
        # Circular orbit — periapsis direction is undefined; use ν = 0
        return (e=e, nu_deg=0.0)
    end
    cos_nu = clamp(dot(e_vec, pos) / (e * r), -1.0, 1.0)
    nu     = acos(cos_nu)
    dot(pos, vel) < 0.0 && (nu = 2π - nu)
    return (e=e, nu_deg=rad2deg(nu))
end

# ── SHA256 reader (from manifest — avoids re-hashing the feather) ─────────────
function read_sha256(dir::String)::Union{String, Nothing}
    manifest = joinpath(dir, "simulation_results.manifest.toml")
    isfile(manifest) || return nothing
    d = TOML.parsefile(manifest)
    # Structure: d["files"]["feather"]["sha256"]
    return get(get(get(d, "files", Dict()), "feather", Dict()), "sha256", nothing)
end

# ── Folder-name builder ───────────────────────────────────────────────────────
# Strip any existing trailing _e{digits.digits} so it is never doubled, then
# append _e{e:.4f}_nu{nu_deg:.4f}.
function new_folder_name(old_name::String, e::Float64, nu_deg::Float64)::String
    # Strip everything from the first _e{digits.digits} onward (handles both
    # bare names, _e-only names, and already-suffixed _e..._nu... names).
    base = replace(old_name, r"_e\d+\.\d+.*$" => "")
    @sprintf("%s_e%.4f_nu%.4f", base, e, nu_deg)
end

# ── Collect all feather files ─────────────────────────────────────────────────
feather_files = String[]
for (root, _, files) in walkdir(BASE_DIR)
    "simulation_results.feather" in files &&
        push!(feather_files, joinpath(root, "simulation_results.feather"))
end
sort!(feather_files)

isempty(feather_files) && (println("No feather files found under $BASE_DIR."); exit(0))

println("Found $(length(feather_files)) feather file(s) under:\n  $BASE_DIR\n")
DRY_RUN && println("DRY RUN — pass --execute to apply renames\n")

n_rename  = 0
n_skip    = 0
n_clash   = 0
n_deleted = 0

for (idx, feather_path) in enumerate(feather_files)
    sched_dir  = dirname(feather_path)      # …/T489261s/gve_inc_e0.0000
    old_name   = basename(sched_dir)
    parent_dir = dirname(sched_dir)         # …/T489261s

    # Read first row only (Arrow.Table is memory-mapped — fast)
    tbl  = Arrow.Table(feather_path)
    pos0 = Float64[tbl.sc1_pos_1[1], tbl.sc1_pos_2[1], tbl.sc1_pos_3[1]]
    vel0 = Float64[tbl.sc1_vel_1[1], tbl.sc1_vel_2[1], tbl.sc1_vel_3[1]]
    oe   = compute_e_and_nu(pos0, vel0)

    new_name = new_folder_name(old_name, oe.e, oe.nu_deg)
    new_dir  = joinpath(parent_dir, new_name)

    # ── Skip legacy bare-name folders unless --all was passed ─────────────────
    if !ALL_FOLDERS && !occursin(r"_e\d+\.\d+", old_name)
        global n_skip += 1
        print("\r[$idx/$(length(feather_files))]")
        continue
    end

    # ── Progress tick ──────────────────────────────────────────────────────────
    print("\r[$idx/$(length(feather_files))]")

    if old_name == new_name
        global n_skip += 1
        continue
    end

    # Relative path for concise output
    rel = relpath(sched_dir, BASE_DIR)
    @printf("\n  %-55s  e=%.6f  ν=%9.4f°\n    → %s\n", rel, oe.e, oe.nu_deg, new_name)

    if isdir(new_dir)
        # ── Clash: destination already exists ─────────────────────────────────
        sha_src = read_sha256(sched_dir)
        sha_dst = read_sha256(new_dir)
        if sha_src !== nothing && sha_src == sha_dst
            @printf("    ✓ identical (SHA256 match) — %s redundant source\n",
                    DRY_RUN ? "would delete" : "deleting")
            !DRY_RUN && rm(sched_dir; recursive=true)
            global n_deleted += 1
        else
            @printf("    ✗ different data — skipped (manual review needed)\n")
            global n_clash += 1
        end
        continue   # skip n_rename increment below
    end

    if !DRY_RUN
        # 1. Rename the directory
        mv(sched_dir, new_dir)

        # 2. Fix the absolute path baked into the manifest
        manifest_path = joinpath(new_dir, "simulation_results.manifest.toml")
        if isfile(manifest_path)
            txt = read(manifest_path, String)
            txt = replace(txt, sched_dir => new_dir)
            write(manifest_path, txt)
        end
    end

    global n_rename += 1
end

println("\n")
if DRY_RUN
    println("Would rename  : $n_rename folder(s)")
    println("Would delete  : $n_deleted folder(s) (identical, redundant)")
    println("Already OK    : $n_skip folder(s)")
    n_clash > 0 && println("Would skip    : $n_clash folder(s) (clash with different data — manual review)")
    println("\nRun with --execute to apply.")
else
    println("Renamed       : $n_rename folder(s)")
    println("Deleted       : $n_deleted folder(s) (identical redundant sources)")
    println("Skipped (ok)  : $n_skip folder(s)")
    n_clash > 0 && println("Skipped clash : $n_clash folder(s) (different data — manual review needed)")
    println("\n⚠  plot_oe_from_feather.jl reads folder names — update its sched_folder")
    println("   logic (or re-run with the new --target-nu-deg flag) if you plan to plot again.")
end
