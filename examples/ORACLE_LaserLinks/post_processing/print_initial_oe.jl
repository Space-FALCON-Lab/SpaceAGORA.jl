#!/usr/bin/env julia

const _PP_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
import Pkg; Pkg.activate(_PP_REPO_ROOT; io=devnull)
# Reads a simulation_results.feather (or its parent directory) and prints the
# starting eccentricity (e) and true anomaly (ν) of sc1.
#
# Usage:
#   julia --project=. ORACLE/post_processing/print_initial_oe.jl \
#       output/paper_plot_mode/h1000km_t850km/ih0.0deg_it1.0deg/N1/T489261s/gve_ecc

using Arrow
using LinearAlgebra
using Printf

const MU = 3.986004418e14   # Earth gravitational parameter [m³/s²]

function compute_e_and_nu(pos::Vector{Float64}, vel::Vector{Float64})
    r     = norm(pos)
    h_vec = cross(pos, vel)
    e_vec = cross(vel, h_vec) / MU - pos / r
    e     = norm(e_vec)

    if e < 1e-10
        # circular orbit — nu undefined, return 0
        return (e=e, nu_deg=0.0)
    end

    cos_nu = clamp(dot(e_vec, pos) / (e * r), -1.0, 1.0)
    nu     = acos(cos_nu)
    # dot(r_vec, v_vec) > 0 means spacecraft is moving away from periapsis (0–π range)
    # dot(r_vec, v_vec) < 0 means spacecraft is moving toward periapsis (π–2π range)
    if dot(pos, vel) < 0.0
        nu = 2π - nu
    end
    return (e=e, nu_deg=rad2deg(nu))
end

function find_feather(path::String)
    if isfile(path) && endswith(path, ".feather")
        return path
    elseif isdir(path)
        candidate = joinpath(path, "simulation_results.feather")
        isfile(candidate) && return candidate
        # search one level deeper (e.g. T*s sub-folder)
        for entry in readdir(path; join=true)
            isdir(entry) || continue
            c = joinpath(entry, "simulation_results.feather")
            isfile(c) && return c
        end
    end
    return nothing
end

if isempty(ARGS)
    println("Usage: julia --project=. ORACLE/post_processing/print_initial_oe.jl <path>")
    println("  <path> can be the feather file itself or a directory containing it.")
    exit(1)
end

input_path = ARGS[1]
# Resolve relative to repo root if not absolute
if !isabspath(input_path)
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    input_path = joinpath(repo_root, input_path)
end

feather_path = find_feather(input_path)
if feather_path === nothing
    println("ERROR: could not find simulation_results.feather under: $input_path")
    exit(1)
end

tbl = Arrow.Table(feather_path)
pos0 = Float64[tbl.sc1_pos_1[1], tbl.sc1_pos_2[1], tbl.sc1_pos_3[1]]
vel0 = Float64[tbl.sc1_vel_1[1], tbl.sc1_vel_2[1], tbl.sc1_vel_3[1]]

result = compute_e_and_nu(pos0, vel0)

println("Feather: $feather_path")
@printf("  Starting e   = %.6f\n",       result.e)
@printf("  Starting ν   = %.4f deg\n",   result.nu_deg)
