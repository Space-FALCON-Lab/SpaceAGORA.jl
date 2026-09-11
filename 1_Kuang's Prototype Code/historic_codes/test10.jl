import GLMakie
GLMakie.activate!()
import GeometryBasics: Point3f, Vec3f
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra, StaticArrays
using Plots, Printf
using Statistics 
using Revise # for automatically track changes to any files included using include() and update them in your current Julia session.
using FileIO  # Optional, for loading textures
using GeometryBasics  # For creating 3D geometries
using Observables  # For reactive programming in Makie
using DelimitedFiles  # For saving CSVs without extra deps

print("\033c")  # Clears the terminal on Windows

#############
# Constants #
#############
const MU       = 3.986004418e14         # Earth μ [m^3/s^2]
const C        = 3.0e8                  # speed of light [m/s]
const R_EARTH  = 6_378_137.0           # Earth mean radius [m]
#const R_EARTH  = 6328136.6             # Earth mean radius [m]
const R_ATMDEF = R_EARTH + 100_000.0    # default "atmosphere radius" (Kármán line) [m]
const Ẑ       = SVector(0.0, 0.0, 1.0)

@inline idx(i, off) = 6*(i-1) + off  # state indexing helper

#############
# Functions #
#############
include("functions/1_LOS_Metrics.jl")  
include("functions/2_Laser_Forces_ver2.jl")
include("functions/3_Dynamics.jl")
include("functions/4_Diagnostics.jl")
include("functions/5_OE_Converters.jl")
include("functions/6_OE_and_dv_in_RTN.jl")
include("functions/7_Plots.jl")
include("functions/8_LoS_time_series.jl")
include("functions/9_Runners.jl")
include("functions/10_Animation_ver2.jl")
include("functions/12_CSV_Write_Read.jl")

############
# Settings #
############
Base.@kwdef struct OracleOptions
    helpers::Int                    = 10
    helper_altitude_km::Float64     = 1050.0
    target_altitude_km::Float64     = 1000.0
    target_inclination_deg::Float64 = 0.0
    helper_inclination_deg::Float64 = 0.0
    target_nu_deg::Float64          = 0.0
    target_ecc::Float64             = 0.0
    orbits::Float64                 = 10.0
    schedule::Symbol                = :naive_next_entering
    laser_range_km::Float64         = 200.0
    laser_power_w::Float64          = 10_000.0
    magnification::Float64          = 100.0
    beta::Float64                   = 1.0
    eta::Float64                    = 1.0
    mass_kg::Float64                = 227.0
    dt_max_s::Float64               = 10.0
    planet::Symbol                  = :earth  # planet symbol passed to make_no_gram_planet
    paper_grid::Bool                = false
    feather_only::Bool              = false
    output_dir::String              = "output"
    timeseries_points::Int          = 1001
    animate::Bool                   = false
end

########
# Main #
########
println("\nSTART SIMULATION:")

opts = OracleOptions()

### 1. Define orbits ###
helper_num = opts.helpers

helper_oe = [
    (a_m = R_EARTH + opts.helper_altitude_km*1e3, e = 0.0,
     i_deg = opts.helper_inclination_deg, Ω_deg = 0.0, ω_deg = 0.0,
     ν_deg = (360.0 / helper_num) * (j - 1)) for j in 1:helper_num
]
target_orbit = (a_m = R_EARTH + opts.target_altitude_km*1e3, e = opts.target_ecc,
                i_deg = opts.target_inclination_deg, Ω_deg = 0.0, ω_deg = 0.0,
                ν_deg = opts.target_nu_deg)
oe = vcat(helper_oe, [target_orbit])

### 2. Laser / cavity parameters ###
Pm  = zeros(length(oe), length(oe))
cav = Dict{Tuple{Int,Int},Dict{Symbol,Any}}()
for j in 1:helper_num
    cav[(j, helper_num + 1)] = Dict(:B => opts.magnification, :Pin => opts.laser_power_w)
end

### 3. Compute simulation duration from orbits ###
a_target  = R_EARTH + opts.target_altitude_km*1e3
T_orbit   = 2π * sqrt(a_target^3 / MU)
T_seconds = opts.orbits * T_orbit

### 4. Run the simulation ###
println(@sprintf("\nCase: target_h=%.0f km, target_i=%.1f deg, target_nu=%.2f deg | helpers=%d",
    opts.target_altitude_km, opts.target_inclination_deg, opts.target_nu_deg, helper_num))

elapsed = @elapsed begin
    sol, p, logs, masses = run_open_cavity_multi(oe;
        mass_kg      = opts.mass_kg,
        Pm           = Pm,
        cavity       = cav,
        use_los      = true,
        min_range    = 0.0,
        max_range    = opts.laser_range_km * 1e3,
        stop_on_dv   = false,
        T_seconds    = 63071, #T_seconds,
        verbose      = false,
        result_plots = true,
        target_only  = true,
        IMG_DIR      = joinpath(@__DIR__, opts.output_dir, "images") * "/",
        helper_num   = helper_num,
        gve_schedule = :none #opts.schedule,  # :none | :naive_next_entering | :positive_along_track | :gve_sma | :gve_ecc | :gve_inc | :gve_raan | :gve_argp
    )
end
println(@sprintf("Simulation runtime: %.2f s (%.2f min)", elapsed, elapsed/60))

### 5. Save results to CSV ###
case_csv_dir = joinpath(@__DIR__, opts.output_dir, @sprintf("CSV/target_h%.0fkm_i%.1fdeg_nu%.2fdeg",
    opts.target_altitude_km, opts.target_inclination_deg, opts.target_nu_deg))
mkpath(case_csv_dir)
save_timeseries_csv(sol, p, helper_oe, target_orbit, csv_dir=case_csv_dir)

### 6. Animation ###
if opts.animate
    fig1, controls1 = animate_all_satellites_3d_smooth_helper_target(sol, p, helper_num; tail=1000, show_earth=true, Δt=0.1)
    GLMakie.display(fig1)
end