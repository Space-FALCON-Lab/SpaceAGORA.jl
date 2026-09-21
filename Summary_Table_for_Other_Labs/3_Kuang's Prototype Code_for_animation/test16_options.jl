# xvfb-run -a julia --startup-file=no --project=2_SpaceAGORA.jl "3_Kuang's Prototype Code_for_animation/test16_options.jl"
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
const Ẑ       = SVector(0.0, 0.0, 1.0)

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
include("functions/13_Feather_Write.jl")
include("extract_feather_encounters.jl")

############
# Settings #
############
Base.@kwdef struct OracleOptions
    helpers::Int                    = 100
    helper_altitude_km::Float64     = 1000.0
    target_altitude_km::Float64     = 1050.0
    target_inclination_deg::Float64 = 5.0
    helper_inclination_deg::Float64 = 0.0
    target_nu_deg::Float64          = 0.0
    target_ecc::Float64             = 0.0
    orbits::Float64                 = 10.0
    schedule::Symbol                = :gve_sma
    laser_range_km::Float64         = 200.0
    laser_power_w::Float64          = 10_000.0
    magnification::Float64          = 100.0
    beta::Float64                   = 1.0
    eta::Float64                    = 1.0
    mass_kg::Float64                = 227.0
    dt_max_s::Float64               = 10.0
    planet::Symbol                  = :earth
    paper_grid::Bool                = false
    feather_only::Bool              = false
    output_dir::String              = "output"
    timeseries_points::Int          = 1001
    use_los::Bool                   = true
    min_range_km::Float64           = 0.0
    use_J2::Bool                    = true
    useDrag::Bool                   = false
    verbose::Bool                   = false
    result_plots::Bool              = true
    target_only::Bool               = true
    animate::Bool                   = true
    show_earth::Bool                = true
    duration_seconds::Float64       = 100.0
    animation_fps::Float64          = 30.0
    helper_trails::Union{Bool,Vector{Int}} = false
    target_trail::Bool              = true
    show_projections::Bool          = true
    helper_projections::Bool        = false
    target_projections::Bool        = true
    radial_exaggeration::Float64    = 40.0
    axis_limit_km::Union{Nothing,Float64} = nothing
    radial_ticks_km::Union{Nothing,Vector{Float64}} = [500.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0]
end

########
# Main #
########
function run_animation_case(opts::OracleOptions=OracleOptions())
println("\nSTART SIMULATION (test16_options):")
opts.planet == :earth || throw(ArgumentError("This prototype supports only planet=:earth"))
!opts.paper_grid || throw(ArgumentError("paper_grid is not supported by this single-case prototype"))
opts.beta == 1.0 && opts.eta == 1.0 || throw(ArgumentError("This prototype supports only beta=eta=1; use laser_power_w and magnification to configure its force model"))
opts.timeseries_points == 1001 || throw(ArgumentError("timeseries_points is a compatibility field, not a sampling control; this prototype saves every 10 seconds plus the endpoint"))
opts.schedule in (:none, :gve_sma, :gve_ecc, :gve_inc, :gve_raan, :gve_argp) ||
    throw(ArgumentError("Unsupported prototype schedule: $(opts.schedule)"))
isfinite(opts.orbits) && opts.orbits > 0 || throw(ArgumentError("orbits must be positive and finite"))
isfinite(opts.dt_max_s) && opts.dt_max_s > 0 || throw(ArgumentError("dt_max_s must be positive and finite"))

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
        use_los      = opts.use_los,
        min_range    = opts.min_range_km * 1e3,
        max_range    = opts.laser_range_km * 1e3,
        stop_on_dv   = false,
        T_seconds    = T_seconds,
        dt_max_s     = opts.dt_max_s,
        use_J2       = opts.use_J2,
        useDrag      = opts.useDrag,
        verbose      = opts.verbose,
        result_plots = opts.result_plots && !opts.feather_only,
        target_only  = opts.target_only,
        IMG_DIR      = joinpath(@__DIR__, opts.output_dir, "images") * "/",
        helper_num   = helper_num,
        gve_schedule = opts.schedule,
        record_events = true
    )
end
println(@sprintf("Simulation runtime: %.2f s (%.2f min)", elapsed, elapsed/60))

### 5. Save results to CSV ###
paths = scenario_paths(joinpath(@__DIR__, opts.output_dir), opts, sol.t[end];
    source="prototype", schedule=opts.schedule, use_J2=opts.use_J2)
case_csv_dir = paths.csv
feather_path = save_timeseries_feather(sol, p; feather_dir=paths.feather)
if !opts.feather_only
    mkpath(case_csv_dir)
    save_timeseries_csv(sol, p, helper_oe, target_orbit, csv_dir=case_csv_dir)
    extract_feather_encounters(feather_path; output_dir=paths.csv)
end

### 6. Animation ###
if opts.animate && !opts.feather_only
    video_path = joinpath(@__DIR__, opts.output_dir, "videos", basename(paths.feather), "animation.mp4")
    animate_all_satellites_3d_smooth_helper_target(sol, p, helper_num;
        show_earth=opts.show_earth, output_file=video_path,
        duration_seconds=opts.duration_seconds, animation_fps=opts.animation_fps,
        helper_trails=opts.helper_trails, target_trail=opts.target_trail,
        show_projections=opts.show_projections,
        helper_projections=opts.helper_projections, target_projections=opts.target_projections,
        radial_exaggeration=opts.radial_exaggeration,
        reference_altitude_km=opts.helper_altitude_km, axis_limit_km=opts.axis_limit_km,
        radial_ticks_km=opts.radial_ticks_km)
end
return (; sol, params=p, paths, feather_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_animation_case()
end
