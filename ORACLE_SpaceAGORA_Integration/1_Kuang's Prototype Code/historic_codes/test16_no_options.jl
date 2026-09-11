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

########
# Main #
########
println("\nSTART SIMULATION (test16_no_options — copied from test11):")

### 1. Define orbits ###
helper_num             = 10
target_altitude_km     = 1000.0
target_inclination_deg = 0.0
nu_deg                 = 0.0

helper_oe = [
    (a_m=R_EARTH+1050e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=(360.0 / helper_num) * (j - 1)) for j in 1:helper_num
]
target_orbit = (a_m=R_EARTH + target_altitude_km*1e3, e=0.0, i_deg=target_inclination_deg, Ω_deg=0.0, ω_deg=0.0, ν_deg=nu_deg)
oe = vcat(helper_oe, [target_orbit])

### 2. Laser / cavity parameters ###
Pm  = zeros(length(oe), length(oe))
cav = Dict{Tuple{Int,Int},Dict{Symbol,Any}}()
for j in 1:helper_num
    cav[(j, helper_num + 1)] = Dict(:B => 100.0, :Pin => 1e4)
end

### 3. Run the simulation ###
println(@sprintf("\nCase: target_h=%.0f km, target_i=%.1f deg, target_nu=%.2f deg | helpers=%d",
    target_altitude_km, target_inclination_deg, nu_deg, helper_num))

elapsed = @elapsed begin
    sol, p, logs, masses = run_open_cavity_multi(oe;
        mass_kg      = 227,
        Pm           = Pm,
        cavity       = cav,
        use_los      = true,
        min_range    = 0.0,
        max_range    = 200e3,
        stop_on_dv   = false,
        T_seconds    = 63071,
        verbose      = false,
        result_plots = true,
        target_only  = true,
        IMG_DIR      = "Kuang's Prototype Code/output/images/",
        helper_num   = helper_num,
        gve_schedule = :none,
    )
end
println(@sprintf("Simulation runtime: %.2f s (%.2f min)", elapsed, elapsed/60))

### 4. Save results to CSV ###
case_csv_dir = @sprintf("Kuang's Prototype Code/output/CSV/target_h%.0fkm_i%.1fdeg_nu%.2fdeg",
    target_altitude_km, target_inclination_deg, nu_deg)
mkpath(case_csv_dir)
save_timeseries_csv(sol, p, helper_oe, target_orbit, csv_dir=case_csv_dir)

### 5. Animation (disabled for output comparison) ###
# fig1, controls1 = animate_all_satellites_3d_smooth_helper_target(sol, p, helper_num; tail=1000, show_earth=true, Δt=0.1)
# GLMakie.display(fig1)
