import GLMakie
import GeometryBasics: Point3f, Vec3f
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra, StaticArrays
using Plots, Printf
using Statistics 
using Revise # for automatically track changes to any files included using include() and update them in your current Julia session.
using FileIO  # Optional, for loading textures
using GeometryBasics  # For creating 3D geometries
using Observables  # For reactive programming in Makie

print("\033c")  # Clears the terminal on Windows

# --------- Constants ---------
const MU       = 3.986004418e14         # Earth μ [m^3/s^2]
const C        = 3.0e8                  # speed of light [m/s]
const R_EARTH  = 6_378_137.0           # Earth mean radius [m]
const R_ATMDEF = R_EARTH + 100_000.0    # default "atmosphere radius" (Kármán line) [m]
const Ẑ       = SVector(0.0, 0.0, 1.0)

@inline idx(i, off) = 6*(i-1) + off  # state indexing helper

# --------- Functions ---------
#include("functions/0_Module_Setup.jl")
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

# --------- Example Run ---------
println("\nSTART SIMULATION:")


helper_num = 1 # Number of helper satellites
# Generate orbital elements for helper satellites
helper_oe = [
    (a_m=R_EARTH+1000e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=(360.0 / helper_num) * (i - 1)) for i in 1:helper_num
]
# Target orbit parameters
target_orbit = (a_m=R_EARTH+1000e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.5)
# Add the target satellite to the orbital elements
oe = vcat(helper_oe, [target_orbit])

Pm = zeros(length(oe), length(oe)); # initialize power matrix with zeros
# # assume all helpers can send single-pass laser to the target
# for i in 1:helper_num
#     Pm[i, helper_num + 1] = 1.0e5; # 100 kW single-pass from helpers to target
# end

cav = Dict{Tuple{Int, Int}, Dict{Symbol, Any}}() # Initialize the cav dictionary
for i in 1:helper_num # Construct cav for helper satellites (1 to helper_num) linking to the target satellite (index N)
    cav[(i, helper_num + 1)] = Dict(:B => 100.0, :Pin => 1e4)
end

min_range = 0.0 # min range limit for laser forces in meters
max_range = 200e3 # max range limit for laser forces in meters

sol, p, logs = run_open_cavity_multi(oe;
                                    mass_kg=227,
                                    Pm=Pm,
                                    cavity=cav,
                                    use_los=true,
                                    min_range=min_range, max_range=max_range, # unit: meters
                                    stop_on_dv=false,
                                    T_seconds=10*3600, # default to automatic based on orbits
                                    verbose=true, result_plots=true,target_only=true,
                                    IMG_DIR="Kuang's Prototype Code/images/",
                                    helper_num = helper_num
                                    )

# ---------- Post-Processing ---------
# println("=============== Link Duty Cycle Analysis ===============")
# for i in 1:length(oe), j in 1:length(oe)
#         if i != j
#                 info = link_duty_and_estimate(sol, p; i=i, j=j)
#                 @printf("Link %d→%d: duty ≈ %.1f%%, on_time ≈ %.1f s, a_nom ≈ %.4g m/s², DV upper bound ≈ %.2f m/s\n",
#                                 i, j, 100*info.duty, info.on_time, info.a_nom, info.dv_upper_bound)
#         end
# end

# Display the first animation (satellite orbits)
fig1, controls1 = animate_all_satellites_3d_smooth_helper_target(sol, p, helper_num; tail=1000, show_earth=true, Δt=0.1)
GLMakie.display(fig1)  # Explicitly display the first figure