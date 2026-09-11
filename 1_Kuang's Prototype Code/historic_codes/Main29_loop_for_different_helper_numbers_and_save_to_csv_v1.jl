import Makie
import WGLMakie
WGLMakie.activate!()
const GLMakie = Makie
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

########
# Main #
########
println("\nSTART SIMULATION:")

for i in vcat(1, 5:15:160)
    ### 1. Define helper satellite orbits and target orbit ###
    helper_num = i # Number of helper satellites
    # Generate orbital elements for helper satellites
    helper_oe = [
        (a_m=R_EARTH+1500e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=(360.0 / helper_num) * (i - 1)) for i in 1:helper_num
    ]
    # Target orbit parameters
    target_orbit = (a_m=R_EARTH+1700e3, e=0.0, i_deg=1.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0)
    # Add the target satellite to the orbital elements
    oe = vcat(helper_oe, [target_orbit])

    ### 2. Define laser power matrix and cavity parameters ###
    Pm = zeros(length(oe), length(oe)); # initialize power matrix with zero
    cav = Dict{Tuple{Int, Int}, Dict{Symbol, Any}}() # Initialize the cav dictionary
    for i in 1:helper_num # Construct cav for helper satellites (1 to helper_num) linking to the target satellite (index N)
        cav[(i, helper_num + 1)] = Dict(:B => 100.0, :Pin => 1e4)
    end
    min_range = 0.0 # min range limit for laser forces in meters
    max_range = 200e3 # max range limit for laser forces in meters

    ### 3. Run the simulation ###
    elapsed = @elapsed begin # tic
    sol, p, logs = run_open_cavity_multi(oe;
                                        mass_kg=227,
                                        Pm=Pm,
                                        cavity=cav,
                                        use_los=true,
                                        min_range=min_range, max_range=max_range, # unit: meters
                                        stop_on_dv=false,
                                        T_seconds=300*3600, # default to automatic based on orbits
                                        verbose=false, result_plots=true, target_only=true,
                                        IMG_DIR="Kuang's Prototype Code/images/",
                                        helper_num = helper_num
                                        )
    end
    println(@sprintf("Simulation runtime: %.2f s (%.2f min)", elapsed, elapsed/60)) # toc

    ### 4. Save results to CSV ###
    save_timeseries_csv(sol, p, helper_oe, target_orbit) 
end