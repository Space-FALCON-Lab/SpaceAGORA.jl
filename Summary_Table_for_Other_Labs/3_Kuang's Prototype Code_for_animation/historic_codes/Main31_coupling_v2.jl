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

# print("\033c")  # Clears the terminal on Windows

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

helper_counts = 1 #vcat(150:50:300) #vcat(50:50:100) # vcat(50:50:300) #vcat(150:50:300) #vcat(1, 50:50:300) # vcat(50:50:200)
target_altitudes_km = 1010 #850 # [950.0, 1000.0, 1050.0, 1150.0] # [850.0, 950.0, 1000.0, 1050.0, 1150.0] # [800]
target_inclinations_deg = 0.0 #[0.0, 0.5, 1.0] #[0.0, 0.5, 1.0, 1.5]
# TEST helper_counts = vcat(50:50:100) for [950.0, 1000.0, 1050.0, 1150.0] and [0.0, 0.5, 1.0, 1.5] !!!
for alt_km in target_altitudes_km
    for inc_deg in target_inclinations_deg
        case_csv_dir = normpath(joinpath(@__DIR__, "output", "CSV",
                            @sprintf("target_h%.0fkm_i%.1fdeg", alt_km, inc_deg)))
        mkpath(case_csv_dir)

        for helper_num in helper_counts
            ### 1. Define helper satellite orbits and target orbit ###
            # Generate orbital elements for helper satellites
            helper_oe = [
                (a_m=R_EARTH+1000e3, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=(360.0 / helper_num) * (j - 1)) for j in 1:helper_num
            ]
            # Target orbit parameters
            target_orbit = (a_m=R_EARTH + alt_km * 1e3, e=0.0, i_deg=inc_deg, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0)
            # Add the target satellite to the orbital elements
            oe = vcat(helper_oe, [target_orbit])

            ### 2. Define laser power matrix and cavity parameters ###
            Pm = zeros(length(oe), length(oe)); # initialize power matrix with zero
            cav = Dict{Tuple{Int, Int}, Dict{Symbol, Any}}() # Initialize the cav dictionary
            for j in 1:helper_num # Construct cav for helper satellites (1 to helper_num) linking to the target satellite (index N)
                cav[(j, helper_num + 1)] = Dict(:B => 100.0, :Pin => 1e4)
            end
            min_range = 0.0 # min range limit for laser forces in meters
            max_range = 200e3 # max range limit for laser forces in meters

            println(@sprintf("\nCase: target_h=%.0f km, target_i=%.1f deg | helpers=%d", alt_km, inc_deg, helper_num))

            ### 3. Run the simulation ###
            img_dir = normpath(joinpath(@__DIR__, "output", "images"))
            mkpath(img_dir)
            elapsed = @elapsed begin # tic
            sol, p, logs = run_open_cavity_multi(oe;
                                                mass_kg=227,
                                                Pm=Pm,
                                                cavity=cav,
                                                use_los=true,
                                                min_range=min_range, max_range=max_range, # unit: meters
                                                stop_on_dv=false,
                                                T_seconds=2*3600, # default to automatic based on orbits
                                                verbose=false, result_plots=false, target_only=true,
                                                IMG_DIR=img_dir,
                                                helper_num = helper_num,
                                                use_J2=true, useDrag=false
                                                )
            end
            println(@sprintf("Simulation runtime: %.2f s (%.2f min)", elapsed, elapsed/60)) # toc

            ### 4. Save results to CSV ###
            save_timeseries_csv(sol, p, helper_oe, target_orbit, csv_dir=case_csv_dir,
                                include_B=true, include_Pin=true, include_Range=true, include_J2=true)
        end
    end
end