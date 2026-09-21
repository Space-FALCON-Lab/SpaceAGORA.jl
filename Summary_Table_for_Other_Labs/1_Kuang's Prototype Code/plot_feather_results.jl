module FeatherResultPlots

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "2_SpaceAGORA.jl"); io=devnull)
using Arrow, LinearAlgebra, StaticArrays, Printf, Statistics, Plots
import GLMakie

const MU = 3.986004418e14
const C = 3.0e8
const R_EARTH = 6_378_137.0
const R_ATMDEF = R_EARTH + 100_000.0
const Ẑ = SVector(0.0, 0.0, 1.0)
idx(satellite, component) = 6*(satellite-1) + component

include("functions/1_LOS_Metrics.jl")
include("functions/2_Laser_Forces_ver2.jl")
include("functions/3_Dynamics.jl")
include("functions/4_Diagnostics.jl")
include("functions/5_OE_Converters.jl")
include("functions/6_OE_and_dv_in_RTN.jl")
include("functions/7_Plots.jl")

function plot_feather_results(directory; laser_power_w, magnification, schedule,
        use_J2=true, output_dir=joinpath(dirname(dirname(abspath(directory))),
            "images", basename(normpath(directory))))
    table = Arrow.Table(joinpath(directory, "trajectory.feather"))
    metadata = Arrow.getmetadata(table)
    metadata["source"] == "Kuang prototype" || error("Expected prototype Feather data")
    spacecraft_ids = sort!([parse(Int, matched.captures[1]) for column in propertynames(table)
        for matched in (match(r"^sc(\d+)_pos_1$", string(column)),) if matched !== nothing])
    count = length(spacecraft_ids)
    spacecraft_ids == collect(1:count) || error("Expected contiguous spacecraft IDs")
    target = parse(Int, metadata["target_id"])
    times = Float64.(table.time)
    states = [zeros(6*count) for time in times]
    for satellite in spacecraft_ids, (field, offset) in (("pos", 0), ("vel", 3)), component in 1:3
        values = getproperty(table, Symbol("sc$(satellite)_$(field)_$(component)"))
        for sample in eachindex(times)
            states[sample][idx(satellite, offset+component)] = values[sample]
        end
    end
    sol = (; t=times, u=states)
    masses = [Float64(first(getproperty(table, Symbol("sc$(satellite)_mass")))) for satellite in spacecraft_ids]
    helpers = filter(!=(target), spacecraft_ids)
    params = Dict{Symbol,Any}(:mu=>MU, :c=>C, :N=>count, :masses=>masses,
        :Pmatrix=>zeros(count, count), :use_los=>true,
        :cavity=>Dict((helper, target)=>Dict(:B=>magnification, :Pin=>laser_power_w) for helper in helpers),
        :R_atm=>R_ATMDEF, :atm_clearance=>0.0, :min_range=>0.0,
        :max_range=>parse(Float64, metadata["maximum_range_m"]),
        :helper_ids=>helpers, :target_ids=>[target], :use_J2=>use_J2,
        :useDrag=>false, :gve_schedule=>String(schedule), :gve_target_idx=>target)
    _, delta_v = delta_v_RTN_time_series(sol, params)
    saved_delta_v = permutedims(hcat(table.dv_r_accumulated, table.dv_t_accumulated, table.dv_n_accumulated))
    difference = maximum(abs, delta_v[target] - saved_delta_v)
    isapprox(delta_v[target], saved_delta_v; rtol=1e-10, atol=1e-12) ||
        error("Reconstructed laser delta-v disagrees with saved diagnostics: $difference m/s")
    mkpath(output_dir)
    for subdirectory in ("r_RTN_sat", "v_RTN_sat", "a_RTN_sat", "dv_from_laser_in_RTN_for_sat",
            "F_from_laser_in_RTN_for_sat", "delta_P_from_laser_in_RTN_for_sat")
        mkpath(joinpath(output_dir, subdirectory))
    end
    plot_orbits(sol; IMG_DIR=output_dir, fn="satellite_orbits.png")
    plot_angmom_two_axes(sol, masses; IMG_DIR=output_dir, fn="angular_momentum.png")
    plot_momentum_two_axes(sol, masses; IMG_DIR=output_dir, fn="linear_momentum.png")
    plot_orbit_energy(sol, params; IMG_DIR=output_dir, fn="orbital_energy.png", subplot=true)
    plot_orbit_energy_individual_satellites(sol, params; IMG_DIR=output_dir, fn="orbit_energy_individual.png")
    plot_orbit_energy_total(sol, params; IMG_DIR=output_dir, fn="orbital_energy_total.png")
    report_and_plot_r_RTN(sol, params; sat=target, IMG_DIR=output_dir, fn_prefix="r_RTN_sat/r_RTN_sat_$target", R_only=true)
    report_and_plot_v_RTN(sol, params; sat=target, IMG_DIR=output_dir, fn_prefix="v_RTN_sat/v_RTN_sat_$target", R_only=true)
    report_and_plot_a_RTN(sol, params; sat=target, IMG_DIR=output_dir, fn_prefix="a_RTN_sat/a_RTN_sat_$target", show_a=false, show_a_gravity=false, show_a_laser=true, R_only=true)
    report_and_plot_dv_RTN(sol, params; sat=target, IMG_DIR=output_dir, fn_prefix="dv_from_laser_in_RTN_for_sat/dv_from_laser_in_RTN_for_sat_$target", RTN_separate=true)
    report_and_plot_F_RTN(sol, params; sat=target, IMG_DIR=output_dir, fn_prefix="F_from_laser_in_RTN_for_sat/F_from_laser_in_RTN_for_sat_$target", RTN_separate=true)
    report_and_plot_delta_P_RTN(sol, params; sat=target, IMG_DIR=output_dir, fn_prefix="delta_P_from_laser_in_RTN_for_sat/delta_P_from_laser_in_RTN_for_sat_$target", R_only=true)
    report_and_plot_OE(sol, MU; sat=target, IMG_DIR=output_dir, fn_prefix="orbital_elements_sat/orbital_elements_sat_$target")
    report_and_plot_rp_ra(sol, MU, R_EARTH; sat=target, IMG_DIR=output_dir, fn_prefix="apogee_perigee")
    report_and_plot_OE_diff(sol, MU; sat1=1, sat2=2, IMG_DIR=output_dir, fn_prefix="orbital_elements_diff")
    open(joinpath(output_dir, "README.md"), "w") do io
        println(io, "# Plots from saved Feather data\n")
        println(io, "Source: `$(abspath(directory))`. No orbit simulation was rerun.\n")
        println(io, "Reuses the target-only result plots from test16_options.jl; target spacecraft is $target. Orbital-element differences retain the reference's spacecraft 1 versus 2 comparison.\n")
        println(io, "Positions, velocities and masses come from saved samples. Laser force, acceleration, delta-v and impulse are reconstructed with the prototype diagnostic routines: power=$laser_power_w W, magnification=$magnification, schedule=$schedule, J2=$use_J2, drag=false. These parameters are supplied explicitly because the Feather metadata does not contain all of them.\n")
        println(io, "Maximum reconstructed versus stored target RTN delta-v difference: $difference m/s. Sampling is $(metadata["output_interval_s"]) seconds plus the endpoint; plots do not restore unsaved solver steps.")
    end
    println("Saved Feather result plots to $output_dir")
    return output_dir
end

end