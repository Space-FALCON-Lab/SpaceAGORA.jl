# Renders the visual record of the aerobraking reconstruction study
# (Figures 1-5 of docs/spaceagora_aerobraking_reconstruction_record.md) from
# the CSVs produced by run_validation.jl and density_compare.jl.
#
#   Figure 1  Odyssey accelerometer-derived / MarsGRAM density ratio
#   Figure 2  Odyssey apoapsis altitude (truth, Tolson, MarsGRAM comparator)
#   Figure 3  Odyssey periapsis altitude (truth, Tolson, MarsGRAM comparator)
#   Figure 4  VEX apoapsis altitude (digitized truth, VenusGRAM reconstruction)
#   Figure 5  VEX periapsis walk-in (digitized truth, VenusGRAM reconstruction)
#
# Figures are only produced for runs whose outputs exist under the results
# root; missing runs are reported and skipped.
#
# Usage:
#   julia --project=. benchmarks/studies/telemetry_validation/plot_results.jl \
#       [--out-root=DIR]

include(joinpath(@__DIR__, "common.jl"))

using CSV
using DataFrames
using Plots
using Printf

const PLOT_OPTS = parse_kv_args(copy(ARGS), ())

function load_errors(root::String, run_name::String)
    path = joinpath(root, run_name, "telemetry_orbit_accuracy_errors.csv")
    isfile(path) || return nothing
    return CSV.read(path, DataFrame)
end

event_slice(df::DataFrame, event::String) = sort(df[df.event .== event, :], :telemetry_axis)

function plot_apsis_overlay(
    truth_source::DataFrame,
    event::String,
    overlays::Vector{<:Tuple},
    title_str::String,
    ylabel_str::String,
    out_path::String
)
    truth = event_slice(truth_source, event)
    plt = plot(
        truth.telemetry_axis,
        truth.telemetry_value_km;
        label="reconstructed trajectory",
        color=:navy,
        linewidth=2,
        xlabel="campaign orbit",
        ylabel=ylabel_str,
        title=title_str,
        legend=:best,
        dpi=200
    )
    for (df, label, style, color) in overlays
        sliced = event_slice(df, event)
        plot!(plt, sliced.telemetry_axis, sliced.sim_interp_value_km;
              label=label, linestyle=style, color=color, linewidth=1.5)
    end
    savefig(plt, out_path)
    println("  saved $out_path")
    return nothing
end

function main()
    root = results_root(get(PLOT_OPTS, "out-root", ""))
    fig_dir = joinpath(root, "figures")
    mkpath(fig_dir)
    println("Rendering record figures from $root")

    tolson = load_errors(root, "odyssey_tolson")
    marsgram = load_errors(root, "odyssey_marsgram")
    vex = load_errors(root, "vex_venusgram")
    density_csv = joinpath(root, "density_compare", "odyssey_accel_vs_marsgram_density.csv")

    # Figure 1: direct density comparison (requires density_compare.jl output).
    if isfile(density_csv)
        ratio_df = CSV.read(density_csv, DataFrame)
        per_pass = combine(groupby(ratio_df, :P), :density_ratio => (r -> exp(sum(log, r) / length(r))) => :ratio)
        sort!(per_pass, :P)
        plt = plot(
            per_pass.P,
            per_pass.ratio;
            label="density ratio",
            color=:purple,
            linewidth=1.2,
            yscale=:log2,
            yticks=([0.25, 0.5, 1.0, 2.0, 4.0], ["0.25", "0.5", "1", "2", "4"]),
            xlabel="campaign orbit",
            ylabel="accelerometer-derived / MarsGRAM density",
            title="Odyssey density: accelerometer / MarsGRAM",
            legend=:topright,
            dpi=200
        )
        hline!(plt, [1.0]; label="parity", color=:gray, linestyle=:dash)
        out1 = joinpath(fig_dir, "figure1_odyssey_density_ratio.png")
        savefig(plt, out1)
        println("  saved $out1")
    else
        println("  figure 1 skipped: $density_csv not found (run density_compare.jl)")
    end

    # Figures 2-3: Odyssey apsis overlays.
    if tolson !== nothing
        overlays = Tuple[(tolson, "Tolson-density reconstruction", :dash, :purple)]
        marsgram !== nothing && pushfirst!(overlays, (marsgram, "MarsGRAM comparator", :dot, :darkorange))
        plot_apsis_overlay(
            tolson, "apo", overlays,
            "Mars Odyssey apoapsis altitude", "apoapsis altitude [km]",
            joinpath(fig_dir, "figure2_odyssey_apoapsis.png")
        )
        plot_apsis_overlay(
            tolson, "peri", overlays,
            "Mars Odyssey periapsis altitude", "periapsis altitude [km]",
            joinpath(fig_dir, "figure3_odyssey_periapsis.png")
        )
    else
        println("  figures 2-3 skipped: odyssey_tolson errors CSV not found (run run_validation.jl)")
    end

    # Figures 4-5: VEX apsis overlays.
    if vex !== nothing
        overlays = Tuple[(vex, "VenusGRAM reconstruction", :dash, :darkorange)]
        plot_apsis_overlay(
            vex, "apo", overlays,
            "Venus Express apoapsis altitude (digitized truth)", "apoapsis altitude [km]",
            joinpath(fig_dir, "figure4_vex_apoapsis.png")
        )
        plot_apsis_overlay(
            vex, "peri", overlays,
            "Venus Express periapsis walk-in (digitized truth)", "periapsis altitude [km]",
            joinpath(fig_dir, "figure5_vex_periapsis.png")
        )
    else
        println("  figures 4-5 skipped: vex_venusgram errors CSV not found (run run_validation.jl)")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
