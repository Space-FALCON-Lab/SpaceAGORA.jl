#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module AbsoluteRNPlots

include("plot_pooled_orbital_elements.jl")
include("plot_paired_dv_RTN.jl")
using .PooledOrbitalElementPlots
using CSV, Plots

const OE = PooledOrbitalElementPlots.OE
const OUTPUT = joinpath(OE.StudyConfig.DEFAULT_OUTPUT, "images")
const IDENTIFIERS = (:scheduler, :target_altitude_km, :delta_i_deg, :helpers,
    :orbits, :time_s, :row, :source_trajectory)
const RTN_FAMILIES = (
    (quantity=:rtn, folder="", csv="cumulative_dv_RTN"),
    (quantity=:absolute, folder="1_absolute_RTN", csv="cumulative_absolute_dv_RTN"),
    (quantity=:split_n, folder="2_signed_RT_positive_negative_N", csv="cumulative_dv_RT_N_branches"),
    (quantity=:branches, folder="3_positive_negative_RTN", csv="cumulative_dv_RTN_branches"),
)
const SELECTED_ALTITUDES = (1000, 1050, 1100)
const SELECTED_FAMILIES = (
    (quantity=:oe, folder="orbital_elements/4_combined_1000_1050_1100km", csv="orbital_element_deltas"),
    (quantity=:absolute, folder="delta_v_RTN/5_absolute_RTN_1000_1050_1100km", csv="cumulative_absolute_dv_RTN"),
    (quantity=:split_n, folder="delta_v_RTN/6_signed_RT_positive_negative_N_1000_1050_1100km", csv="cumulative_dv_RT_N_branches"),
    (quantity=:branches, folder="delta_v_RTN/7_positive_negative_RTN_1000_1050_1100km", csv="cumulative_dv_RTN_branches"),
)

select_fields(row, fields) = NamedTuple{fields}(Tuple(getproperty(row, field) for field in fields))

function write_endpoints(path, rows, window, fields; altitudes=OE.ALTITUDES)
    selected = filter(row -> row.orbits == window && row.target_altitude_km in altitudes, rows)
    expected = Set((policy, altitude, inclination, helpers) for policy in OE.POLICIES
        for altitude in altitudes for inclination in OE.INCLINATIONS for helpers in OE.COUNTS)
    actual = Set((row.scheduler, row.target_altitude_km, row.delta_i_deg, row.helpers) for row in selected)
    length(selected) == length(expected) && actual == expected ||
        error("Expected the complete $(length(expected))-case grid for $window orbits")
    CSV.write(path, [select_fields(row, (IDENTIFIERS..., fields...)) for row in selected])
end

function render_selected_altitudes(rows, thresholds, windows, output)
    selected_thresholds = filter(t -> t.target_altitude_km in SELECTED_ALTITUDES, thresholds)
    expected = Set((h, i) for h in SELECTED_ALTITUDES for i in OE.INCLINATIONS)
    length(selected_thresholds) == length(expected) &&
        Set((t.target_altitude_km, t.delta_i_deg) for t in selected_thresholds) == expected ||
        error("Expected nine continuous-access thresholds for the selected altitudes")
    gr()
    for family in SELECTED_FAMILIES
        directory = joinpath(output, family.folder)
        mkpath(directory)
        CSV.write(joinpath(directory, "continuous_access_thresholds.csv"), selected_thresholds)
        fields = OE.quantity_spec(family.quantity).fields
        family.quantity == :oe && (fields = (fields..., :raan_defined))
        for window in windows
            write_endpoints(joinpath(directory, "$(family.csv)_$(window)orbits.csv"),
                rows, window, fields; altitudes=SELECTED_ALTITUDES)
            OE.make_plot(rows, selected_thresholds, "combined", window, directory;
                quantity=family.quantity, altitudes=SELECTED_ALTITUDES)
        end
    end
    nothing
end

function render_all(rows, thresholds, windows, output)
    oe_output = joinpath(output, "orbital_elements")
    separate = joinpath(oe_output, "1_seperate_plots_by_sceduler_type")
    pooled = joinpath(oe_output, "2_separate_plots_by_N-sceduler_h-i")
    combined = joinpath(oe_output, "3_combined_plots")
    for directory in (separate, pooled, combined)
        mkpath(directory)
    end
    CSV.write(joinpath(oe_output, "continuous_access_thresholds.csv"), thresholds)
    gr()
    for window in windows
        write_endpoints(joinpath(oe_output, "orbital_element_deltas_$(window)orbits.csv"),
            rows, window, (OE.quantity_spec(:oe).fields..., :raan_defined))
        for policy in OE.POLICIES
            OE.make_plot(rows, thresholds, policy, window, separate)
        end
        OE.make_plot(rows, thresholds, "combined", window, combined)
        for xfield in (:helpers, :target_altitude_km)
            figure = PooledOrbitalElementPlots.make_plot(rows, window, xfield)
            stem = xfield == :helpers ? "OE_vs_helpers_by_scheduler" : "OE_vs_target_altitude_by_inclination"
            for extension in ("png", "pdf")
                path = joinpath(pooled, "$(stem)_$(window)orbits.$extension")
                savefig(figure, path)
                println("Saved $path")
            end
        end
    end
    for family in RTN_FAMILIES
        directory = joinpath(output, "delta_v_RTN", family.folder)
        mkpath(directory)
        CSV.write(joinpath(directory, "continuous_access_thresholds.csv"), thresholds)
        for window in windows
            write_endpoints(joinpath(directory, "$(family.csv)_$(window)orbits.csv"),
                rows, window, OE.quantity_spec(family.quantity).fields)
            OE.make_plot(rows, thresholds, "combined", window, directory; quantity=family.quantity)
        end
    end
    nothing
end

function main(args=ARGS)
    options = OE.parse_options(vcat(["--output", OUTPUT], args))
    rows, thresholds = OE.load_grid(options.input, options.windows; include_rtn=true)
    render_all(rows, thresholds, options.windows, options.output)
    render_selected_altitudes(rows, thresholds, options.windows, options.output)
    PairedRTNPlots.render_all(rows, thresholds, options.windows,
        joinpath(options.output, "delta_v_RTN", "8_paired_RTN_by_scheduler"))
    println("Validated $(length(rows)) endpoints from $(OE.CASE_COUNT) schema-3, checksum-matched trajectories.")
    println("Saved $(15 * length(options.windows)) figures in both PNG and PDF, with endpoint and threshold CSVs.")
    nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    AbsoluteRNPlots.main()
end
