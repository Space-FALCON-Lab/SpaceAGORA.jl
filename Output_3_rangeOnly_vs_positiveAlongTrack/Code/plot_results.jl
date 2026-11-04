if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "9_SpaceAGORA.jl-main_full integration")))
end

module ComparisonPlots

using CSV, DataFrames, Printf, SHA, TOML
get!(ENV, "GKSwstype", "100")
using Plots

include("layout.jl")
using .ComparisonLayout

const POLICIES = ("naive_next_entering", "positive_along_track")
const LABELS = ("Range-only (next entering)", "Positive along-track")
const COLORS = ("#0072B2", "#D55E00")
const ELEMENTS = (
    (:a_km, "Semimajor axis [km]"),
    (:eccentricity, "Eccentricity"),
    (:inclination_deg, "Inclination [deg]"),
    (:raan_deg, "RAAN [deg, unwrapped]"),
)

filehash(path) = open(io -> bytes2hex(sha256(io)), path)

function find_case(root, policy)
    paper_root = joinpath(root, "paper_plot_mode")
    isdir(paper_root) || error("Missing paper_plot_mode directory: $paper_root")
    candidates = String[]
    for (directory, _, files) in walkdir(paper_root)
        "run_info.toml" in files || continue
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        info["policy"] == policy && push!(candidates, directory)
    end
    length(candidates) == 1 ||
        error("Expected exactly one $policy case under $paper_root; found $(length(candidates))")
    only(candidates)
end

function read_result(root, policy)
    directory = find_case(root, policy)
    info = TOML.parsefile(joinpath(directory, "run_info.toml"))
    info["status"] == "verified" || error("Run is not verified: $directory")
    path = joinpath(directory, "target_orbital_elements.csv")
    filehash(path) == info["checksums"]["target_orbital_elements.csv"] ||
        error("Target CSV checksum mismatch: $path")
    table = CSV.read(path, DataFrame)
    all(name -> all(isfinite, table[!, name]), names(table)) ||
        error("Nonfinite target data: $path")
    nrow(table) == info["trajectory_rows"] || error("Incorrect target row count")
    first(table.time_s) == 0 && last(table.time_s) == info["settings"]["duration_s"] ||
        error("Incomplete target time span")
    first(table.target_orbits) == 0 &&
        isapprox(last(table.target_orbits), info["settings"]["orbits"]; atol=1e-12, rtol=0) ||
        error("Incorrect orbit-count axis")
    table.delta_a_km == table.a_km .- first(table.a_km) ||
        error("SMA changes must equal a(t) - a(0): $path")
    (; table, info)
end

function plot_results(root=DEFAULT_OUTPUT)
    results = [read_result(root, policy) for policy in POLICIES]
    settings = first(results).info["settings"]
    all(result -> result.info["settings"] == settings, results) ||
        error("The two runs do not share identical physical/numerical settings")
    results[1].table.time_s == results[2].table.time_s ||
        error("The two runs do not share a sampling grid")
    output = joinpath(root, "images", "orbital_elements")
    mkpath(output)
    title = @sprintf("Target %.0f km | %d helpers at %.0f km | relative inclination %.1f deg | J2",
        settings["target_altitude_km"], settings["helpers"], settings["helper_altitude_km"],
        abs(settings["target_inclination_deg"] - settings["helper_inclination_deg"]))
    files = String[]
    for (axis, xlabel, suffix) in (
        (:time_days, "Elapsed time [days]", "time"),
        (:target_orbits, "Initial target orbital periods (t / T0)", "orbits"),
    )
        for (elements, layout, dimensions, prefix) in (
            (ELEMENTS, (2, 2), (1500, 950), "orbital_elements"),
            ((first(ELEMENTS),), (1, 1), (1200, 750), "sma"),
            (((:delta_a_km, "Change in semimajor axis [km]"),),
                (1, 1), (1200, 750), "delta_sma"),
        )
            panels = map(elements) do (column, ylabel)
                panel = plot(; xlabel, ylabel, legend=:best, grid=true,
                    framestyle=:box, tickfontsize=10, guidefontsize=12, legendfontsize=10)
                for index in eachindex(results)
                    table = results[index].table
                    plot!(panel, table[!, axis], table[!, column]; label=LABELS[index],
                        color=COLORS[index],
                        linestyle=(length(elements) == 1 || index == 1) ? :solid : :dash,
                        linewidth=1.35)
                end
                panel
            end
            figure = plot(panels...; layout, size=dimensions, dpi=180,
                plot_title=title, plot_titlefontsize=13, margin=6Plots.mm)
            for extension in ("png", "pdf")
                name = "$(prefix)_vs_$suffix.$extension"
                savefig(figure, joinpath(output, name))
                filesize(joinpath(output, name)) > 1000 || error("Empty plot: $name")
                push!(files, name)
            end
        end
    end
    summary = DataFrame([merge((policy=POLICIES[index],),
        NamedTuple(results[index].table[end, :])) for index in eachindex(results)])
    CSV.write(joinpath(root, "comparison_summary.csv"), summary)
    open(joinpath(output, "plot_info.toml"), "w") do io
        TOML.print(io, Dict("input_status" => "verified", "rows_per_policy" => nrow(results[1].table),
            "plot_source_sha256" => filehash(@__FILE__),
            "figure_groups" => ["four_panel", "sma_only", "sma_change"],
            "sma_change_definition" => "a(t) - a(0) in km, using each policy's first saved SMA",
            "raan" => "unwrapped across every saved sample",
            "elements" => "osculating, not orbit-averaged",
            "files_sha256" => Dict(name => filehash(joinpath(output, name)) for name in files)))
    end
    println("Saved four-panel, SMA, and SMA-change comparisons versus time and orbits (PNG and PDF) in $output")
    output
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: julia plot_results.jl [output_directory]")
    ComparisonPlots.plot_results(isempty(ARGS) ? ComparisonPlots.DEFAULT_OUTPUT : abspath(only(ARGS)))
end
