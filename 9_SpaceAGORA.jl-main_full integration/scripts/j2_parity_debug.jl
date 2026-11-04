using DataFrames
using Statistics

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

# Optional controls:
# - SPACEAGORA_DEBUG_COMPARE_J2=1
# - SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT=file_c20|planet_j2
# - SPACEAGORA_TELEMETRY_J2_SOURCE_PLANET_SCENARIOS=name1,name2

include(joinpath(@__DIR__, "tb_matrix_debug_defs_only.jl"))

function _j2_combined_xyz_rmse(summary::DataFrame, scenario::String)::Float64
    rows = summary[(summary.scenario .== scenario) .& in.(summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
    return sqrt(sum(Float64.(rows.rmse_km) .^ 2))
end

# See docs/quality/scenario_matrix_debug.md for references and invocation.
function run_j2_parity_debug(;
    runner=ScenarioMatrixDebugSupport._run_basilisk_scenario_matrix_result_once,
    check_inputs=ScenarioMatrixDebugSupport.require_matrix_inputs
)
    check_inputs()
    result = runner()
    summary = result.summary

    j2_rows = summary[in.(summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])) .& occursin.("_j2_", String.(summary.scenario)), :]
    scenarios = sort(unique(String.(j2_rows.scenario)))

    println("J2 scenario parity summary")
    println("  J2 source default: ", get(ENV, "SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT", "file_c20"))
    println("  Planet override list: ", get(ENV, "SPACEAGORA_TELEMETRY_J2_SOURCE_PLANET_SCENARIOS", ""))
    println("  Compare J2 analytic/generic: ", get(ENV, "SPACEAGORA_DEBUG_COMPARE_J2", "0"))
    println()

    table = DataFrame(
        scenario=String[],
        rmse_x_km=Float64[],
        rmse_y_km=Float64[],
        rmse_z_km=Float64[],
        rmse_xyz_norm_km=Float64[]
    )

    for scenario in scenarios
        rows = summary[(summary.scenario .== scenario) .& in.(summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
        x = rows[rows.event .== "state_x_time", :].rmse_km[1]
        y = rows[rows.event .== "state_y_time", :].rmse_km[1]
        z = rows[rows.event .== "state_z_time", :].rmse_km[1]
        push!(table, (scenario, x, y, z, _j2_combined_xyz_rmse(summary, scenario)))
    end

    sort!(table, :rmse_xyz_norm_km, rev=true)
    show(stdout, MIME("text/plain"), table)
    println()

    println("\nWorst J2 cases:")
    nshow = min(5, nrow(table))
    for i in 1:nshow
        r = table[i, :]
        println("  ", r.scenario, ": xyz_norm_rmse=", r.rmse_xyz_norm_km, " km")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_j2_parity_debug()
end
