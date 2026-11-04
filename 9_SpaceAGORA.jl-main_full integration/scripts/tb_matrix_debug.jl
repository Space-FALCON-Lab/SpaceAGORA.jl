using DataFrames
using Statistics

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "tb_matrix_debug_defs_only.jl"))

function _combined_xyz_rmse(summary::DataFrame, scenario::String)::Union{Float64, Nothing}
    rows = summary[(summary.scenario .== scenario) .& in.(summary.event, Ref(["state_x_time", "state_y_time", "state_z_time"])), :]
    if nrow(rows) != 3
        return nothing
    end
    return sqrt(sum(Float64.(rows.rmse_km) .^ 2))
end

# The legacy GMAT report/CSV label maps to Basilisk_Examples_Full references
# using reference_target=:basilisk; see docs/quality/scenario_matrix_debug.md.
function run_tb_matrix_debug(;
    basilisk_runner=ScenarioMatrixDebugSupport._run_basilisk_scenario_matrix_result_once,
    stk_runner=ScenarioMatrixDebugSupport._run_stk_scenario_matrix_result_once,
    check_inputs=() -> ScenarioMatrixDebugSupport.require_matrix_inputs(; stk=true),
    csv_path::AbstractString=joinpath(REPO_ROOT, "scripts", "tb_matrix_rmse_m.csv")
)
    # Check both reference sets before either potentially expensive solve.
    check_inputs()
    bodies = ["earth", "mars", "venus", "moon"]
    gravity_tags = ["j0", "j2", "j50"]
    tb_tags = ["tbfalse", "tbtrue"]

    println("Running GMAT-target scenario matrix (this may take a while)...")
    gmat_result = basilisk_runner()
    gmat_summary = gmat_result.summary
    println("GMAT-target matrix done. Summary rows: ", nrow(gmat_summary))

    println("Running STK-target scenario matrix (this may take a while)...")
    stk_result = stk_runner()
    stk_summary = stk_result.summary
    println("STK-target matrix done. Summary rows: ", nrow(stk_summary))

    missing_pairs = Tuple{String, String}[]  # (scenario, target)
    values_m = Dict{Tuple{String, String}, Float64}()  # (scenario, target) => meters

    for body in bodies, gtag in gravity_tags, tbtag in tb_tags
        scenario = "$(body)_$(gtag)_$(tbtag)"
        gmat_rmse_km = _combined_xyz_rmse(gmat_summary, scenario)
        if gmat_rmse_km === nothing
            push!(missing_pairs, (scenario, "GMAT"))
        else
            values_m[(scenario, "GMAT")] = gmat_rmse_km * 1000.0
        end

        stk_rmse_km = _combined_xyz_rmse(stk_summary, scenario)
        if stk_rmse_km === nothing
            push!(missing_pairs, (scenario, "STK"))
        else
            values_m[(scenario, "STK")] = stk_rmse_km * 1000.0
        end
    end

    fmt(x) = x === nothing ? "MISSING" : string(round(x, digits=3))

    # Machine-readable CSV, in meters, full precision, for downstream reuse
    # (paper table + any future re-derivation) without rerunning the matrix.
    open(csv_path, "w") do io
        println(io, "planet,gravity,target,third_body,rmse_m")
        for body in bodies, gtag in gravity_tags, tbtag in tb_tags, target in ["GMAT", "STK"]
            scenario = "$(body)_$(gtag)_$(tbtag)"
            v = get(values_m, (scenario, target), nothing)
            println(io, body, ",", gtag, ",", target, ",", tbtag, ",", v === nothing ? "" : v)
        end
    end
    println("\nWrote machine-readable results to: ", csv_path)

    println()
    println("=== Combined-axis position RMSE (meters) ===")
    println(rpad("Planet", 8), " | ", rpad("Gravity", 8), " | ", rpad("GMAT tbfalse", 14), " | ", rpad("GMAT tbtrue", 14), " | ", rpad("STK tbfalse", 14), " | ", rpad("STK tbtrue", 14))
    println(repeat("-", 90))

    gravity_label = Dict("j0" => "L=0", "j2" => "J2", "j50" => "L=50")

    for body in bodies, gtag in gravity_tags
        scen_false = "$(body)_$(gtag)_tbfalse"
        scen_true = "$(body)_$(gtag)_tbtrue"
        gmat_f = get(values_m, (scen_false, "GMAT"), nothing)
        gmat_t = get(values_m, (scen_true, "GMAT"), nothing)
        stk_f = get(values_m, (scen_false, "STK"), nothing)
        stk_t = get(values_m, (scen_true, "STK"), nothing)
        println(
            rpad(body, 8), " | ", rpad(gravity_label[gtag], 8), " | ",
            rpad(fmt(gmat_f), 14), " | ", rpad(fmt(gmat_t), 14), " | ",
            rpad(fmt(stk_f), 14), " | ", rpad(fmt(stk_t), 14)
        )
    end

    println()
    if isempty(missing_pairs)
        println("No missing scenario/target pairs.")
    else
        println("MISSING scenario/target pairs:")
        for (scenario, target) in missing_pairs
            println("  ", scenario, " (", target, ")")
        end
    end

    println()
    println("Potential no-op third-body cases (tbfalse == tbtrue exactly):")
    any_dup = false
    for body in bodies, gtag in gravity_tags, target in ["GMAT", "STK"]
        scen_false = "$(body)_$(gtag)_tbfalse"
        scen_true = "$(body)_$(gtag)_tbtrue"
        vf = get(values_m, (scen_false, target), nothing)
        vt = get(values_m, (scen_true, target), nothing)
        if vf !== nothing && vt !== nothing && vf == vt
            println("  ", body, "_", gtag, " (", target, "): both = ", vf, " m")
            any_dup = true
        end
    end
    if !any_dup
        println("  none")
    end

    println()
    println("=== SCRIPT DONE ===")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_tb_matrix_debug()
end
