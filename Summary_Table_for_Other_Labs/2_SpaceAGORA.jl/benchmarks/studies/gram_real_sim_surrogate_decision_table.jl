const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using CSV
using DataFrames

const PLANETS = ("mars", "venus", "earth", "titan")
const SCENARIOS = ("drag_passage", "entry", "orbit")

const THRESH_RHO_REL_P95 = Dict(
    "drag_passage" => 0.50,
    "entry" => 0.30,
    "orbit" => 0.05
)
const THRESH_TEMP_ABS_MAX_K = Dict(
    "drag_passage" => 8.0,
    "entry" => 7.0,
    "orbit" => 15.0
)
const THRESH_WIND_ABS_P95_MS = Dict(
    "drag_passage" => 120.0,
    "entry" => 100.0,
    "orbit" => 200.0
)
const THRESH_RUNTIME_SPEEDUP = 0.90

function _summary_path(planet::String)
    return joinpath(REPO_ROOT, "output", "gram_real_sim_surrogate_matrix_summary_$(planet).csv")
end

function _runtime_path(planet::String)
    return joinpath(REPO_ROOT, "output", "gram_real_sim_surrogate_matrix_runtime_$(planet).csv")
end

function _load_summary(planet::String)
    path = _summary_path(planet)
    isfile(path) || return DataFrame()
    df = CSV.read(path, DataFrame)
    df[!, :planet] .= planet
    return df
end

function _load_runtime(planet::String)
    path = _runtime_path(planet)
    isfile(path) || return DataFrame()
    df = CSV.read(path, DataFrame)
    df[!, :planet] .= planet
    return df
end

function build_table()
    summary = vcat([_load_summary(p) for p in PLANETS]...; cols=:union)
    runtime = vcat([_load_runtime(p) for p in PLANETS]...; cols=:union)

    rows = NamedTuple[]
    for planet in PLANETS
        runtime_p = runtime[runtime.planet .== planet, :]
        summary_p = summary[summary.planet .== planet, :]
        for scenario in SCENARIOS
            idx = findfirst((summary_p.scenario .== scenario))
            point_ok = any((runtime_p.scenario .== scenario) .& (runtime_p.mode .== "point_to_point") .& (runtime_p.retcode .== "Success"))
            surr_ok = any((runtime_p.scenario .== scenario) .& (runtime_p.mode .== "offline_surrogate") .& (runtime_p.retcode .== "Success"))
            if idx === nothing
                status = if point_ok && surr_ok
                    "MISSING_SUMMARY"
                elseif xor(point_ok, surr_ok)
                    "PARTIAL"
                else
                    "NO_DATA"
                end
                push!(rows, (
                    planet=planet,
                    scenario=scenario,
                    status=status,
                    rho_rel_p95=NaN,
                    temp_abs_max=NaN,
                    wind_abs_p95=NaN,
                    runtime_speedup=NaN,
                    rho_rel_p95_threshold=THRESH_RHO_REL_P95[scenario],
                    temp_abs_max_threshold=THRESH_TEMP_ABS_MAX_K[scenario],
                    wind_abs_p95_threshold=THRESH_WIND_ABS_P95_MS[scenario],
                    runtime_speedup_threshold=THRESH_RUNTIME_SPEEDUP,
                    pass_rho=false,
                    pass_temp=false,
                    pass_wind=false,
                    pass_runtime=false,
                    pass_all=false
                ))
            else
                r = summary_p[idx, :]
                pass_rho = r.rho_rel_p95 <= THRESH_RHO_REL_P95[scenario]
                pass_temp = r.temp_abs_max <= THRESH_TEMP_ABS_MAX_K[scenario]
                pass_wind = r.wind_abs_p95 <= THRESH_WIND_ABS_P95_MS[scenario]
                pass_runtime = r.runtime_speedup >= THRESH_RUNTIME_SPEEDUP
                pass_all = point_ok && surr_ok && pass_rho && pass_temp && pass_wind && pass_runtime
                push!(rows, (
                    planet=planet,
                    scenario=scenario,
                    status=pass_all ? "PASS" : "FAIL",
                    rho_rel_p95=r.rho_rel_p95,
                    temp_abs_max=r.temp_abs_max,
                    wind_abs_p95=r.wind_abs_p95,
                    runtime_speedup=r.runtime_speedup,
                    rho_rel_p95_threshold=THRESH_RHO_REL_P95[scenario],
                    temp_abs_max_threshold=THRESH_TEMP_ABS_MAX_K[scenario],
                    wind_abs_p95_threshold=THRESH_WIND_ABS_P95_MS[scenario],
                    runtime_speedup_threshold=THRESH_RUNTIME_SPEEDUP,
                    pass_rho=pass_rho,
                    pass_temp=pass_temp,
                    pass_wind=pass_wind,
                    pass_runtime=pass_runtime,
                    pass_all=pass_all
                ))
            end
        end
    end

    df = DataFrame(rows)
    sort!(df, [:planet, :scenario])
    return df
end

function run()
    table = build_table()
    out = joinpath(REPO_ROOT, "output", "gram_real_sim_surrogate_decision_table.csv")
    CSV.write(out, table)
    println("Saved decision table: $out")
    show(table, allrows=true, allcols=true)
    println()
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end
