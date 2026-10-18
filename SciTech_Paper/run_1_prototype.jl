using DelimitedFiles
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using Printf
using StaticArrays
using Statistics
using TOML

module KuangPrototype
using LinearAlgebra
using StaticArrays
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using Printf
using Statistics
using DelimitedFiles

const MU = 3.98600436233e14
const C = 299_792_458.0
const R_EARTH = 6_378_136.6
const R_ATMDEF = R_EARTH + 100_000.0
const Ẑ = SVector(0.0, 0.0, 1.0)
idx(i, off) = 6 * (i - 1) + off

const PROTOTYPE_ROOT = normpath(joinpath(@__DIR__, "..", "1_Kuang's Prototype Code"))
const FUNCTIONS = joinpath(PROTOTYPE_ROOT, "functions")
include(joinpath(FUNCTIONS, "1_LOS_Metrics.jl"))
include(joinpath(FUNCTIONS, "2_Laser_Forces_ver2.jl"))
include(joinpath(FUNCTIONS, "3_Dynamics.jl"))
include(joinpath(FUNCTIONS, "4_Diagnostics.jl"))
include(joinpath(FUNCTIONS, "5_OE_Converters.jl"))
include(joinpath(FUNCTIONS, "6_OE_and_dv_in_RTN.jl"))
include(joinpath(FUNCTIONS, "7_Plots.jl"))
include(joinpath(FUNCTIONS, "8_LoS_time_series.jl"))
include(joinpath(FUNCTIONS, "9_Runners.jl"))
include(joinpath(FUNCTIONS, "12_CSV_Write_Read.jl"))
end

const PAPER_DIR = @__DIR__
length(ARGS) in (0, 2, 4) ||
    error("Usage: run_1_prototype.jl [helper_count orbit_count [schedule target_inclination_deg]]")
const N_HELPERS = isempty(ARGS) ? 20 : parse(Int, ARGS[1])
const ORBITS = isempty(ARGS) ? 5 : parse(Int, ARGS[2])
const SCHEDULE = length(ARGS) == 4 ? Symbol(ARGS[3]) : :gve_sma
const TARGET_INC_DEG = length(ARGS) == 4 ? parse(Float64, ARGS[4]) : 0.0
N_HELPERS > 0 && ORBITS > 0 || error("Helper and orbit counts must be positive.")
SCHEDULE in (:gve_sma, :gve_inc) || error("Supported comparison schedules: gve_sma, gve_inc.")
isfinite(TARGET_INC_DEG) && 0 <= TARGET_INC_DEG <= 180 ||
    error("Target inclination must be finite and between 0 and 180 degrees.")
const CASE_DIR = isempty(ARGS) ? "" : "$(N_HELPERS)_helpers_$(ORBITS)_orbits" *
    (length(ARGS) == 4 ? "_$(SCHEDULE)_target_inc_$(TARGET_INC_DEG)_deg" : "")
const MASS_KG = 227.0
const POWER_W = 10_000.0
const MAGNIFICATION = 100.0
const RANGE_M = 200e3
const TARGET_INDEX = N_HELPERS + 1

function target_sma(u)
    position = SVector{3, Float64}(u[KuangPrototype.idx(TARGET_INDEX, 1):KuangPrototype.idx(TARGET_INDEX, 3)])
    velocity = SVector{3, Float64}(u[KuangPrototype.idx(TARGET_INDEX, 4):KuangPrototype.idx(TARGET_INDEX, 6)])
    return inv(2 / norm(position) - dot(velocity, velocity) / KuangPrototype.MU)
end

function target_inclination_deg(u)
    position = SVector{3, Float64}(u[KuangPrototype.idx(TARGET_INDEX, 1):KuangPrototype.idx(TARGET_INDEX, 3)])
    velocity = SVector{3, Float64}(u[KuangPrototype.idx(TARGET_INDEX, 4):KuangPrototype.idx(TARGET_INDEX, 6)])
    h = cross(position, velocity)
    return rad2deg(atan(hypot(h[1], h[2]), h[3]))
end

function selected_helper(u, cavity_map, target_radius)
    positions = Matrix{Float64}(undef, 3, TARGET_INDEX)
    for satellite in 1:TARGET_INDEX
        positions[:, satellite] .= u[KuangPrototype.idx(satellite, 1):KuangPrototype.idx(satellite, 3)]
    end
    selected = KuangPrototype._gve_select_best_cavity(cavity_map, TARGET_INDEX,
        u, positions, SCHEDULE, true, KuangPrototype.R_ATMDEF, 0.0,
        0.0, RANGE_M, KuangPrototype.MU)
    return isempty(selected) ? 0 : only(keys(selected))[1]
end

function run_case(; output_root=PAPER_DIR)
    helper_radius = KuangPrototype.R_EARTH + 1000e3
    target_radius = KuangPrototype.R_EARTH + 1050e3
    target_period = 2pi * sqrt(target_radius^3 / KuangPrototype.MU)
    duration = ORBITS * target_period
    helper_orbits = [
        (a_m=helper_radius, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0,
            ν_deg=360.0 * (helper - 1) / N_HELPERS)
        for helper in 1:N_HELPERS]
    target_orbit = (a_m=target_radius, e=0.0, i_deg=TARGET_INC_DEG, Ω_deg=0.0,
        ω_deg=0.0, ν_deg=0.0)
    orbit_list = vcat(helper_orbits, [target_orbit])
    cavity_map = Dict{Tuple{Int, Int}, Dict{Symbol, Any}}(
        (helper, TARGET_INDEX) => Dict(:B => MAGNIFICATION, :Pin => POWER_W)
        for helper in 1:N_HELPERS)
    plot_cache = joinpath(output_root, "runs", CASE_DIR, "prototype", "plot_cache")
    mkpath(plot_cache)
    elapsed = @elapsed solution, parameters, _, _ = KuangPrototype.run_open_cavity_multi(
        orbit_list;
        mass_kg=MASS_KG,
        cavity=cavity_map,
        use_los=true,
        R_atm=KuangPrototype.R_ATMDEF,
        atm_clearance=0.0,
        min_range=0.0,
        max_range=RANGE_M,
        stop_on_dv=false,
        T_seconds=duration,
        verbose=false,
        result_plots=false,
        target_only=true,
        IMG_DIR=plot_cache,
        helper_num=N_HELPERS,
        use_J2=false,
        useDrag=false,
        gve_schedule=SCHEDULE,
        gve_target_idx=TARGET_INDEX)
    string(solution.retcode) == "Success" || error("Simulation failed: $(solution.retcode)")
    solution.t[end] == duration || error("Simulation ended before the requested duration.")
    initial_inc = target_inclination_deg(first(solution.u))
    final_inc = target_inclination_deg(last(solution.u))
    isapprox(initial_inc, TARGET_INC_DEG; atol=1e-10) ||
        error("Initial target inclination does not match the requested value.")

    sampled_helpers = [selected_helper(u, cavity_map, target_radius) for u in solution.u]
    changes = count(index -> sampled_helpers[index] != sampled_helpers[index - 1],
        2:length(sampled_helpers))
    active_samples = count(>(0), sampled_helpers)
    final_elements = KuangPrototype.rv2coe(
        SVector{3, Float64}(solution.u[end][KuangPrototype.idx(TARGET_INDEX, 1):KuangPrototype.idx(TARGET_INDEX, 3)]),
        SVector{3, Float64}(solution.u[end][KuangPrototype.idx(TARGET_INDEX, 4):KuangPrototype.idx(TARGET_INDEX, 6)]),
        KuangPrototype.MU)
    result = (
        code="1_Kuang's Prototype Code",
        helper_count=N_HELPERS,
        schedule=string(SCHEDULE),
        initial_target_altitude_m=1050e3,
        initial_helper_altitude_m=1000e3,
        initial_helper_inclination_deg=0.0,
        initial_inclination_deg=initial_inc,
        final_inclination_deg=final_inc,
        delta_inclination_deg=final_inc - initial_inc,
        requested_target_orbits=ORBITS,
        duration_s=duration,
        target_period_s=target_period,
        initial_sma_m=target_radius,
        final_sma_m=final_elements.a,
        delta_sma_m=final_elements.a - target_radius,
        final_eccentricity=final_elements.e,
        selection_samples=length(solution.t),
        link_activations=count(index -> sampled_helpers[index - 1] == 0 &&
            sampled_helpers[index] > 0, 2:length(sampled_helpers)),
        sampled_selection_changes=changes,
        active_updates=active_samples,
        final_helper=last(sampled_helpers),
        runtime_s=elapsed,
        retcode=string(solution.retcode),
    )

    result_dir = joinpath(output_root, "results", CASE_DIR, "prototype")
    mkpath(result_dir)
    open(joinpath(result_dir, "summary.toml"), "w") do io
        TOML.print(io, Dict{String, Any}(string(k) => v for (k, v) in pairs(result)))
    end
    open(joinpath(result_dir, "schedule.csv"), "w") do io
        println(io, "time_s,target_orbits,active_helper")
        for (time, helper) in zip(solution.t, sampled_helpers)
            println(io, time, ",", time / target_period, ",", helper)
        end
    end
    println(result)
    return result
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_case()
