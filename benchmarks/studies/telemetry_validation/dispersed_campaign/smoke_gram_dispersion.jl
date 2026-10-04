# Smoke test for a dispersed Odyssey MarsGRAM campaign: do per-sample
# `gram_seed` / `gram_perturbation_scales` overrides reach the density the
# propagator actually uses?
#
# Two checks, three variants of the odyssey_marsgram record manifest
# (manifests/odyssey_marsgram.toml), each with its own GRAM instance:
#   V0  seed 1001, scales (0,0,0,0)   -- the record configuration
#   V1  seed 11,   scales (1,1,1,1)   -- GRAM's documented nominal scales
#   V2  seed 22,   scales (1,1,1,1)
#
# A. Point queries. Each variant's density model is built by the same
#    TelemetryVerification builder the runner uses and walked along one fixed
#    synthetic path. At every point it records the density SpaceAGORA's RHS
#    hook returns (EnvironmentModels._gram_point_density) and, from the same
#    GRAM instance right after that call, GRAM's own perturbedDensity
#    (get_density_state). If the hook density is identical across variants
#    while perturbedDensity differs, the seed and scales are reaching GRAM's
#    random walk but not the trajectory.
# B. One-orbit solves. The record scenario truncated to ONE orbit, run through
#    the runner's own _run_simulation_dataframe for each variant; the saved
#    position/velocity histories are compared.
#
# Usage: julia --project=. benchmarks/studies/telemetry_validation/dispersed_campaign/smoke_gram_dispersion.jl [--out=DIR]

include(joinpath(@__DIR__, "..", "common.jl"))
load_gramsuite!()

using SpaceAGORA
using SpaceAGORA.TelemetryVerification
using DataFrames
using CSV
using TOML
using Printf

const TV = SpaceAGORA.TelemetryVerification
const EM = SpaceAGORA.SimulationModel.EnvironmentModels

const OPTS = parse_kv_args(copy(ARGS))
const OUT = abspath(get(OPTS, "out", joinpath(REPO_ROOT, "results", "odyssey_dispersion_smoke")))
mkpath(OUT)

const RECORD_MANIFEST = joinpath(MANIFEST_DIR, "odyssey_marsgram.toml")
const VARIANTS = [
    (tag="V0_seed1001_scales0", seed=1001, scales=[0.0, 0.0, 0.0, 0.0]),
    (tag="V1_seed11_scales1",   seed=11,   scales=[1.0, 1.0, 1.0, 1.0]),
    (tag="V2_seed22_scales1",   seed=22,   scales=[1.0, 1.0, 1.0, 1.0]),
]

function variant_cfg(seed::Int, scales::Vector{Float64})
    doc = TOML.parsefile(RECORD_MANIFEST)
    sc = doc["scenarios"][1]
    sc["atmosphere_truth"]["gram_seed"] = seed
    sc["atmosphere_truth"]["gram_perturbation_scales"] = scales
    path = joinpath(mktempdir(), "manifest.toml")
    # Relative telemetry paths in the manifest resolve against the repo root.
    cd(REPO_ROOT) do
        open(path, "w") do io
            TOML.print(io, doc)
        end
    end
    cfgs = cd(() -> TV._load_scenarios_from_manifest(path), REPO_ROOT)
    return only(cfgs)
end

# ── A. point queries ─────────────────────────────────────────────────────────
const N_POINTS = 400
path_h(i) = 1000.0 * (100.0 + 60.0 * abs(1.0 - 2.0 * (i - 1) / (N_POINTS - 1)))  # 160 -> 100 -> 160 km
path_lat(i) = deg2rad(-60.0 + 120.0 * (i - 1) / (N_POINTS - 1))
path_lon(i) = deg2rad(mod(10.0 + 0.9 * (i - 1), 360.0))
path_t(i) = 5.0 * (i - 1)

point_rows = NamedTuple[]
# The runner furnishes SPICE kernels when it builds the planet for a solve
# (_make_orbit_args -> _planet_from_name); the point queries below bypass that,
# so do it here first or GRAM's ephemeris lookup has no leapseconds kernel.
TV._planet_from_name("mars")
for v in VARIANTS
    cfg = variant_cfg(v.seed, v.scales)
    @assert cfg.atmosphere_truth.gram_seed == v.seed
    @assert collect(cfg.atmosphere_truth.gram_perturbation_scales) == v.scales
    model = TV._scenario_density_model(cfg)
    core = model.core
    get_density_state = Base.invokelatest(getproperty, core.gram, :get_density_state)
    for i in 1:N_POINTS
        rho, T, w = EM._gram_point_density(model, path_h(i), path_lat(i), path_lon(i), path_t(i), true)
        ds = Base.invokelatest(get_density_state, core.gram_atmosphere)
        push!(point_rows, (variant=v.tag, seed=v.seed, scale=v.scales[1], i=i,
                           h_m=path_h(i), lat_deg=rad2deg(path_lat(i)), lon_deg=rad2deg(path_lon(i)),
                           t_s=path_t(i), rho_rhs_kgm3=rho,
                           gram_mean_density_kgm3=Float64(ds.density),
                           gram_perturbed_density_kgm3=Float64(ds.perturbedDensity),
                           gram_density_sigma_pct=Float64(ds.densityStandardDeviation)))
    end
end
points = DataFrame(point_rows)
CSV.write(joinpath(OUT, "point_queries.csv"), points)

println("\n=== A. point queries ($(N_POINTS) points per variant)")
base = points[points.variant .== VARIANTS[1].tag, :]
for v in VARIANTS
    d = points[points.variant .== v.tag, :]
    rel_rhs = maximum(abs.(d.rho_rhs_kgm3 .- base.rho_rhs_kgm3) ./ base.rho_rhs_kgm3)
    rel_pert = maximum(abs.(d.gram_perturbed_density_kgm3 .- d.gram_mean_density_kgm3) ./ d.gram_mean_density_kgm3)
    rhs_eq_mean = maximum(abs.(d.rho_rhs_kgm3 .- d.gram_mean_density_kgm3) ./ d.gram_mean_density_kgm3)
    @printf("%-22s max|rho_rhs - rho_rhs(V0)|/rho = %.6e   max|perturbed - mean|/mean = %.6e   max|rho_rhs - mean|/mean = %.6e\n",
            v.tag, rel_rhs, rel_pert, rhs_eq_mean)
end
let a = points[points.variant .== VARIANTS[2].tag, :], b = points[points.variant .== VARIANTS[3].tag, :]
    @printf("V1 vs V2: max rel diff rho_rhs = %.6e, perturbedDensity = %.6e\n",
            maximum(abs.(a.rho_rhs_kgm3 .- b.rho_rhs_kgm3) ./ a.rho_rhs_kgm3),
            maximum(abs.(a.gram_perturbed_density_kgm3 .- b.gram_perturbed_density_kgm3) ./ a.gram_mean_density_kgm3))
end

# ── B. one-orbit solves ──────────────────────────────────────────────────────
println("\n=== B. one-orbit solves")
finals = NamedTuple[]
histories = Dict{String, DataFrame}()
for v in VARIANTS
    cfg = variant_cfg(v.seed, v.scales)
    args = TV._make_orbit_args(cfg, 1)
    args = TV._with_study_settings(args; quick=false)
    run = cd(() -> TV._run_simulation_dataframe(args, cfg.name, cfg.atmosphere_truth, :full;
                                                 extra_callbacks=TV._scenario_extra_callbacks(cfg)),
             REPO_ROOT)
    df = run.results_df
    histories[v.tag] = df
    CSV.write(joinpath(OUT, "one_orbit_$(v.tag).csv"), df)
    last_row = df[end, :]
    @printf("%-22s rows=%d solve=%.1f s retcode=%s\n", v.tag, nrow(df), run.elapsed_s, run.solver_info.solver_retcode)
    push!(finals, (variant=v.tag, rows=nrow(df), solve_s=run.elapsed_s, retcode=run.solver_info.solver_retcode))
end
CSV.write(joinpath(OUT, "one_orbit_summary.csv"), DataFrame(finals))

numeric_cols(df) = [n for n in names(df) if eltype(df[!, n]) <: Union{Missing, Real}]
ref = histories[VARIANTS[1].tag]
for v in VARIANTS[2:end]
    df = histories[v.tag]
    if nrow(df) != nrow(ref)
        println("$(v.tag): row count differs from V0 ($(nrow(df)) vs $(nrow(ref))) -- trajectories differ")
        continue
    end
    maxdiff = 0.0; worst = ""
    for n in intersect(numeric_cols(ref), numeric_cols(df))
        d = maximum(abs.(coalesce.(df[!, n], 0.0) .- coalesce.(ref[!, n], 0.0)))
        if d > maxdiff
            maxdiff = d; worst = n
        end
    end
    @printf("%-22s vs V0: max abs difference over all numeric columns = %.6e (column %s)\n", v.tag, maxdiff, worst)
end
println("\nSaved to $OUT")
