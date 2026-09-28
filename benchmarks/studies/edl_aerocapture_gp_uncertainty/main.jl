const STUDY_ROOT = @__DIR__
const STUDY_PROJECT = joinpath(STUDY_ROOT, "Project.toml")

if something(Base.active_project(), "") != STUDY_PROJECT
    import Pkg
    Pkg.activate(STUDY_ROOT; io=devnull)
end

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using Random
using Statistics

using SpaceAGORA
import GRAMSuite

const SM = SpaceAGORA.SimulationModel
const RuntimeServices = SpaceAGORA.RuntimeServices

include(joinpath(STUDY_ROOT, "corridor.jl"))
include(joinpath(STUDY_ROOT, "gram_prior.jl"))
include(joinpath(STUDY_ROOT, "gram_correlation.jl"))
include(joinpath(STUDY_ROOT, "prior_sources.jl"))
include(joinpath(STUDY_ROOT, "merra2.jl"))
include(joinpath(STUDY_ROOT, "merra2_native.jl"))
include(joinpath(STUDY_ROOT, "truth_sources.jl"))
include(joinpath(STUDY_ROOT, "gp_models.jl"))
include(joinpath(STUDY_ROOT, "scoring.jl"))

struct StudySpec
    out_dir::String
    seed::Int
    truth_seed::Int
    prior_seed::Int
    n_dispersion::Int
    case_limit::Int
    save_pointwise::Bool
    lat_scale_deg::Float64
    lon_scale_deg::Float64
    alt_scale_km::Float64
    mean_bases::Vector{String}
    truth_source::TruthSource
    aerocapture_entry_alt_km::Float64
    aerocapture_exit_alt_km::Float64
    kernels::Vector{String}
    gram_ref_alt_km::Float64
    kernel_time_scale_s::Float64
    amplitude_mode::String
    lambda_bounds::Union{Nothing, Tuple{Float64, Float64}}
    prior_source::PriorSource
end

function _usage()
    return """
    Usage:
      julia --project=benchmarks/studies/edl_aerocapture_gp_uncertainty benchmarks/studies/edl_aerocapture_gp_uncertainty/main.jl [options]

    Options:
      --out-dir output/edl_aerocapture_gp_uncertainty
      --seed 42
      --truth-seed 10042
      --prior-seed 20042
      --n-dispersion 24
      --case-limit 0             # 0 runs all cases
      --save-pointwise true
      --report-only false                     # rebuild RESULTS.md from CSVs already in --out-dir
      --lat-scale-deg 20.5                    # GRAM's Lh of 2279 km at 35 km, as degrees of latitude
      --lon-scale-deg 21.8                    # the same 2279 km as degrees of longitude at 20N
      --alt-scale-km 18.33                    # GRAM's Lz at 35 km
      --mean-basis none,constant,linear_alt   # comma list; residual-GP mean function
      --kernels squared_exponential,matern32,matern52,gram_exponential
      --gram-scale-ref-alt-km 35.0            # altitude at which gram_exponential reads GRAM's scales
      --kernel-time-scale-hours auto          # auto uses GRAM's Lt at that altitude; off drops the time axis
      --prior-source gram                     # gram | nrlmsise00; nrlmsise00 breaks the MERRA-2 lineage
      --prior-sigma-rel auto                  # auto borrows GRAM's ensemble spread; or a constant relative sigma
      --msis-f107a 150.0  --msis-f107 150.0  --msis-ap 4.0
      --lambda-bounds auto                    # 'lo,hi', or auto: floor 1.0 with no mean basis, 0.1 with one
      --amplitude-mode prior_scaled           # prior_scaled folds the GRAM ensemble sigma into the kernel;
                                              # stationary reproduces the legacy global amplitude plus variance fusion

      --aerocapture-entry-alt-km 125.0        # highest altitude the sensing pass reaches
      --aerocapture-exit-alt-km 60.0          # lowest altitude the sensing pass reaches

    Truth source (see truth_sources.jl):
      --truth-source merra2                   # one of: $(join(TRUTH_SOURCE_NAMES, ", "))
      --merra2-hour-code 0                    # 0 = from epoch hour, 1..8 = 00Z..21Z, 9 = all-hours mean
      --merra2-dispersion 1.0                 # multiples of MERRA-2's own interannual sigma; 0 = raw field
      --merra2-blend-width-km 10.0            # taper above the 0.1 mb ceiling back to GRAM nominal
      --merra2-time-decorrelation true        # evolve the truth anomaly as AR(1) with GRAM's Lt (merra2 only)
      --merra2-native-dir <path>              # day-specific model-level granules for truth-source merra2_native
      --merra2-native-anomaly-top-km 52.0     # cap the anomaly below GRAM's own model-blending seam
      --truth-epoch-shift-hours 36.0          # gram_epoch_shift only
      --field-bias 0.08                       # synthetic_field only, planted log-density offset
      --field-amplitude 0.10                  # synthetic_field only, log-density anomaly std
      --field-lat-scale-deg 6.0
      --field-lon-scale-deg 12.0
      --field-alt-scale-km 25.0

    `gram_walk` reproduces the original configuration and is unsound: GRAM's
    perturbedDensity is a call-history random walk, not a function of position,
    so its aerocapture and EDL residuals are uncorrelated (measured: -0.13) and
    no position-indexed estimator can beat the prior on it.
    """
end

function _parse_cli(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            println(_usage())
            exit(0)
        end
        startswith(arg, "--") || throw(ArgumentError("Unsupported argument '$arg'."))
        body = arg[3:end]
        if occursin("=", body)
            key, value = split(body, "="; limit=2)
            opts[key] = value
        else
            if i == length(args)
                opts[body] = "true"
            else
                i += 1
                opts[body] = args[i]
            end
        end
        i += 1
    end
    return opts
end

@inline _get(opts::Dict{String, String}, key::String, default::String) = get(opts, key, default)
@inline _get_int(opts::Dict{String, String}, key::String, default::Int) = parse(Int, get(opts, key, string(default)))
@inline _get_float(opts::Dict{String, String}, key::String, default::Float64) = parse(Float64, get(opts, key, string(default)))
@inline _get_bool(opts::Dict{String, String}, key::String, default::Bool) = lowercase(strip(get(opts, key, string(default)))) in ("1", "true", "yes", "on")

function _parse_mean_bases(raw::AbstractString)::Vector{String}
    bases = String[strip(x) for x in split(raw, ",") if !isempty(strip(x))]
    isempty(bases) && throw(ArgumentError("--mean-basis needs at least one entry."))
    for b in bases
        b in MEAN_BASIS_NAMES ||
            throw(ArgumentError("Unsupported mean basis '$b'. Use one of: $(join(MEAN_BASIS_NAMES, ", "))."))
    end
    return bases
end

function _parse_lambda_bounds(raw::AbstractString)::Union{Nothing, Tuple{Float64, Float64}}
    lowercase(strip(raw)) == "auto" && return nothing
    parts = split(raw, ",")
    length(parts) == 2 || throw(ArgumentError("--lambda-bounds takes 'lo,hi'."))
    lo = parse(Float64, strip(parts[1])); hi = parse(Float64, strip(parts[2]))
    (0.0 < lo <= hi) || throw(ArgumentError("--lambda-bounds must satisfy 0 < lo <= hi."))
    return (lo, hi)
end

function _parse_amplitude_mode(raw::AbstractString)::String
    mode = strip(raw)
    mode in AMPLITUDE_MODE_NAMES ||
        throw(ArgumentError("Unsupported amplitude mode '$mode'. Use one of: $(join(AMPLITUDE_MODE_NAMES, ", "))."))
    return mode
end

function _parse_time_scale(raw::AbstractString, gram_ref_alt_km::Float64)::Float64
    token = lowercase(strip(raw))
    token == "auto" && return gram_time_scale_s(gram_ref_alt_km)
    token in ("off", "none", "inf") && return Inf
    hours = tryparse(Float64, token)
    hours === nothing && throw(ArgumentError("--kernel-time-scale-hours must be auto, off, or a number of hours."))
    hours > 0.0 || throw(ArgumentError("--kernel-time-scale-hours must be positive."))
    return 3600.0 * hours
end

function _parse_kernels(raw::AbstractString)::Vector{String}
    names = String[strip(x) for x in split(raw, ",") if !isempty(strip(x))]
    isempty(names) && throw(ArgumentError("--kernels needs at least one entry."))
    for k in names
        k in KERNEL_NAMES ||
            throw(ArgumentError("Unsupported kernel '$k'. Use one of: $(join(KERNEL_NAMES, ", "))."))
    end
    return names
end

function spec_from_args(args::Vector{String})::StudySpec
    opts = _parse_cli(args)
    seed = _get_int(opts, "seed", 42)
    truth_seed = _get_int(opts, "truth-seed", seed + 10_000)
    gram_ref_alt_km = _get_float(opts, "gram-scale-ref-alt-km", DEFAULT_GRAM_REF_ALT_KM)
    truth_source = truth_source_from_name(
        _get(opts, "truth-source", "merra2");
        seed=truth_seed,
        epoch_shift_s=3600.0 * _get_float(opts, "truth-epoch-shift-hours", 36.0),
        merra2_hour_code=_get_int(opts, "merra2-hour-code", 0),
        merra2_dispersion=_get_float(opts, "merra2-dispersion", 1.0),
        merra2_blend_width_km=_get_float(opts, "merra2-blend-width-km", 10.0),
        merra2_time_decorrelation=_get_bool(opts, "merra2-time-decorrelation", true),
        merra2_native_dir=_get(opts, "merra2-native-dir", ""),
        merra2_native_anomaly_top_km=_get_float(opts, "merra2-native-anomaly-top-km", 52.0),
        field_bias=_get_float(opts, "field-bias", 0.08),
        field_amplitude=_get_float(opts, "field-amplitude", 0.10),
        field_lat_scale_deg=_get_float(opts, "field-lat-scale-deg", 6.0),
        field_lon_scale_deg=_get_float(opts, "field-lon-scale-deg", 12.0),
        field_alt_scale_km=_get_float(opts, "field-alt-scale-km", 25.0),
    )
    return StudySpec(
        _get(opts, "out-dir", joinpath(STUDY_ROOT, "output")),
        seed,
        truth_seed,
        _get_int(opts, "prior-seed", seed + 20_000),
        _get_int(opts, "n-dispersion", 24),
        _get_int(opts, "case-limit", 0),
        _get_bool(opts, "save-pointwise", true),
        _get_float(opts, "lat-scale-deg", 20.5),
        _get_float(opts, "lon-scale-deg", 21.8),
        _get_float(opts, "alt-scale-km", 18.33),
        _parse_mean_bases(_get(opts, "mean-basis", "none,constant,linear_alt")),
        truth_source,
        _get_float(opts, "aerocapture-entry-alt-km", 125.0),
        _get_float(opts, "aerocapture-exit-alt-km", 60.0),
        _parse_kernels(_get(opts, "kernels", join(KERNEL_NAMES, ","))),
        gram_ref_alt_km,
        _parse_time_scale(_get(opts, "kernel-time-scale-hours", "auto"), gram_ref_alt_km),
        _parse_amplitude_mode(_get(opts, "amplitude-mode", "prior_scaled")),
        _parse_lambda_bounds(_get(opts, "lambda-bounds", "auto")),
        prior_source_from_name(
            _get(opts, "prior-source", "gram");
            f107a=_get_float(opts, "msis-f107a", 150.0),
            f107=_get_float(opts, "msis-f107", 150.0),
            ap=_get_float(opts, "msis-ap", 4.0),
            sigma_rel=lowercase(strip(_get(opts, "prior-sigma-rel", "auto"))) == "auto" ?
                NaN : _get_float(opts, "prior-sigma-rel", NaN),
        )
    )
end

function _kernel_comparison(summary_df::DataFrame)::DataFrame
    gp_rows = filter(:kernel => !=("gram_baseline"), summary_df)
    return combine(
        groupby(gp_rows, [:parameterization, :kernel, :mean_basis]),
        :rmse => mean => :mean_rmse,
        :weighted_rmse => mean => :mean_weighted_rmse,
        :rmse_log => mean => :mean_rmse_log,
        :weighted_rmse_log => mean => :mean_weighted_rmse_log,
        :nlpd => mean => :mean_nlpd,
        :weighted_nlpd => mean => :mean_weighted_nlpd,
        :coverage_1sigma => mean => :mean_coverage_1sigma,
        :coverage_2sigma => mean => :mean_coverage_2sigma,
    )
end

@inline function _clip_positive(x::Float64)::Float64
    return max(x, 1.0e-18)
end

function _noisy_measurements(rng::AbstractRNG, truth::Vector{Float64})::Vector{Float64}
    meas = Vector{Float64}(undef, length(truth))
    @inbounds for i in eachindex(truth)
        sigma = 0.05 * truth[i]
        meas[i] = _clip_positive(truth[i] + sigma * randn(rng))
    end
    return meas
end

function _baseline_vectors(samples::Vector{GramSample})
    mu = Float64[s.mean_density for s in samples]
    sigma = Float64[max(s.std_density, 1.0e-18) for s in samples]
    return mu, sigma
end

function _metric_row(
    case::StudyCase,
    parameterization::String,
    kernel::String,
    mean_basis::String,
    metrics::MetricSummary
)
    return (
        case_id=case.case_id,
        anchor_time_utc=string(case.anchor_time),
        gap_s=case.gap_s,
        lat_offset_deg=case.lat_offset_deg,
        lon_offset_deg=case.lon_offset_deg,
        parameterization=parameterization,
        kernel=kernel,
        mean_basis=mean_basis,
        rmse=metrics.rmse,
        weighted_rmse=metrics.weighted_rmse,
        rmse_log=metrics.rmse_log,
        weighted_rmse_log=metrics.weighted_rmse_log,
        nlpd=metrics.nlpd,
        weighted_nlpd=metrics.weighted_nlpd,
        coverage_1sigma=metrics.coverage_1sigma,
        coverage_2sigma=metrics.coverage_2sigma,
        weighted_coverage_1sigma=metrics.weighted_coverage_1sigma,
        weighted_coverage_2sigma=metrics.weighted_coverage_2sigma
    )
end

function _pointwise_rows(
    case::StudyCase,
    parameterization::String,
    kernel::String,
    mean_basis::String,
    points::Vector{TrajectoryPoint},
    truth::Vector{Float64},
    mu::Vector{Float64},
    sigma::Vector{Float64},
    weights::Vector{Float64}
)
    rows = Vector{NamedTuple}(undef, length(points))
    @inbounds for i in eachindex(points)
        p = points[i]
        rows[i] = (
            case_id=case.case_id,
            parameterization=parameterization,
            kernel=kernel,
            mean_basis=mean_basis,
            dt_utc=string(p.dt),
            elapsed_s=p.elapsed_s,
            lat_deg=p.lat_deg,
            lon_deg=p.lon_deg,
            alt_m=p.alt_m,
            q_weight=weights[i],
            truth_density=truth[i],
            pred_mean_density=mu[i],
            pred_std_density=sigma[i]
        )
    end
    return rows
end

"""
    run_config_row(spec) -> NamedTuple

Everything needed to reproduce a run and to label its results. Written to
`run_config.csv` so `--report-only` describes the run that produced the numbers
rather than whatever defaults happen to be on the command line at report time.
"""
function run_config_row(spec::StudySpec)
    lh, lz = gram_correlation_scales(spec.gram_ref_alt_km)
    return (
        truth_source=truth_source_name(spec.truth_source),
        truth_position_indexed=is_position_indexed(spec.truth_source),
        prior_source=prior_source_name(spec.prior_source),
        prior_shares_merra2_lineage=shares_lineage_with_merra2(spec.prior_source),
        aerocapture_entry_alt_km=spec.aerocapture_entry_alt_km,
        aerocapture_exit_alt_km=spec.aerocapture_exit_alt_km,
        amplitude_mode=spec.amplitude_mode,
        mean_bases=join(spec.mean_bases, "|"),
        kernels=join(spec.kernels, "|"),
        gram_ref_alt_km=spec.gram_ref_alt_km,
        gram_lh_km=lh,
        gram_lz_km=lz,
        kernel_time_scale_hr=isfinite(spec.kernel_time_scale_s) ? spec.kernel_time_scale_s / 3600.0 : NaN,
        lambda_bounds=spec.lambda_bounds === nothing ? "auto" : "$(spec.lambda_bounds[1]),$(spec.lambda_bounds[2])",
        lat_scale_deg=spec.lat_scale_deg,
        lon_scale_deg=spec.lon_scale_deg,
        alt_scale_km=spec.alt_scale_km,
        n_dispersion=spec.n_dispersion,
        seed=spec.seed,
        truth_seed=spec.truth_seed,
        prior_seed=spec.prior_seed,
    )
end

function _write_results_markdown(
    spec::StudySpec, summary_df::DataFrame, mean_fn_df::DataFrame;
    config=run_config_row(spec),
)
    out_dir = spec.out_dir
    results_path = joinpath(STUDY_ROOT, "RESULTS.md")
    io = IOBuffer()
    println(io, "# Results")
    println(io)
    println(io, "Generated on $(Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS")) UTC.")
    println(io)
    # Provenance comes from the recorded configuration, never from the current
    # command line: under --report-only those differ.
    println(io, "Truth source: `$(config.truth_source)` (position-indexed: $(config.truth_position_indexed)).")
    println(io, "Prior source: `$(config.prior_source)`.")
    if config.prior_shares_merra2_lineage && startswith(String(config.truth_source), "merra2")
        println(io)
        println(io, "> Prior and truth share a lineage: GRAM's nominal below about `65 km` is")
        println(io, "> the MERRA-2 climatology, so the residual here is a specific day minus")
        println(io, "> its own monthly mean. Run `--prior-source nrlmsise00` for an")
        println(io, "> independent prior.")
    end
    println(io, "Aerocapture sensing band: `$(config.aerocapture_entry_alt_km)` to `$(config.aerocapture_exit_alt_km) km`.")
    let t = isfinite(config.kernel_time_scale_hr) ? "`$(round(config.kernel_time_scale_hr, digits=2)) hr`" : "disabled"
        println(io, "GRAM correlation scales at `$(config.gram_ref_alt_km) km`: horizontal `$(round(config.gram_lh_km, digits=1)) km`, vertical `$(round(config.gram_lz_km, digits=2)) km`, temporal $t.")
    end
    println(io, "Kernels: $(join(("`" .* split(config.kernels, "|") .* "`"), ", ")).")
    println(io, "Mean bases: $(join(("`" .* split(config.mean_bases, "|") .* "`"), ", ")).")
    println(io, "Amplitude mode: `$(config.amplitude_mode)`, lambda bounds `$(config.lambda_bounds)`.")
    println(io)
    if !is_position_indexed(spec.truth_source)
        println(io, "> This truth source is a call-history random walk, not a spatial field.")
        println(io, "> No position-indexed estimator can beat the prior on it; these numbers")
        println(io, "> measure nothing about the estimator.")
        println(io)
    end

    baseline = filter(row -> row.kernel == "gram_baseline", summary_df)
    candidates = filter(row -> row.kernel != "gram_baseline", summary_df)

    # Aggregate across cases *before* ranking. Ranking raw rows would pick the
    # single luckiest (case, kernel, parameterization) combination and report it
    # as the model's performance — on the 36-case MERRA-2 run that reads as
    # 73.7% against GRAM where the model's actual mean is 49.2%.
    if isempty(candidates)
        println(io, "No GP results were generated.")
    else
        n_cases = length(unique(summary_df.case_id))
        base_log = mean(baseline.weighted_rmse_log)
        base_nlpd = mean(baseline.weighted_nlpd)
        base_c1 = mean(baseline.weighted_coverage_1sigma)
        base_c2 = mean(baseline.weighted_coverage_2sigma)

        models = combine(
            groupby(candidates, [:parameterization, :kernel, :mean_basis]),
            :weighted_rmse_log => mean => :wrmse_log,
            :weighted_nlpd => mean => :wnlpd,
            :weighted_coverage_1sigma => mean => :cov1,
            :weighted_coverage_2sigma => mean => :cov2,
            nrow => :n,
        )
        models.improvement = 100.0 .* (base_log .- models.wrmse_log) ./ max(base_log, 1.0e-18)
        sort!(models, :wrmse_log)

        println(io, "Means over $n_cases cases. Nominal coverage is `0.683` and `0.954`.")
        println(io)
        println(io, "## GRAM Baseline")
        println(io)
        println(io, "| weighted log RMSE | weighted NLPD | coverage `1-sigma` | coverage `2-sigma` |")
        println(io, "|---|---|---|---|")
        println(io, "| $(round(base_log, sigdigits=5)) | $(round(base_nlpd, sigdigits=5)) | $(round(base_c1, digits=3)) | $(round(base_c2, digits=3)) |")
        println(io)
        if base_c1 < 0.55
            println(io, "The baseline is well under-dispersed: GRAM's ensemble is roughly half as")
            println(io, "wide as the truth field requires. Every model below inherits that.")
            println(io)
        end

        println(io, "## Models, Ranked by Mean Weighted Log-Density RMSE")
        println(io)
        println(io, "| parameterization | kernel | mean basis | vs GRAM | weighted NLPD | cov `1-sigma` | cov `2-sigma` |")
        println(io, "|---|---|---|---|---|---|---|")
        for row in eachrow(first(models, 12))
            println(io, "| `$(row.parameterization)` | `$(row.kernel)` | `$(row.mean_basis)` | " *
                        "$(round(row.improvement, digits=2))% | $(round(row.wnlpd, digits=2)) | " *
                        "$(round(row.cov1, digits=3)) | $(round(row.cov2, digits=3)) |")
        end
        println(io)
        println(io, "Accuracy and calibration do not rank together: read both columns before")
        println(io, "picking a configuration.")
        println(io)

        println(io, "## By Mean Basis")
        println(io)
        println(io, "| mean basis | vs GRAM | weighted NLPD | cov `1-sigma` | cov `2-sigma` |")
        println(io, "|---|---|---|---|---|")
        for basis in split(config.mean_bases, "|")
            rows = filter(r -> r.mean_basis == basis, models)
            isempty(rows) && continue
            println(io, "| `$basis` | $(round(mean(rows.improvement), digits=2))% | " *
                        "$(round(mean(rows.wnlpd), digits=2)) | $(round(mean(rows.cov1), digits=3)) | " *
                        "$(round(mean(rows.cov2), digits=3)) |")
        end
        println(io)

        if :gap_s in propertynames(candidates) && length(unique(candidates.gap_s)) > 1
            println(io, "## By Aerocapture-to-EDL Gap")
            println(io)
            println(io, "| gap | vs GRAM |")
            println(io, "|---|---|")
            for gap in sort(unique(candidates.gap_s))
                rows = filter(r -> r.gap_s == gap, candidates)
                b = filter(r -> r.gap_s == gap, baseline)
                isempty(rows) || isempty(b) && continue
                bl = mean(b.weighted_rmse_log)
                println(io, "| $(round(gap / 3600.0, digits=1)) hr | " *
                            "$(round(100.0 * (bl - mean(rows.weighted_rmse_log)) / bl, digits=2))% |")
            end
            println(io)
        end
    end

    if !isempty(mean_fn_df)
        fits = filter(row -> row.mean_basis != "none", mean_fn_df)
        if !isempty(fits)
            println(io, "## Fitted Mean Function")
            println(io)
            println(io, "`beta_constant` is the bulk log-density offset the aerocapture pass")
            println(io, "inferred relative to the GRAM prior. Unlike the zero-mean GP term it")
            println(io, "extrapolates into the unsensed band.")
            println(io)
            println(io, "| parameterization | mean basis | beta_constant | beta_alt_slope | lambda | signal fraction |")
            println(io, "|---|---|---|---|---|---|")
            grouped = combine(
                groupby(fits, [:parameterization, :mean_basis]),
                :beta_constant => mean => :beta_constant,
                :beta_alt_slope => mean => :beta_alt_slope,
                :gp_amplitude => mean => :gp_amplitude,
                :signal_fraction => mean => :signal_fraction,
                :lambda_at_bound => (v -> any(v)) => :any_at_bound,
            )
            for row in eachrow(grouped)
                flag = row.any_at_bound ? " (some at bound)" : ""
                println(io, "| `$(row.parameterization)` | `$(row.mean_basis)` | $(round(row.beta_constant, sigdigits=4)) | $(round(row.beta_alt_slope, sigdigits=4)) | $(round(row.gp_amplitude, sigdigits=4))$flag | $(round(row.signal_fraction, digits=3)) |")
            end
            println(io)
            println(io, "`signal fraction` is the share of training-residual variance above the")
            println(io, "`5%` measurement-noise floor. Values near zero mean the pass carried")
            println(io, "little information and `lambda` is barely identifiable.")
            println(io)
        end
    end

    println(io, "## Output Files")
    println(io)
    println(io, "- `$(joinpath(out_dir, "summary_metrics.csv"))`")
    println(io, "- `$(joinpath(out_dir, "kernel_comparison.csv"))`")
    println(io, "- `$(joinpath(out_dir, "mean_function_fits.csv"))`")
    println(io, "- `$(joinpath(out_dir, "cases.csv"))`")
    println(io, "- `$(joinpath(out_dir, "pointwise_predictions.csv"))` when `--save-pointwise=true`")
    write(results_path, String(take!(io)))
    return results_path
end

function run_study(spec::StudySpec=spec_from_args(collect(ARGS)))
    mkpath(spec.out_dir)
    rng = MersenneTwister(spec.seed)
    gram_seed_plan(spec.truth_seed, spec.prior_seed, spec.n_dispersion)
    cases = default_study_cases()
    spec.case_limit >= 0 || throw(ArgumentError("--case-limit must be nonnegative."))
    if spec.case_limit > 0
        cases = first(cases, min(spec.case_limit, length(cases)))
    end

    metric_rows = NamedTuple[]
    point_rows = NamedTuple[]
    case_rows = NamedTuple[]
    mean_fn_rows = NamedTuple[]

    is_position_indexed(spec.truth_source) || @warn """
    Truth source '$(truth_source_name(spec.truth_source))' is not a function of position.
    GRAM's perturbedDensity is a call-history random walk, so the aerocapture and EDL
    residuals are uncorrelated and no position-indexed estimator can beat the prior.
    Results from this configuration measure nothing about the estimator.\
    """

    gp_template = GPConfig(
        "squared_exponential",
        spec.lat_scale_deg,
        spec.lon_scale_deg,
        spec.alt_scale_km,
        1.0e-9
    )

    for case in cases
        println("Running $(case.case_id) at $(case.anchor_time), gap=$(round(case.gap_s / 3600.0, digits=2)) hr")
        pair = build_trajectory_pair(case;
            aerocapture_entry_alt_m=1.0e3 * spec.aerocapture_entry_alt_km,
            aerocapture_exit_alt_m=1.0e3 * spec.aerocapture_exit_alt_km)
        aero_truth, edl_truth = truth_profiles(spec.truth_source, pair, first(pair.aerocapture).dt)
        measurements = _noisy_measurements(rng, aero_truth)

        aero_gram = prior_samples(
            spec.prior_source, pair.aerocapture, first(pair.aerocapture).dt;
            prior_seed=spec.prior_seed, n_dispersion=spec.n_dispersion,
        )
        edl_gram = prior_samples(
            spec.prior_source, pair.edl, first(pair.aerocapture).dt;
            prior_seed=spec.prior_seed, n_dispersion=spec.n_dispersion,
        )

        push!(
            case_rows,
            (
                case_id=case.case_id,
                anchor_time_utc=string(case.anchor_time),
                gap_s=case.gap_s,
                lat_offset_deg=case.lat_offset_deg,
                lon_offset_deg=case.lon_offset_deg,
                truth_seed=spec.truth_seed,
                prior_seed=spec.prior_seed,
                truth_source=truth_source_name(spec.truth_source),
                truth_position_indexed=is_position_indexed(spec.truth_source),
                aerocapture_entry_alt_km=spec.aerocapture_entry_alt_km,
                aerocapture_exit_alt_km=spec.aerocapture_exit_alt_km,
                aerocapture_points=length(pair.aerocapture),
                edl_points=length(pair.edl)
            )
        )

        base_mu, base_sigma = _baseline_vectors(edl_gram)
        base_metrics = score_predictions(edl_truth, base_mu, base_sigma, pair.q_weights)
        push!(metric_rows, _metric_row(case, "gram_prior", "gram_baseline", "none", base_metrics))
        if spec.save_pointwise
            append!(point_rows, _pointwise_rows(case, "gram_prior", "gram_baseline", "none", pair.edl, edl_truth, base_mu, base_sigma, pair.q_weights))
        end

        for par in PARAMETERIZATIONS
            train_targets = Vector{Float64}(undef, length(pair.aerocapture))
            train_noise = Vector{Float64}(undef, length(pair.aerocapture))
            train_prior_sigma = Vector{Float64}(undef, length(pair.aerocapture))
            for i in eachindex(pair.aerocapture)
                yt, nv = _measurement_target(par, measurements[i], aero_gram[i])
                pm, pv = _prior_target(par, aero_gram[i])
                train_targets[i] = yt - pm
                train_noise[i] = nv
                train_prior_sigma[i] = sqrt(pv)
            end

            for kernel_name in spec.kernels, mean_basis in spec.mean_bases
                cfg = GPConfig(
                    kernel_name,
                    gp_template.lat_scale_deg,
                    gp_template.lon_scale_deg,
                    gp_template.alt_scale_km,
                    gp_template.jitter,
                    mean_basis,
                    spec.gram_ref_alt_km,
                    spec.kernel_time_scale_s,
                    spec.amplitude_mode
                )
                fitted = fit_gp_residual(
                    pair.aerocapture, train_targets, train_noise, cfg;
                    prior_sigma=train_prior_sigma, lambda_bounds=spec.lambda_bounds
                )
                pred_mu, pred_sigma = predict_density_profile(
                    par,
                    fitted,
                    pair.edl,
                    edl_gram
                )
                metrics = score_predictions(edl_truth, pred_mu, pred_sigma, pair.q_weights)
                push!(metric_rows, _metric_row(case, par.name, kernel_name, mean_basis, metrics))
                report = mean_function_report(fitted)
                push!(
                    mean_fn_rows,
                    (
                        case_id=case.case_id,
                        parameterization=par.name,
                        kernel=kernel_name,
                        mean_basis=mean_basis,
                        beta_constant=report.beta_constant,
                        beta_alt_slope=report.beta_alt_slope,
                        center_alt_km=report.center_alt_km,
                        gp_amplitude=sqrt(report.amplitude2),
                        amplitude_mode=report.amplitude_mode,
                        lambda=report.lambda,
                        log_marginal=report.log_marginal,
                        signal_fraction=report.signal_fraction,
                        lambda_at_bound=report.lambda_at_bound,
                    )
                )
                if spec.save_pointwise
                    append!(point_rows, _pointwise_rows(case, par.name, kernel_name, mean_basis, pair.edl, edl_truth, pred_mu, pred_sigma, pair.q_weights))
                end
            end
        end
    end

    summary_df = DataFrame(metric_rows)
    cases_df = DataFrame(case_rows)
    kernel_df = _kernel_comparison(summary_df)
    CSV.write(joinpath(spec.out_dir, "summary_metrics.csv"), summary_df)
    CSV.write(joinpath(spec.out_dir, "cases.csv"), cases_df)
    CSV.write(joinpath(spec.out_dir, "kernel_comparison.csv"), kernel_df)
    CSV.write(joinpath(spec.out_dir, "mean_function_fits.csv"), DataFrame(mean_fn_rows))
    CSV.write(joinpath(spec.out_dir, "run_config.csv"), DataFrame([run_config_row(spec)]))
    if spec.save_pointwise
        CSV.write(joinpath(spec.out_dir, "pointwise_predictions.csv"), DataFrame(point_rows))
    end
    results_path = _write_results_markdown(spec, summary_df, DataFrame(mean_fn_rows))
    println("Wrote results summary to $results_path")
    return (summary=summary_df, cases=cases_df, results_path=results_path)
end

"""
    regenerate_report(spec)

Rewrite `RESULTS.md` from CSVs already in `spec.out_dir`, without re-running the
study. Invoked by `--report-only true`.
"""
function regenerate_report(spec::StudySpec)
    summary_path = joinpath(spec.out_dir, "summary_metrics.csv")
    isfile(summary_path) || throw(ErrorException("No summary_metrics.csv in $(spec.out_dir)."))
    mean_fn_path = joinpath(spec.out_dir, "mean_function_fits.csv")
    summary_df = CSV.read(summary_path, DataFrame)
    mean_fn_df = isfile(mean_fn_path) ? CSV.read(mean_fn_path, DataFrame) : DataFrame()
    config_path = joinpath(spec.out_dir, "run_config.csv")
    config = if isfile(config_path)
        NamedTuple(first(CSV.read(config_path, DataFrame)))
    else
        @warn """
        No run_config.csv in $(spec.out_dir); the header will describe the current
        command line rather than the run that produced these numbers. Re-run the
        study to record it.\
        """
        run_config_row(spec)
    end
    path = _write_results_markdown(spec, summary_df, mean_fn_df; config)
    println("Rewrote $path from $(spec.out_dir)")
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    local_spec = spec_from_args(collect(ARGS))
    if _get_bool(_parse_cli(collect(ARGS)), "report-only", false)
        regenerate_report(local_spec)
    else
        run_study(local_spec)
    end
end
