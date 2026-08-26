const STUDY_ROOT = @__DIR__
const STUDY_PROJECT = joinpath(STUDY_ROOT, "Project.toml")

if something(Base.active_project(), "") != STUDY_PROJECT
    import Pkg
    Pkg.activate(STUDY_ROOT; io=devnull)
end

# main.jl pulls in the whole study (corridor, priors, truth sources, GP,
# scoring) and defines StudySpec/spec_from_args. Its run guard is keyed on
# PROGRAM_FILE, so including it here does not launch a study.
include(joinpath(STUDY_ROOT, "main.jl"))
include(joinpath(STUDY_ROOT, "accel_filter.jl"))

# Compares the aerocapture-trained GP against an in-situ density-ratio filter of
# the kind used in Tracy/Falcone/Manchester, and against the two combined.
#
# Two metrics, and the second is the one that matters for guidance:
#
#  - "at vehicle": density error where the vehicle currently is. The filter is
#    strong here by construction; it is measuring that point.
#  - "forward": density error over the *remaining* trajectory as predicted at
#    each guidance call. A predictor-corrector re-simulates to parachute deploy
#    every call, so this is what it consumes. A reactive filter can only hold its
#    current ratio forward; a spatial model can vary it with position.

const VEHICLE = EntryVehicle()
const FILTER_CFG = DensityFilterConfig()

function run_case(case::StudyCase, spec::StudySpec, rng::AbstractRNG)
    pair = build_trajectory_pair(case;
        aerocapture_entry_alt_m=1.0e3 * spec.aerocapture_entry_alt_km,
        aerocapture_exit_alt_m=1.0e3 * spec.aerocapture_exit_alt_km)
    t0 = first(pair.aerocapture).dt
    aero_truth, edl_truth = truth_profiles(spec.truth_source, pair, t0)

    aero_prior = prior_samples(spec.prior_source, pair.aerocapture, t0;
        prior_seed=spec.prior_seed, n_dispersion=spec.n_dispersion)
    edl_prior = prior_samples(spec.prior_source, pair.edl, t0;
        prior_seed=spec.prior_seed, n_dispersion=spec.n_dispersion)

    # --- aerocapture-trained GP, best configuration from the study -----------
    par = Parameterization("log_density")
    cfg = GPConfig("gram_exponential", spec.lat_scale_deg, spec.lon_scale_deg,
                   spec.alt_scale_km, 1.0e-9, "constant", spec.gram_ref_alt_km,
                   spec.kernel_time_scale_s, "prior_scaled")
    meas = [max(aero_truth[i] + 0.05 * aero_truth[i] * randn(rng), 1.0e-18) for i in eachindex(aero_truth)]
    targets = Vector{Float64}(undef, length(pair.aerocapture))
    noises = similar(targets)
    psig = similar(targets)
    for i in eachindex(pair.aerocapture)
        yt, nv = _measurement_target(par, meas[i], aero_prior[i])
        pm, pv = _prior_target(par, aero_prior[i])
        targets[i] = yt - pm; noises[i] = nv; psig[i] = sqrt(pv)
    end
    gp = fit_gp_residual(pair.aerocapture, targets, noises, cfg; prior_sigma=psig)
    rho_gp, _ = predict_density_profile(par, gp, pair.edl, edl_prior)

    rho_prior = Float64[s.mean_density for s in edl_prior]
    speeds = corridor_speeds(pair.edl)

    # --- in-situ filter, against the plain prior and against the GP ----------
    k_prior, _, obs = run_density_filter(FILTER_CFG, VEHICLE, pair.edl, speeds, edl_truth, rho_prior; rng)
    k_gp, _, _ = run_density_filter(FILTER_CFG, VEHICLE, pair.edl, speeds, edl_truth, rho_gp; rng)

    w = pair.q_weights
    at_vehicle(p) = sqrt(sum(w .* (log.(p ./ edl_truth)) .^ 2) / sum(w))

    return (
        case_id=case.case_id,
        gap_hr=case.gap_s / 3600.0,
        observable_frac=count(obs) / length(obs),
        mean_speed_kms=mean(speeds) * 1e-3,
        # at the vehicle
        v_prior=at_vehicle(rho_prior),
        v_gp=at_vehicle(rho_gp),
        v_filter=at_vehicle(k_prior .* rho_prior),
        v_gp_filter=at_vehicle(k_gp .* rho_gp),
        # forward-predicted over the remaining trajectory
        f_prior=forward_prediction_error((j, i) -> rho_prior[i], edl_truth, w),
        f_gp=forward_prediction_error((j, i) -> rho_gp[i], edl_truth, w),
        f_filter=forward_prediction_error((j, i) -> k_prior[j] * rho_prior[i], edl_truth, w),
        f_gp_filter=forward_prediction_error((j, i) -> k_gp[j] * rho_gp[i], edl_truth, w),
        # Correlation-decayed coupling. Multiplying the whole forward profile by
        # the filter's current scalar re-introduces exactly the defect that makes
        # a reactive filter poor ahead of the vehicle: it smears a local ratio
        # over a region where the true ratio varies. Instead the filter's
        # correction is trusted at the vehicle and relaxes into the GP's own
        # spatial prediction with distance, using GRAM's vertical correlation
        # scale as the relaxation length.
        f_gp_filter_decay=forward_prediction_error(
            (j, i) -> begin
                dz_km = abs(pair.edl[i].alt_m - pair.edl[j].alt_m) * 1e-3
                _, lz = gram_correlation_scales(pair.edl[j].alt_m * 1e-3)
                rho_gp[i] * exp(log(k_gp[j]) * exp(-dz_km / lz))
            end, edl_truth, w),
        f_filter_decay=forward_prediction_error(
            (j, i) -> begin
                dz_km = abs(pair.edl[i].alt_m - pair.edl[j].alt_m) * 1e-3
                _, lz = gram_correlation_scales(pair.edl[j].alt_m * 1e-3)
                rho_prior[i] * exp(log(k_prior[j]) * exp(-dz_km / lz))
            end, edl_truth, w),
    )
end

function main()
    spec = spec_from_args(collect(ARGS))
    rng = MersenneTwister(spec.seed)
    cases = default_study_cases()
    spec.case_limit > 0 && (cases = first(cases, min(spec.case_limit, length(cases))))

    rows = [run_case(c, spec, rng) for c in cases]
    df = DataFrame(rows)
    mkpath(spec.out_dir)
    CSV.write(joinpath(spec.out_dir, "filter_comparison.csv"), df)

    @printf("\n%d cases | truth %s | prior %s | sensing %.0f-%.0f km | vehicle beta = %.0f kg/m^2\n",
            nrow(df), truth_source_name(spec.truth_source), prior_source_name(spec.prior_source),
            spec.aerocapture_entry_alt_km, spec.aerocapture_exit_alt_km, ballistic_coefficient(VEHICLE))
    @printf("mean corridor speed %.2f km/s | drag accel above noise floor on %.0f%% of the corridor\n\n",
            mean(df.mean_speed_kms), 100 * mean(df.observable_frac))

    base_v = mean(df.v_prior); base_f = mean(df.f_prior)
    println("Weighted RMS log-density error (lower is better), means over cases")
    println()
    @printf("%-28s %12s %10s %12s %10s\n", "", "at vehicle", "vs prior", "forward", "vs prior")
    for (lab, v, f) in (
        ("prior alone", df.v_prior, df.f_prior),
        ("in-situ filter only", df.v_filter, df.f_filter),
        ("aerocapture GP only", df.v_gp, df.f_gp),
        ("GP + filter (naive)", df.v_gp_filter, df.f_gp_filter),
        ("filter, decayed forward", df.v_filter, df.f_filter_decay),
        ("GP + filter (decayed)", df.v_gp_filter, df.f_gp_filter_decay),
    )
        @printf("%-28s %12.5f %9.1f%% %12.5f %9.1f%%\n", lab, mean(v),
                100 * (1 - mean(v) / base_v), mean(f), 100 * (1 - mean(f) / base_f))
    end
    println()
    println("by gap (forward-prediction error):")
    for g in sort(unique(df.gap_hr))
        sub = filter(r -> r.gap_hr == g, df)
        @printf("  %.0f hr:  filter %.5f   GP %.5f   GP+filter(decayed) %.5f\n",
                g, mean(sub.f_filter), mean(sub.f_gp), mean(sub.f_gp_filter_decay))
    end
    println("\nWrote $(joinpath(spec.out_dir, "filter_comparison.csv"))")
    return df
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
