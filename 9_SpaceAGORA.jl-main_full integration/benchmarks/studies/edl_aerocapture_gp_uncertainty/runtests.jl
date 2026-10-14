using Dates
using LinearAlgebra
using Printf
using Random
using Statistics
using Test
using SpaceAGORA
import GRAMSuite

const STUDY_ROOT = @__DIR__
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

@testset "Independent GRAM seed plan" begin
    plan = gram_seed_plan(10_042, 20_042, 4)
    @test plan.truth_seed == 10_042
    @test plan.prior_nominal_seed == 20_042
    @test plan.prior_dispersion_seeds == [20_043, 20_044, 20_045, 20_046]
    @test_throws ArgumentError gram_seed_plan(10, 10, 4)
    @test_throws ArgumentError gram_seed_plan(10, 20, 1)
end

@testset "GP kernel estimation smoke" begin
    t0 = DateTime(2024, 3, 20, 18)
    points = [
        TrajectoryPoint(t0 + Second(2 * i), 2.0 * i, 20.0 + 0.1 * i, -60.0 + 0.1 * i, 60.0e3 + 1.0e3 * i)
        for i in 0:4
    ]
    residuals = [0.2, -0.1, 0.05, 0.15, -0.05]
    noise = fill(0.01, length(points))
    prior = fill(GramSample(1.0e-4, 1.0e-5), length(points))
    for kernel in KERNEL_NAMES
        model = fit_gp_residual(points, residuals, noise, GPConfig(kernel, 1.0, 1.0, 10.0, 1.0e-9))
        mu, variance = predict_gp_residual(model, points[3])
        means, stds = predict_density_profile(
            Parameterization("density_scale_factor"),
            model,
            points,
            prior,
        )
        @test isfinite(mu)
        @test variance > 0.0
        @test all(>(0.0), means)
        @test all(>(0.0), stds)
    end
end

@testset "Seeded GRAM truth and prior smoke" begin
    initial_dt = DateTime(2024, 3, 20, 18)
    point = TrajectoryPoint(initial_dt, 0.0, 20.0, -60.0, 35.0e3)
    pair = TrajectoryPair([point], [point], [1.0])
    truth_a, _ = gram_truth_profiles(pair, initial_dt; seed=10_042)
    truth_b, _ = gram_truth_profiles(pair, initial_dt; seed=10_043)
    prior = gram_prior_samples([point], initial_dt; prior_seed=20_042, n_dispersion=2)

    @test truth_a[1] > 0.0
    @test truth_a[1] != truth_b[1]
    @test prior[1].mean_density > 0.0
    @test prior[1].std_density > 0.0
end

@testset "Truth sources are position-indexed" begin
    @test !is_position_indexed(GramWalkTruth())
    @test is_position_indexed(GramEpochShiftTruth())
    @test is_position_indexed(SyntheticFieldTruth())
    @test truth_source_name(truth_source_from_name("synthetic_field")) == "synthetic_field"
    @test_throws ArgumentError truth_source_from_name("not_a_source")

    # The property GRAM's perturbedDensity lacks: the anomaly must be a function
    # of position alone, so the same query always returns the same value and the
    # aerocapture and EDL passes see the same field.
    src = SyntheticFieldTruth(; seed=7, bias=0.05, amplitude=0.1)
    t0 = DateTime(2024, 3, 20, 18)
    a = TrajectoryPoint(t0, 0.0, 20.0, -60.0, 35.0e3)
    b = TrajectoryPoint(t0 + Hour(6), 21_600.0, 20.0, -60.0, 35.0e3)
    @test synthetic_log_anomaly(src, a) == synthetic_log_anomaly(src, b)
    @test synthetic_log_anomaly(src, a) == synthetic_log_anomaly(src, a)
    @test synthetic_log_anomaly(src, a) != synthetic_log_anomaly(src, TrajectoryPoint(t0, 0.0, 40.0, -60.0, 35.0e3))

    # The field prior has the requested moments. Sampled across realizations at a
    # fixed point rather than across a single realization, because one draw over
    # this domain spans only a few correlation lengths per axis.
    draws = [
        synthetic_log_anomaly(SyntheticFieldTruth(; seed=s, bias=0.05, amplitude=0.1), a)
        for s in 1:600
    ]
    @test isapprox(mean(draws), 0.05; atol=0.015)
    @test isapprox(std(draws), 0.1; rtol=0.15)
end

@testset "Mean basis extrapolates a planted bias" begin
    t0 = DateTime(2024, 3, 20, 18)
    # Training points confined to 60-125 km, prediction far below at 35 km --
    # the geometry that made the zero-mean GP revert to the prior.
    train = [
        TrajectoryPoint(t0 + Second(2 * i), 2.0 * i, 20.0, -60.0, 1.0e3 * (125.0 - 0.5 * i))
        for i in 0:130
    ]
    far = TrajectoryPoint(t0 + Hour(3), 10_800.0, 20.0, -60.0, 35.0e3)
    bias = 0.30
    rng = MersenneTwister(11)
    residuals = [bias + 0.02 * randn(rng) for _ in train]
    noise = fill(0.05^2, length(train))

    zero_mean = fit_gp_residual(train, residuals, noise, GPConfig("matern32", 1.5, 1.5, 12.0, 1.0e-9, "none"))
    with_mean = fit_gp_residual(train, residuals, noise, GPConfig("matern32", 1.5, 1.5, 12.0, 1.0e-9, "constant"))

    mu_zero, var_zero = predict_gp_residual(zero_mean, far)
    mu_mean, var_mean = predict_gp_residual(with_mean, far)

    # Two length scales below the lowest training point the zero-mean GP has shed
    # most of the bias; the parametric mean carries all of it.
    @test abs(mu_zero) < 0.15 * bias
    @test isapprox(mu_mean, bias; rtol=0.05)
    @test abs(mu_mean - bias) < abs(mu_zero - bias) / 5.0
    @test var_mean > 0.0 && isfinite(var_mean)
    @test isempty(zero_mean.beta)
    @test length(with_mean.beta) == 1
    @test isapprox(mean_function_report(with_mean).beta_constant, bias; rtol=0.05)

    # Both still interpolate the training data.
    for model in (zero_mean, with_mean)
        mu_in, _ = predict_gp_residual(model, train[65])
        @test isapprox(mu_in, bias; atol=0.05)
    end

    # linear_alt recovers an altitude trend that the constant basis cannot.
    slope_per_km = -0.004
    sloped = [bias + slope_per_km * (p.alt_m * 1e-3 - 90.0) + 0.01 * randn(rng) for p in train]
    lin = fit_gp_residual(train, sloped, noise, GPConfig("matern32", 1.5, 1.5, 12.0, 1.0e-9, "linear_alt"))
    rep = mean_function_report(lin)
    @test isapprox(rep.beta_alt_slope / rep.center_alt_km * 0.0 + rep.beta_alt_slope / 12.0, slope_per_km; rtol=0.1)
    mu_lin, var_lin = predict_gp_residual(lin, far)
    @test isapprox(mu_lin, bias + slope_per_km * (35.0 - 90.0); rtol=0.15)
    # Extrapolating the basis must cost variance relative to a point inside it.
    _, var_inside = predict_gp_residual(lin, train[65])
    @test var_lin > var_inside
end

@testset "Mean basis validation" begin
    @test _basis_dim("none") == 0
    @test _basis_dim("constant") == 1
    @test _basis_dim("linear_alt") == 2
    @test_throws ArgumentError _basis_dim("quadratic")
    # Default GPConfig constructor stays zero-mean.
    @test GPConfig("matern32", 1.0, 1.0, 10.0, 1.0e-9).mean_basis == "none"
end

@testset "MERRA-2 reader" begin
    @test merra2_hour_code(0) == 1
    @test merra2_hour_code(18) == 7
    @test merra2_hour_code(21) == 8
    @test merra2_hour_code(23) == 1          # rolls over, matching buildDataFileName
    @test endswith(merra2_file_path(3, 7), joinpath("18Z", "MERRA2_3hr_18Z_03.bin"))
    @test endswith(merra2_file_path(10, 9), joinpath("All Mean", "MERRA2All_10.bin"))
    @test_throws ArgumentError merra2_file_path(13, 9)
    @test_throws ArgumentError merra2_file_path(3, 11)

    grid = load_merra2_grid(3, 9)
    @test (grid.n_pres, grid.n_lat, grid.n_lon) == (42, 90, 180)
    @test load_merra2_grid(3, 9) === grid    # cached

    # Ideal gas closes on the vendored grid, which pins the block order: if the
    # dens block were misidentified this would not hold.
    lat, lon = 20.0, -72.0
    rho, rel = merra2_density(grid, lat, lon, 35.0e3)
    @test isfinite(rho) && rho > 0.0
    @test 0.0 < rel < 0.1
    @test 5.0e-3 < rho < 2.0e-2             # ~8e-3 kg/m^3 at 35 km

    # Monotone decreasing with altitude through the domain.
    profile = [first(merra2_density(grid, lat, lon, 1.0e3 * z)) for z in 10.0:5.0:60.0]
    @test all(isfinite, profile)
    @test all(diff(profile) .< 0.0)

    # Vertical domain: about 0.1 km to 64 km.
    ceiling = merra2_ceiling_m(grid, lat, lon)
    @test 60.0e3 < ceiling < 70.0e3
    @test !isfinite(first(merra2_density(grid, lat, lon, ceiling + 5.0e3)))
    @test !isfinite(first(merra2_density(grid, lat, lon, -2.0e3)))

    # Longitude wraps rather than erroring at the seam.
    @test first(merra2_density(grid, lat, 359.5, 35.0e3)) ≈
          first(merra2_density(grid, lat, -0.5, 35.0e3))

    # Time-of-day slots differ from the all-hours mean, but only slightly: this
    # is why Merra2Truth disperses by MERRA-2's own sigma instead of relying on
    # the slot anomaly alone.
    slot = load_merra2_grid(3, 7)
    rho_slot, _ = merra2_density(slot, lat, lon, 35.0e3)
    @test rho_slot != rho
    @test abs(log(rho_slot / rho)) < 0.02
end

@testset "MERRA-2 truth source" begin
    src = Merra2Truth(; dispersion=1.0, seed=5)
    @test truth_source_name(src) == "merra2"
    @test is_position_indexed(src)
    @test truth_source_name(truth_source_from_name("merra2")) == "merra2"
    @test_throws ArgumentError Merra2Truth(; hour_code=12)
    @test_throws ArgumentError Merra2Truth(; dispersion=-1.0)
    @test_throws ArgumentError Merra2Truth(; blend_width_km=0.0)

    initial_dt = DateTime(2024, 3, 20, 18)
    alts = [125.0e3, 70.0e3, 60.0e3, 35.0e3, 10.0e3]
    points = [TrajectoryPoint(initial_dt, 0.0, 20.0, -72.0, a) for a in alts]
    pair = TrajectoryPair(points, points, ones(length(points)))
    truth, edl_truth = truth_profiles(src, pair, initial_dt)
    @test truth == edl_truth
    @test all(isfinite, truth)
    @test all(>(0.0), truth)
    @test all(diff(truth) .> 0.0)   # alts run 125 km down to 10 km

    # Deterministic in position: the property GRAM's walk lacks.
    again, _ = truth_profiles(src, pair, initial_dt)
    @test truth == again

    # Above the ceiling the anomaly tapers out and truth returns to GRAM nominal.
    grid = load_merra2_grid(3, merra2_hour_code(18))
    field = _build_fourier_field(src)
    model = _build_gram_model(_to_initial_time(initial_dt), "earth", 1)
    _, taper_high, inside_high = merra2_anomaly(src, grid, field, model, points[1], initial_dt)
    anom_low, taper_low, inside_low = merra2_anomaly(src, grid, field, model, points[4], initial_dt)
    @test inside_low && taper_low == 1.0
    @test !inside_high && taper_high == 0.0
    @test abs(anom_low) > 0.0

    # dispersion=0 leaves the raw MERRA-2 anomaly, which is small because GRAM's
    # nominal is built from the same reanalysis.
    raw = Merra2Truth(; dispersion=0.0)
    anom_raw, _, _ = merra2_anomaly(raw, grid, _build_fourier_field(raw), model, points[4], initial_dt)
    @test abs(anom_raw) < 0.05
    @test anom_raw != anom_low
end

@testset "GRAM correlation scales" begin
    lh35, lz35 = gram_correlation_scales(35.0)
    @test lh35 == 2279.4 && lz35 == 18.33
    lh, lz = gram_correlation_scales(32.5)          # midpoint of 30 and 35
    @test lh ≈ 0.5 * (2240.3 + 2279.4)
    @test lz ≈ 0.5 * (17.44 + 18.33)
    @test gram_correlation_scales(-10.0) == (117.5, 2.65)      # clamped low
    @test gram_correlation_scales(500.0) == (5040.0, 48.06)    # clamped high
    @test length(GRAM_RS_ALT_KM) == length(GRAM_RS_LH_KM) == length(GRAM_RS_LZ_KM)
    @test issorted(GRAM_RS_ALT_KM)

    @test chordal_distance_km(20.0, -72.0, 20.0, -72.0) == 0.0
    # 1 degree of latitude is about 111 km.
    @test isapprox(chordal_distance_km(20.0, -72.0, 21.0, -72.0), 111.2; rtol=0.01)
    # Chordal tracks great-circle closely at the separations that matter here.
    gc = 6371.0 * acos(sind(20.0) * sind(20.0) + cosd(20.0) * cosd(20.0) * cosd(25.0))
    @test isapprox(chordal_distance_km(20.0, 0.0, 20.0, 25.0), gc; rtol=0.02)
end

@testset "gram_exponential kernel" begin
    t0 = DateTime(2024, 3, 20, 18)
    cfg = GPConfig("gram_exponential", 1.5, 1.5, 12.0, 1.0e-9, "none", 35.0)
    a = TrajectoryPoint(t0, 0.0, 20.0, -72.0, 60.0e3)
    b = TrajectoryPoint(t0, 0.0, 20.0, -72.0, 35.0e3)   # 25 km below, same column

    train = [TrajectoryPoint(t0 + Second(2i), 2.0i, 20.0, -72.0, 1.0e3 * (125.0 - 0.5i)) for i in 0:130]
    resid = fill(0.2, length(train))
    noise = fill(0.05^2, length(train))
    model = fit_gp_residual(train, resid, noise, cfg)
    mu, var = predict_gp_residual(model, b)
    @test isfinite(mu) && var > 0.0

    # GRAM's vertical scale is 18.33 km at 35 km, so 25 km of extrapolation
    # retains far more correlation than the study default of 12 km with matern32.
    k_gram = _kernel(cfg, _feature_row(a), _feature_row(b), 1.0)
    k_study = _kernel(GPConfig("matern32", 1.5, 1.5, 12.0, 1.0e-9), _feature_row(a), _feature_row(b), 1.0)
    @test k_gram ≈ exp(-25.0 / 18.33)
    @test isapprox(k_gram / k_study, 2.05; rtol=0.02)

    # Separable: horizontal and vertical factors multiply.
    c = TrajectoryPoint(t0, 0.0, 25.0, -72.0, 35.0e3)
    lh, lz = gram_correlation_scales(35.0)
    dx = chordal_distance_km(20.0, -72.0, 25.0, -72.0)
    @test _kernel(cfg, _feature_row(a), _feature_row(c), 1.0) ≈ exp(-dx / lh) * exp(-25.0 / lz)

    # Kernel matrix stays positive definite on the real corridor geometry.
    pair = build_trajectory_pair(first(default_study_cases()); aerocapture_exit_alt_m=35.0e3)
    m2 = fit_gp_residual(pair.aerocapture, randn(MersenneTwister(2), length(pair.aerocapture)),
                         fill(0.05^2, length(pair.aerocapture)), cfg)
    @test all(>(0.0), diag(m2.chol.U))
    for p in pair.edl[1:40:end]
        _, v = predict_gp_residual(m2, p)
        @test v > 0.0
    end
end

@testset "GRAM time scale" begin
    # Lt = max(3 hr, 0.735 day * h_km^0.116), from getCorrelationCoefficients.
    @test gram_time_scale_s(35.0) ≈ 86400.0 * 0.735 * 35.0^0.116
    @test isapprox(gram_time_scale_s(35.0) / 3600.0, 26.6; rtol=0.01)
    @test isapprox(gram_time_scale_s(125.0) / 3600.0, 30.9; rtol=0.01)
    @test gram_time_scale_s(0.0) == 10800.0            # floor binds only at the ground
    @test gram_time_scale_s(125.0) > gram_time_scale_s(10.0)
end

@testset "Truth time decorrelation" begin
    initial_dt = DateTime(2024, 3, 20, 18)
    pt(lat, alt, offset_s) = TrajectoryPoint(
        initial_dt + Millisecond(round(Int, 1000 * offset_s)), offset_s, lat, -72.0, alt
    )

    # Same places, sampled now and six hours later.
    now_pts = [pt(18.0 + 0.5i, 1.0e3 * (55.0 - 2.0i), 2.0i) for i in 0:20]
    later_pts = [TrajectoryPoint(p.dt + Hour(6), p.elapsed_s, p.lat_deg, p.lon_deg, p.alt_m) for p in now_pts]
    pair = TrajectoryPair(now_pts, later_pts, ones(length(now_pts)))

    frozen = Merra2Truth(; dispersion=1.0, seed=3, time_decorrelation=false)
    evolving = Merra2Truth(; dispersion=1.0, seed=3, time_decorrelation=true)

    a_frozen, b_frozen = truth_profiles(frozen, pair, initial_dt)
    a_eq, b_eq = truth_profiles(evolving, pair, initial_dt)

    @test all(isfinite, b_eq) && all(>(0.0), b_eq)
    # Reproducible: the AR(1) innovation is seeded, not redrawn per call.
    @test truth_profiles(evolving, pair, initial_dt)[2] == b_eq

    # Frozen truth gives the same field at both epochs apart from GRAM nominal's
    # own mild elapsed-time drift; the evolving one genuinely moves.
    drift_frozen = maximum(abs.(log.(b_frozen ./ a_frozen)))
    drift_evolving = maximum(abs.(log.(b_eq ./ a_eq)))
    @test drift_evolving > 3.0 * drift_frozen

    # The aerocapture leg starts at the anchor, so it is essentially undecayed.
    @test maximum(abs.(log.(a_eq ./ a_frozen))) < 0.02 * drift_evolving + 1.0e-3

    # Truth must not depend on where the sensing pass flew, or sweeping
    # --aerocapture-exit-alt-km stops being a controlled comparison. The AR(1)
    # innovation is scaled and referenced off the prediction corridor only.
    case = first(default_study_cases())
    edl_of(exit_km) = truth_profiles(
        evolving, build_trajectory_pair(case; aerocapture_exit_alt_m=1.0e3 * exit_km),
        first(build_trajectory_pair(case).aerocapture).dt,
    )[2]
    @test edl_of(60.0) == edl_of(45.0)
    @test edl_of(60.0) == edl_of(25.0)

    # rho_t follows GRAM: 6 hr at 35 km retains about 0.80.
    rho6 = exp(-6 * 3600.0 / gram_time_scale_s(35.0))
    @test isapprox(rho6, 0.80; atol=0.01)
    @test isapprox(sqrt(1 - rho6^2), 0.60; atol=0.02)
end

@testset "Kernel time axis" begin
    t0 = DateTime(2024, 3, 20, 18)
    a = TrajectoryPoint(t0, 0.0, 20.0, -72.0, 40.0e3)
    b = TrajectoryPoint(t0 + Hour(6), 0.0, 20.0, -72.0, 40.0e3)  # same place, 6 hr later

    lt = gram_time_scale_s(35.0)
    timed = GPConfig("gram_exponential", 20.5, 21.8, 18.33, 1.0e-9, "none", 35.0, lt)
    untimed = GPConfig("gram_exponential", 20.5, 21.8, 18.33, 1.0e-9, "none", 35.0, Inf)

    # Co-located points separated only in time.
    @test _kernel(untimed, _feature_row(a), _feature_row(b), 1.0) ≈ 1.0
    @test _kernel(timed, _feature_row(a), _feature_row(b), 1.0) ≈ exp(-6 * 3600.0 / lt)
    @test _kernel(timed, _feature_row(a), _feature_row(a), 1.0) ≈ 1.0

    # The default constructor wires GRAM's Lt at the reference altitude.
    @test GPConfig("gram_exponential", 20.5, 21.8, 18.33, 1.0e-9, "none", 35.0).time_scale_s ≈ lt
    @test GPConfig("matern32", 1.0, 1.0, 10.0, 1.0e-9).time_scale_s ≈ gram_time_scale_s(DEFAULT_GRAM_REF_ALT_KM)

    # Euclidean kernels pick the time axis up too.
    m_timed = GPConfig("matern32", 20.5, 21.8, 18.33, 1.0e-9, "none", 35.0, lt)
    m_untimed = GPConfig("matern32", 20.5, 21.8, 18.33, 1.0e-9, "none", 35.0, Inf)
    @test _kernel(m_timed, _feature_row(a), _feature_row(b), 1.0) <
          _kernel(m_untimed, _feature_row(a), _feature_row(b), 1.0)

    # Features carry a shared absolute time origin, not per-corridor elapsed_s.
    pair = build_trajectory_pair(first(default_study_cases()); aerocapture_exit_alt_m=45.0e3)
    @test _feature_row(first(pair.edl))[4] > _feature_row(last(pair.aerocapture))[4]

    # Still positive definite on the real corridor with the time axis live.
    model = fit_gp_residual(pair.aerocapture, randn(MersenneTwister(4), length(pair.aerocapture)),
                            fill(0.05^2, length(pair.aerocapture)), timed)
    @test all(>(0.0), diag(model.chol.U))
    for p in pair.edl[1:40:end]
        _, v = predict_gp_residual(model, p)
        @test v > 0.0
    end
end

@testset "Prior-scaled amplitude" begin
    t0 = DateTime(2024, 3, 20, 18)
    train = [TrajectoryPoint(t0 + Second(2i), 2.0i, 20.0, -72.0, 1.0e3 * (125.0 - 0.5i)) for i in 0:130]
    far = TrajectoryPoint(t0 + Hour(3), 10_800.0, 20.0, -72.0, 35.0e3)

    # Prior sigma varies 30x across the corridor, exactly the range a single
    # global amplitude cannot represent.
    s_of(alt_km) = 0.01 + 0.10 * clamp((alt_km - 55.0) / 70.0, 0.0, 1.0)
    train_sigma = [s_of(p.alt_m * 1e-3) for p in train]
    rng = MersenneTwister(21)
    resid = [train_sigma[i] * randn(rng) for i in eachindex(train)]
    noise = fill(0.05^2, length(train))

    base = GPConfig("gram_exponential", 20.5, 21.8, 18.33, 1.0e-9, "none", 35.0)
    scaled = GPConfig("gram_exponential", 20.5, 21.8, 18.33, 1.0e-9, "none", 35.0,
                      gram_time_scale_s(35.0), "prior_scaled")

    @test_throws ArgumentError fit_gp_residual(train, resid, noise, scaled)   # needs prior sigma
    model = fit_gp_residual(train, resid, noise, scaled; prior_sigma=train_sigma)
    @test length(model.prior_sigma) == length(train)
    @test isfinite(model.log_marginal)
    lambda = sqrt(model.amplitude2)
    @test 0.3 < lambda < 4.0
    @test mean_function_report(model).lambda ≈ lambda

    @test_throws ArgumentError predict_gp_residual(model, far)               # needs prior sigma

    # Far from data the posterior returns to the prior: lambda^2 s(x)^2.
    s_far = s_of(35.0)
    isolated = TrajectoryPoint(t0, 0.0, -60.0, 100.0, 35.0e3)   # other side of the planet
    _, var_iso = predict_gp_residual(model, isolated, s_far)
    @test isapprox(var_iso, lambda^2 * s_far^2; rtol=0.02)

    # Near data it tightens, and never below zero.
    _, var_in = predict_gp_residual(model, train[65], train_sigma[65])
    @test 0.0 < var_in < lambda^2 * train_sigma[65]^2

    # The predictive sigma tracks the local prior sigma rather than a constant.
    _, v_low = predict_gp_residual(model, isolated, s_of(15.0))
    _, v_high = predict_gp_residual(model, TrajectoryPoint(t0, 0.0, -60.0, 100.0, 120.0e3), s_of(120.0))
    @test isapprox(sqrt(v_high / v_low), s_of(120.0) / s_of(15.0); rtol=0.05)

    # The legacy path caps the reported variance at the kernel amplitude
    # everywhere; the new one does not.
    legacy = fit_gp_residual(train, resid, noise, base)
    gram_hi = GramSample(1.0e-6, 1.0e-6 * s_of(120.0))
    par = Parameterization("log_density")
    _, prior_var_hi = _prior_target(par, gram_hi)
    _, dv = predict_gp_residual(legacy, isolated)
    fused = 1.0 / (1.0 / prior_var_hi + 1.0 / dv)
    @test fused < legacy.amplitude2
    @test fused < prior_var_hi                       # shrinks a prior it learned nothing about
    _, new_var = predict_gp_residual(model, isolated, sqrt(prior_var_hi))
    @test isapprox(new_var, lambda^2 * prior_var_hi; rtol=0.02)

    # lambda responds to genuine ensemble mis-dispersion: residuals three times
    # wider than the stated prior should pull it up.
    wide = [3.0 * r for r in resid]
    inflated = fit_gp_residual(train, wide, noise, scaled; prior_sigma=train_sigma)
    @test sqrt(inflated.amplitude2) > 1.5 * lambda

    # predict_density_profile applies no fusion in prior_scaled mode.
    gram = [GramSample(1.0e-3, 1.0e-3 * s_of(p.alt_m * 1e-3)) for p in train]
    mu, sd = predict_density_profile(Parameterization("density_scale_factor"), model, train, gram)
    @test all(isfinite, mu) && all(>(0.0), sd)
end

@testset "Stationary mode is unchanged" begin
    t0 = DateTime(2024, 3, 20, 18)
    pts = [TrajectoryPoint(t0 + Second(2i), 2.0i, 20.0 + 0.1i, -60.0 + 0.1i, 60.0e3 + 1.0e3 * i) for i in 0:4]
    r = [0.2, -0.1, 0.05, 0.15, -0.05]
    nv = fill(0.01, length(pts))
    cfg = GPConfig("matern32", 1.0, 1.0, 10.0, 1.0e-9)
    @test cfg.amplitude_mode == "stationary"
    m = fit_gp_residual(pts, r, nv, cfg)
    @test m.amplitude2 ≈ max(var(r), 1.0e-6)      # unchanged global amplitude
    @test isempty(m.prior_sigma)
    mu, v = predict_gp_residual(m, pts[3])
    @test isfinite(mu) && v > 0.0
end

# Builds a granule with exactly the schema the MERRA-2 File Specification gives
# for M2I3NVASM: `tzyx` variables that NCDatasets presents as
# (lon, lat, lev, time), 72 model levels indexed top-down, fill 1e15.
#
# The planted atmosphere has log-density linear in latitude, longitude, time and
# geometric altitude, so every interpolation the reader performs is exact on it
# and the recovered values can be checked to machine precision. That pins the
# level flip, the geopotential conversion, the fill handling and the time axis
# without needing a 2.1 GB download.
function _write_synthetic_native_granule(path::String, day::Date; corrupt_levels::Bool=false)
    nlev = corrupt_levels ? 40 : MERRA2_NATIVE_LEVELS
    lats = collect(15.0:0.5:25.0)
    lons = collect(-75.0:0.625:-55.0)
    alts_km = collect(range(0.5, 80.0; length=nlev))          # ascending geometric km
    hours = collect(0.0:3.0:21.0)

    rho0, scale_km = 1.2, 7.5
    c_lat, c_lon, c_t = 0.011, -0.004, 0.006
    plant(z_km, lat, lon, hr) = rho0 * exp(-z_km / scale_km) *
        exp(c_lat * (lat - 20.0) + c_lon * (lon + 65.0) + c_t * hr)
    geopotential(z_m) = EARTH_RADIUS_M * z_m / (EARTH_RADIUS_M + z_m)

    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        defDim(ds, "lon", length(lons)); defDim(ds, "lat", length(lats))
        defDim(ds, "lev", nlev);         defDim(ds, "time", length(hours))
        defVar(ds, "lon", Float64, ("lon",))[:] = lons
        defVar(ds, "lat", Float64, ("lat",))[:] = lats
        defVar(ds, "lev", Float64, ("lev",))[:] = collect(1.0:nlev)
        tv = defVar(ds, "time", Int32, ("time",); attrib=[
            "units" => "minutes since $(Dates.format(day, dateformat"yyyy-mm-dd")) 00:00:00",
        ])
        tv[:] = Int32.(round.(hours .* 60))

        dims = ("lon", "lat", "lev", "time")
        pl = defVar(ds, "PL", Float64, dims; attrib=["_FillValue" => 1.0e15])
        ta = defVar(ds, "T", Float64, dims; attrib=["_FillValue" => 1.0e15])
        hh = defVar(ds, "H", Float64, dims; attrib=["_FillValue" => 1.0e15])
        qv = defVar(ds, "QV", Float64, dims; attrib=["_FillValue" => 1.0e15])

        T0 = 250.0
        for (it, hr) in enumerate(hours), k in 1:nlev
            kk = nlev - k + 1                                   # write top-down
            z_km = alts_km[kk]
            for (j, lat) in enumerate(lats), (i, lon) in enumerate(lons)
                # Bottom two layers below terrain in one corner: exercise fill.
                if kk <= 2 && lat > 24.0 && lon > -56.0
                    pl[i, j, k, it] = 1.0e15; ta[i, j, k, it] = 1.0e15
                    hh[i, j, k, it] = 1.0e15; qv[i, j, k, it] = 1.0e15
                    continue
                end
                rho = plant(z_km, lat, lon, hr)
                ta[i, j, k, it] = T0
                qv[i, j, k, it] = 0.0
                pl[i, j, k, it] = rho * MERRA2_NATIVE_R_DRY * T0
                hh[i, j, k, it] = geopotential(1.0e3 * z_km)
            end
        end
    end
    return (; lats, lons, alts_km, hours, plant)
end

@testset "MERRA-2 native reader" begin
    @test merra2_native_stream(1985) == "100"
    @test merra2_native_stream(1995) == "200"
    @test merra2_native_stream(2005) == "300"
    @test merra2_native_stream(2024) == "400"
    @test_throws ArgumentError merra2_native_stream(1970)
    @test merra2_native_filename(Date(2024, 3, 20)) == "MERRA2_400.inst3_3d_asm_Nv.20240320.nc4"

    dir = mktempdir()
    day = Date(2024, 3, 20)
    path = joinpath(dir, merra2_native_filename(day))
    truth = _write_synthetic_native_granule(path, day)

    t0 = DateTime(2024, 3, 20, 1); t1 = DateTime(2024, 3, 20, 20)
    w = load_merra2_native([path], t0, t1)
    @test size(w.density, 1) == MERRA2_NATIVE_LEVELS
    @test issorted(w.lats) && issorted(w.lons)
    @test length(w.times) >= 7

    # Levels come back ascending in height despite the top-down file order.
    col = w.height[:, 5, 5, 1]
    @test issorted(filter(isfinite, col))
    @test isapprox(minimum(filter(isfinite, col)), 500.0; rtol=0.02)
    @test isapprox(maximum(filter(isfinite, col)), 80_000.0; rtol=0.02)

    # Every interpolation is exact on the planted field.
    for (lat, lon, z_km, hr) in (
        (20.0, -65.0, 35.0, 6.0), (18.25, -61.5625, 12.0, 6.0),
        (21.7, -58.3, 55.0, 7.5), (17.1, -70.9, 3.0, 16.25),
    )
        got = merra2_native_density(w, lat, lon, 1.0e3 * z_km, DateTime(2024, 3, 20) + Millisecond(round(Int, 3.6e6 * hr)))
        @test isapprox(got, truth.plant(z_km, lat, lon, hr); rtol=1.0e-6)
    end

    # Time really is interpolated, not snapped to the nearest analysis.
    a = merra2_native_density(w, 20.0, -65.0, 35.0e3, DateTime(2024, 3, 20, 6))
    b = merra2_native_density(w, 20.0, -65.0, 35.0e3, DateTime(2024, 3, 20, 9))
    mid = merra2_native_density(w, 20.0, -65.0, 35.0e3, DateTime(2024, 3, 20, 7, 30))
    @test a != b
    @test isapprox(mid, sqrt(a * b); rtol=1.0e-6)     # log-linear in time

    # Ceiling is the native ~80 km, not the 64 km of the pressure-level grid.
    ceiling = merra2_native_ceiling_m(w, 20.0, -65.0)
    @test isapprox(ceiling, 80_000.0; rtol=0.02)
    @test ceiling > merra2_ceiling_m(load_merra2_grid(3, 9), 20.0, -65.0) + 10_000.0
    @test !isfinite(merra2_native_density(w, 20.0, -65.0, ceiling + 5.0e3, DateTime(2024, 3, 20, 6)))

    # Outside the loaded window, and outside the loaded time span.
    @test !isfinite(merra2_native_density(w, 40.0, -65.0, 35.0e3, DateTime(2024, 3, 20, 6)))
    @test !isfinite(merra2_native_density(w, 20.0, -65.0, 35.0e3, DateTime(2024, 3, 21, 12)))

    # Filled (below-terrain) columns do not leak into a neighbour's answer.
    @test isfinite(merra2_native_density(w, 24.5, -55.5, 40.0e3, DateTime(2024, 3, 20, 6)))
    @test !isfinite(merra2_native_density(w, 24.5, -55.5, 600.0, DateTime(2024, 3, 20, 6)))

    # A granule with the wrong vertical size is rejected rather than transposed.
    bad = joinpath(dir, "bad.nc4")
    _write_synthetic_native_granule(bad, day; corrupt_levels=true)
    @test_throws ErrorException load_merra2_native([bad], t0, t1)

    # A missing granule says how to get one.
    err = try
        load_merra2_native([joinpath(dir, "nope.nc4")], t0, t1); nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing && occursin("fetch_merra2_native", err)
end

@testset "MERRA-2 native truth source" begin
    src = Merra2NativeTruth()
    @test truth_source_name(src) == "merra2_native"
    @test is_position_indexed(src)
    @test truth_source_name(truth_source_from_name("merra2_native")) == "merra2_native"
    @test_throws ArgumentError Merra2NativeTruth(; blend_width_km=0.0)
    @test "merra2_native" in TRUTH_SOURCE_NAMES

    # Without granules it fails with an actionable message rather than silently
    # falling back to the climatology.
    case = first(default_study_cases())
    pair = build_trajectory_pair(case)
    empty_dir = mktempdir()
    err = try
        truth_profiles(Merra2NativeTruth(; data_dir=empty_dir), pair, first(pair.aerocapture).dt); nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing && occursin("Earthdata", err)
end

@testset "Prior sources" begin
    @test prior_source_name(GramPrior()) == "gram"
    @test prior_source_name(Nrlmsise00Prior()) == "nrlmsise00"
    @test shares_lineage_with_merra2(GramPrior())
    @test !shares_lineage_with_merra2(Nrlmsise00Prior())
    @test_throws ArgumentError prior_source_from_name("jacchia")
    @test_throws ArgumentError Nrlmsise00Prior(; f107=-1.0)
    @test_throws ArgumentError Nrlmsise00Prior(; sigma_rel=0.0)

    t0 = DateTime(2024, 3, 20, 19, 30)
    pt = TrajectoryPoint(t0, 0.0, 20.0, -65.0, 35.0e3)
    src = Nrlmsise00Prior()
    rho = nrlmsise00_density(src, pt)
    # Radians, not degrees: the radian call reproduces MERRA-2's 8.4e-3 at 35 km.
    @test isapprox(rho, 8.4e-3; rtol=0.05)
    @test rho > nrlmsise00_density(src, TrajectoryPoint(t0, 0.0, 20.0, -65.0, 45.0e3))

    pts = [TrajectoryPoint(t0, 0.0, 20.0, -65.0, 1.0e3 * z) for z in (15.0, 35.0, 55.0)]
    fixed = prior_samples(Nrlmsise00Prior(; sigma_rel=0.05), pts, t0)
    @test all(s -> s.mean_density > 0.0, fixed)
    @test all(s -> isapprox(s.std_density / s.mean_density, 0.05; rtol=1.0e-6), fixed)
    @test [s.mean_density for s in fixed] == [nrlmsise00_density(src, p) for p in pts]
end
