# Truth sources for the study.
#
# The estimator can only learn a residual field that is actually a function of
# position. GRAM's `perturbedDensity` is not: `PerturbedAtmosphere::setPosition`
# measures displacement since the previous call on that instance and
# `updatePerturbations` accumulates it, so the returned density depends on the
# instance's call history (see `data/GRAMSuite.jl/src/GRAMSuite.jl` on
# `_gram_wind_mode`). Measured on this study's own corridors, walking one
# instance along the same path twice differs by up to 44.5%, and the
# aerocapture and EDL residuals at matched altitude correlate at -0.13.
#
# GRAM *nominal* density is a deterministic function of (lat, lon, alt, t)
# (repeatability measured at exactly 0), so every sound truth source here is
# built on the nominal field.

abstract type TruthSource end

"""
    GramWalkTruth(seed, planet_name)

Legacy truth: GRAM's `perturbedDensity` walked along the trajectory. Retained
only to reproduce the pre-fix results. The residual it produces is not a
spatial field and is not learnable by any position-indexed estimator.
"""
struct GramWalkTruth <: TruthSource
    seed::Int
    planet_name::String
end
GramWalkTruth(; seed::Int=10_042, planet_name::String="earth") = GramWalkTruth(seed, planet_name)

"""
    GramEpochShiftTruth(epoch_shift_s, planet_name)

Truth is the GRAM nominal field at an epoch shifted by `epoch_shift_s` from the
epoch the onboard prior assumes. Deterministic in position, and the
truth-minus-prior residual is genuine day-to-day atmospheric variability that
the prior has no way to know. This is the realistic same-model case.
"""
struct GramEpochShiftTruth <: TruthSource
    epoch_shift_s::Float64
    planet_name::String
end
GramEpochShiftTruth(; epoch_shift_s::Float64=36.0 * 3600.0, planet_name::String="earth") =
    GramEpochShiftTruth(epoch_shift_s, planet_name)

"""
    SyntheticFieldTruth(...)

Truth is GRAM nominal multiplied by `exp(bias + amplitude * f(lat, lon, alt))`,
where `f` is an exact draw from a zero-mean, unit-variance squared-exponential
GP with the given length scales, realized through random Fourier features so it
evaluates as a closed-form function of position.

This is the validation configuration: the residual field is a known,
learnable, spatially smooth signal with known hyperparameters, so a correct
estimator must recover it. `bias` plants a constant log-density offset — the
component a stationary zero-mean GP structurally cannot extrapolate and that
the mean-function basis in `gp_models.jl` exists to capture.
"""
struct SyntheticFieldTruth <: TruthSource
    seed::Int
    bias::Float64
    amplitude::Float64
    lat_scale_deg::Float64
    lon_scale_deg::Float64
    alt_scale_km::Float64
    n_features::Int
    planet_name::String
end

function SyntheticFieldTruth(;
    seed::Int=10_042,
    bias::Float64=0.08,
    amplitude::Float64=0.10,
    lat_scale_deg::Float64=6.0,
    lon_scale_deg::Float64=12.0,
    alt_scale_km::Float64=25.0,
    n_features::Int=256,
    planet_name::String="earth",
)
    n_features >= 1 || throw(ArgumentError("SyntheticFieldTruth needs at least one Fourier feature."))
    amplitude >= 0.0 || throw(ArgumentError("SyntheticFieldTruth amplitude must be nonnegative."))
    return SyntheticFieldTruth(
        seed, bias, amplitude, lat_scale_deg, lon_scale_deg, alt_scale_km, n_features, planet_name
    )
end

"""
    Merra2Truth(...)

Truth is the MERRA-2 field vendored with Earth-GRAM, read directly by
[`merra2.jl`](merra2.jl) and applied as a multiplicative anomaly on the GRAM
nominal profile:

    rho_truth = rho_gram * exp(taper * (log(rho_m2 / rho_gram) + dispersion * sigma_m2 * g))

`rho_m2` is the monthly climatology for the epoch's month and time-of-day slot,
`sigma_m2` is MERRA-2's own relative interannual density standard deviation at
that point, and `g` is a smooth unit-variance random field. The anomaly is
deterministic in position; only its `dispersion` component is random, and that
is drawn once per `seed` rather than per call.

`dispersion = 0` gives the raw MERRA-2 field, which measures how GRAM's nominal
departs from the reanalysis it is built from. `dispersion = 1` scales a smooth
anomaly by MERRA-2's observed year-to-year spread, which is the realistic
"specific day" case an EDL study needs: the vendored grids are climatologies,
so a single slot is not an independent weather day. Measured on the default
corridor, the 18Z-minus-all-hours anomaly alone is only `0.05-0.5%` in
log-density, well under the `5%` measurement noise, while `sigma_m2` is
`0.7-1.6%`.

MERRA-2 spans 1000 mb to 0.1 mb, roughly `0.1 km` to `64 km`. Above the local
ceiling the anomaly is held at its ceiling value and tapered to zero over
`blend_width_km`, so truth reverts continuously to GRAM nominal.

With `time_decorrelation` the anomaly evolves as a stationary AR(1) in time
using GRAM's own temporal scale:

    a(x, t) = rho_t * a(x, 0) + sqrt(1 - rho_t^2) * sigma_a * g_innov(x),
    rho_t   = exp(-t / Lt(h)),   Lt from `gram_time_scale_s`

so `Corr(a(x,0), a(x,t)) = rho_t` and the marginal variance is preserved. The
innovation field `g_innov` is an independent unit-variance draw carrying GRAM's
*spatial* scales, `sigma_a` is the RMS of the base anomaly over the corridor, and
`t` runs from the anchor epoch, so the aerocapture pass samples an essentially
undecayed field while the EDL pass samples it `gap_s` later.

Without this the anomaly is frozen: the field the estimator learns during
aerocapture is bit-identical to the field it predicts hours later, which is
information a real flight does not have. `Lt` is about `26.6 hr` at `35 km`, so
a `6 hr` gap retains `rho_t = 0.80` and injects `0.61` of an independent field.
"""
struct Merra2Truth <: TruthSource
    hour_code::Int          # 0 derives the slot from the epoch hour; 9 is the all-hours mean
    dispersion::Float64
    seed::Int
    lat_scale_deg::Float64
    lon_scale_deg::Float64
    alt_scale_km::Float64
    blend_width_km::Float64
    time_decorrelation::Bool
    planet_name::String
end

function Merra2Truth(;
    hour_code::Int=0,
    dispersion::Float64=1.0,
    seed::Int=10_042,
    lat_scale_deg::Float64=6.0,
    lon_scale_deg::Float64=12.0,
    alt_scale_km::Float64=15.0,
    blend_width_km::Float64=10.0,
    time_decorrelation::Bool=true,
    planet_name::String="earth",
)
    (hour_code == 0 || 1 <= hour_code <= 9) ||
        throw(ArgumentError("MERRA-2 hour code must be 0 (auto) or 1..9, got $hour_code."))
    dispersion >= 0.0 || throw(ArgumentError("Merra2Truth dispersion must be nonnegative."))
    blend_width_km > 0.0 || throw(ArgumentError("Merra2Truth blend width must be positive."))
    return Merra2Truth(
        hour_code, dispersion, seed, lat_scale_deg, lon_scale_deg, alt_scale_km,
        blend_width_km, time_decorrelation, planet_name,
    )
end

"""
    Merra2NativeTruth(...)

Truth from day-specific MERRA-2 native model-level reanalysis, applied as a
multiplicative anomaly on GRAM nominal exactly as [`Merra2Truth`](@ref) does:

    rho_truth = rho_gram * exp(taper * log(rho_m2 / rho_gram))

The difference is what `rho_m2` is. [`Merra2Truth`](@ref) reads the monthly
time-of-day climatologies vendored with GRAM, which top out near `64 km` and
contain no weather, so it needs a synthetic dispersion term for realistic
amplitude and an AR(1) surrogate for time evolution. This source reads the real
`3-hourly` analysis on the `72` native model levels: the anomaly is genuine
weather, it evolves genuinely across the aerocapture-to-EDL gap, and the ceiling
moves from about `64 km` to about `80 km`, which puts the GRAM/MERRA-2 blending
seam well above the scored band.

Neither `dispersion` nor `time_decorrelation` exists here because neither is
needed.

Granules are not vendored. Fetch them with `scripts/fetch_merra2_native.sh`
(a NASA Earthdata login is required) or point `SPACEAGORA_MERRA2_NATIVE_PATH` at
an existing collection.
"""
struct Merra2NativeTruth <: TruthSource
    data_dir::String
    anomaly_top_km::Float64
    blend_width_km::Float64
    lat_margin_deg::Float64
    lon_margin_deg::Float64
    planet_name::String
end

function Merra2NativeTruth(;
    data_dir::String="",
    anomaly_top_km::Float64=52.0,
    blend_width_km::Float64=10.0,
    lat_margin_deg::Float64=6.0,
    lon_margin_deg::Float64=6.0,
    planet_name::String="earth",
)
    blend_width_km > 0.0 || throw(ArgumentError("Merra2NativeTruth blend width must be positive."))
    anomaly_top_km > 0.0 || throw(ArgumentError("Merra2NativeTruth anomaly top must be positive."))
    return Merra2NativeTruth(data_dir, anomaly_top_km, blend_width_km, lat_margin_deg, lon_margin_deg, planet_name)
end

const TRUTH_SOURCE_NAMES = ("merra2_native", "merra2", "gram_walk", "gram_epoch_shift", "synthetic_field")

"""
    truth_source_from_name(name; kwargs...) -> TruthSource

Build a truth source from its CLI name.
"""
function truth_source_from_name(
    name::AbstractString;
    seed::Int=10_042,
    planet_name::String="earth",
    epoch_shift_s::Float64=36.0 * 3600.0,
    merra2_hour_code::Int=0,
    merra2_dispersion::Float64=1.0,
    merra2_blend_width_km::Float64=10.0,
    merra2_time_decorrelation::Bool=true,
    merra2_native_dir::String="",
    merra2_native_anomaly_top_km::Float64=52.0,
    field_bias::Float64=0.08,
    field_amplitude::Float64=0.10,
    field_lat_scale_deg::Float64=6.0,
    field_lon_scale_deg::Float64=12.0,
    field_alt_scale_km::Float64=25.0,
)::TruthSource
    if name == "merra2_native"
        return Merra2NativeTruth(;
            data_dir=merra2_native_dir,
            anomaly_top_km=merra2_native_anomaly_top_km,
            blend_width_km=merra2_blend_width_km,
            planet_name,
        )
    elseif name == "merra2"
        return Merra2Truth(;
            hour_code=merra2_hour_code,
            dispersion=merra2_dispersion,
            seed,
            lat_scale_deg=field_lat_scale_deg,
            lon_scale_deg=field_lon_scale_deg,
            alt_scale_km=field_alt_scale_km,
            blend_width_km=merra2_blend_width_km,
            time_decorrelation=merra2_time_decorrelation,
            planet_name,
        )
    elseif name == "gram_walk"
        return GramWalkTruth(; seed, planet_name)
    elseif name == "gram_epoch_shift"
        return GramEpochShiftTruth(; epoch_shift_s, planet_name)
    elseif name == "synthetic_field"
        return SyntheticFieldTruth(;
            seed,
            bias=field_bias,
            amplitude=field_amplitude,
            lat_scale_deg=field_lat_scale_deg,
            lon_scale_deg=field_lon_scale_deg,
            alt_scale_km=field_alt_scale_km,
            planet_name,
        )
    end
    throw(ArgumentError("Unsupported truth source '$name'. Use one of: $(join(TRUTH_SOURCE_NAMES, ", "))."))
end

truth_source_name(::Merra2NativeTruth) = "merra2_native"
truth_source_name(::Merra2Truth) = "merra2"
truth_source_name(::GramWalkTruth) = "gram_walk"
truth_source_name(::GramEpochShiftTruth) = "gram_epoch_shift"
truth_source_name(::SyntheticFieldTruth) = "synthetic_field"

"""
    is_position_indexed(src) -> Bool

Whether the truth residual this source produces is a function of position, and
therefore learnable in principle. False only for [`GramWalkTruth`](@ref).
"""
is_position_indexed(::Merra2NativeTruth) = true
is_position_indexed(::Merra2Truth) = true
is_position_indexed(::GramWalkTruth) = false
is_position_indexed(::GramEpochShiftTruth) = true
is_position_indexed(::SyntheticFieldTruth) = true

# --- Random Fourier feature field -------------------------------------------

struct FourierField
    freqs::Matrix{Float64}   # 3 x M, rows are lat/lon/alt in scaled units
    phases::Vector{Float64}  # M
    scale::NTuple{3, Float64}
    norm::Float64
end

function _build_fourier_field(
    seed::Int, n_features::Int, lat_scale_deg::Float64, lon_scale_deg::Float64, alt_scale_km::Float64
)::FourierField
    rng = MersenneTwister(seed)
    return FourierField(
        randn(rng, 3, n_features),
        2.0 * pi .* rand(rng, n_features),
        (lat_scale_deg, lon_scale_deg, alt_scale_km),
        sqrt(2.0 / n_features),
    )
end

_build_fourier_field(src::SyntheticFieldTruth) = _build_fourier_field(
    src.seed, src.n_features, src.lat_scale_deg, src.lon_scale_deg, src.alt_scale_km
)

_build_fourier_field(src::Merra2Truth) = _build_fourier_field(
    src.seed, 256, src.lat_scale_deg, src.lon_scale_deg, src.alt_scale_km
)

const KM_PER_DEG_LAT = 111.2

"""
    _build_innovation_field(src, ref_lat_deg, ref_alt_km) -> FourierField

The independent field mixed in by the AR(1) time decorrelation. It carries
GRAM's *spatial* correlation scales at `ref_alt_km`, converted to degrees at
`ref_lat_deg`, so truth uses GRAM's space and time correlation model
consistently. Seeded off `src.seed` by a fixed offset so it is independent of
the dispersion field but still reproducible.
"""
function _build_innovation_field(src::Merra2Truth, ref_lat_deg::Float64, ref_alt_km::Float64)::FourierField
    lh_km, lz_km = gram_correlation_scales(ref_alt_km)
    lat_scale = lh_km / KM_PER_DEG_LAT
    lon_scale = lh_km / (KM_PER_DEG_LAT * max(cosd(ref_lat_deg), 0.1))
    return _build_fourier_field(src.seed + 7919, 256, lat_scale, lon_scale, lz_km)
end

function _evaluate_field(field::FourierField, point::TrajectoryPoint)::Float64
    u1 = point.lat_deg / field.scale[1]
    u2 = point.lon_deg / field.scale[2]
    u3 = (point.alt_m * 1e-3) / field.scale[3]
    acc = 0.0
    @inbounds for j in eachindex(field.phases)
        acc += cos(field.freqs[1, j] * u1 + field.freqs[2, j] * u2 + field.freqs[3, j] * u3 + field.phases[j])
    end
    return field.norm * acc
end

"""
    synthetic_log_anomaly(src, point) -> Float64

The planted log-density anomaly at `point`, for tests and for reporting the
signal a run is asking the estimator to recover.
"""
function synthetic_log_anomaly(src::SyntheticFieldTruth, point::TrajectoryPoint)::Float64
    return src.bias + src.amplitude * _evaluate_field(_build_fourier_field(src), point)
end

# --- Truth evaluation --------------------------------------------------------

"""
    truth_profiles(src, pair, initial_dt) -> (aerocapture_truth, edl_truth)

Truth density along both corridors, in kg/m^3.
"""
function truth_profiles(
    src::GramWalkTruth,
    pair::TrajectoryPair,
    initial_dt::DateTime,
)::Tuple{Vector{Float64}, Vector{Float64}}
    return gram_truth_profiles(pair, initial_dt; planet_name=src.planet_name, seed=src.seed)
end

function truth_profiles(
    src::GramEpochShiftTruth,
    pair::TrajectoryPair,
    initial_dt::DateTime,
)::Tuple{Vector{Float64}, Vector{Float64}}
    shifted_dt = initial_dt + Millisecond(round(Int, 1000 * src.epoch_shift_s))
    model = _build_gram_model(_to_initial_time(shifted_dt), src.planet_name, 1)
    evaluate(points) = [
        _gram_density_at_point(model, p; elapsed_time_s=_elapsed_from_initial_s(p, initial_dt))
        for p in points
    ]
    return evaluate(pair.aerocapture), evaluate(pair.edl)
end

function truth_profiles(
    src::SyntheticFieldTruth,
    pair::TrajectoryPair,
    initial_dt::DateTime,
)::Tuple{Vector{Float64}, Vector{Float64}}
    field = _build_fourier_field(src)
    model = _build_gram_model(_to_initial_time(initial_dt), src.planet_name, 1)
    function evaluate(points)
        out = Vector{Float64}(undef, length(points))
        @inbounds for i in eachindex(points)
            rho = _gram_density_at_point(model, points[i]; elapsed_time_s=_elapsed_from_initial_s(points[i], initial_dt))
            out[i] = rho * exp(src.bias + src.amplitude * _evaluate_field(field, points[i]))
        end
        return out
    end
    return evaluate(pair.aerocapture), evaluate(pair.edl)
end

@inline function _smoothstep(x::Float64)::Float64
    t = clamp(x, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)
end

"""
    merra2_anomaly(src, grid, model, point, initial_dt) -> (anomaly, taper, inside)

The log-density anomaly MERRA-2 implies at `point`, the taper applied to it, and
whether the point sits inside MERRA-2's vertical domain. Split out so tests and
diagnostics can inspect the field without running a whole case.
"""
function merra2_anomaly(
    src::Merra2Truth,
    grid::Merra2Grid,
    field::FourierField,
    gram_model,
    point::TrajectoryPoint,
    initial_dt::DateTime,
)::Tuple{Float64, Float64, Bool}
    ceiling_m = merra2_ceiling_m(grid, point.lat_deg, point.lon_deg)
    isfinite(ceiling_m) || return (0.0, 0.0, false)
    top_m = ceiling_m - 1.0e3
    inside = point.alt_m <= top_m
    alt_eval = min(point.alt_m, top_m)

    rho_m2, rel_sigma = merra2_density(grid, point.lat_deg, point.lon_deg, alt_eval)
    isfinite(rho_m2) || return (0.0, 0.0, false)

    eval_point = inside ? point : TrajectoryPoint(
        point.dt, point.elapsed_s, point.lat_deg, point.lon_deg, alt_eval
    )
    rho_gram = _gram_density_at_point(
        gram_model, eval_point; elapsed_time_s=_elapsed_from_initial_s(eval_point, initial_dt)
    )
    rho_gram > 0.0 || return (0.0, 0.0, false)

    anomaly = log(rho_m2 / rho_gram)
    if src.dispersion > 0.0
        anomaly += src.dispersion * rel_sigma * _evaluate_field(field, eval_point)
    end
    taper = inside ? 1.0 :
        1.0 - _smoothstep((point.alt_m - top_m) * 1.0e-3 / src.blend_width_km)
    return anomaly, taper, inside
end

function truth_profiles(
    src::Merra2Truth,
    pair::TrajectoryPair,
    initial_dt::DateTime,
)::Tuple{Vector{Float64}, Vector{Float64}}
    hour_code = src.hour_code == 0 ? merra2_hour_code(Dates.hour(initial_dt)) : src.hour_code
    grid = load_merra2_grid(Dates.month(initial_dt), hour_code)
    field = _build_fourier_field(src)
    model = _build_gram_model(_to_initial_time(initial_dt), src.planet_name, 1)

    points = vcat(pair.aerocapture, pair.edl)
    n = length(points)
    n_aero = length(pair.aerocapture)
    rho_gram = Vector{Float64}(undef, n)
    anomaly = Vector{Float64}(undef, n)
    taper = Vector{Float64}(undef, n)
    inside = falses(n)
    @inbounds for i in 1:n
        p = points[i]
        rho_gram[i] = _gram_density_at_point(
            model, p; elapsed_time_s=_elapsed_from_initial_s(p, initial_dt)
        )
        anomaly[i], taper[i], inside[i] = merra2_anomaly(src, grid, field, model, p, initial_dt)
    end

    if src.time_decorrelation
        # Scale the innovation to the base anomaly's own spread, measured where
        # MERRA-2 actually has data so the tapered upper corridor does not
        # deflate it. Measured over the *prediction* corridor only: sampling the
        # aerocapture leg as well would let the sensing geometry change the truth
        # field, so sweeping `--aerocapture-exit-alt-km` would no longer be a
        # controlled comparison.
        edl_range = (n_aero + 1):n
        sampled = anomaly[edl_range][inside[edl_range]]
        sigma_a = isempty(sampled) ? 0.0 : sqrt(mean(abs2, sampled))
        ref_lat = mean(p.lat_deg for p in pair.edl)
        ref_alt = mean(p.alt_m for p in pair.edl) * 1e-3
        innovation = _build_innovation_field(src, ref_lat, ref_alt)
        @inbounds for i in 1:n
            p = points[i]
            dt_s = _elapsed_from_initial_s(p, initial_dt)
            rho_t = exp(-dt_s / gram_time_scale_s(p.alt_m * 1e-3))
            anomaly[i] = rho_t * anomaly[i] +
                sqrt(max(1.0 - rho_t^2, 0.0)) * sigma_a * _evaluate_field(innovation, p)
        end
    end

    out = @. rho_gram * exp(taper * anomaly)
    return out[1:n_aero], out[(n_aero + 1):end]
end

"""
    truth_profiles(src::Merra2NativeTruth, pair, initial_dt)

One window covering both corridors and the gap between them is loaded once, then
every point is evaluated at its own timestamp so the anomaly evolves with the
real analysis.
"""
function truth_profiles(
    src::Merra2NativeTruth,
    pair::TrajectoryPair,
    initial_dt::DateTime,
)::Tuple{Vector{Float64}, Vector{Float64}}
    points = vcat(pair.aerocapture, pair.edl)
    n_aero = length(pair.aerocapture)
    lats = [p.lat_deg for p in points]
    lons = [p.lon_deg for p in points]
    times = [p.dt for p in points]

    window = load_merra2_native_span(
        minimum(times), maximum(times);
        dir=isempty(src.data_dir) ? merra2_native_dir() : src.data_dir,
        lat_min=minimum(lats) - src.lat_margin_deg,
        lat_max=maximum(lats) + src.lat_margin_deg,
        lon_min=minimum(lons) - src.lon_margin_deg,
        lon_max=maximum(lons) + src.lon_margin_deg,
    )
    model = _build_gram_model(_to_initial_time(initial_dt), src.planet_name, 1)

    out = Vector{Float64}(undef, length(points))
    @inbounds for i in eachindex(points)
        p = points[i]
        rho_gram = _gram_density_at_point(
            model, p; elapsed_time_s=_elapsed_from_initial_s(p, initial_dt)
        )
        ceiling_m = merra2_native_ceiling_m(window, p.lat_deg, p.lon_deg)
        if !isfinite(ceiling_m)
            out[i] = rho_gram
            continue
        end
        # Cap the anomaly below GRAM's own model-blending transition rather than
        # at MERRA-2's data ceiling. Measured on 2024-03-20 over the corridor,
        # log(rho_m2 / rho_gram) sits at +0.03..+0.06 through 20-55 km and then
        # breaks to -0.12..-0.14 by 65-75 km: that step is GRAM handing off from
        # its MERRA-2-based lower atmosphere to its upper-atmosphere model, not
        # weather, and it is what dominated the weighted metrics. Truth therefore
        # carries the reanalysis anomaly only where the two models are genuinely
        # comparable and reverts to GRAM nominal above.
        top_m = min(ceiling_m - 1.0e3, 1.0e3 * src.anomaly_top_km)
        inside = p.alt_m <= top_m
        alt_eval = min(p.alt_m, top_m)
        rho_m2 = merra2_native_density(window, p.lat_deg, p.lon_deg, alt_eval, p.dt)
        if !isfinite(rho_m2) || rho_gram <= 0.0
            out[i] = rho_gram
            continue
        end

        anomaly = if inside
            log(rho_m2 / rho_gram)
        else
            eval_point = TrajectoryPoint(p.dt, p.elapsed_s, p.lat_deg, p.lon_deg, alt_eval)
            rho_gram_eval = _gram_density_at_point(
                model, eval_point; elapsed_time_s=_elapsed_from_initial_s(eval_point, initial_dt)
            )
            rho_gram_eval > 0.0 ? log(rho_m2 / rho_gram_eval) : 0.0
        end
        taper = inside ? 1.0 :
            1.0 - _smoothstep((p.alt_m - top_m) * 1.0e-3 / src.blend_width_km)
        out[i] = rho_gram * exp(taper * anomaly)
    end
    return out[1:n_aero], out[(n_aero + 1):end]
end
