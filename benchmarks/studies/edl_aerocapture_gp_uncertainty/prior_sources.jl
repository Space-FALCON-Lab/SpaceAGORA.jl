# Onboard prior sources.
#
# The prior is a mean profile plus a per-point spread; the GP learns a residual
# against the mean and `prior_scaled` folds the spread into its kernel.
#
# Why this is swappable. GRAM's nominal below about 65 km *is* MERRA-2 — the
# vendored climatology is the MERRA-2 pressure-level product monthly-averaged.
# With `merra2_native` truth, prior and truth therefore share a lineage, and the
# residual the estimator learns is "a specific MERRA-2 day minus the MERRA-2
# monthly climatology". That is a real, meaningful signal, but it is not an
# independent model-versus-observation comparison, and the under-dispersion
# result it produces is a self-consistency statement about GRAM rather than an
# external validation.
#
# NRLMSISE-00 is fitted to a different observational record and is not
# MERRA-2-derived, so swapping it in breaks that lineage. Measured against
# MERRA-2 on 2024-03-20 over this corridor, it is not the weak stratospheric
# prior one might expect: RMS log-anomaly 0.0724 against GRAM's 0.0701 overall,
# and *smaller* through the scored band (0.011-0.031 at 25-45 km against GRAM's
# 0.032-0.048).
#
# Note this is NRLMSISE-00 (2000), not NRLMSIS 2.0/2.1. The 2.x line specifically
# overhauled the 0-100 km region; if a Julia binding for it appears, it is the
# better choice here and drops in behind this same interface.

using SatelliteToolboxAtmosphericModels
const ATMOS_MODELS = SatelliteToolboxAtmosphericModels.AtmosphericModels

abstract type PriorSource end

"""
    GramPrior()

GRAM nominal as the mean, and the standard deviation of `n_dispersion` seeded
GRAM `perturbedDensity` realizations as the spread. The original behaviour.
"""
struct GramPrior <: PriorSource end

"""
    Nrlmsise00Prior(; f107a, f107, ap, sigma_rel)

NRLMSISE-00 as the prior mean.

`sigma_rel` sets the relative spread. `NaN`, the default, borrows GRAM's
dispersed-ensemble sigma so that swapping the prior changes *only the mean* and
the comparison against a `GramPrior` run isolates the lineage question. Give it
a number to use a constant relative sigma instead, which decouples the spread
from GRAM entirely at the cost of no longer being an ensemble.

Solar indices default to moderate fixed values. Below about 60 km the solar
terms are negligible, so this matters little for the scored band; set them from
real indices for a run that reaches into the thermosphere.
"""
struct Nrlmsise00Prior <: PriorSource
    f107a::Float64
    f107::Float64
    ap::Float64
    sigma_rel::Float64
end

function Nrlmsise00Prior(; f107a::Float64=150.0, f107::Float64=150.0, ap::Float64=4.0, sigma_rel::Float64=NaN)
    (f107a > 0.0 && f107 > 0.0 && ap >= 0.0) ||
        throw(ArgumentError("NRLMSISE-00 solar indices must be positive (ap nonnegative)."))
    (isnan(sigma_rel) || sigma_rel > 0.0) ||
        throw(ArgumentError("prior sigma_rel must be positive, or NaN to borrow GRAM's ensemble."))
    return Nrlmsise00Prior(f107a, f107, ap, sigma_rel)
end

const PRIOR_SOURCE_NAMES = ("gram", "nrlmsise00")

function prior_source_from_name(
    name::AbstractString;
    f107a::Float64=150.0, f107::Float64=150.0, ap::Float64=4.0, sigma_rel::Float64=NaN,
)::PriorSource
    name == "gram" && return GramPrior()
    name == "nrlmsise00" && return Nrlmsise00Prior(; f107a, f107, ap, sigma_rel)
    throw(ArgumentError("Unsupported prior source '$name'. Use one of: $(join(PRIOR_SOURCE_NAMES, ", "))."))
end

prior_source_name(::GramPrior) = "gram"
prior_source_name(::Nrlmsise00Prior) = "nrlmsise00"

"""
    shares_lineage_with_merra2(src) -> Bool

Whether this prior is derived from MERRA-2, and therefore not independent of a
`merra2` or `merra2_native` truth source.
"""
shares_lineage_with_merra2(::GramPrior) = true
shares_lineage_with_merra2(::Nrlmsise00Prior) = false

"""
    nrlmsise00_density(src, point) -> Float64

NRLMSISE-00 total mass density in kg/m^3. Latitude and longitude go in as
**radians** — verified against MERRA-2, which the radian call reproduces to
`8.4055e-3` against `8.4e-3` at `35 km` where the degree call gives `7.72e-3`.
"""
@inline function nrlmsise00_density(src::Nrlmsise00Prior, point::TrajectoryPoint)::Float64
    state = ATMOS_MODELS.nrlmsise00(
        point.dt, point.alt_m, deg2rad(point.lat_deg), deg2rad(point.lon_deg),
        src.f107a, src.f107, src.ap,
    )
    return Float64(state.total_density)
end

"""
    prior_samples(src, points, initial_dt; planet_name, prior_seed, n_dispersion) -> Vector{GramSample}

Prior mean and spread at each point.
"""
function prior_samples(
    ::GramPrior, points::Vector{TrajectoryPoint}, initial_dt::DateTime;
    planet_name::String="earth", prior_seed::Int=20_042, n_dispersion::Int=24,
)::Vector{GramSample}
    return gram_prior_samples(points, initial_dt; planet_name, prior_seed, n_dispersion)
end

function prior_samples(
    src::Nrlmsise00Prior, points::Vector{TrajectoryPoint}, initial_dt::DateTime;
    planet_name::String="earth", prior_seed::Int=20_042, n_dispersion::Int=24,
)::Vector{GramSample}
    # Borrowing GRAM's ensemble spread costs a full dispersed sweep, but keeps
    # the only difference from a GramPrior run in the mean.
    borrowed = isnan(src.sigma_rel) ?
        gram_prior_samples(points, initial_dt; planet_name, prior_seed, n_dispersion) :
        GramSample[]

    out = Vector{GramSample}(undef, length(points))
    @inbounds for i in eachindex(points)
        mu = nrlmsise00_density(src, points[i])
        mu > 0.0 && isfinite(mu) ||
            throw(ErrorException("NRLMSISE-00 returned $mu at $(points[i].alt_m * 1e-3) km."))
        rel = isnan(src.sigma_rel) ?
            borrowed[i].std_density / max(borrowed[i].mean_density, 1.0e-18) : src.sigma_rel
        out[i] = GramSample(mu, max(mu * rel, mu * 1.0e-6))
    end
    return out
end
