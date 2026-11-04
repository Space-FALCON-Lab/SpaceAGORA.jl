const STUDY_ROOT = @__DIR__
const STUDY_PROJECT = joinpath(STUDY_ROOT, "Project.toml")

if something(Base.active_project(), "") != STUDY_PROJECT
    import Pkg
    Pkg.activate(STUDY_ROOT; io=devnull)
end

using Dates
using Printf
using Random
using Statistics

include(joinpath(STUDY_ROOT, "corridor.jl"))
include(joinpath(STUDY_ROOT, "gram_correlation.jl"))
include(joinpath(STUDY_ROOT, "merra2.jl"))

# Measures the empirical spatial correlation of the MERRA-2 log-density anomaly
# and compares it with the scales Earth-GRAM assumes, to check whether the
# correlation structure the estimator relies on is actually present in the data.
#
# Two ensembles, because a single one cannot answer both questions:
#
# - Diurnal: the eight time-of-day slots per month as departures from that
#   month's all-hours mean, 96 fields. Good for horizontal structure. Not usable
#   vertically: the diurnal signal is dominated by vertically propagating tides
#   with a wavelength of roughly 20-30 km, so the vertical correlation oscillates
#   in sign rather than decaying.
# - Seasonal: the twelve monthly all-hours means as departures from the annual
#   mean, 12 fields. Deep and coherent with no tidal wavelength, so it is the one
#   to use vertically, at the cost of a much smaller ensemble.

const BAND_LAT_MIN = 0.0     # the study corridor sits near 20N
const BAND_LAT_MAX = 45.0
const PROBE_LEVELS = (20, 25, 30, 33, 35, 38, 41)

function _band_indices(grid::Merra2Grid)
    grid_lat = (grid.n_lat - 1) / 180.0
    lo = max(1, floor(Int, (BAND_LAT_MIN + 90.0) * grid_lat) + 1)
    hi = min(grid.n_lat, floor(Int, (BAND_LAT_MAX + 90.0) * grid_lat) + 1)
    return lo:hi
end

_grid_lat_deg(grid::Merra2Grid, i::Int) = (i - 1) * 180.0 / (grid.n_lat - 1) - 90.0
_grid_lon_deg(grid::Merra2Grid, j::Int) = (j - 1) * 360.0 / grid.n_lon

"""
    anomaly_ensemble(:diurnal | :seasonal) -> (members, band, mean_alt_km)

Log-density anomaly fields over the latitude band, with the ensemble mean
removed so each member is a zero-mean departure.
"""
function anomaly_ensemble(kind::Symbol)
    reference = load_merra2_grid(3, 9)
    band = _band_indices(reference)
    members = Array{Float64, 3}[]

    if kind === :diurnal
        for month in 1:12
            base = log.(load_merra2_grid(month, 9).dens[:, band, :])
            for slot in 1:8
                push!(members, log.(load_merra2_grid(month, slot).dens[:, band, :]) .- base)
            end
        end
    elseif kind === :seasonal
        monthly = [log.(load_merra2_grid(m, 9).dens[:, band, :]) for m in 1:12]
        annual = sum(monthly) ./ 12.0
        members = [m .- annual for m in monthly]
    else
        throw(ArgumentError("Unsupported ensemble $kind. Use :diurnal or :seasonal."))
    end

    ensemble_mean = sum(members) ./ length(members)
    members = [m .- ensemble_mean for m in members]
    mean_alt_km = [mean(filter(isfinite, reference.hgt[k, band, :])) * 1e-3 for k in 1:reference.n_pres]
    return members, band, mean_alt_km
end

@inline function _standardize(values::Vector{Float64})
    all(isfinite, values) || return nothing
    s = std(values)
    s > 0.0 || return nothing
    m = mean(values)
    return (values .- m) ./ s
end

"""
    horizontal_efolding_km(members, grid, band, level; n_base, max_km, bin_km) -> Float64

The `e`-folding scale of an exponential fitted to the binned correlation of the
anomaly against chordal separation at one pressure level.
"""
function horizontal_efolding_km(
    members::Vector{Array{Float64, 3}},
    grid::Merra2Grid,
    band,
    level::Int;
    n_base::Int=200,
    max_km::Float64=6000.0,
    bin_km::Float64=250.0,
    seed::Int=0,
)::Float64
    n_lat = length(band)
    n_lon = grid.n_lon
    z = Array{Union{Nothing, Vector{Float64}}, 2}(nothing, n_lat, n_lon)
    for i in 1:n_lat, j in 1:n_lon
        z[i, j] = _standardize([m[level, i, j] for m in members])
    end
    usable = [(i, j) for i in 1:n_lat, j in 1:n_lon if z[i, j] !== nothing]
    isempty(usable) && return NaN

    n_bins = ceil(Int, max_km / bin_km)
    num = zeros(n_bins)
    den = zeros(Int, n_bins)
    rng = MersenneTwister(seed)
    for _ in 1:n_base
        i0, j0 = usable[rand(rng, 1:length(usable))]
        lat0 = _grid_lat_deg(grid, band[i0])
        lon0 = _grid_lon_deg(grid, j0)
        z0 = z[i0, j0]
        for (i, j) in usable
            d = chordal_distance_km(lat0, lon0, _grid_lat_deg(grid, band[i]), _grid_lon_deg(grid, j))
            b = floor(Int, d / bin_km) + 1
            (1 <= b <= n_bins) || continue
            num[b] += mean(z0 .* z[i, j])
            den[b] += 1
        end
    end

    centers = Float64[bin_km * (b - 0.5) for b in 1:n_bins]
    corr = [den[b] > 0 ? num[b] / den[b] : NaN for b in 1:n_bins]
    keep = findall(b -> isfinite(corr[b]) && corr[b] > 0.05 && centers[b] < 4000.0, 1:n_bins)
    length(keep) >= 4 || return NaN
    x = centers[keep]
    y = log.(corr[keep])
    slope = (mean(x .* y) - mean(x) * mean(y)) / (mean(x .^ 2) - mean(x)^2)
    return slope < 0.0 ? -1.0 / slope : NaN
end

"""
    vertical_correlation(members, level) -> Vector{Float64}

Correlation of the anomaly at every pressure level against a reference level, in
the same grid column, averaged over the band.
"""
function vertical_correlation(members::Vector{Array{Float64, 3}}, level::Int)::Vector{Float64}
    n_pres, n_lat, n_lon = size(first(members))
    out = zeros(n_pres)
    counts = zeros(Int, n_pres)
    for i in 1:n_lat, j in 1:n_lon
        zref = _standardize([m[level, i, j] for m in members])
        zref === nothing && continue
        for k in 1:n_pres
            zk = _standardize([m[k, i, j] for m in members])
            zk === nothing && continue
            out[k] += mean(zref .* zk)
            counts[k] += 1
        end
    end
    return [counts[k] > 0 ? out[k] / counts[k] : NaN for k in 1:n_pres]
end

function main()
    grid = load_merra2_grid(3, 9)
    diurnal, band, alt_km = anomaly_ensemble(:diurnal)
    seasonal, _, _ = anomaly_ensemble(:seasonal)
    @printf("band %.0fN-%.0fN, diurnal ensemble %d members, seasonal %d\n\n",
            BAND_LAT_MIN, BAND_LAT_MAX, length(diurnal), length(seasonal))

    println("=== Horizontal: MERRA-2 e-folding vs GRAM's tabulated Lh ===")
    println(" alt_km    MERRA-2_km    GRAM_km    ratio")
    for k in PROBE_LEVELS
        measured = horizontal_efolding_km(diurnal, grid, band, k)
        lh, _ = gram_correlation_scales(alt_km[k])
        @printf("  %5.1f    %9.0f    %7.0f    %5.2f\n", alt_km[k], measured, lh, measured / lh)
    end

    println("\n=== Vertical: seasonal ensemble about 30.8 km ===")
    corr = vertical_correlation(seasonal, 30)
    _, lz = gram_correlation_scales(alt_km[30])
    println(" dz_km   MERRA-2   GRAM exp(-|dz|/Lz)")
    for k in 20:2:40
        dz = alt_km[k] - alt_km[30]
        @printf("  %+5.0f    %6.2f    %6.2f\n", dz, corr[k], exp(-abs(dz) / lz))
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
