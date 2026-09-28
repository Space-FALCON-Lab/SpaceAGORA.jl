# Earth-GRAM's own tabulated spatial correlation scales.
#
# Transcribed from `Earth/source/EarthAtmosphereData.cpp::rsData`, whose columns
# are `xlbar, xsigl, xlmin, xscale, zlbar, zsigl, zlmin, zscale, wr`;
# `initializeData` maps `xscale -> xScaleSmall` (horizontal, km) and
# `zscale -> zScaleSmall` (vertical, km).
#
# `EarthAtmosphere::getCorrelationCoefficients` combines them as
#
#     rho = exp(-( dx / Lh + dz / Lz + dt / Lt ))
#
# a *separable exponential* on an L1 metric, not the Euclidean squared-exponential
# or Matern form the study's other kernels use.
#
# Measured against the vendored MERRA-2 grids over 0-45N (see README), the
# horizontal scales are confirmed to within a few percent between 18 km and
# 53 km: MERRA-2 gives 2270 km at 30.8 km where GRAM tabulates 2240 km.

const GRAM_RS_ALT_KM = [
    0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,
    65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0,
    140.0, 160.0, 180.0, 200.0,
]

const GRAM_RS_LH_KM = [
    117.5, 819.9, 820.8, 1397.2, 2040.4, 2127.5, 2240.3, 2279.4, 2308.6, 2318.4,
    2300.4, 2240.5, 2164.0, 2107.6, 2038.9, 1951.2, 2043.0, 2140.8, 2239.2,
    2515.5, 2791.8, 3068.8, 3345.9, 3622.9, 3900.0, 4189.2, 4482.0, 4770.0, 5040.0,
]

const GRAM_RS_LZ_KM = [
    2.65, 10.56, 12.32, 14.36, 16.01, 16.16, 17.44, 18.33, 18.78, 19.09, 19.32,
    19.39, 19.32, 19.08, 18.63, 18.00, 19.04, 20.08, 22.04, 23.99, 25.94, 27.44,
    28.95, 30.45, 31.96, 35.92, 40.00, 44.04, 48.06,
]

const EARTH_MEAN_RADIUS_KM = 6371.0

"""
    gram_time_scale_s(alt_km) -> Float64

GRAM's temporal correlation scale, transcribed from
`EarthAtmosphere::getCorrelationCoefficients`:

    Lt = max(3 hr, 0.735 day * h_km^0.116)

About `23 hr` at `10 km`, `26.6 hr` at `35 km`, and `30.9 hr` at `125 km`, so
the `3 hr` floor never binds over this study's corridor.
"""
function gram_time_scale_s(alt_km::Float64)::Float64
    return max(10800.0, 86400.0 * 0.735 * abs(alt_km)^0.116)
end

"""
    gram_correlation_scales(alt_km) -> (horizontal_km, vertical_km)

GRAM's assumed correlation scales at an altitude, linearly interpolated and
clamped to the table's ends.
"""
function gram_correlation_scales(alt_km::Float64)::Tuple{Float64, Float64}
    n = length(GRAM_RS_ALT_KM)
    alt_km <= GRAM_RS_ALT_KM[1] && return (GRAM_RS_LH_KM[1], GRAM_RS_LZ_KM[1])
    alt_km >= GRAM_RS_ALT_KM[n] && return (GRAM_RS_LH_KM[n], GRAM_RS_LZ_KM[n])
    j = searchsortedfirst(GRAM_RS_ALT_KM, alt_km)
    t = (alt_km - GRAM_RS_ALT_KM[j - 1]) / (GRAM_RS_ALT_KM[j] - GRAM_RS_ALT_KM[j - 1])
    return (
        muladd(t, GRAM_RS_LH_KM[j] - GRAM_RS_LH_KM[j - 1], GRAM_RS_LH_KM[j - 1]),
        muladd(t, GRAM_RS_LZ_KM[j] - GRAM_RS_LZ_KM[j - 1], GRAM_RS_LZ_KM[j - 1]),
    )
end

"""
    chordal_distance_km(lat1_deg, lon1_deg, lat2_deg, lon2_deg) -> Float64

Straight-line separation between two points on a sphere of mean Earth radius.
Chordal rather than great-circle because `exp(-d / L)` on the chordal metric is
positive definite on the sphere by embedding in R^3, while the geodesic form is
not guaranteed to be. The two agree to under `1.5%` out to `3000 km`, well
inside the scales involved here.
"""
function chordal_distance_km(lat1_deg::Float64, lon1_deg::Float64, lat2_deg::Float64, lon2_deg::Float64)::Float64
    phi1 = deg2rad(lat1_deg); phi2 = deg2rad(lat2_deg)
    lam1 = deg2rad(lon1_deg); lam2 = deg2rad(lon2_deg)
    x1 = cos(phi1) * cos(lam1); y1 = cos(phi1) * sin(lam1); z1 = sin(phi1)
    x2 = cos(phi2) * cos(lam2); y2 = cos(phi2) * sin(lam2); z2 = sin(phi2)
    return EARTH_MEAN_RADIUS_KM * sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
end
