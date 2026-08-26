"""
    WAMIPETruthSource

Study-local configuration for WAM-IPE neutral-density truth. Coordinates passed
to this adapter are latitude and longitude in degrees and altitude in metres,
which matches the corridor representation used by this study.
"""
struct WAMIPETruthSource
    interpolator::WamIPE.Density.WAMInterpolator
end

function WAMIPETruthSource(; product::String="wfs", interpolation::Symbol=:sciml)
    return WAMIPETruthSource(WamIPE.Density.WAMInterpolator(
        product=product,
        interpolation=interpolation,
    ))
end

"""
    sample_wam_ipe(source, times, lats_deg, lons_deg, alts_m)

Retrieve neutral-density truth for corresponding points. This deliberately uses
WamIPE's optimized trajectory API: it groups points sharing a NetCDF file pair
and keeps each pair open while its queries are evaluated.
"""
function sample_wam_ipe(
    source::WAMIPETruthSource,
    times::AbstractVector{<:DateTime},
    lats_deg::AbstractVector{<:Real},
    lons_deg::AbstractVector{<:Real},
    alts_m::AbstractVector{<:Real},
)::Vector{Float64}
    n = length(times)
    (length(lats_deg) == n && length(lons_deg) == n && length(alts_m) == n) ||
        throw(ArgumentError("WAM-IPE query vectors must share length."))
    all(>(0), alts_m) || throw(ArgumentError("WAM-IPE requires altitudes strictly above 0 m."))
    return WamIPE.Density.get_density_trajectory_optimised(
        source.interpolator,
        times,
        lats_deg,
        lons_deg,
        alts_m;
        angles_in_deg=true,
    )
end
