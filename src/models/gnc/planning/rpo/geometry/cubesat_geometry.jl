"""
    RPOCubeSatGeometry

Chaser geometry for RPO clearance queries: a box with body-frame half extents
`half_extents_body` in metres. Build one with `RPOCubeSatGeometry(; dims_m)`.
"""
struct RPOCubeSatGeometry
    half_extents_body::SVector{3, Float64}
end

"""
    RPOCubeSatGeometry(; dims_m=(0.1, 0.1, 0.3))

Build a box chaser from its full body-frame dimensions in metres, each finite and
positive. The default is a 3U CubeSat.
"""
function RPOCubeSatGeometry(; dims_m=(0.1, 0.1, 0.3))
    dims = SVector{3, Float64}(dims_m)
    any(x -> !isfinite(x) || x <= 0.0, dims) &&
        throw(ArgumentError("RPOCubeSatGeometry dimensions must be finite and positive."))
    return RPOCubeSatGeometry(0.5 .* dims)
end

