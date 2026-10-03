"""Simple chaser geometry model for RPO clearance calculations."""
struct RPOCubeSatGeometry
    half_extents_body::SVector{3, Float64}
end

"""Simple chaser geometry model for RPO clearance calculations."""
function RPOCubeSatGeometry(; dims_m=(0.1, 0.1, 0.3))
    dims = SVector{3, Float64}(dims_m)
    any(x -> !isfinite(x) || x <= 0.0, dims) &&
        throw(ArgumentError("RPOCubeSatGeometry dimensions must be finite and positive."))
    return RPOCubeSatGeometry(0.5 .* dims)
end

