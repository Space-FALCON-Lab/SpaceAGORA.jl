module OrbitalGeometry
using StaticArrays

# Inertial axes to perifocal axes for radian orbital angles. Keep the original
# component expressions and their evaluation order; callers own vector scaling
# and any direction normalization. This is an internal, non-exported owner.
@inline function _perifocal_from_inertial(i::Float64, ω::Float64, Ω::Float64)::SMatrix{3, 3, Float64, 9}
    return @SMatrix [
        -sin(Ω)*cos(i)*sin(ω)+cos(Ω)*cos(ω)  cos(Ω)*cos(i)*sin(ω)+sin(Ω)*cos(ω)  sin(i)*sin(ω)
        -sin(Ω)*cos(i)*cos(ω)-cos(Ω)*sin(ω)  cos(Ω)*cos(i)*cos(ω)-sin(Ω)*sin(ω)  sin(i)*cos(ω)
         sin(Ω)*sin(i)                       -cos(Ω)*sin(i)                       cos(i)
    ]
end
end
