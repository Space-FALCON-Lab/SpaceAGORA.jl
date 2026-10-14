# Absolute RTN error components for logarithmic mission plots.
# Inputs are equal-length Cartesian error and reference-state sample vectors.
# Outputs retain the error input units and sample order, with an eps(Float64)
# floor. This owner needs no mission setup, plotting backend or SPICE bindings.
using LinearAlgebra: norm, cross, dot
using StaticArrays: SVector

function _rtn_error_components(
    err_x::AbstractVector{<:Real},
    err_y::AbstractVector{<:Real},
    err_z::AbstractVector{<:Real},
    ref_x::AbstractVector{<:Real},
    ref_y::AbstractVector{<:Real},
    ref_z::AbstractVector{<:Real},
    ref_vx::AbstractVector{<:Real},
    ref_vy::AbstractVector{<:Real},
    ref_vz::AbstractVector{<:Real}
)
    n = length(err_x)
    err_r = Vector{Float64}(undef, n)
    err_t = Vector{Float64}(undef, n)
    err_n = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        r_ref = SVector{3, Float64}(ref_x[i], ref_y[i], ref_z[i])
        v_ref = SVector{3, Float64}(ref_vx[i], ref_vy[i], ref_vz[i])
        err = SVector{3, Float64}(err_x[i], err_y[i], err_z[i])
        r_mag = norm(r_ref)
        h_vec = cross(r_ref, v_ref)
        h_mag = norm(h_vec)
        if r_mag <= eps(Float64) || h_mag <= eps(Float64)
            throw(ArgumentError("Cannot construct RTN frame at sample $i: degenerate reference position/angular momentum."))
        end
        r_hat = r_ref / r_mag
        n_hat = h_vec / h_mag
        t_hat = cross(n_hat, r_hat)
        err_r[i] = max(abs(dot(err, r_hat)), eps(Float64))
        err_t[i] = max(abs(dot(err, t_hat)), eps(Float64))
        err_n[i] = max(abs(dot(err, n_hat)), eps(Float64))
    end
    return err_r, err_t, err_n
end
