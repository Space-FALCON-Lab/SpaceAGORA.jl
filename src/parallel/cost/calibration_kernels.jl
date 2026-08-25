# Synthetic kernels for machine-constant calibration.
#
# These exist because the real harmonics kernel cannot identify its own
# constants. Measured across every configuration it can produce -- degree and
# order crossed -- `simd_terms` and `coeff_touches` stay 0.998 correlated
# (exactly 1.000 for a zonal field, where both are linear in L), because the
# arithmetic is what consumes the table and the two cannot be varied apart.
#
# That collinearity is fatal here rather than merely inconvenient, because the
# two routing candidates weight `coeff_touches` differently: the flat SIMD batch
# walks the coefficient table once per worker, `satellite_batch` walks it once
# per satellite. Two coefficient splits that predict flat identically can
# predict batch 49% apart (measured at N=1024, a=12, L=50 full field: 38.1x
# versus 19.6x for the batch/flat ratio). A fit that cannot separate the two
# terms therefore matches whatever it was shown and mis-ranks the candidate it
# was not.
#
# The synthetic kernel breaks the tie by taking the touch count and the
# arithmetic-per-touch count as independent arguments, so the calibration design
# matrix can be made orthogonal by construction rather than merely
# better-conditioned. Its access pattern deliberately mirrors the real one --
# column-major table, fixed row, varying column, so each read strides by
# `rows * 8` bytes and lands on its own cache line -- because a constant
# measured against a different access pattern will not transfer, and the
# cold-prediction check against the real kernel is what catches it if it does not.

"""
    calib_batch_kernel!(out, table, row, n_touch, simd_per_touch) -> Float64

Synthetic stand-in for the harmonics phase-4 accumulation loop.

Performs `n_touch` strided reads from `table` (fixed `row`, varying column) and,
after each, `simd_per_touch` vectorized passes over the whole of `out`. So

    coeff_touches = n_touch
    simd_terms    = n_touch * simd_per_touch      (per lane of `out`)

and the two are independently controllable, which is the entire point.
"""
function calib_batch_kernel!(
    out::Vector{Float64},
    table::Matrix{Float64},
    row::Int,
    n_touch::Int,
    simd_per_touch::Int,
)::Float64
    B = length(out)
    @inbounds for t in 1:n_touch
        c = table[row, t]
        for _ in 1:simd_per_touch
            @simd for b in 1:B
                out[b] = muladd(c, 1.0000001, out[b])
            end
        end
    end
    return @inbounds out[1]
end

"""
    calib_touch_kernel!(table, row, n_touch) -> Float64

Pure strided table walk with no arithmetic beyond the accumulation needed to
keep it alive. Isolates the cost of a coefficient touch at whatever cache level
`table` occupies.
"""
function calib_touch_kernel!(table::Matrix{Float64}, row::Int, n_touch::Int)::Float64
    acc = 0.0
    @inbounds for t in 1:n_touch
        acc += table[row, t]
    end
    return acc
end

"""
    calib_table(rows, cols) -> Matrix{Float64}

Coefficient-table stand-in. `rows` sets the stride between consecutive touches
(`rows * 8` bytes), matching a harmonics model of degree `rows - 2`; `rows * cols
* 8` is the footprint that decides which cache level the walk hits.
"""
function calib_table(rows::Int, cols::Int)::Matrix{Float64}
    tab = Matrix{Float64}(undef, rows, cols)
    @inbounds for j in 1:cols, i in 1:rows
        tab[i, j] = 1.0 + (i + j) * 1e-9
    end
    return tab
end

"""
    calib_scalar_kernel!(state, n_items) -> Float64

Per-satellite scalar setup/teardown stand-in: non-vectorized work proportional
to satellite count, mirroring phases 1 and 5 of the harmonics kernel.
"""
function calib_scalar_kernel!(state::Vector{Float64}, n_items::Int)::Float64
    acc = 0.0
    @inbounds for i in 1:n_items
        idx = ((i - 1) % length(state)) + 1
        v = state[idx]
        acc += sqrt(v) * 1.0000001 + v * v
    end
    return acc
end
