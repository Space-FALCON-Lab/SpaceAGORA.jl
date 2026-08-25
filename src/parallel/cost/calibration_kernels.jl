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
    ws::Matrix{Float64},
    table::Matrix{Float64},
    row::Int,
    n_touch::Int,
    simd_per_touch::Int,
)::Float64
    # `ws` is B x depth, mirroring the harmonics kernel's A workspace
    # (batch x (L+3) x (M+2)), and successive SIMD passes walk successive depth
    # slices rather than hammering one vector.
    #
    # The depth dimension is not decoration. At batch width 1024 with L=20 the
    # real A workspace is 1024*23*22 doubles = 4.1 MB, far outside any cache,
    # while a B-element vector stand-in is 8 KB and L1-resident. Calibrating
    # against the 8 KB version measured pure arithmetic throughput and missed
    # the memory traffic that dominates the real kernel at large batch -- which
    # is why the model under-predicted the serial case by ~2x and kept choosing
    # it.
    B, depth = size(ws)
    @inbounds for t in 1:n_touch
        c = table[row, t]
        for s in 1:simd_per_touch
            d = ((t * simd_per_touch + s) % depth) + 1
            @simd for b in 1:B
                ws[b, d] = muladd(c, 1.0000001, ws[b, d])
            end
        end
    end
    return @inbounds ws[1, 1]
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
const CALIB_SCALAR_FLOPS = 6

function calib_scalar_kernel!(out::Vector{Float64}, state::Vector{Float64}, n_items::Int)::Float64
    # `state` must be at least `n_items` long. An earlier version wrapped the
    # index with `((i-1) % length(state)) + 1` to be safe, and that integer
    # modulo costs 20-40 cycles against the six flops it was supposed to be
    # measuring -- so ns_per_scalar_item was reporting the cost of a division
    # instruction, roughly an order of magnitude high, and the model duly
    # concluded that the unvectorized phases dominated everything.
    # Results go to independent slots, NOT into a running accumulator.
    #
    # An accumulator creates a loop-carried dependency and makes the kernel
    # latency-bound, which measured 0.274 ns/flop -- 24x the SIMD rate. The
    # phases this is meant to represent (harmonics gather and back-transform)
    # are per-satellite independent and therefore throughput-bound, so a
    # latency-bound stand-in overstates them by roughly that factor and led the
    # model to believe the unvectorized phases dominated a zonal solve.
    @inbounds for i in 1:n_items
        v = state[i]
        # Exactly CALIB_SCALAR_FLOPS flops, and no sqrt: a transcendental would
        # make the per-flop rate depend on how many of them the counted kernel
        # happens to contain, which is not something the count can express.
        out[i] = muladd(v, 1.0000001, v) + v * v - v * 0.5
    end
    return @inbounds out[1]
end
