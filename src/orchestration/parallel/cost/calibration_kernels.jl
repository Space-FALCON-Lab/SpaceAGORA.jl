# Synthetic kernels for the machine calibration's speedup ladders.
#
# Each kernel is dominated by one thing, so a ladder run over it attributes what
# it measures to that thing rather than to an unlabelled mixture:
#
#   calib_batch_kernel!   arithmetic on cache-resident data, shaped like the
#                         harmonics batch loop. Its ladder gives the dispatch
#                         primitive's own serial fraction and contention, the
#                         baseline every workload pays.
#   calib_alloc_kernel!   the allocator and the garbage collector, both
#                         process-wide. This is the mechanism behind the
#                         measured contrast between a constellation that
#                         allocates nothing per RHS pass and scales, and one
#                         that allocates tens of GiB over a solve and inverts.
#
# Neither body may be elidable. Both return a value derived from the elements
# they touch: an empty loop is optimised away, and the resulting measurement is
# wrong without looking wrong.

"""
    calib_batch_kernel!(out, table, row, n_touch, simd_per_touch) -> Float64

Synthetic stand-in for the harmonics phase-4 accumulation loop.

Performs `n_touch` strided reads from `table` (fixed `row`, varying column,
starting again at the first column after the last) and, after each,
`simd_per_touch` vectorized passes over a slice of `out`. The access pattern
mirrors the real one -- column-major table, fixed row, varying column, so each
read strides by `rows * 8` bytes and lands on its own cache line.
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
    B, depth = size(ws)
    checkbounds(table, row, 1)
    checkbounds(ws, 1, 1)
    # The ladders ask for more touches than the table has columns (up to 4096
    # of 64), so the walk wraps round the row, as the harmonics kernel walks its
    # table again for every satellite. Widening the table to the touch count
    # instead would move it out of the cache level it was sized for. A counter,
    # not `mod1`, so the wrap adds no integer division.
    cols = size(table, 2)
    col = 0
    @inbounds for t in 1:n_touch
        col = col == cols ? 1 : col + 1
        c = table[row, col]
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
    calib_alloc_kernel!(n_arrays, n_elems) -> Float64

Allocate `n_arrays` fresh vectors of `n_elems` elements, touch each, and return
a checksum. The allocation is the point, so the arrays are not reused and not
preallocated.
"""
function calib_alloc_kernel!(n_arrays::Int, n_elems::Int)::Float64
    acc = 0.0
    for _ in 1:n_arrays
        v = Vector{Float64}(undef, n_elems)
        @inbounds for i in 1:n_elems
            v[i] = i * 1.0000001
        end
        @inbounds acc += v[1] + v[n_elems]
    end
    return acc
end
