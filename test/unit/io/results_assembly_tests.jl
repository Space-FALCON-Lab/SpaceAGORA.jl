using Test
using SpaceAGORA
using DataFrames
using CSV
using Arrow
using StaticArrays

const IOO = SpaceAGORA.SimulationModel.IOOutputs
const _SaveData = SpaceAGORA.SimulationModel.SaveData

# `_direct_assembly_columns!` is a fast path under `_append_save_field_columns!`
# and must be indistinguishable from the generic comprehension path it bypasses
# -- not merely equal in value, but identical in column order, column element
# type, and therefore in the CSV and Arrow bytes, since the feather schema is
# part of the published output. This is the oracle: the generic path written out
# exactly as `_append_save_field_columns!` runs it.
function _generic_per_satellite_columns!(df, name::Symbol, prefix::String, saved_data, num_sats::Int)
    field_series = [snapshot[name] for snapshot in saved_data]
    for sat_idx in 1:num_sats
        sat_series = [value[sat_idx] for value in field_series]
        IOO._append_series_columns!(df, "sc$(sat_idx)_$(prefix)", sat_series)
    end
    return nothing
end

struct _SaveFieldStub
    name::Symbol
    per_satellite::Bool
    column_prefix::String
end

function _assert_paths_agree(label, name, prefix, make_row, num_sats, n_rows)
    saved_data = [begin
        snapshot = _SaveData()
        snapshot[name] = make_row(row, num_sats)
        snapshot
    end for row in 1:n_rows]
    times = collect(1.0:n_rows)
    direct = DataFrame(time=times)
    generic = DataFrame(time=times)

    @testset "$label" begin
        @test IOO._direct_assembly_columns!(direct, _SaveFieldStub(name, true, prefix), saved_data, num_sats)
        _generic_per_satellite_columns!(generic, name, prefix, saved_data, num_sats)

        @test names(direct) == names(generic)
        @test [eltype(direct[!, c]) for c in names(direct)] ==
              [eltype(generic[!, c]) for c in names(generic)]
        @test isequal(direct, generic)

        direct_csv, generic_csv = tempname(), tempname()
        CSV.write(direct_csv, direct); CSV.write(generic_csv, generic)
        @test read(direct_csv) == read(generic_csv)

        direct_arrow, generic_arrow = tempname(), tempname()
        Arrow.write(direct_arrow, direct); Arrow.write(generic_arrow, generic)
        @test read(direct_arrow) == read(generic_arrow)

        rm.((direct_csv, generic_csv, direct_arrow, generic_arrow); force=true)
    end
end

_assert_paths_agree("3-vector field, several satellites", :position, "pos",
    (row, n) -> [SVector{3, Float64}(row + i, row - i, row * i) for i in 1:n], 7, 5)
_assert_paths_agree("scalar field, several satellites", :altitude, "altitude",
    (row, n) -> [Float64(100row + i) for i in 1:n], 7, 5)
_assert_paths_agree("one satellite", :position, "pos",
    (row, n) -> [SVector{3, Float64}(row, -row, 2row) for _ in 1:n], 1, 4)
_assert_paths_agree("one row", :position, "pos",
    (row, n) -> [SVector{3, Float64}(row, -row, 2row) for _ in 1:n], 3, 1)

# Shapes the fast path must decline, leaving them to the generic path. A `false`
# return has to mean "nothing was written", or the generic path would then
# append a second copy of the same columns.
@testset "declines shapes it does not own" begin
    function _declines(name, rows, num_sats)
        df = DataFrame()
        took = IOO._direct_assembly_columns!(df, _SaveFieldStub(name, true, "x"), rows, num_sats)
        @test !took
        @test isempty(names(df))
    end

    # Ragged rows: no fixed column count, so the width is not knowable up front.
    _declines(:x, [begin s = _SaveData(); s[:x] = [[1.0, 2.0], [3.0]]; s end for _ in 1:2], 2)
    # Non-numeric leaves still belong to the recursive generic path.
    _declines(:x, [begin s = _SaveData(); s[:x] = [(a = 1.0, b = 2.0)]; s end for _ in 1:2], 1)
    # No rows, and a row shorter than the constellation.
    _declines(:x, _SaveData[], 1)
    _declines(:x, [begin s = _SaveData(); s[:x] = [1.0, 2.0]; s end,
                   begin s = _SaveData(); s[:x] = [1.0]; s end], 2)
end
