using Test
using SpaceAGORA
using DataFrames
using CSV
using Arrow
using StaticArrays
using Random

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

_assert_paths_agree("uniform ranges, all entries the same width", :x, "x",
    (row, n) -> [(row + 10i):(row + 10i + 2) for i in 1:n], 3, 4)
_assert_paths_agree("uniform plain vectors", :x, "x",
    (row, n) -> [[Float64(row + i), Float64(row - i)] for i in 1:n], 3, 4)

# Ragged per-satellite entries: the component count is not a property of the
# field, so the fast path must decline and let the generic path name the columns
# per entry. Both directions were live defects when the width was inferred from
# `isbitstype(S)` -- true for `UnitRange{Int}`, whose length is runtime data.
@testset "ragged entries fall back, and the fallback is correct" begin
    # First entry short: the fast path wrote four columns where the generic path
    # writes five, silently losing sc2_x_3.
    short_first = [begin s = _SaveData(); s[:x] = r; s end for r in ([1:2, 10:12], [2:3, 20:22])]
    direct = DataFrame()
    @test !IOO._direct_assembly_columns!(direct, _SaveFieldStub(:x, true, "x"), short_first, 2)
    @test isempty(names(direct))

    # The real entry point must still produce the generic path's five columns.
    full = DataFrame()
    IOO._append_save_field_columns!(full, _SaveFieldStub(:x, true, "x"), short_first, 2)
    @test names(full) == ["sc1_x_1", "sc1_x_2", "sc2_x_1", "sc2_x_2", "sc2_x_3"]
    @test full[!, "sc2_x_3"] == [12, 22]

    # First entry long: worse than a dropped column. The write loop is
    # `@inbounds`, so reading past the end of `10:11` fabricated a third
    # component (12) that was never in the input.
    long_first = [begin s = _SaveData(); s[:x] = [1:3, 10:11]; s end]
    direct_long = DataFrame()
    @test !IOO._direct_assembly_columns!(direct_long, _SaveFieldStub(:x, true, "x"), long_first, 2)
    @test isempty(names(direct_long))

    full_long = DataFrame()
    IOO._append_save_field_columns!(full_long, _SaveFieldStub(:x, true, "x"), long_first, 2)
    @test names(full_long) == ["sc1_x_1", "sc1_x_2", "sc1_x_3", "sc2_x_1", "sc2_x_2"]
end

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

# The cases above are hand-picked, which is their weakness: they encode what the
# author thought the fast path had to handle, so they cannot falsify the
# author's own assumptions. The `isbitstype` defect survived them for exactly
# that reason — no hand-written case fed the fast path a bitstype whose length
# was runtime data.
#
# This generates the shapes instead and compares the real entry point against
# the generic path for every one of them. Element kind, satellite count, row
# count and per-satellite component counts are all drawn, with ragged trials
# included deliberately, since raggedness is the axis the hand-written cases
# under-covered. Seeded, so a failure is reproducible.
#
# Validated by running it against the defective implementation: 31 of 344 cases
# diverged there, and none do here.
@testset "generated shapes agree with the generic path" begin
    function generated_entry(rng, kind, len)
        kind === :svector && return SVector{len, Float64}(ntuple(_ -> rand(rng), len))
        kind === :vector  && return [rand(rng) for _ in 1:len]
        start = rand(rng, 1:50)
        return start:(start + len - 1)
    end

    rng = MersenneTwister(20260911)
    checked = 0
    for _ in 1:400
        num_sats = rand(rng, 1:4)
        n_rows   = rand(rng, 1:4)
        kind     = rand(rng, (:svector, :vector, :range, :scalar))
        ragged   = rand(rng, Bool)
        # An SVector carries its length in its type, so a ragged column of them
        # is not a shape this code can be handed.
        kind === :svector && ragged && continue
        widths = ragged ? [rand(rng, 1:4) for _ in 1:num_sats] :
                          fill(rand(rng, 1:4), num_sats)

        saved_data = [begin
            snapshot = _SaveData()
            snapshot[:x] = kind === :scalar ?
                [rand(rng) for _ in 1:num_sats] :
                [generated_entry(rng, kind, widths[i]) for i in 1:num_sats]
            snapshot
        end for _ in 1:n_rows]

        got, want = DataFrame(), DataFrame()
        IOO._append_save_field_columns!(got, _SaveFieldStub(:x, true, "x"), saved_data, num_sats)
        _generic_per_satellite_columns!(want, :x, "x", saved_data, num_sats)

        @test names(got) == names(want)
        @test [eltype(got[!, c]) for c in names(got)] == [eltype(want[!, c]) for c in names(want)]
        @test isequal(got, want)

        got_arrow, want_arrow = tempname(), tempname()
        Arrow.write(got_arrow, got); Arrow.write(want_arrow, want)
        @test read(got_arrow) == read(want_arrow)
        rm.((got_arrow, want_arrow); force=true)
        checked += 1
    end
    @test checked > 300
end
