module TerrainModelsTests

using Test
using JSON
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const REPO_ROOT = dirname(dirname(realpath(pathof(SpaceAGORA))))
include(joinpath(REPO_ROOT, "docs", "public_api_symbols.jl"))
const TERRAIN_SYMBOLS = (
    :AbstractTerrainModel, :NoTerrainModel, :DEMGrid, :DEMTerrainModel,
    :terrain_height, :terrain_radius, :load_dem_grid, :load_site_terrain,
    :dem_grid_covers,
)
const RADIUS = 3_396_200.0
const HEIGHTS = Float32[100 110 120 130; 200 210 220 230; 300 310 320 330]

struct ConstantTerrain <: AbstractTerrainModel
    height::Float64
end
SpaceAGORA.terrain_height(t::ConstantTerrain, ::Real, ::Real) = t.height

_plane(lat, lon) = 245.0 - 100.0 * lat + 10.0 * lon
_clamped_plane(lat, lon) = _plane(clamp(lat, 0.5, 2.5), clamp(lon, 10.5, 13.5))
_grid(; kwargs...) = DEMGrid(HEIGHTS, 0.0, 3.0, 10.0, 14.0; kwargs...)
_model(g=_grid(); kwargs...) = DEMTerrainModel([g]; reference_radius_m=RADIUS, kwargs...)

function _header()
    return Dict{String,Any}(
        "rows" => 3, "cols" => 4, "lat_min" => 0.0, "lat_max" => 3.0,
        "lon_min" => 10.0, "lon_max" => 14.0,
        "reference_radius_m" => RADIUS, "source" => "synthetic test", "units" => "m",
        "layout" => "row-major, north to south, west to east, little-endian float32",
    )
end

# Explicit bytes make this fixture independent of the host's byte order and
# Julia's column-major array storage. No native assets or downloads are used.
function _little_endian_rows(heights)
    bytes = UInt8[]
    for row in axes(heights, 1), col in axes(heights, 2)
        bits = reinterpret(UInt32, Float32(heights[row, col]))
        for shift in (0, 8, 16, 24)
            push!(bytes, UInt8((bits >> shift) & 0xff))
        end
    end
    return bytes
end

function _write_json(path, value)
    open(path, "w") do io
        JSON.print(io, value)
    end
    return path
end

function _write_grid(dir, stem; meta=_header(), bytes=_little_endian_rows(HEIGHTS), ext=".json")
    path = joinpath(dir, stem * ext)
    _write_json(path, meta)
    write(joinpath(dir, stem * ".f32"), bytes)
    return path
end

function _site(entries=[Dict("name" => "coarse", "reference_radius_m" => RADIUS)])
    return Dict{String,Any}(
        "site" => Dict("lat_deg" => 1.5, "lon_deg" => 12.0, "name" => "synthetic site"),
        "dem" => entries,
    )
end

@testset "terrain public API and documentation" begin
    specs = PublicAPISymbols.public_api_specs(SpaceAGORA)
    for symbol in TERRAIN_SYMBOLS
        @test symbol in names(SpaceAGORA)
        @test getproperty(SpaceAGORA, symbol) === getproperty(SM, symbol)
        @test Base.Docs.doc(Base.Docs.Binding(SpaceAGORA, symbol)) !== nothing
        @test count(spec -> spec.symbol === symbol, specs) == 1
    end
end

@testset "regional cell-centred interpolation and priority" begin
    g = _grid(name="coarse")
    model = _model(g; fallback_height_m=-5.0)
    @test size(g) == (3, 4)
    @test g.heights isa Matrix{Float32}
    # These analytic values cover the original PR121 interpolation intentions,
    # with a nonsquare matrix that exposes transposed north/east axes.
    for lat in (0.5, 1.5, 2.5), lon in (10.5, 11.5, 12.5, 13.5)
        @test terrain_height(model, lat, lon) ≈ _plane(lat, lon) atol=1e-10
    end
    for (lat, lon) in ((1.5, 12.0), (2.0, 10.5), (1.25, 11.75), (2.25, 12.25))
        @test terrain_height(model, lat, lon) ≈ _plane(lat, lon) atol=1e-10
    end
    # Bounds are cell edges: clamp to the nearest cell centre inside coverage.
    for lat in (0.0, 0.25, 1.5, 2.75, 3.0), lon in (10.0, 10.25, 12.0, 13.75, 14.0)
        @test dem_grid_covers(g, lat, lon)
        @test terrain_height(model, lat, lon) ≈ _clamped_plane(lat, lon) atol=1e-10
    end
    for (lat, lon) in ((prevfloat(0.0), 12.0), (nextfloat(3.0), 12.0),
                       (1.5, prevfloat(10.0)), (1.5, nextfloat(14.0)))
        @test !dem_grid_covers(g, lat, lon)
        @test terrain_height(model, lat, lon) == -5.0
    end
    for shift in (-720.0, -360.0, 0.0, 360.0, 720.0)
        @test terrain_height(model, 2.5, 10.5 + shift) == 100.0
    end
    fine = DEMGrid(fill(7.0f0, 2, 2), 1.0, 2.0, 11.0, 12.0; name="fine")
    preferred = DEMTerrainModel([fine, g]; reference_radius_m=RADIUS)
    reversed = DEMTerrainModel([g, fine]; reference_radius_m=RADIUS)
    @test terrain_height(preferred, 1.5, 11.5) == 7.0
    @test terrain_height(preferred, 2.5, 10.5) == 100.0
    @test terrain_height(reversed, 1.5, 11.5) == 210.0
end

@testset "longitude wrapping terminates and preserves regional coverage" begin
    g = DEMGrid(HEIGHTS, 0.0, 3.0, 350.0, 370.0)
    model = _model(g; fallback_height_m=-5.0)
    for lon in (-5.0, 355.0, 715.0, -365.0)
        @test dem_grid_covers(g, 1.5, lon)
        @test terrain_height(model, 1.5, lon) == 205.0
    end
    @test terrain_height(model, 1.5, 0.0) == terrain_height(model, 1.5, 360.0) == 215.0
    @test dem_grid_covers(g, 1.5, 0.0)
    @test dem_grid_covers(g, 1.5, 350.0)
    @test terrain_height(model, 1.5, 350.0) == 200.0
    @test dem_grid_covers(g, 1.5, 10.0)
    @test terrain_height(model, 1.5, 10.0) == 230.0
    # Compare the unshifted seam boundary: adding 360 can round this just-outside
    # query back onto the regional edge and wrongly suppress the fallback.
    for lon in (nextfloat(10.0), prevfloat(350.0))
        @test !dem_grid_covers(g, 1.5, lon)
        @test terrain_height(model, 1.5, lon) == -5.0
    end
    @test !dem_grid_covers(g, 1.5, 180.0)
    @test terrain_height(model, 1.5, 180.0) == -5.0
    for west in (-10.0, 350.0, 710.0)
        shifted = DEMGrid(HEIGHTS, 0.0, 3.0, west, west + 20.0)
        @test shifted.lon_min == 350.0
        @test shifted.lon_max == 370.0
        @test terrain_height(_model(shifted), 1.5, -5.0) == 205.0
    end
    # A nonzero western origin catches loss of that offset when a huge query
    # is reduced. The oracle reduces the query independently with mod(lon,360).
    wide = DEMGrid(Float32[0 10 20 30; 0 10 20 30], -1.0, 1.0, 1.0, 359.0)
    wide_model = _model(wide; fallback_height_m=-99.0)
    for lon in (1e18, -1e18, 1e100, -1e100, 1e300, -1e300, floatmax(Float64))
        wrapped = mod(lon, 360.0)
        covered = 1.0 <= wrapped <= 359.0
        expected = covered ? 10.0 * (clamp(wrapped, 45.75, 314.25) - 45.75) / 89.5 : -99.0
        @test dem_grid_covers(wide, 0.0, lon) == covered
        @test terrain_height(wide_model, 0.0, lon) ≈ expected atol=1e-12
    end
end

@testset "constructors own data and reject invalid terrain" begin
    supplied = copy(HEIGHTS)
    g = DEMGrid(supplied, 0.0, 3.0, 10.0, 14.0)
    supplied[1, 1] = -800
    @test g.heights[1, 1] == 100
    input_grids = [g]
    a = DEMTerrainModel(input_grids; reference_radius_m=RADIUS)
    b = DEMTerrainModel(input_grids; reference_radius_m=RADIUS)
    @test a.grids[1].reference_radius_m == RADIUS
    @test_throws ArgumentError DEMTerrainModel(a.grids; reference_radius_m=RADIUS + 1)
    empty!(input_grids)
    g.heights[1, 1] = -900
    @test terrain_height(a, 2.5, 10.5) == 100
    @test terrain_height(b, 2.5, 10.5) == 100
    a.grids[1].heights[1, 1] = -700
    @test terrain_height(b, 2.5, 10.5) == 100
    @test g.reference_radius_m === nothing
    @test_throws ArgumentError DEMTerrainModel(DEMGrid[]; reference_radius_m=RADIUS)
    for dims in ((1, 3), (3, 1), (0, 2), (2, 0))
        @test_throws ArgumentError DEMGrid(zeros(dims...), 0.0, 3.0, 10.0, 14.0)
    end
    for bad in (NaN, Inf, -Inf, floatmax(Float64))
        h = Float64.(HEIGHTS); h[2, 2] = bad
        @test_throws ArgumentError DEMGrid(h, 0.0, 3.0, 10.0, 14.0)
    end
    for bounds in ((3.0, 0.0, 10.0, 14.0), (0.0, 0.0, 10.0, 14.0),
                   (-91.0, 3.0, 10.0, 14.0), (0.0, 91.0, 10.0, 14.0),
                   (0.0, 3.0, 14.0, 10.0), (0.0, 3.0, 10.0, 10.0),
                   (0.0, 3.0, 0.0, 360.0), (0.0, 3.0, 350.0, 711.0))
        @test_throws ArgumentError DEMGrid(HEIGHTS, bounds...)
    end
    for index in 1:4, bad in (NaN, Inf, -Inf)
        bounds = [0.0, 3.0, 10.0, 14.0]; bounds[index] = bad
        @test_throws ArgumentError DEMGrid(HEIGHTS, bounds...)
    end
    # One ULP per cell can still round adjacent half-cell centres onto the same
    # longitude. Such a grid cannot represent four distinct interpolation sites.
    @test_throws ArgumentError DEMGrid(fill(1f0, 2, 4), 0.0, 2.0,
                                       350.0, 350.0 + 4 * eps(350.0))
    for bad in (0.0, -1.0, Inf, -Inf, NaN)
        @test_throws ArgumentError _grid(reference_radius_m=bad)
        @test_throws ArgumentError DEMTerrainModel([_grid()]; reference_radius_m=bad)
    end
    for bad in (NaN, Inf, -Inf)
        @test_throws ArgumentError _model(; fallback_height_m=bad)
    end
    overflow = big"1e10000"
    @test_throws ArgumentError DEMGrid(HEIGHTS, big"0", overflow, 10, 14)
    @test_throws ArgumentError DEMTerrainModel([_grid()]; reference_radius_m=overflow)
    @test_throws ArgumentError _model(; fallback_height_m=overflow)
end

@testset "query validation and explicit reference sphere" begin
    g = _grid(reference_radius_m=RADIUS)
    model = _model(g)
    @test terrain_radius(model, 2.5, 10.5) == RADIUS + 100
    @test terrain_radius(model, 2.5, 10.5, RADIUS) == terrain_radius(model, 2.5, 10.5)
    @test terrain_height(NoTerrainModel(), 1.0, 2.0) == 0.0
    @test terrain_radius(NoTerrainModel(), 1.0, 2.0, RADIUS) == RADIUS
    @test terrain_radius(ConstantTerrain(12.5), 1.0, 2.0, RADIUS) == RADIUS + 12.5
    @test_throws ArgumentError DEMTerrainModel([g]; reference_radius_m=RADIUS + 1)
    @test_throws ArgumentError terrain_radius(model, 1.5, 12.0, RADIUS + 1)
    @test_throws ArgumentError DEMTerrainModel([g, _grid(reference_radius_m=RADIUS + 1)]; reference_radius_m=RADIUS)
    for (lat, lon) in ((91.0, 12.0), (-91.0, 12.0), (NaN, 12.0), (Inf, 12.0),
                       (-Inf, 12.0), (1.0, NaN), (1.0, Inf), (1.0, -Inf))
        @test_throws ArgumentError dem_grid_covers(g, lat, lon)
        @test_throws ArgumentError terrain_height(model, lat, lon)
        @test_throws ArgumentError terrain_height(NoTerrainModel(), lat, lon)
        @test_throws ArgumentError terrain_radius(ConstantTerrain(12.5), lat, lon, RADIUS)
    end
    for bad in (0.0, -1.0, NaN, Inf, -Inf)
        @test_throws ArgumentError terrain_radius(NoTerrainModel(), 1.0, 2.0, bad)
    end
    @test_throws ArgumentError terrain_radius(ConstantTerrain(NaN), 1.0, 2.0, RADIUS)
end

@testset "local DEM binary and header contract" begin
    mktempdir() do dir
        bytes = _little_endian_rows(HEIGHTS)
        @test bytes[1:4] == UInt8[0x00, 0x00, 0xc8, 0x42] # 100.0f0, little-endian
        @test length(bytes) == 48
        path = _write_grid(dir, "non.square"; ext=".JSON")
        loaded = load_dem_grid(path)
        @test loaded.heights == HEIGHTS
        @test loaded.reference_radius_m == RADIUS
        @test terrain_height(_model(loaded), 1.5, 12.0) == 215.0
        minimal = _header(); delete!(minimal, "units"); delete!(minimal, "layout")
        @test load_dem_grid(_write_grid(dir, "minimal"; meta=minimal)).heights == HEIGHTS
        producer = _header()
        producer["reference_radius_m"] = 1_737_400.0
        producer["units"] = "m above the 1737.4 km sphere"
        producer["format"] = "little-endian float32"
        @test load_dem_grid(_write_grid(dir, "producer"; meta=producer)).reference_radius_m == 1_737_400.0
        producer["reference_radius_m"] = RADIUS
        @test_throws ArgumentError load_dem_grid(_write_grid(dir, "wrong_units_radius"; meta=producer))
        for bytes_bad in (UInt8[], bytes[1:end-1], vcat(bytes, UInt8[0]), vcat(bytes, bytes[1:4]))
            @test_throws ArgumentError load_dem_grid(_write_grid(dir, "wrong_size"; bytes=bytes_bad))
        end
        missing = _write_grid(dir, "missing")
        rm(joinpath(dir, "missing.f32"))
        @test_throws ArgumentError load_dem_grid(missing)
        for bad in (NaN32, Inf32, -Inf32)
            heights = copy(HEIGHTS); heights[2, 3] = bad
            @test_throws ArgumentError load_dem_grid(_write_grid(dir, "invalid_height"; bytes=_little_endian_rows(heights)))
        end
        for key in ("rows", "cols", "lat_min", "lat_max", "lon_min", "lon_max", "reference_radius_m")
            meta = _header(); delete!(meta, key)
            @test_throws ArgumentError load_dem_grid(_write_grid(dir, "missing_key"; meta=meta))
        end
        for key in ("rows", "cols"), bad in (0, 1, -2, 3.5, 3.0, true, "3", nothing, typemax(Int))
            meta = _header(); meta[key] = bad
            @test_throws ArgumentError load_dem_grid(_write_grid(dir, "bad_dimension"; meta=meta))
        end
        for (key, bad) in (("lat_min", -91), ("lat_max", 91), ("lat_max", 0),
                           ("lon_max", 10), ("lon_max", 370), ("lon_max", 371),
                           ("reference_radius_m", 0), ("reference_radius_m", -1),
                           ("reference_radius_m", nothing), ("units", "km"),
                           ("lat_min", nothing), ("lon_max", "14"),
                           ("format", "big-endian float32"),
                           ("layout", "column-major, big-endian float64"))
            meta = _header(); meta[key] = bad
            @test_throws ArgumentError load_dem_grid(_write_grid(dir, "bad_metadata"; meta=meta))
        end
        @test_throws ArgumentError load_dem_grid(_write_grid(dir, "wrong_extension"; ext=".txt"))
        @test_throws ArgumentError load_dem_grid(_write_grid(dir, "nonobject"; meta=Any[]))
    end
end

@testset "site files preserve grid priority and enforce one known radius" begin
    mktempdir() do dir
        _write_grid(dir, "coarse")
        fine_meta = _header()
        merge!(fine_meta, Dict{String,Any}("rows" => 2, "cols" => 2, "lat_min" => 1.0,
                                          "lat_max" => 2.0, "lon_min" => 11.0, "lon_max" => 12.0))
        _write_grid(dir, "fine"; meta=fine_meta, bytes=_little_endian_rows(fill(7f0, 2, 2)))
        entries = [Dict("name" => "fine", "reference_radius_m" => RADIUS),
                   Dict("name" => "coarse", "reference_radius_m" => RADIUS)]
        meta = _site(entries)
        site_path = _write_json(joinpath(dir, "site.json"), meta)
        model, site = load_site_terrain(site_path)
        @test site.name == "synthetic site"
        @test site.height_m == 7.0
        @test site.lat_deg == 1.5 && site.lon_deg == 12.0
        @test site.reference_radius_m == model.reference_radius_m == RADIUS
        @test terrain_height(model, 2.5, 10.5) == 100.0
        @test [g.name for g in model.grids] == ["fine", "coarse"]
        for altered in (
            _site(Any[]),
            _site([Dict("name" => "coarse")]),
            _site([Dict("name" => "coarse", "reference_radius_m" => RADIUS + 1)]),
            _site([entries[1], Dict("name" => "coarse", "reference_radius_m" => RADIUS + 1)]),
        )
            @test_throws ArgumentError load_site_terrain(_write_json(site_path, altered))
        end
        # A new top-level field cannot stand in for the missing per-DEM datum.
        top_level_only = _site([Dict("name" => "coarse")])
        top_level_only["reference_radius_m"] = RADIUS
        @test_throws ArgumentError load_site_terrain(_write_json(site_path, top_level_only))
        for key in ("lat_deg", "lon_deg")
            invalid = _site(); invalid["site"][key] = nothing
            @test_throws ArgumentError load_site_terrain(_write_json(site_path, invalid))
        end
        invalid_lat = _site(); invalid_lat["site"]["lat_deg"] = 91.0
        @test_throws ArgumentError load_site_terrain(_write_json(site_path, invalid_lat))
        different_header = _header(); different_header["reference_radius_m"] = RADIUS + 1
        _write_grid(dir, "coarse"; meta=different_header)
        @test_throws ArgumentError load_site_terrain(_write_json(site_path, _site(entries)))
        mixed_matching_entries = [entries[1], Dict("name" => "coarse", "reference_radius_m" => RADIUS + 1)]
        @test_throws ArgumentError load_site_terrain(_write_json(site_path, _site(mixed_matching_entries)))
        unknown_header = _header(); delete!(unknown_header, "reference_radius_m")
        _write_grid(dir, "coarse"; meta=unknown_header)
        @test_throws ArgumentError load_site_terrain(_write_json(site_path, _site()))
    end
end

end # module TerrainModelsTests
