using Test
using SpaceAGORA, GRAMSuite

# SpaceAGORA's grids moved to data/planets/gram_surrogate; the pinned wrapper
# still probes the former adjacent data/GRAM_surrogate location. The extension
# supplies the relocated directory as a search root only when the wrapper's own
# package-local grid is absent, so the selection order is unchanged:
# package-local, SpaceAGORA, bundled GRAM Suite grids.
const EXT = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt)
const PLANET_KEY = GRAMSuite._gram_normalize_planet_name("mars")
const GRID_NAME = "$(PLANET_KEY)_surrogate.jls"

function _with_roots(f; package_grid::Bool, spaceagora_grid::Bool)
    mktempdir() do tmp
        package_root = joinpath(tmp, "wrapper")
        spaceagora_dir = joinpath(tmp, "relocated")
        mkpath(joinpath(package_root, "data", "GRAM_surrogate"))
        mkpath(spaceagora_dir)
        package_grid && write(joinpath(package_root, "data", "GRAM_surrogate", GRID_NAME), "grid")
        spaceagora_grid && write(joinpath(spaceagora_dir, GRID_NAME), "grid")
        f(package_root, spaceagora_dir)
    end
end

@testset "Relocated SpaceAGORA grid directory keeps the wrapper's selection order" begin
    @test EXT._grid_planet_filename("Mars") == GRID_NAME
    @test endswith(EXT._SPACEAGORA_GRID_DIR, joinpath("data", "planets", "gram_surrogate"))

    _with_roots(; package_grid=false, spaceagora_grid=true) do package_root, spaceagora_dir
        roots = EXT._grid_search_roots("mars", String[]; package_root=package_root, spaceagora_dir=spaceagora_dir)
        @test roots == [spaceagora_dir]
        # End to end through the wrapper: the relocated grid is the one selected.
        resolved = GRAMSuite.resolve_gram_grid_file("mars"; search_roots=roots)
        @test resolved == normpath(joinpath(spaceagora_dir, GRID_NAME))
        # Control for the original defect: without the extension's root the
        # wrapper no longer knows the relocated directory at all.
        bundled = try
            GRAMSuite.resolve_gram_grid_file("mars")
        catch err
            err isa ArgumentError ? nothing : rethrow()
        end
        @test bundled != normpath(joinpath(spaceagora_dir, GRID_NAME))
    end

    _with_roots(; package_grid=true, spaceagora_grid=true) do package_root, spaceagora_dir
        # The wrapper's package-local grid keeps priority: no root is added.
        @test EXT._grid_search_roots("mars", String[]; package_root=package_root, spaceagora_dir=spaceagora_dir) == String[]
    end

    _with_roots(; package_grid=false, spaceagora_grid=false) do package_root, spaceagora_dir
        # Nothing to add when the relocated grid is not present either.
        @test EXT._grid_search_roots("mars", String[]; package_root=package_root, spaceagora_dir=spaceagora_dir) == String[]
    end

    _with_roots(; package_grid=false, spaceagora_grid=true) do package_root, spaceagora_dir
        # Caller-supplied roots and non-vector values pass through untouched.
        callers = [joinpath(package_root, "elsewhere")]
        @test EXT._grid_search_roots("mars", callers; package_root=package_root, spaceagora_dir=spaceagora_dir) === callers
        @test EXT._grid_search_roots("mars", ("tuple",); package_root=package_root, spaceagora_dir=spaceagora_dir) === ("tuple",)
    end

    # Explicit files bypass discovery entirely: the wrapper's own validation
    # still rejects a missing explicit path.
    @test_throws ArgumentError SpaceAGORA.GRAMGridAtmosphereModel(planet="mars", surrogate_file=joinpath(tempdir(), "absent_surrogate.jls"))
end
