using Test
using Random
using SpaceAGORA

module RPOSTLDecoderBoundaryFixture
    function write_binary(path, triangles)
        open(path, "w") do io
            # A binary STL may start with "solid"; its exact size identifies it.
            write(io, vcat(Vector{UInt8}("solid binary fixture"), zeros(UInt8, 60)))
            write(io, UInt32(size(triangles, 2) ÷ 3))
            for i in 1:3:size(triangles, 2)
                write(io, Float32[0, 0, 1])
                write(io, Float32.(triangles[:, i:i+2]))
                write(io, UInt16(19))
            end
        end
    end
end

@testset "RPO CAD reuses binary decoding and retains sampling semantics" begin
    rpo = SpaceAGORA.RPOStationAssets
    structure = SpaceAGORA.SimulationModel.Structure
    @test rpo._stl_binary_triangles === structure._stl_binary_triangles
    @test SpaceAGORA.load_rpo_station_cad_triangles === rpo.load_rpo_station_cad_triangles
    @test SpaceAGORA.load_rpo_station_cad_pointcloud === rpo.load_rpo_station_cad_pointcloud
    mktempdir() do dir
        triangles = [-0.0 1.25 0.0 2.0 3.0 2.0; 0.0 0.0 2.5 0.0 0.0 4.0; 0.5 0.5 0.5 1.0 1.0 1.0]
        binary = joinpath(dir, "two.stl")
        RPOSTLDecoderBoundaryFixture.write_binary(binary, triangles)
        @test rpo._stl_is_binary(binary)
        @test isequal(rpo._load_stl_triangles(binary), triangles)
        @test isequal(rpo._load_stl_triangles(binary), structure._load_stl_triangles(binary))
        # Recorded from the original reader with the same two-triangle fixture.
        expected = [0.5644852510983267 -1.4331001406462835 0.643590858127711;
            -0.08862943573234006 -1.2422674342655142 -1.1650890427272573;
            0.25 -0.25 0.25]
        @test isequal(rpo._stl_pointcloud(binary; n_points=3, rng=MersenneTwister(42)), expected)
        for seed in (0, 42), center in (false, true), scale in (1.0, 0.001, 2.5)
            sampled_rng = MersenneTwister(seed)
            reference_rng = MersenneTwister(seed)
            points = rpo._stl_pointcloud(binary; n_points=7, rng=sampled_rng, center, scale)
            @test size(points) == (3, 7)
            @test isequal(points, rpo._stl_pointcloud(binary; n_points=7, rng=MersenneTwister(seed), center, scale))
            # One draw selects a triangle, and two select a point within it.
            for _ in 1:21
                rand(reference_rng)
            end
            @test rand(sampled_rng) === rand(reference_rng)
        end
        @test_throws ArgumentError rpo._stl_pointcloud(binary; n_points=0)
        ascii = joinpath(dir, "ascii.stl")
        write(ascii, "solid fixture\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nvertex 2 3 4\nendsolid\n")
        @test !rpo._stl_is_binary(ascii)
        @test size(rpo._load_stl_triangles(ascii)) == (3, 4)
        @test size(structure._load_stl_triangles(ascii)) == (3, 3)
        write(ascii, "solid fixture\nvertex 1 2\nendsolid\n")
        @test_throws BoundsError rpo._load_stl_triangles(ascii)
        @test size(structure._load_stl_triangles(ascii)) == (3, 0)
    end
end
