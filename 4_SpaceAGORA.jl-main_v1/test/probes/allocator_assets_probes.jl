using Test
using LinearAlgebra
using Random
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

if !isdefined(@__MODULE__, :RPOStationAssets)
    include(joinpath(REPO_ROOT, "src", "assets", "rpo_station_assets.jl"))
end

const SM = SimulationModel
const CH = SM.ControlHooks
const ASSETS = RPOStationAssets

"Fake integrator state carrying only the attitude fields the allocator reads."
_probe_att_state(q::SVector{4, Float64}, ω::SVector{3, Float64}) = (sc=[(q=q, ω=ω)],)

"Minimal ASCII STL with a single unit right triangle in the z=0 plane."
const _PROBE_ASCII_STL = """
solid probe
facet normal 0.0 0.0 1.0
 outer loop
  vertex 0.0 0.0 0.0
  vertex 1.0 0.0 0.0
  vertex 0.0 1.0 0.0
 endloop
endfacet
endsolid probe
"""

"Write a 1-triangle binary STL (84-byte header block + one 50-byte record)."
function _probe_write_binary_stl(path::AbstractString, v1, v2, v3)
    open(path, "w") do io
        write(io, zeros(UInt8, 80))          # header
        write(io, UInt32(1))                 # triangle count
        write(io, Float32.([0.0, 0.0, 1.0])) # normal
        write(io, Float32.(v1))
        write(io, Float32.(v2))
        write(io, Float32.(v3))
        write(io, UInt16(0))                 # attribute byte count
    end
    return path
end

@testset "coverage debt allocator and station asset probes" begin

    @testset "reaction wheel allocator" begin
        q_id = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        ω0 = SVector{3, Float64}(0.0, 0.0, 0.0)

        # Both gains zero -> early-return zero torque regardless of state.
        passive = CH.RPOMPCControlModel()
        @test passive.attitude_kp == 0.0 && passive.rate_kd == 0.0
        τ_off = CH.rpo_reaction_wheel_torque_command(
            passive,
            _probe_att_state(SVector{4, Float64}(0.5, -0.5, 0.5, -0.5), SVector{3, Float64}(9.0, -9.0, 9.0)),
            1,
        )
        @test τ_off === SVector{3, Float64}(0.0, 0.0, 0.0)

        # Full PD path with an unsaturated (infinite) torque budget.
        pd = CH.RPOMPCControlModel(attitude_kp=2.0, rate_kd=0.5)
        @test !isfinite(pd.max_rw_torque_nm)
        q = SVector{4, Float64}(0.1, -0.2, 0.05, 0.97)
        ω = SVector{3, Float64}(0.01, -0.02, 0.03)
        τ = CH.rpo_reaction_wheel_torque_command(pd, _probe_att_state(q, ω), 1)
        @test τ ≈ -2.0 .* SVector{3, Float64}(0.1, -0.2, 0.05) .- 0.5 .* ω
        # PD torque opposes the attitude error and the body rate.
        @test dot(τ, q[1:3]) < 0.0

        # Quaternion double cover: -q is the same physical attitude, so the
        # shadow-quaternion branch (q4 < 0) must command the identical torque.
        τ_shadow = CH.rpo_reaction_wheel_torque_command(pd, _probe_att_state(-q, ω), 1)
        @test τ_shadow ≈ τ
        @test (-q)[4] < 0.0  # confirms the negated-scalar branch was exercised

        # Equilibrium (zero error, zero rate) commands exactly zero torque.
        τ_eq = CH.rpo_reaction_wheel_torque_command(pd, _probe_att_state(q_id, ω0), 1)
        @test τ_eq == SVector{3, Float64}(0.0, 0.0, 0.0)

        # Rate-damping-only model must not take the zero-gain early return.
        damper = CH.RPOMPCControlModel(attitude_kp=0.0, rate_kd=2.0)
        ω_spin = SVector{3, Float64}(0.4, -0.1, 0.25)
        τ_damp = CH.rpo_reaction_wheel_torque_command(damper, _probe_att_state(q_id, ω_spin), 1)
        @test τ_damp ≈ -2.0 .* ω_spin
        @test dot(τ_damp, ω_spin) < 0.0

        # Per-axis saturation: axes over budget clamp to +-limit, the axis
        # inside the budget passes through untouched.
        saturating = CH.RPOMPCControlModel(attitude_kp=10.0, rate_kd=0.0, max_rw_torque_nm=0.5)
        q_big = SVector{4, Float64}(0.9, 0.01, -0.9, 0.1)
        τ_sat = CH.rpo_reaction_wheel_torque_command(saturating, _probe_att_state(q_big, ω0), 1)
        @test τ_sat == SVector{3, Float64}(-0.5, -0.1, 0.5)
        @test norm(τ_sat, Inf) <= saturating.max_rw_torque_nm + eps()

        # Finite budget that is not active leaves the command unchanged.
        roomy = CH.RPOMPCControlModel(attitude_kp=2.0, rate_kd=0.5, max_rw_torque_nm=100.0)
        τ_roomy = CH.rpo_reaction_wheel_torque_command(roomy, _probe_att_state(q, ω), 1)
        @test τ_roomy ≈ τ

        # sat_idx selects which spacecraft's attitude state is used.
        two_sat = (sc=[(q=q_id, ω=ω0), (q=q, ω=ω)],)
        τ_target = CH.rpo_reaction_wheel_torque_command(pd, two_sat, 2)
        @test τ_target ≈ τ
        @test CH.rpo_reaction_wheel_torque_command(pd, two_sat, 1) == SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    @testset "station asset paths" begin
        demo_path = ASSETS.station_geometry_path()
        @test demo_path == joinpath(REPO_ROOT, "data", "rpo", "station_geometry", "demo", "station_pointcloud.csv")
        @test ASSETS.station_geometry_path(:demo) == demo_path

        gw_path = ASSETS.station_cad_path()
        @test gw_path == joinpath(REPO_ROOT, "data", "rpo", "station_geometry", "gateway", "Gateway_Core.stl")
        @test ASSETS.station_cad_path(:gateway) == gw_path
        @test ASSETS.station_cad_path(:gateway_core) == gw_path
        # CAD-backed kinds route the geometry-path helper to the CAD path.
        @test ASSETS.station_geometry_path(:gateway) == gw_path
        @test ASSETS.station_geometry_path(:gateway_core) == gw_path

        @test_throws ArgumentError ASSETS.station_geometry_path(:nonexistent_kind)
        @test_throws ArgumentError ASSETS.station_cad_path(:demo)
        @test_throws ArgumentError ASSETS.load_rpo_station_pointcloud(:gateway)
        @test_throws ArgumentError ASSETS.load_rpo_station_cad_triangles(:nonexistent_kind)
        @test_throws ArgumentError ASSETS.load_rpo_station_cad_pointcloud(:nonexistent_kind)
    end

    @testset "station demo pointcloud csv" begin
        demo_path = ASSETS.station_geometry_path(:demo)
        if isfile(demo_path)
            # Repo ships the demo asset: exercise the parser on the real file.
            pc = ASSETS.load_rpo_station_pointcloud()
            @test size(pc, 1) == 3
            @test size(pc, 2) >= 1
            @test all(isfinite, pc)
        else
            # Missing-file error path first, then drive the CSV parser with a
            # temporary fixture that is removed again afterwards.
            @test_throws ArgumentError ASSETS.load_rpo_station_pointcloud()
            demo_dir = dirname(demo_path)
            made_dir = !isdir(demo_dir)
            made_dir && mkpath(demo_dir)
            try
                write(demo_path, "# station demo probe fixture\n\n1.0,2.0,3.0\n4.0,5.0,6.0\n")
                pc = ASSETS.load_rpo_station_pointcloud(:demo)
                @test pc == [1.0 4.0; 2.0 5.0; 3.0 6.0]

                # Rows without exactly x,y,z entries are rejected.
                write(demo_path, "1.0,2.0\n")
                @test_throws ArgumentError ASSETS.load_rpo_station_pointcloud()
            finally
                rm(demo_path; force=true)
                made_dir && isdir(demo_dir) && isempty(readdir(demo_dir)) && rm(demo_dir)
            end
            @test !isfile(demo_path)  # fixture cleaned up
        end
    end

    @testset "STL parsing helpers" begin
        mktempdir() do dir
            # ASCII STL: header sniffing must classify it as non-binary.
            ascii_path = joinpath(dir, "probe_ascii.stl")
            write(ascii_path, _PROBE_ASCII_STL)
            @test !ASSETS._stl_is_binary(ascii_path)

            # Files shorter than the 84-byte binary header are never binary.
            tiny_path = joinpath(dir, "probe_tiny.stl")
            write(tiny_path, "solid t\n")
            @test !ASSETS._stl_is_binary(tiny_path)

            # ASCII vertex parsing preserves the triangle exactly.
            tri_ascii = ASSETS._load_stl_triangles(ascii_path)
            @test size(tri_ascii) == (3, 3)
            @test tri_ascii[:, 1] == [0.0, 0.0, 0.0]
            @test tri_ascii[:, 2] == [1.0, 0.0, 0.0]
            @test tri_ascii[:, 3] == [0.0, 1.0, 0.0]

            # Binary STL: exact size match on the 50-byte record layout.
            bin_path = _probe_write_binary_stl(
                joinpath(dir, "probe_bin.stl"),
                [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0],
            )
            @test filesize(bin_path) == 84 + 50
            @test ASSETS._stl_is_binary(bin_path)
            tri_bin = ASSETS._load_stl_triangles(bin_path)
            @test size(tri_bin) == (3, 3)
            @test tri_bin[:, 2] ≈ [2.0, 0.0, 0.0]
            @test tri_bin[:, 3] ≈ [0.0, 2.0, 0.0]

            # Centering shifts the bounding-box midpoint onto the origin.
            centered = ASSETS._center_triangles!(copy(tri_bin))
            mids = (vec(minimum(centered; dims=2)) .+ vec(maximum(centered; dims=2))) ./ 2.0
            @test all(abs.(mids) .< 1.0e-12)

            # Point sampling stays on the (uncentered) triangle surface.
            rng = MersenneTwister(11)
            pts = ASSETS._stl_pointcloud(ascii_path; n_points=200, rng=rng, center=false)
            @test size(pts) == (3, 200)
            @test all(pts[3, :] .== 0.0)
            @test all(pts[1, :] .>= -1.0e-12)
            @test all(pts[2, :] .>= -1.0e-12)
            @test all(pts[1, :] .+ pts[2, :] .<= 1.0 + 1.0e-12)  # barycentric fold u+v<=1

            # Scale is applied before centering/sampling.
            pts_scaled = ASSETS._stl_pointcloud(ascii_path; n_points=50, rng=MersenneTwister(3), center=false, scale=2.0)
            @test all(pts_scaled[1, :] .+ pts_scaled[2, :] .<= 2.0 + 1.0e-12)
            @test maximum(pts_scaled[1, :] .+ pts_scaled[2, :]) > 1.0  # actually uses the scaled triangle

            # Error paths of the sampler.
            @test_throws ArgumentError ASSETS._stl_pointcloud(ascii_path; n_points=0)

            empty_path = joinpath(dir, "probe_empty.stl")
            write(empty_path, "solid empty\nendsolid empty\n")
            @test_throws ArgumentError ASSETS._stl_pointcloud(empty_path)  # no triangles

            degen_path = joinpath(dir, "probe_degenerate.stl")
            write(degen_path, """
            solid degen
            facet normal 0.0 0.0 1.0
             outer loop
              vertex 0.0 0.0 0.0
              vertex 1.0 0.0 0.0
              vertex 2.0 0.0 0.0
             endloop
            endfacet
            endsolid degen
            """)
            @test_throws ArgumentError ASSETS._stl_pointcloud(degen_path)  # zero total area
        end
    end

    @testset "gateway CAD loaders" begin
        gw_path = ASSETS.station_cad_path(:gateway)
        if isfile(gw_path)
            tri = ASSETS.load_rpo_station_cad_triangles()
            @test size(tri, 1) == 3
            @test size(tri, 2) % 3 == 0
            @test size(tri, 2) > 0
            # Default centering puts the bounding-box midpoint at the origin.
            mids = (vec(minimum(tri; dims=2)) .+ vec(maximum(tri; dims=2))) ./ 2.0
            @test all(abs.(mids) .< 1.0e-6)

            raw = ASSETS.load_rpo_station_cad_triangles(:gateway_core; center=false)
            scaled = ASSETS.load_rpo_station_cad_triangles(:gateway; center=false, scale=2.0)
            @test scaled ≈ 2.0 .* raw

            pc = ASSETS.load_rpo_station_cad_pointcloud(; n_points=64)
            @test size(pc) == (3, 64)
            # Default rng is a fixed seed: repeat calls are deterministic.
            @test pc == ASSETS.load_rpo_station_cad_pointcloud(; n_points=64)
            # Sampled surface points live inside the centered CAD bounding box.
            lo = vec(minimum(tri; dims=2)) .- 1.0e-9
            hi = vec(maximum(tri; dims=2)) .+ 1.0e-9
            @test all(lo .<= vec(minimum(pc; dims=2)))
            @test all(vec(maximum(pc; dims=2)) .<= hi)

            pc_alt = ASSETS.load_rpo_station_cad_pointcloud(:gateway_core; n_points=8, rng=MersenneTwister(5), center=false, scale=0.5)
            @test size(pc_alt) == (3, 8)
            lo_raw = 0.5 .* vec(minimum(raw; dims=2)) .- 1.0e-9
            hi_raw = 0.5 .* vec(maximum(raw; dims=2)) .+ 1.0e-9
            @test all(lo_raw .<= vec(minimum(pc_alt; dims=2)))
            @test all(vec(maximum(pc_alt; dims=2)) .<= hi_raw)
        else
            @test_throws ArgumentError ASSETS.load_rpo_station_cad_triangles()
            @test_throws ArgumentError ASSETS.load_rpo_station_cad_pointcloud()
        end
    end
end
