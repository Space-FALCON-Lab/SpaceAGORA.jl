include("plot_firing_latitude.jl")
using .FiringLatitude
using CSV, DataFrames, LinearAlgebra, StaticArrays, TOML, Test

const FL = FiringLatitude

function synthetic_state(u_deg, relative_deg; helper_normal_m=0.0)
    u, relative, inclination = deg2rad(u_deg), deg2rad(relative_deg), deg2rad(1.0)
    radial = SVector(cos(u), sin(u) * cos(inclination), sin(u) * sin(inclination))
    transverse = SVector(-sin(u), cos(u) * cos(inclination), cos(u) * sin(inclination))
    normal = cross(radial, transverse)
    target = 7_478_136.6radial
    helper = 7_378_136.6 * (cos(relative) * radial + sin(relative) * transverse) +
        helper_normal_m * normal
    vcat(helper, zeros(3), target, 7300transverse, zeros(11))
end

struct AnalyticTrajectory
    t::Vector{Float64}
end

function (::AnalyticTrajectory)(time)
    state = synthetic_state(359 + 2time, -0.15 + 0.2time)
    state[13:23] .= collect(1:11) .* time
    state
end

@testset "Target argument of latitude and projected helper phase" begin
    for angle in (0, 45, 90, 180, 270, 359)
        g = FL.geometry(synthetic_state(angle, 0.4; helper_normal_m=100))
        @test abs(FL.CR.OE.signed_angle(g.u - deg2rad(angle))) < 1e-12
        @test rad2deg(g.relative_phase) ≈ 0.4 atol=1e-12
        @test g.helper_normal_m ≈ 100 atol=1e-8
    end
    g = FL.geometry(synthetic_state(90, 0.4))
    @test rad2deg(g.u) ≈ 90
    @test rad2deg(asin(g.r[3] / norm(g.r))) ≈ 1
    @test rad2deg(FL.geometry(synthetic_state(30, -0.4)).relative_phase) ≈ -0.4
    state = synthetic_state(120, -0.4; helper_normal_m=100)
    node_angle = deg2rad(70)
    rotation = [cos(node_angle) -sin(node_angle) 0; sin(node_angle) cos(node_angle) 0; 0 0 1]
    for indices in (1:3, 4:6, 7:9, 10:12)
        state[indices] = rotation * state[indices]
    end
    rotated = FL.geometry(state)
    @test rad2deg(rotated.u) ≈ 120
    @test rad2deg(rotated.relative_phase) ≈ -0.4
    @test rotated.helper_normal_m ≈ 100 atol=1e-8
end

@testset "Mutual helper-target orbital nodes" begin
    for angle in (-179, -30, 0, 30, 90, 179)
        state = synthetic_state(45, 0.4)
        g = FL.geometry(state)
        expected_node = cosd(angle) * g.node + sind(angle) * cross(g.normal, g.node)
        helper_normal = cosd(1) * g.normal + sind(1) * cross(g.normal, expected_node)
        state[1:3] = 7_378_136.6expected_node
        state[4:6] = 7350cross(helper_normal, expected_node)
        actual_u = FL.mutual_node_u(state)
        @test rad2deg(actual_u) ≈ angle atol=1e-10
        actual_node = cos(actual_u) * g.node + sin(actual_u) * cross(g.normal, g.node)
        @test dot(actual_node, helper_normal) ≈ 0 atol=1e-12
        @test dot(actual_node, g.normal) ≈ 0 atol=1e-12
        @test dot(cross(g.normal, actual_node), helper_normal) > 0
        @test dot(cross(g.normal, -actual_node), helper_normal) < 0
        state[4:6] .= 0
        @test_throws ErrorException FL.mutual_node_u(state)
        state[4:6] = 7350cross(g.normal, expected_node)
        @test_throws ErrorException FL.mutual_node_u(state)
    end
end

if "--saved" in ARGS
    @testset "Saved firing distributions and provenance" begin
        directory = FL.OUTPUT
        manifest = TOML.parsefile(joinpath(directory, "validation.toml"))
        @test manifest["passed"]
        @test manifest["helper_counts"] == [300, 800]
        @test manifest["orbit_node_markers_deg"] == [0, 180, 360]
        @test manifest["orbit_node_marker_tolerance_deg"] == FL.NODE_LABEL_TOLERANCE_DEG
        for extension in ("png", "pdf"), suffix in ("", "_clean")
            @test haskey(manifest["figure_checksums"],
                "laser_on_time_and_inclination_vs_argument_of_latitude$suffix.$extension")
        end
        for group in ("csv_checksums", "figure_checksums", "source_code_checksums")
            for (name, checksum) in manifest[group]
                path = isabspath(name) ? name : joinpath(directory, name)
                @test FL.SC.filehash(path) == checksum
            end
        end
        latitude = CSV.read(joinpath(directory, "latitude_bins.csv"), DataFrame)
        joint = CSV.read(joinpath(directory, "joint_latitude_relative_phase.csv"), DataFrame)
        summary = CSV.read(joinpath(directory, "case_summary.csv"), DataFrame)
        @test nrow(latitude) == 72
        @test nrow(joint) == 2880
        @test summary.helpers == [300, 800]
        for result in eachrow(summary)
            @test -FL.NODE_LABEL_TOLERANCE_DEG < result.mutual_ascending_node_min_deg <=
                result.mutual_ascending_node_max_deg < FL.NODE_LABEL_TOLERANCE_DEG
            @test manifest["convergence"][string(result.helpers)]["max_node_extrema_change_deg"] < 1e-4
            bins = latitude[latitude.helpers .== result.helpers, :]
            @test nrow(bins) == 36
            @test bins.latitude_start_deg ≈ collect(0:10:350)
            @test bins.latitude_end_deg ≈ collect(10:10:360)
            @test all(>(0), bins.laser_on_s)
            @test bins.laser_on_hours ≈ bins.laser_on_s ./ 3600
            @test sum(bins.laser_on_s) ≈ result.firing_s atol=1e-6 rtol=0
            @test sum(bins.fraction_of_all_firing) ≈ 1
            @test sum(bins.laser_di_deg) ≈ result.laser_di_deg atol=1e-12 rtol=0
            @test sum(bins.normal_dv_mps) ≈ result.normal_dv_mps atol=1e-10 rtol=0
            @test bins.positive_laser_di_deg + bins.negative_laser_di_deg ≈
                bins.laser_di_deg atol=1e-12
            @test bins.mean_inclination_rate_deg_per_s .* bins.laser_on_s ≈ bins.laser_di_deg
            @test result.laser_di_deg + result.J2_di_deg ≈ result.total_di_deg atol=1e-10
            @test result.positive_laser_di_deg + result.negative_laser_di_deg ≈ result.laser_di_deg
            @test FL.SC.filehash(result.source_trajectory) == result.source_trajectory_sha256
            @test FL.SC.filehash(result.source_firing) == result.source_firing_sha256
            for row in eachrow(bins)
                cells = joint[(joint.helpers .== result.helpers) .&
                    (joint.latitude_start_deg .== row.latitude_start_deg), :]
                @test nrow(cells) == 40
                @test cells.relative_phase_start_deg ≈ collect(-2:0.1:1.9)
                @test cells.relative_phase_end_deg ≈ collect(-1.9:0.1:2)
                @test all(>=(0), cells.laser_on_s)
                @test sum(cells.laser_on_s) ≈ row.laser_on_s atol=1e-6 rtol=0
                @test sum(cells.laser_di_deg) ≈ row.laser_di_deg atol=1e-12 rtol=0
            end
        end
    end
end

@testset "Dense bin crossings, wraparound and time conservation" begin
    solution = AnalyticTrajectory([0.0, 1.0])
    times = FL.breakpoints(t -> FL.geometry(solution(t)), 0.0, 1.0)
    @test times ≈ [0.0, 0.25, 0.5, 0.75, 1.0] atol=1e-8
    bins = zeros(36, 12)
    joint, joint_di = zeros(36, 40), zeros(36, 40)
    FL.accumulate!(bins, joint, joint_di, solution)
    @test sum(bins[:, 1]) ≈ 1.0
    @test bins[36, 1] ≈ 0.5 atol=1e-8
    @test bins[1, 1] ≈ 0.5 atol=1e-8
    @test all(iszero, bins[2:35, :])
    @test vec(sum(bins[:, 2:end]; dims=1)) ≈ collect(1:11)
    @test sum(joint) ≈ 1.0
    @test sum(joint_di) ≈ 1.0
    @test count(>(1e-8), joint) == 4
    @test FL.bin_index(FL.U_EDGES, 0.0) == 1
    @test FL.bin_index(FL.U_EDGES, 2pi) == 36
    @test_throws ErrorException FL.bin_index(FL.U_EDGES, -0.01)
    @test_throws ErrorException FL.bin_index(FL.RELATIVE_EDGES, deg2rad(3.0))
    @test_throws ErrorException FL.breakpoints(
        t -> (u=deg2rad(10 - t), relative_phase=0.0), 0.0, 1.0)
end

@testset "Normal-force geometry and Gauss inclination sign" begin
    p = (; planet=(μ=FL.SC.MU, J2=FL.SC.J2, Rp_e=FL.SC.EARTH_RADIUS_M),
        mass=FL.SC.MASS_KG, force_N=100 * 10_000 / FL.SC.LIGHT_SPEED)
    for angle in (0, 90, 180, 270)
        state = synthetic_state(angle, 0.4; helper_normal_m=100)
        g = FL.geometry(state)
        derivative = zeros(length(state))
        FL.held_link_dynamics!(derivative, state, p, 0.0)
        expected_normal = -p.force_N / p.mass * g.helper_normal_m / g.distance
        @test derivative[17] ≈ expected_normal atol=1e-16
        @test derivative[13] ≈ norm(g.r) / norm(g.h) * cos(g.u) * expected_normal atol=1e-18
        @test derivative[15] + derivative[16] ≈ derivative[13]
        @test derivative[18] + derivative[19] ≈ derivative[17]
        angle == 0 && @test derivative[13] < 0
        angle == 180 && @test derivative[13] > 0
        angle in (90, 270) && @test abs(derivative[13]) < 1e-18
    end
end
