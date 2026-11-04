using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using DataFrames
using Arrow
using TOML

# Load the actual study and backfill entry points without running a campaign.
const STUDY_DIR = normpath(joinpath(@__DIR__, "..", "..", "benchmarks", "studies", "aerobraking_perturbation_mc"))
include(joinpath(STUDY_DIR, "backfill_density_history.jl"))
include(joinpath(STUDY_DIR, "analytical_perturbation_models.jl"))
const SM = SpaceAGORA.SimulationModel

struct StudyProbeFrame <: SM.AbstractTypes.AbstractEphemeridesModel
    rotation::SMatrix{3, 3, Float64, 9}
end
SM.ephemerides_time_seconds(initial_time, ::StudyProbeFrame) = 0.0
SM.planet_frame_lpi(planet, et::Float64, model::StudyProbeFrame) = model.rotation

struct StudyAnalyticalPlanet
    Rp_e::Float64
    ω::SVector{3, Float64}
    rotation::SMatrix{3, 3, Float64, 9}
end
SM.planet_frame_lpi(planet::StudyAnalyticalPlanet, et::Float64, ::SM.SpiceEphemeridesModel) = planet.rotation

struct StudyProbeAtmosphere <: SM.AbstractTypes.AbstractDensityModel
    wind::SVector{3, Float64}
end
SpaceAGORA.getDensity(model::StudyProbeAtmosphere, h::Float64, lat::Float64,
    lon::Float64, t::Float64, wind::Bool, p) = begin
    @assert wind "Study diagnostics must preserve the runtime's existing true wind request"
    (2.0e-7, 250.0, model.wind)
end

function study_probe_args(rotation, wind)
    return (
        environment_model=(planet=SM.make_no_gram_planet(:earth),
            ephemerides_model=StudyProbeFrame(rotation),
            density_model=StudyProbeAtmosphere(wind), EI=500.0, wind=false),
        mission_configuration=(keplerian=false,),
        initial_time=SM.InitialTime(year=2020),
    )
end

function study_probe_table(planet)
    r = planet.Rp_e + 150e3
    # Off-axis samples expose both spin transport and local-wind conventions.
    positions = [r * SVector(0.8, 0.6, 0.0), r * normalize(SVector(2.0, -1.0, 3.0)),
        r * normalize(SVector(-1.0, 2.0, 1.0))]
    velocities = [SVector(-4200.0, 5600.0, 1500.0), SVector(6200.0, 3300.0, -1300.0),
        SVector(3100.0, -2100.0, 6200.0)]
    return DataFrame(time=[0.0, 10.0, 20.0], sc1_pos_1=getindex.(positions, 1),
        sc1_pos_2=getindex.(positions, 2), sc1_pos_3=getindex.(positions, 3),
        sc1_vel_1=getindex.(velocities, 1), sc1_vel_2=getindex.(velocities, 2),
        sc1_vel_3=getindex.(velocities, 3), sc1_mass=fill(450.0, 3))
end

function independent_pressure(args, row)
    pos = SVector(row.sc1_pos_1, row.sc1_pos_2, row.sc1_pos_3)
    vel = SVector(row.sc1_vel_1, row.sc1_vel_2, row.sc1_vel_3)
    planet = args.environment_model.planet
    rotation = args.environment_model.ephemerides_model.rotation
    rp = rotation * pos
    # Write out the cross product to avoid sharing the conversion under test.
    wx, wy, wz = planet.ω
    transport = SVector(wy * rp[3] - wz * rp[2], wz * rp[1] - wx * rp[3], wx * rp[2] - wy * rp[1])
    vp = rotation * vel - transport
    alt, lat, lon = SM.SimulationCallbacks.rtolatlong(rp, planet)
    east = SVector(-sin(lon), cos(lon), 0.0)
    north = SVector(-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat))
    up = SVector(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
    wE, wN, wU = args.environment_model.density_model.wind
    airflow = vp - (wE * east + wN * north + wU * up)
    return (; position=rp, velocity=vp, altitude=alt, q=0.5 * 2.0e-7 * dot(airflow, airflow))
end

@testset "Aerobraking study frame and density backfill" begin
    @testset "analytical helper keeps its position and spherical-altitude convention" begin
        AP = AnalyticalPerturbationModels
        rotation = @SMatrix [0.8 0.0 0.6; 0.0 1.0 0.0; -0.6 0.0 0.8]
        planet = StudyAnalyticalPlanet(6378e3, SVector(0.0, 0.0, 7.292115e-5), rotation)
        ctx = AP.CaseContext("earth", planet, 1.0, 6500e3, 7500e3, SVector(1.0, 0.0, 0.0),
            7000e3, 10.0, 450.0, 50.0, nothing, nothing, nothing, nothing, nothing, 0.0)
        pos, vel = SVector(5000e3, 4000e3, 1000e3), SVector(-4000.0, 5000.0, 2000.0)
        frame = AP._planet_frame(ctx, pos, vel, 10.0)
        rp = rotation * pos
        expected_velocity = rotation * vel - SVector(-planet.ω[3] * rp[2], planet.ω[3] * rp[1], 0.0)
        @test frame.pos_pp == rp
        @test frame.vel_pp ≈ expected_velocity rtol=1e-14
        @test frame.alt_m == norm(rp) - planet.Rp_e
        @test frame.lat_rad == asin(clamp(rp[3] / norm(rp), -1.0, 1.0))
        @test frame.lon_rad == atan(rp[2], rp[1])
    end

    @testset "streamed and backfilled pressure use current frame and local wind" begin
        identity_rotation = SMatrix{3, 3, Float64}(I)
        theta = deg2rad(37.0)
        tilted_rotation = @SMatrix [cos(theta) 0.0 sin(theta); 0.0 1.0 0.0; -sin(theta) 0.0 cos(theta)]
        for rotation in (identity_rotation, tilted_rotation), wind in (SVector(0.0, 0.0, 0.0), SVector(110.0, -70.0, 35.0))
            args = study_probe_args(rotation, wind)
            table = study_probe_table(args.environment_model.planet)
            history = _sample_density_history(args, table)
            writer = MC.OrbitChunkWriter(args, (norbits=1,), "unused", 6000.0)
            params = (args=args, shared_buffers=(density_models=SM.AbstractTypes.AbstractDensityModel[],))
            expectations = [independent_pressure(args, row) for row in eachrow(table)]
            for (i, row) in enumerate(eachrow(table))
                state, _, _ = _state_from_row(row)
                frame = MC._planet_frame_sample(args, state, Float64(row.time))
                reference = expectations[i]
                @test frame.pos_pp ≈ reference.position rtol=1e-14
                @test frame.vel_pp ≈ reference.velocity rtol=1e-14
                @test history.qdyn[i] ≈ reference.q rtol=1e-13
                sampled = MC._stream_density_update!(writer, state, params, Float64(row.time))
                @test sampled.q ≈ reference.q rtol=1e-13
                @test sampled.q == history.qdyn[i]
                @test sampled.alt == history.altitude[i]
            end
            @test writer.max_dynamic_pressure ≈ maximum(x.q for x in expectations) rtol=1e-13
            expected_integral = sum(5.0 * (expectations[i-1].q + expectations[i].q) for i in 2:3)
            @test writer.integrated_dynamic_pressure ≈ expected_integral rtol=1e-13
        end
    end

    @testset "real Arrow repair preserves states and dry-run preserves bytes" begin
        rotation = @SMatrix [0.8 0.0 0.6; 0.0 1.0 0.0; -0.6 0.0 0.8]
        args = study_probe_args(rotation, SVector(110.0, -70.0, 35.0))
        table = study_probe_table(args.environment_model.planet)
        expected = [independent_pressure(args, row).q for row in eachrow(table)]
        mktempdir() do dir
            path = joinpath(dir, "trajectory_with_active_force.feather")
            Arrow.write(path, table)
            before = read(path)
            preview = _write_density_columns!(path, args; dry_run=true, force=true)
            @test read(path) == before
            @test preview.dataframe.dynamic_pressure_pa ≈ expected rtol=1e-13
            repaired = _write_density_columns!(path, args; force=true)
            after = DataFrame(Arrow.Table(path))
            @test after.dynamic_pressure_pa ≈ expected rtol=1e-13
            @test select(after, names(table)) == table
            @test repaired.max_dynamic_pressure ≈ maximum(expected) rtol=1e-13
        end
    end

    @testset "backfill follows producer frame and recorded epoch" begin
        info = _parse_case_dir("case_001_earth_nominal_apo_05000km_two_body_density_none_ms1p00_inc093_aop080")
        planet = SM.make_no_gram_planet(:earth)
        planet_cache = Dict{Symbol, Any}(:earth => planet)
        density_cache = Dict((info.planet, info.dynamics_case, info.density_scale) => SM.NoAtmosphereModel())
        historical = _backfill_args(info, planet_cache, density_cache)
        @test historical.environment_model.ephemerides_model isa SM.SpiceEphemeridesModel
        @test historical.initial_time.year == 2020
        @test historical.initial_time.month == 1
        @test historical.initial_time.day == 1
        mktempdir() do dir
            legacy = @test_logs (:warn, r"Legacy study output") _backfill_frame_configuration(dir)
            @test legacy.ephemerides_model isa SM.SpiceEphemeridesModel
            metadata = MC._study_frame_metadata()
            @test metadata["ephemerides_model"] == "SpiceEphemeridesModel"
            metadata["initial_time"]["year"] = 2023
            metadata["initial_time"]["second"] = 12.5
            path = joinpath(dir, "manifest.toml")
            open(io -> TOML.print(io, Dict("frame" => metadata)), path, "w")
            recorded = _backfill_frame_configuration(dir)
            restored = _backfill_args(info, planet_cache, density_cache; frame=recorded)
            @test restored.initial_time.year == 2023
            @test restored.initial_time.second == 12.5
            @test restored.environment_model.ephemerides_model isa SM.SpiceEphemeridesModel
            metadata["ephemerides_model"] = "unknown"
            open(io -> TOML.print(io, Dict("frame" => metadata)), path, "w")
            before = read(path)
            @test_throws ArgumentError backfill_run!(dir; force=true)
            @test read(path) == before
            metadata["ephemerides_model"] = "SpiceEphemeridesModel"
            delete!(metadata["initial_time"], "hour")
            open(io -> TOML.print(io, Dict("frame" => metadata)), path, "w")
            @test_throws ArgumentError _backfill_frame_configuration(dir)
        end
    end
end
