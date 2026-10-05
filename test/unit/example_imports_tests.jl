using Test
using SpaceAGORA

@testset "shared example imports preserve caller module aliases" begin
    common_path = normpath(joinpath(@__DIR__, "..", "..", "examples", "common.jl"))
    aliases = quote
        const SimulationEngine = SpaceAGORA.SimulationEngine
        const SimulationModel = SpaceAGORA.SimulationModel
        const RuntimeServices = SpaceAGORA.RuntimeServices
        const SM = SimulationModel
    end

    for predeclared in (false, true)
        @testset "aliases declared before include: $predeclared" begin
            # Separate caller modules match examples and analytical studies.
            # Including common.jl must not load GRAM or start a simulation.
            probe_module = Module(gensym(:ExampleImports))
            Core.eval(probe_module, :(using SpaceAGORA))
            predeclared && Core.eval(probe_module, aliases)
            Base.include(probe_module, common_path)

            # AnalyticalPerturbationModels declares SM after its common include.
            # An imported binding rejects this otherwise valid local constant.
            Core.eval(probe_module, aliases)
            @test probe_module.SimulationEngine === SpaceAGORA.SimulationEngine
            @test probe_module.SimulationModel === SpaceAGORA.SimulationModel
            @test probe_module.RuntimeServices === SpaceAGORA.RuntimeServices
            @test probe_module.SM === SpaceAGORA.SimulationModel
            @test probe_module.run_simulation === SpaceAGORA.run_simulation
            @test probe_module.quat_mult === SpaceAGORA.SimulationModel.quat_mult
            @test probe_module.make_example_config === SpaceAGORA.TelemetryVerification.make_example_config
            @test probe_module.make_three_body_spacecraft === SpaceAGORA.TelemetryVerification.make_three_body_spacecraft
            @test probe_module.run_and_report === SpaceAGORA.TelemetryVerification.run_and_report
        end
    end
end

module ExampleResultTableBoundaryTests
using Test
using SpaceAGORA: SimulationModel
const SimulationConfiguration = SimulationModel.SimulationConfiguration

# The saved-table owner can load independently of common.jl and the mission
# plotting helper. Keep the caller's project and working directory intact.
const PROJECT_BEFORE_INCLUDE = Base.active_project()
const DIRECTORY_BEFORE_INCLUDE = pwd()
include(joinpath(@__DIR__, "..", "..", "examples", "support", "aerobraking_result_tables.jl"))

@testset "saved-result access without mission or plotting setup" begin
    @test Base.active_project() == PROJECT_BEFORE_INCLUDE
    @test pwd() == DIRECTORY_BEFORE_INCLUDE
    for binding in (:Plots, :PlotlyJS, :SPICE, :RuntimeServices, :REPO_ROOT,
                    :AerobrakingMissionSpiceConfig)
        @test !isdefined(@__MODULE__, binding)
    end

    mktempdir() do directory
        csv_path = joinpath(directory, "saved.csv")
        write(csv_path, "other,time\n3,2\n4,-1\n")
        df = _read_simulation_results(csv_path)
        @test propertynames(df) == [:other, :time]
        @test _require_float_column(df, :time) == [2.0, -1.0]
        @test _require_float_column(df, :time) isa Vector{Float64}
        @test_throws ArgumentError _require_float_column(df, :absent)
        @test_throws ArgumentError _read_simulation_results(joinpath(directory, "missing.csv"))
    end
end
end # module ExampleResultTableBoundaryTests

module ExamplePlotResultsTests
using Test

# Load the same definitions as the examples, without running a mission, plotting,
# furnishing kernels, or loading a native atmosphere model.
include(joinpath(@__DIR__, "..", "..", "examples", "common.jl"))
include(joinpath(@__DIR__, "..", "..", "examples", "aerobraking_mission_plot_utils.jl"))

function _fixture_config(results_directory::String)
    craft = make_three_body_spacecraft(
        bus_dims=(1.0, 1.0, 1.0), panel_dims=(0.01, 1.0, 0.5),
        bus_mass=100.0, panel_mass_each=1.0, panel_offset_y=1.0,
        ic=SM.InitialCondition(ra=4.5e6, rp=3.8e6, i=0.0, ω=0.0, Ω=0.0, ν=0.0),
        prop_mass=0.0, id=1
    )
    return make_example_config(
        planet=SpaceAGORA.make_no_gram_planet(:mars), spacecraft=craft, mission_time=40.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.ConstantGravityModel(),), density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(), EI_km=2.0,
        orientation_sim=false, keplerian=true, verbose=false, results=true,
        results_directory=results_directory
    )
end

_capture_error(f) = try
    f()
    nothing
catch err
    err
end

@testset "aerobraking saved-result readers" begin
    mktempdir() do directory
        results_dir = mkdir(joinpath(directory, "results"))
        args = _fixture_config(results_dir)
        csv_path = joinpath(results_dir, "simulation_results.csv")

        @testset "missing results fail before any plotting" begin
            consumers = (
                () -> _derive_orbit_extrema_from_results(csv_path, args.environment_model.planet),
                () -> _save_drag_along_velocity_plot(args),
                () -> _save_aero_sideways_components_plot(args),
                () -> _simulation_position_samples(args),
                () -> _simulation_velocity_samples(args),
                () -> _trajectory_marker_times(args)
            )
            for consume in consumers
                err = _capture_error(consume)
                @test err isa ArgumentError
                @test sprint(showerror, err) ==
                    "ArgumentError: Simulation results CSV not found at $(abspath(csv_path))."
            end
        end

        @testset "independent position and velocity columns, in kilometers" begin
            # Reordered columns and distinct signs/scales expose swapped axes,
            # accidental position/velocity reuse, and changes to time units.
            write(csv_path, "sc1_vel_3,sc1_pos_2,time,sc1_pos_3,sc1_vel_1,sc1_pos_1,sc1_vel_2\n" *
                "-900,-2000,10,3000,100,1000,200\n" *
                "1200,5000,25,-6000,-400,-4000,500\n")
            @test _simulation_position_samples(args) ==
                ([10.0, 25.0], [1.0, -4.0], [-2.0, 5.0], [3.0, -6.0])
            @test _simulation_velocity_samples(args) ==
                ([10.0, 25.0], [0.1, -0.4], [0.2, 0.5], [-0.9, 1.2])

            # Each reader requires only its own vector's columns.
            write(csv_path, "time,sc1_pos_1,sc1_pos_2,sc1_pos_3\n7,2000,-3000,4000\n")
            @test _simulation_position_samples(args) == ([7.0], [2.0], [-3.0], [4.0])
            velocity_error = _capture_error(() -> _simulation_velocity_samples(args))
            @test velocity_error isa ArgumentError
            @test occursin("sc1_vel_1", sprint(showerror, velocity_error))
            write(csv_path, "time,sc1_vel_1,sc1_vel_2,sc1_vel_3\n8,-200,300,-400\n")
            @test _simulation_velocity_samples(args) == ([8.0], [-0.2], [0.3], [-0.4])
            position_error = _capture_error(() -> _simulation_position_samples(args))
            @test position_error isa ArgumentError
            @test occursin("sc1_pos_1", sprint(showerror, position_error))
        end

        @testset "header-only results" begin
            write(csv_path, "time,sc1_pos_1,sc1_pos_2,sc1_pos_3,sc1_vel_1,sc1_vel_2,sc1_vel_3,sc1_altitude\n")
            for sample in (_simulation_position_samples(args), _simulation_velocity_samples(args))
                @test sample == (Float64[], Float64[], Float64[], Float64[])
                @test all(values -> values isa Vector{Float64}, sample)
            end
            @test _trajectory_marker_times(args) == (
                apoapsis_s=Float64[], periapsis_s=Float64[],
                atmosphere_entry_s=Float64[], atmosphere_exit_s=Float64[]
            )
            extrema_error = _capture_error(() -> _derive_orbit_extrema_from_results(csv_path, args.environment_model.planet))
            @test extrema_error isa ArgumentError
            @test sprint(showerror, extrema_error) ==
                "ArgumentError: Need at least 3 saved samples to derive periapsis/apoapsis extrema."
        end

        @testset "trajectory markers and event-table precedence" begin
            write(csv_path, "time,sc1_altitude\n0,3000\n10,1000\n20,3000\n30,1000\n40,3000\n")
            markers = (
                apoapsis_s=[20.0], periapsis_s=[10.0, 30.0],
                atmosphere_entry_s=[5.0, 25.0], atmosphere_exit_s=[15.0, 35.0]
            )
            @test _trajectory_marker_times(args) == markers
            event_path = joinpath(results_dir, "periapsis_events.csv")
            write(event_path, "time_s,altitude_km\n9.5,1\n31.25,2\n")
            @test _trajectory_marker_times(args) == merge(markers, (periapsis_s=[9.5, 31.25],))
            # An existing table without usable times suppresses the fallback.
            for event_table in ("time_s\n", "orbit,altitude_km\n1,1\n")
                write(event_path, event_table)
                @test _trajectory_marker_times(args) == merge(markers, (periapsis_s=Float64[],))
            end
        end

        @testset "orbit extrema retain correction-reset filtering" begin
            write(csv_path, "time,sc1_pos_1,sc1_pos_2,sc1_pos_3,sc1_vel_1,sc1_vel_2,sc1_vel_3,sc1_altitude,sc1_latitude_deg,sc1_longitude_deg\n" *
                "0,4000000,0,0,-100,0,0,3000,10,100\n" *
                "10,4001000,0,0,100,0,0,1000,11,101\n" *
                "20,4002000,0,0,100,0,0,4000,12,102\n" *
                "30,4003000,0,0,-100,0,0,2000,13,103\n" *
                "40,4004000,0,0,100,0,0,5000,14,104\n")
            extrema = _derive_orbit_extrema_from_results(csv_path, args.environment_model.planet)
            @test extrema.peri == (
                orbit=[1, 2], altitude_km=[1.0, 2.0],
                latitude_deg=[11.0, 13.0], longitude_deg=[101.0, 103.0]
            )
            @test extrema.apo == (orbit=[1, 2], altitude_km=[3.0, 3.0], radius_km=[4000.0, 4002.5])
            write(joinpath(directory, "leg_summary.csv"),
                "reset_time_s,correction_applied\n5,true\n25,false\nNaN,true\n")
            corrected = _derive_orbit_extrema_from_results(csv_path, args.environment_model.planet)
            @test corrected.peri == (
                orbit=[1], altitude_km=[2.0], latitude_deg=[13.0], longitude_deg=[103.0]
            )
            @test corrected.apo == extrema.apo
        end
    end
end
end # module ExamplePlotResultsTests
