using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using CSV
using DataFrames

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine

# Opt-in results thinning (SPACEAGORA_RESULTS_THIN_STRIDE /
# SPACEAGORA_RESULTS_FINAL_ONLY). See docs/architecture/heap_contention.md.
#
# What this pins: with both settings at their default (off), the written
# results output is byte-identical to before the feature existed
# (`_thin_results_segment` is the identity). With thinning on, the written
# set is a STRICT SUBSET of the default set, and every retained row is
# byte-identical to what the default (unthinned) run would have written for
# that same saved time -- not a tolerance match, the same DataFrame row.

@testset "_thin_results_segment: pure-function contract" begin
    times = collect(0.0:1.0:9.0)              # 10 states, indices 1..10
    data = SM.SaveData[SM.SaveData() for _ in times]

    withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => nothing, "SPACEAGORA_RESULTS_FINAL_ONLY" => nothing) do
        @test SE._results_thin_stride() == 1
        @test SE._results_final_only() == false
        out_t, out_d = SE._thin_results_segment(times, data)
        @test out_t === times                  # identity: same object, not a copy
        @test out_d === data
    end

    withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => "1") do
        out_t, _ = SE._thin_results_segment(times, data)
        @test out_t === times
    end

    withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => "3") do
        @test SE._results_thin_stride() == 3
        out_t, out_d = SE._thin_results_segment(times, data)
        # Every 3rd state (1-based: 1,4,7,10) plus the final state, which is
        # already index 10 here, so no extra append happens in this case.
        @test out_t == [times[1], times[4], times[7], times[10]]
        @test out_d == [data[1], data[4], data[7], data[10]]
        @test issubset(Set(out_t), Set(times))
        @test length(out_t) < length(times)
    end

    withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => "4") do
        out_t, _ = SE._thin_results_segment(times, data)
        # Stride 4 over 10 states (1-based: 1,5,9) does NOT land on the final
        # state (index 10) -- it must still be appended.
        @test out_t == [times[1], times[5], times[9], times[10]]
    end

    withenv("SPACEAGORA_RESULTS_FINAL_ONLY" => "1", "SPACEAGORA_RESULTS_THIN_STRIDE" => "3") do
        @test SE._results_final_only() == true
        out_t, out_d = SE._thin_results_segment(times, data)
        # final_only takes priority over stride when both are set.
        @test out_t == [times[end]]
        @test out_d == [data[end]]
    end

    withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => "0", "SPACEAGORA_RESULTS_FINAL_ONLY" => nothing) do
        @test SE._results_thin_stride() == 1     # <= 1 falls back to no thinning
    end
    withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => "not-a-number") do
        @test SE._results_thin_stride() == 1     # unparseable falls back to no thinning
    end

    withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => nothing, "SPACEAGORA_RESULTS_FINAL_ONLY" => nothing) do
        empty_t, empty_d = SE._thin_results_segment(Float64[], SM.SaveData[])
        @test isempty(empty_t)
        @test isempty(empty_d)
    end
end

function _storage_test_config(results_directory::String; mission_s::Float64=60.0, data_rate::Float64=5.0)
    planet = SM.Earth()
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(ra=planet.Rp_e + 550_000.0, rp=planet.Rp_e + 550_000.0,
                           i=53.0, ω=0.0, Ω=10.0, ν=0.0)
    spacecraft = [SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, 1)]
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true, verbose=false, generate_plots=false, normalize=false,
            save_csv=true, results_directory=results_directory
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_s, orientation_sim=false, num_steps_to_save=20,
            data_rate=data_rate
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=300.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
            ephemerides_model=SimpleEphemeridesModel()
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

@testset "results thinning: default is byte-identical, thinned is a proven subset" begin
    mktempdir() do base_dir
        default_dir = joinpath(base_dir, "default")
        thinned_dir = joinpath(base_dir, "thinned")
        final_dir = joinpath(base_dir, "final")
        mkpath(default_dir); mkpath(thinned_dir); mkpath(final_dir)

        args_default = _storage_test_config(default_dir)
        withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => nothing, "SPACEAGORA_RESULTS_FINAL_ONLY" => nothing) do
            run_simulation(args_default; isolate_state=false)
        end
        default_csv = joinpath(default_dir, "simulation_results.csv")
        @test isfile(default_csv)
        default_df = CSV.read(default_csv, DataFrame)
        @test nrow(default_df) > 6   # enough rows for stride=3 to actually thin something

        args_thinned = _storage_test_config(thinned_dir)
        withenv("SPACEAGORA_RESULTS_THIN_STRIDE" => "3") do
            run_simulation(args_thinned; isolate_state=false)
        end
        thinned_csv = joinpath(thinned_dir, "simulation_results.csv")
        @test isfile(thinned_csv)
        thinned_df = CSV.read(thinned_csv, DataFrame)

        # Strict subset: fewer rows, and every thinned row's full record
        # matches some default row exactly (same column values, not just the
        # timestamp) -- i.e. no retained state differs from what the default,
        # unthinned run wrote for that same saved time.
        @test nrow(thinned_df) < nrow(default_df)
        @test names(thinned_df) == names(default_df)
        default_rows = Set(eachrow(default_df))
        for row in eachrow(thinned_df)
            @test row in default_rows
        end
        # The mission-end time must be present even though stride=3 need not
        # land on it exactly.
        @test maximum(thinned_df.time) == maximum(default_df.time)

        args_final = _storage_test_config(final_dir)
        withenv("SPACEAGORA_RESULTS_FINAL_ONLY" => "1") do
            run_simulation(args_final; isolate_state=false)
        end
        final_csv = joinpath(final_dir, "simulation_results.csv")
        @test isfile(final_csv)
        final_df = CSV.read(final_csv, DataFrame)
        @test nrow(final_df) == 1
        @test only(eachrow(final_df)) in default_rows
    end
end
