using Test
using DiffEqCallbacks: SavedValues
using SpaceAGORA
using SpaceAGORA.SimulationModel

const _TR_SM = SpaceAGORA.SimulationModel
const _TR_SC = SpaceAGORA.SimulationModel.SimulationCallbacks

function _trajectory_recorder_test_config(; n_sats::Int=3, mission_s::Float64=60.0, data_rate::Float64=10.0)
    planet = _TR_SM.Earth()
    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link(root=true, m=500.0, ref_area=12.0)
        ic = InitialCondition(
            ra=planet.Rp_e + 550_000.0 + 100.0 * i,
            rp=planet.Rp_e + 550_000.0 + 100.0 * i,
            i=53.0,
            ω=0.0,
            Ω=10.0,
            ν=360.0 * (i - 1) / n_sats,
        )
        push!(
            spacecraft,
            SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, i),
        )
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false,
        ),
        mission_configuration=MissionConfiguration(
            MissionTime,
            true,
            1,
            mission_s,
            false,
            20,
            data_rate,
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=300.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
            ephemerides_model=SimpleEphemeridesModel(),
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            dt_max_orbit=2.0,
        ),
    )
end

function _test_saved_value_matches(got, expected)
    if got isa AbstractArray && expected isa AbstractArray
        @test axes(got) == axes(expected)
        @test got ≈ expected
    else
        @test got == expected
    end
end

@testset "preallocated trajectory recorder" begin
    args = _trajectory_recorder_test_config()
    save_fields = default_save_fields(args)

    reference_saved_values = SavedValues(Float64, _TR_SM.SaveData)
    reference_callback = _TR_SC.get_data_saving_callback(
        length(args.dynamics_model.spacecraft),
        args,
        save_fields,
        reference_saved_values,
    )
    reference_metadata = run_simulation(
        args;
        isolate_state=false,
        return_solver_metadata=true,
        extra_callbacks=(reference_callback,),
    )
    @test reference_metadata.solution === nothing
    @test reference_metadata.retcode == "Success"

    recorder = TrajectoryRecorder(args; capacity=1)
    callback = get_trajectory_recorder_callback(recorder)

    metadata = run_simulation(
        args;
        isolate_state=false,
        return_solver_metadata=true,
        extra_callbacks=(callback,),
    )

    @test metadata.solution === nothing
    @test metadata.retcode == "Success"
    @test recorder.count > 1
    @test collect(trajectory_times(recorder)) ≈ reference_saved_values.t
    @test size(trajectory_positions(recorder)) == (3, 3, recorder.count)
    @test size(trajectory_velocities(recorder)) == (3, 3, recorder.count)
    @test size(trajectory_masses(recorder)) == (3, recorder.count)
    @test all(isfinite, trajectory_positions(recorder))
    @test all(isfinite, trajectory_velocities(recorder))
    @test all(trajectory_masses(recorder) .≈ 500.0)

    recorded_snapshots = trajectory_save_data(recorder)
    @test length(recorded_snapshots) == length(reference_saved_values.saveval)
    for sample_idx in eachindex(recorded_snapshots)
        got = recorded_snapshots[sample_idx]
        expected = reference_saved_values.saveval[sample_idx]
        @test sort!(collect(keys(got)); by=string) == sort!(collect(keys(expected)); by=string)
        for field in save_fields
            _test_saved_value_matches(got[field.name], expected[field.name])
        end
    end

    reset_trajectory_recorder!(recorder)
    @test recorder.count == 0
end
