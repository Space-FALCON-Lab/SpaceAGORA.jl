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


@testset "trajectory recorder validates explicit-capacity cadence" begin
    args = _trajectory_recorder_test_config(n_sats=1, mission_s=6.0, data_rate=2.0)
    for rate in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError TrajectoryRecorder(args; data_rate=rate, capacity=1)
    end
end

@testset "trajectory recorder preserves custom getters and additional defaults" begin
    args = _trajectory_recorder_test_config(n_sats=2, mission_s=6.0, data_rate=2.0)
    fields = SaveField[
        SaveField(:total_mass, (u, t, integrator) -> sum(sc.mass for sc in u.sc)),
        SaveField(:position, (u, t, integrator) -> fill(42.0, length(u.sc)); per_satellite=true),
    ]
    custom = TrajectoryRecorder(args; capacity=1, save_fields=fields)
    defaults = TrajectoryRecorder(args; capacity=1)
    # Model a newly registered default that has no specialized array filler.
    push!(defaults.save_fields, SaveField(:extension_value, (u, t, integrator) -> Float64(t) + 100.0))
    defaults.fallback[:extension_value] = Vector{Any}(undef, length(defaults.t))
    metadata = run_simulation(
        args;
        return_solver_metadata=true,
        extra_callbacks=(get_trajectory_recorder_callback(custom), get_trajectory_recorder_callback(defaults)),
    )
    @test metadata.retcode == "Success"
    @test custom.count > 1
    @test collect(trajectory_times(custom)) == collect(trajectory_times(defaults))
    @test all(value -> value == 1000.0, trajectory_field(custom, :total_mass))
    @test all(value -> value == [42.0, 42.0], trajectory_field(custom, :position))
    @test collect(trajectory_field(defaults, :extension_value)) == collect(trajectory_times(defaults)) .+ 100.0
    snapshots = trajectory_save_data(defaults)
    @test [snapshot[:extension_value] for snapshot in snapshots] == collect(trajectory_times(defaults)) .+ 100.0
    reset_trajectory_recorder!(defaults)
    @test isempty(trajectory_times(defaults))
    @test isempty(trajectory_field(defaults, :extension_value))
end
