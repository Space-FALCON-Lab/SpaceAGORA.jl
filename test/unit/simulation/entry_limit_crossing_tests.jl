module EntryLimitCrossingTests

using Test, Serialization, StaticArrays, SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const CB = SM.SimulationCallbacks

# Zero-density eccentric orbits cross the 500 km interface. Drag enables the
# entry limit, but exerts no force. Different phases exercise successive roots;
# identical phases exercise one simultaneous mask. No SPICE or native GRAM.
function configuration(directory; duration=30000.0, resume=false, simultaneous=false,
                       equal_tolerances=false)
    planet = make_no_gram_planet(:earth)
    spacecraft = map(1:3) do k
        root = Link(root=true, m=100.0, ref_area=2.0)
        ic = InitialCondition(ra=planet.Rp_e + 2000e3, rp=planet.Rp_e + 400e3,
            i=30.0, Ω=40.0k, ω=0.0, ν=simultaneous ? 40.0 : 10.0 + 30.0k)
        SpacecraftModel(links=[root], root=root, initial_condition=ic, id=100+k,
            prop_mass=50.0, inertia_tensor=root.inertia)
    end
    tolerances = equal_tolerances ? IntegrationTolerances(
        reltol_orbit=1e-10, abstol_orbit=1e-10, dt_max_orbit=30.0,
        reltol_atmosphere=1e-10, abstol_atmosphere=1e-10, dt_max_atmosphere=30.0) : IntegrationTolerances()
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, checkpoint_enabled=true,
            checkpoint_interval_s=7000.0, checkpoint_directory=directory,
            resume_from_checkpoint=resume),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, keplerian=true),
        environment_model=EnvironmentModel(planet=planet, EI=500.0,
            density_model=NoAtmosphereModel(), ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=DynamicsModel(spacecraft,
            (InverseSquaredGravityModel(), SM.AerodynamicCoefficientConstant())),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1),
        integration_tolerances=tolerances,
        solver_config=SM.SolverConfig(solver_mode=:tsit5, maxiters=20000))
end

payload(config) = deserialize(SE._checkpoint_paths(config).data)
counts(config) = payload(config).runtime_state.entries.entry_counter

@testset "Entry limits and atmosphere phases share crossing events" begin
    withenv("SPACEAGORA_ENTRY_TARGET_COUNT" => "3") do
        whole_config = configuration(mktempdir())
        whole = run_simulation(whole_config; return_solution=true)
        @test 18000.0 < whole.t[end] < 20000.0
        @test counts(whole_config) == [3, 3, 3]

        # The old two-root path stalls at t≈13009 s, even without a resume.
        # The bounded maxiters makes that regression fail instead of hanging.
        split_dir = mktempdir()
        first_leg = configuration(split_dir; duration=7000.0)
        run_simulation(first_leg)
        @test counts(first_leg) == [1, 1, 1]
        old_payload = payload(first_leg)
        resumed_config = configuration(split_dir; resume=true)
        resumed = run_simulation(resumed_config; return_solution=true)
        @test resumed.t[end] == whole.t[end]
        @test isequal(collect(resumed.u[end]), collect(whole.u[end]))
        @test counts(resumed_config) == [3, 3, 3]

        # A checkpoint without the new count keeps its old zero-count start.
        # It completes three further entries instead of stalling at an exit.
        runtime = Base.structdiff(old_payload.runtime_state, NamedTuple{(:entries,)})
        serialize(SE._checkpoint_paths(first_leg).data, merge(old_payload, (runtime_state=runtime,)))
        legacy = run_simulation(resumed_config; return_solution=true)
        @test whole.t[end] < legacy.t[end] < 30000.0
        @test counts(resumed_config) == [3, 3, 3]

        # Equal phase settings need only the standalone entry-limit callback.
        for equal_tolerances in (false, true)
            config = configuration(mktempdir(); simultaneous=true, equal_tolerances)
            sol = run_simulation(config; return_solution=true)
            @test sol.t[end] < 20000.0
            @test counts(config) == [3, 3, 3]
        end

        # Inactive spacecraft must neither accrue entries nor prevent stopping.
        config = configuration(mktempdir())
        deactivate = CB.PresetTimeCallback([1.0], i -> (i.p.is_active[3] = false);
            save_positions=(false, false))
        sol = run_simulation(config; return_solution=true, extra_callbacks=(deactivate,))
        @test sol.t[end] < 20000.0
        @test counts(config) == [3, 3, 0]
    end
end

end
