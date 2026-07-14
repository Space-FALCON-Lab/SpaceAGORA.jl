# Requires GRAMSuite (guarded at the call site via HAS_GRAMSUITE, set up by
# test/helpers/bootstrap.jl). Modeled directly on the working end-to-end Mars
# GRAM scenario at
# benchmarks/studies/gram_mars_fix_and_constellation_scaling/mars_150km_gram_scaling.jl,
# which exists specifically to exercise the GRAM Mars ephemeris-state bypass
# documented there and in memory [[project_gram_mars_isolated_cspice_bug]].
function build_gram_mars_config()
    planet = Mars("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "Mars50c.csv"), planet
    )

    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=planet.Rp_e + 150e3,
        rp=planet.Rp_e + 150e3,
        i=53.0,
        ω=0.0,
        Ω=10.0,
        ν=0.0
    )
    sc = SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, 1)

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=600.0,
            orientation_sim=false,
            num_steps_to_save=20
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=300.0,
            density_model=Base.invokelatest(GRAMAtmosphereModel; planet_name="mars"),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel([sc], (harmonics, AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end
