# Shared scenario for gram_wind_determinism_tests.jl and the child processes it
# starts. Defines functions only; the including module supplies `SM`, `EM` and
# `WD_SPICE_PATH`.

# Low-perigee LEO spacecraft with drag and winds on, so GRAM's perturbed winds
# reach the trajectory through the relative velocity.
function wind_case_args(density_model; n::Int=2, mission_time::Float64=60.0, wind::Bool=true)
    planet = SM.Earth("", WD_SPICE_PATH)
    spacecraft = SM.SpacecraftModel[]
    for i in 1:n
        root = SM.Link(root=true, m=500.0, ref_area=12.0)
        ic = SM.InitialCondition(ra=planet.Rp_e + 400e3, rp=planet.Rp_e + 130e3,
            i=35.0, ω=40.0, Ω=360.0 * (i - 1) / n, ν=-5.0)
        push!(spacecraft, SM.SpacecraftModel(SM.Joint[], SM.Link[root], root, true,
            root.m, 0.0, root.inertia, 0, 0, ic, i))
    end
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=mission_time,
            orientation_sim=false, num_steps_to_save=10, data_rate=10.0),
        environment_model=SM.EnvironmentModel(planet=planet, EI=600.0,
            density_model=density_model,
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=wind),
        dynamics_model=SM.DynamicsModel(spacecraft,
            (SM.InverseSquaredGravityModel(), SM.AerodynamicCoefficientfM())),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9,
            abstol_orbit=1e-9, dt_max_orbit=5.0))
end

wind_case_model() = EM.GRAMAtmosphereModel(planet_name="earth",
    initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))

function wind_case_final(args; isolate_state::Bool=true)
    sol = SpaceAGORA.run_simulation(args; isolate_state=isolate_state, return_solution=true)
    return (t=copy(sol.t), u=collect(sol.u[end]))
end

# Above Earth-GRAM's lower-atmosphere fairing, where the first-atmosphere
# defect lived: the first model updated there in a process got zero wind
# perturbations and a NaN north/south component.
const WIND_CASE_POINT = (132_840.0, deg2rad(19.42), deg2rad(-70.0), 0.0)

# Raw native winds of the first query on a fresh model, NaN not replaced.
function wind_case_first_native_winds(model)
    core = model.core
    GRAMSuite.point_density_state(core, WIND_CASE_POINT..., true)
    w = Base.invokelatest(getfield(core.gram, :get_winds_state), core.gram_atmosphere)
    return (w.perturbedEWWind, w.perturbedNSWind, w.perturbedVerticalWind)
end

# Run once in a fresh child after Earth is loaded. The Mars system SPK can
# change Earth's ephemeris too, so an Earth cache entry at the same epoch must
# be refreshed before the next atmosphere query. Keeping this in the existing
# children makes the load order independent of earlier test files.
function wind_case_kernel_refresh()
    ext = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt)
    planets = SM.Planets
    epoch = SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
    return withenv("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on") do
        ephemeris() = ext._gram_spice_ephemeris_state(
            "earth", epoch, 0.0, rad2deg(0.3), rad2deg(-1.2))
        before = ephemeris()
        furnished_before = copy(planets._FURNISHED_KERNELS)
        SM.Mars("", WD_SPICE_PATH)
        added_kernels = setdiff(planets._FURNISHED_KERNELS, furnished_before)
        mars_kernel_added = any(path -> basename(path) in
            ("mar097_GRAM.bsp", "mar097.bsp"), added_kernels)
        automatic = ephemeris()
        uncached = withenv(ephemeris, "SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "off")

        used = wind_case_model()
        path = [(150.0e3 - 50.0k, 0.3 + 1e-3k, -1.2 + 2e-3k, 0.5k) for k in 0:20]
        query(model) = [EM._gram_point_density(model, point..., true) for point in path]
        reference = query(wind_case_model())
        query(used); query(used)
        reset = EM.reset_density_model_history!(used)
        return (mars_kernel_added=mars_kernel_added, before=before,
            automatic=automatic, uncached=uncached, reset=reset,
            reference=reference, reseeded=query(used))
    end
end
