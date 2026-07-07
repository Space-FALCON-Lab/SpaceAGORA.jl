################################################################ DEFINE SPACECRAFT HERE ################################################################
# This function creates a SpacecraftModel with the given id, mass, and initial condition.
function _spacecraft(id::Int, mass_kg::Float64, ic::InitialCondition)
# Takes 3 inputs: a numeric ID, the mass, and an initial orbital condition (position/velocity)

    bus = Link(root=true, m=mass_kg, ref_area=1.0)
    # Creates the spacecraft's single structural body — the "bus" (main body).
    # root=true means it's the top of the body tree (no parent).
    # ref_area=1.0 is the reference area (m²) used for drag/pressure calculations.

    return SpacecraftModel(
        Joint[],       # No joints — this is a rigid single-body spacecraft, no articulated parts
        [bus],         # List of all Links (bodies); just the one bus
        bus,           # The root link (the bus itself)
        true,          # keplerian=true — use Keplerian orbital dynamics
        mass_kg,       # Total mass
        0.0,           # Initial drag coefficient or similar scalar — set to zero here
        bus.inertia,   # Moment of inertia tensor, taken from the bus body
        0,             # Number of reaction wheels = 0
        0,             # Number of thrusters = 0
        ic,            # Initial orbital condition (altitude, inclination, true anomaly, etc.)
        id,            # Unique integer ID for this spacecraft (1 = target, 2..N = helpers)
    )
end
########################################################################### END SPACECRAFT DEFINITION ###############################################################

# Compute cumulative orbit count for the TARGET satellite using the instantaneous
# semi-major axis a(t) at each saved time step. Correctly accounts for orbital period
# changes caused by the laser (a grows → T grows → fewer orbits per second).
# orbit_count(t) = ∫₀ᵗ dt'/T(t')  ≈  cumsum(Δt / T_mid(tₖ))
function _orbit_count_from_sol(sol, mu::Float64)
    T_series = [begin
        sc = u.sc[1]
        r  = SVector{3, Float64}(sc.pos)
        v  = SVector{3, Float64}(sc.vel)
        a  = _rv_to_elements(r, v, mu).a
        2pi * sqrt(a^3 / mu)
    end for u in sol.u]
    dt    = diff(Float64.(sol.t))
    T_mid = (T_series[1:end-1] .+ T_series[2:end]) ./ 2
    return cumsum([0.0; dt ./ T_mid])  # length == length(sol.t)
end

# Builds the configuration for the case 2 of the Oracle simulation.
function build_case_config(opts::OracleCase2Options)
    # 1. Sanity-checks every number the user passed in (altitudes > 0, mass > 0, etc.) 
    _validate_options(opts)
    
    # 2. Planet + radii
    planet = make_no_gram_planet(:earth)
    target_radius_m = planet.Rp_e + opts.target_altitude_km * 1e3
    helper_radius_m = planet.Rp_e + opts.helper_altitude_km * 1e3

    # 3. Build the spacecraft array
    spacecraft = SpacecraftModel[]
    push!(spacecraft, _spacecraft(1, opts.mass_kg, InitialCondition(
        target_radius_m, 0.0, opts.target_inclination_deg, 0.0, 0.0, 0.0
    ))) # 3.1. Creates one target satellite (ID 1)

    helper_inclination_deg = opts.helper_inclination_deg
    for helper in 1:opts.helpers
        nu_deg = 360.0 * (helper - 1) / opts.helpers # evenly space them around the orbit
        push!(spacecraft, _spacecraft(helper + 1, opts.mass_kg, InitialCondition(
            helper_radius_m, 0.0, helper_inclination_deg, 0.0, 0.0, nu_deg
        )))
    end # 3.2. Create N helper satellites (IDs 2…N+1)

    # 4. Build the laser model
    laser_model = OpenCavityLaserLinkModel(
        1,  # target is spacecraft #1
        collect(2:(opts.helpers + 1)); # helpers are #2..N+1
        range_m=opts.laser_range_km * 1e3,
        power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        schedule=opts.schedule,
    )

    # 5. Mission time duration in seconds
    target_period_s = 2pi * sqrt(target_radius_m^3 / planet.μ)  # initial period [s] (for reference)
    mission_time_s = opts.orbits * target_period_s

    # 6. The big constructor — wires everything together into one config object
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
        ), # 6.1. Turns off results files, verbose logging, and plots — this script handles output itself

        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=false,
            num_steps_to_save=1000,
            data_rate=max(10.0, mission_time_s / 1000.0),
        ), # 6.2. Sets the stop condition (MissionTime), how many solution snapshots to save (~1000), and the data output rate

        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
        ), # 6.3. Earth gravity + no atmosphere + simple Sun/Moon ephemerides + thermal model

        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(), laser_model)), # 6.4. forces acting on all spacecraft: J2 gravity + the laser link model

        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]), 
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]), # 6.5, no attitude control, no navigation

        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0), # 6.6. Calendar start date (Jan 1 2026) used for ephemerides
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-12,
            abstol_orbit=1e-12,
            dt_max_orbit=opts.dt_max_s,
        ), # 6.7. Very tight ODE tolerances (1e-12) and a max step size so the integrator doesn't skip over fast events

        solver_config=SolverConfig(solver_mode=:tsit5), # 6.8. 	Use the Tsitouras 5th-order Runge-Kutta solver (:tsit5)
    )
    return args, laser_model, target_period_s
end

# Like build_case_config, but enables SpaceAGORA's native output pipeline.
# results_directory should be the per-scenario subfolder (caller computes the path).
# SpaceAGORA will write simulation_results.feather + .csv + .manifest.toml there automatically.
function build_case_config_native(opts::OracleCase2Options, results_directory::String)
    _validate_options(opts)

    planet = make_no_gram_planet(:earth)
    target_radius_m = planet.Rp_e + opts.target_altitude_km * 1e3
    helper_radius_m = planet.Rp_e + opts.helper_altitude_km * 1e3

    spacecraft = SpacecraftModel[]
    push!(spacecraft, _spacecraft(1, opts.mass_kg, InitialCondition(
        target_radius_m, 0.0, opts.target_inclination_deg, 0.0, 0.0, 0.0
    )))
    helper_inclination_deg = opts.helper_inclination_deg
    for helper in 1:opts.helpers
        nu_deg = 360.0 * (helper - 1) / opts.helpers
        push!(spacecraft, _spacecraft(helper + 1, opts.mass_kg, InitialCondition(
            helper_radius_m, 0.0, helper_inclination_deg, 0.0, 0.0, nu_deg
        )))
    end

    laser_model = OpenCavityLaserLinkModel(
        1,
        collect(2:(opts.helpers + 1));
        range_m=opts.laser_range_km * 1e3,
        power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        schedule=opts.schedule,
    )

    target_period_s = 2pi * sqrt(target_radius_m^3 / planet.μ)
    mission_time_s  = opts.orbits * target_period_s

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true,                            # ← let SpaceAGORA write output files
            verbose=false,
            results_directory=results_directory,     # ← per-scenario subfolder
            generate_plots=false,
            save_csv=!opts.feather_only,              # ← false = feather+manifest only; true = also write .csv
            normalize=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=false,
            num_steps_to_save=1000,
            data_rate=max(10.0, mission_time_s / 1000.0),
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(), laser_model)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-12,
            abstol_orbit=1e-12,
            dt_max_orbit=opts.dt_max_s,
        ),
        solver_config=SolverConfig(solver_mode=:tsit5),
    )
    return args, laser_model, target_period_s
end
