function _entry_orbital_elements_from_gamma_v_h(
    planet::AbstractPlanet;
    gamma_deg::Float64,
    v_mps::Float64,
    h_m::Float64
)::NamedTuple{(:a, :e, :ν_deg), Tuple{Float64, Float64, Float64}}
    v_mps > 0.0 || throw(ArgumentError("Entry speed must be > 0 m/s (got $v_mps)."))
    h_m >= 0.0 || throw(ArgumentError("Entry altitude h must be >= 0 m (got $h_m)."))

    r = planet.Rp_e + h_m
    γ = deg2rad(gamma_deg)
    μ = planet.μ
    specific_energy = 0.5 * v_mps^2 - μ / r
    abs(specific_energy) > 1e-12 || throw(ArgumentError("Parabolic entry state is unsupported (specific energy ~ 0)."))

    a = -μ / (2.0 * specific_energy)
    h_spec = r * v_mps * cos(γ)
    p = h_spec^2 / μ
    e_sq = 1.0 - p / a
    e_sq >= -1e-10 || throw(ArgumentError("Invalid entry state (computed e^2=$e_sq < 0). Check gamma/v/h."))
    e = sqrt(max(0.0, e_sq))

    ν_deg = if e <= 1e-10
        180.0
    else
        cos_ν = clamp((p / r - 1.0) / e, -1.0, 1.0)
        rdot = v_mps * sin(γ)
        sin_ν = clamp((rdot * h_spec) / (μ * e), -1.0, 1.0)
        rad2deg(mod(atan(sin_ν, cos_ν), 2π))
    end

    return (a=a, e=e, ν_deg=ν_deg)
end

function make_spacecraft(
    planet::AbstractPlanet;
    id::Int=1,
    ra_alt_m::Float64=550e3,
    rp_alt_m::Float64=500e3,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=170.0,
    with_panel::Bool=true,
    panel_count::Int=1,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing,
    root_mass::Float64=500.0,
    root_area::Float64=12.0,
    panel_mass::Float64=30.0,
    panel_area::Float64=6.0,
    panel_offset_y::Float64=1.2,
    prop_mass::Float64=0.0
)::SpacecraftModel
    root = Link{0}(root=true, m=root_mass, ref_area=root_area)
    links = Link[root]
    if with_panel
        panel_count >= 1 || throw(ArgumentError("panel_count must be >= 1 when with_panel=true (got $panel_count)."))
        if panel_count == 1
            panel = Link{0}(root=false, m=panel_mass, ref_area=panel_area, r=MVector{3, Float64}(0.0, panel_offset_y, 0.0))
            push!(links, panel)
        else
            # Spread appended bodies around the root to keep a balanced 5-body orientation benchmark.
            for panel_idx in 1:panel_count
                θ = 2π * (panel_idx - 1) / panel_count
                panel = Link{0}(
                    root=false,
                    m=panel_mass,
                    ref_area=panel_area,
                    r=MVector{3, Float64}(panel_offset_y * cos(θ), panel_offset_y * sin(θ), 0.0)
                )
                push!(links, panel)
            end
        end
    end

    ra = planet.Rp_e + ra_alt_m
    rp = planet.Rp_e + rp_alt_m

    ic = if isnothing(orientation_state)
        InitialCondition(ra=ra, rp=rp, i=i_deg, ω=ω_deg, Ω=Ω_deg, ν=ν_deg)
    else
        q0, w0 = orientation_state
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = sum(link.m for link in links)
    return SpacecraftModel(
        Joint[],
        links,
        root,
        true,
        dry_mass,
        prop_mass,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end

function make_blunted_cone_entry_spacecraft(
    planet::AbstractPlanet;
    id::Int=1,
    gamma_deg::Float64=-11.5,
    v_mps::Float64=5500.0,
    h_m::Float64=130e3,
    i_deg::Float64=51.6,
    ω_deg::Float64=30.0,
    Ω_deg::Float64=25.0,
    root_mass::Float64=320.0,
    base_radius_m::Float64=0.89,
    body_length_m::Float64=1.2,
    reflection_coefficient::Float64=1.0,
    prop_mass::Float64=0.0
)::SpacecraftModel
    base_radius_m > 0.0 || throw(ArgumentError("base_radius_m must be > 0 (got $base_radius_m)."))
    body_length_m > 0.0 || throw(ArgumentError("body_length_m must be > 0 (got $body_length_m)."))

    entry_oe = _entry_orbital_elements_from_gamma_v_h(
        planet;
        gamma_deg=gamma_deg,
        v_mps=v_mps,
        h_m=h_m
    )

    root = Link{0}(
        root=true,
        m=root_mass,
        ref_area=π * base_radius_m^2,
        dims=MVector{3, Float64}(body_length_m, 2.0 * base_radius_m, 2.0 * base_radius_m),
        reflection_coefficient=reflection_coefficient
    )
    links = Link[root]
    ic = InitialCondition(entry_oe.a, entry_oe.e, i_deg, ω_deg, Ω_deg, entry_oe.ν_deg)
    dry_mass = sum(link.m for link in links)

    return SpacecraftModel(
        Joint[],
        links,
        root,
        true,
        dry_mass,
        prop_mass,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end

function make_constellation(
    planet::AbstractPlanet,
    n::Int;
    with_panel::Bool=false,
    panel_count::Int=1
)::Vector{SpacecraftModel}
    sats = SpacecraftModel[]
    for i in 1:n
        ra_alt = 540e3 + 20e3 * (i - 1)
        rp_alt = 470e3 + 15e3 * (i - 1)
        if rp_alt >= ra_alt
            ra_alt = rp_alt + 8e3
        end
        ν = 140.0 + 180.0 * (i - 1) / n
        push!(
            sats,
            make_spacecraft(
                planet;
                id=i,
                ra_alt_m=ra_alt,
                rp_alt_m=rp_alt,
                ν_deg=ν,
                with_panel=with_panel,
                panel_count=panel_count
            )
        )
    end
    return sats
end

function build_config(;
    planet::AbstractPlanet,
    spacecraft::Vector{SpacecraftModel},
    mission_time_s::Float64,
    orientation_sim::Bool,
    dynamic_effectors::Tuple,
    mission_type::MissionType=MissionTime,
    mission_keplerian::Bool=true,
    mission_orbits::Int=1,
    density_model=NoAtmosphereModel(),
    guidance_effectors::Tuple=(),
    guidance_rates::Vector{Float64}=Float64[],
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    dt_max_orbit::Float64=1.0,
    reltol_orbit::Float64=1e-9,
    abstol_orbit::Float64=1e-9
)::SimulationConfiguration
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=mission_type,
            keplerian=mission_keplerian,
            number_of_orbits=mission_orbits,
            mission_time=mission_time_s,
            orientation_sim=orientation_sim,
            num_steps_to_save=400
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=guidance_effectors, guidance_rates=guidance_rates),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=reltol_orbit,
            abstol_orbit=abstol_orbit,
            dt_max_orbit=dt_max_orbit
        )
    )
end

Base.@kwdef struct _GRAMOfflineSurrogateFileBase
    planet_name::String = "earth"
end

function SimulationModel.EnvironmentModels._gram_point_density(
    model::_GRAMOfflineSurrogateFileBase,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    # Offline benchmark path: keep density lookup file-backed and avoid native point GRAM calls.
    return 0.0, 200.0, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function _build_earth_gram_surrogate_density()
    isfile(EARTH_GRAM_SURROGATE_FILE) || throw(ArgumentError("GRAM surrogate file not found: $(EARTH_GRAM_SURROGATE_FILE)"))
    base_model = _GRAMOfflineSurrogateFileBase(planet_name="earth")
    return GRAMAtmosphereModelSurrogate(base_model, EARTH_GRAM_SURROGATE_FILE, nothing)
end

function _build_earth_gram_point_density()
    gram_root = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0")
    # GRAM Python bindings can be sensitive to world-age in long-lived Julia sessions.
    return Base.invokelatest(
        GRAMAtmosphereModel;
        gram_directory=gram_root,
        gram_data_directory=gram_root,
        spice_directory=SPICE_PATH,
        planet_name="earth",
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
    )
end

function _build_mars_gram_point_density()
    gram_root = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0")
    # GRAM Python bindings can be sensitive to world-age in long-lived Julia sessions.
    return Base.invokelatest(
        GRAMAtmosphereModel;
        gram_directory=gram_root,
        gram_data_directory=gram_root,
        spice_directory=SPICE_PATH,
        planet_name="mars",
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
    )
end

function make_montecarlo_config(seed::Int, planet::Earth, mission_time_s::Float64)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 530e3 + randn(rng) * 20e3
    rp_alt = 490e3 + randn(rng) * 20e3
    rp_alt = max(rp_alt, 120e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 8e3
    end

    sc = make_spacecraft(
        planet;
        id=1,
        with_panel=false,
        root_mass=160.0,
        root_area=1.5,
        prop_mass=20.0,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=28.0 + randn(rng) * 0.8,
        ω_deg=15.0 + randn(rng) * 2.0,
        Ω_deg=20.0 + randn(rng) * 2.0,
        ν_deg=160.0 + randn(rng) * 8.0
    )

    thruster = BaseThrusterModel(
        thrust=[0.6 + rand(rng) * 0.5],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[0.0],
        stop_burn_time=[180.0 + rand(rng) * 40.0],
        Isp=[280.0 + rand(rng) * 40.0]
    )

    return build_config(
        planet=planet,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster,),
        control_rates=[1.0],
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

@inline function _montecarlo_scenario_catalog()::Vector{NamedTuple}
    return NamedTuple[
        (
            variant=:mars_aerobraking,
            name="montecarlo_mars_aerobraking",
            description="Mars aerobraking mission randomized per seed (aero + MarsGRAM point-to-point)"
        ),
        (
            variant=:multi_sat,
            name="montecarlo_multi_sat",
            description="4-spacecraft constellation randomized per seed"
        ),
        (
            variant=:high_accuracy,
            name="montecarlo_high_accuracy",
            description="High-accuracy single-spacecraft mission randomized per seed (L50 harmonics)"
        )
    ]
end

@inline function _active_montecarlo_scenarios()::Vector{NamedTuple}
    return _montecarlo_scenario_catalog()
end

@inline function _montecarlo_batch_mission_time_s(spec::ProfileSpec, variant::Symbol)::Float64
    if variant == :mars_aerobraking
        # Keep Mars aerobraking seeds long enough to include drag-passage behavior.
        return max(spec.montecarlo_mission_s, spec.mission_long_s)
    end
    return spec.montecarlo_mission_s
end

function make_montecarlo_mars_aerobraking_config(
    seed::Int,
    mars::Mars,
    mission_time_s::Float64
)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 4500e3 + randn(rng) * 220e3
    rp_alt = clamp(120e3 + randn(rng) * 18e3, 95e3, 180e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 120e3
    end

    sc = make_spacecraft(
        mars;
        id=1,
        with_panel=false,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=93.0 + randn(rng) * 0.8,
        ω_deg=80.0 + randn(rng) * 3.0,
        Ω_deg=30.0 + randn(rng) * 3.0,
        ν_deg=180.0 + randn(rng) * 10.0
    )

    return build_config(
        planet=mars,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        mission_keplerian=false,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        density_model=_perf_mars_gram_point_density_model(),
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

function make_montecarlo_multi_sat_config(
    seed::Int,
    planet::Earth,
    mission_time_s::Float64
)::SimulationConfiguration
    rng = MersenneTwister(seed)
    spacecraft = SpacecraftModel[]
    for i in 1:4
        ra_alt = 540e3 + 20e3 * (i - 1) + randn(rng) * 8e3
        rp_alt = 470e3 + 15e3 * (i - 1) + randn(rng) * 8e3
        rp_alt = max(rp_alt, 120e3)
        if rp_alt >= ra_alt
            ra_alt = rp_alt + 8e3
        end
        ν = 140.0 + 180.0 * (i - 1) / 4 + randn(rng) * 5.0
        push!(
            spacecraft,
            make_spacecraft(
                planet;
                id=i,
                with_panel=false,
                ra_alt_m=ra_alt,
                rp_alt_m=rp_alt,
                i_deg=35.0 + randn(rng) * 0.4,
                ω_deg=40.0 + randn(rng) * 1.0,
                Ω_deg=10.0 + randn(rng) * 1.0,
                ν_deg=ν
            )
        )
    end
    harmonics20 = _perf_harmonics20_model(planet)
    return build_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(), harmonics20),
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

function make_montecarlo_high_accuracy_config(
    seed::Int,
    planet::Earth,
    mission_time_s::Float64
)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 550e3 + randn(rng) * 12e3
    rp_alt = 500e3 + randn(rng) * 12e3
    rp_alt = max(rp_alt, 140e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 8e3
    end

    sc = make_spacecraft(
        planet;
        id=1,
        with_panel=false,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=35.0 + randn(rng) * 0.3,
        ω_deg=40.0 + randn(rng) * 0.8,
        Ω_deg=10.0 + randn(rng) * 0.8,
        ν_deg=170.0 + randn(rng) * 3.0
    )
    harmonics50 = _perf_harmonics50_model(planet)
    return build_config(
        planet=planet,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(), harmonics50),
        dt_max_orbit=0.5,
        reltol_orbit=1e-10,
        abstol_orbit=1e-10
    )
end

function make_montecarlo_case(
    seed::Int,
    mission_time_s::Float64,
    variant::Symbol,
    planet::Earth;
    mars::Union{Nothing, Mars}=nothing
)::BenchmarkCase
    catalog = _active_montecarlo_scenarios()
    scenario_idx = findfirst(s -> s.variant == variant, catalog)
    scenario_idx === nothing && throw(ArgumentError("Unsupported Monte Carlo scenario variant '$variant'."))
    scenario_meta = catalog[scenario_idx]

    args_template = if variant == :mars_aerobraking
        mars_planet = isnothing(mars) ? perf_worker_mars() : mars
        make_montecarlo_mars_aerobraking_config(seed, mars_planet, mission_time_s)
    elseif variant == :multi_sat
        make_montecarlo_multi_sat_config(seed, planet, mission_time_s)
    elseif variant == :high_accuracy
        make_montecarlo_high_accuracy_config(seed, planet, mission_time_s)
    else
        throw(ArgumentError("Unsupported Monte Carlo scenario variant '$variant'."))
    end

    return BenchmarkCase(
        name=scenario_meta.name,
        category="montecarlo",
        description=scenario_meta.description,
        args_template=args_template,
        run_in_quick=true
    )
end

function build_cases(spec::ProfileSpec, planet::Earth)::Vector{BenchmarkCase}
    harmonics20 = GravitationalHarmonicsModel(20, 20, EARTH_HARMONICS_FILE, planet)
    harmonics50 = GravitationalHarmonicsModel(50, 50, EARTH_HARMONICS_FILE, planet)
    nbody_sun_moon = NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet)
    mars = Mars("", SPICE_PATH)

    q0 = normalize(SVector{4, Float64}(0.15, -0.05, 0.2, 0.96))
    w0 = SVector{3, Float64}(0.01, -0.02, 0.015)
    q1 = normalize(SVector{4, Float64}(0.12, -0.08, 0.24, 0.96))
    w1 = SVector{3, Float64}(0.012, -0.018, 0.02)

    sc_baseline = [make_spacecraft(planet; id=1, with_panel=false)]
    sc_orientation = [make_spacecraft(planet; id=1, with_panel=true, panel_count=4, orientation_state=(q0, w0))]
    sc_entry_shallow = [make_blunted_cone_entry_spacecraft(
        planet;
        id=1,
        gamma_deg=-8.5,
        v_mps=5200.0,
        h_m=135e3,
        root_mass=320.0,
        base_radius_m=0.89,
        body_length_m=1.2,
        prop_mass=0.0,
        i_deg=51.6,
        ω_deg=30.0,
        Ω_deg=25.0
    )]
    sc_entry_nominal = [make_blunted_cone_entry_spacecraft(
        planet;
        id=1,
        gamma_deg=-11.5,
        v_mps=5500.0,
        h_m=130e3,
        root_mass=320.0,
        base_radius_m=0.89,
        body_length_m=1.2,
        prop_mass=0.0,
        i_deg=51.6,
        ω_deg=30.0,
        Ω_deg=25.0
    )]
    sc_entry_steep = [make_blunted_cone_entry_spacecraft(
        planet;
        id=1,
        gamma_deg=-14.5,
        v_mps=5800.0,
        h_m=125e3,
        root_mass=320.0,
        base_radius_m=0.89,
        body_length_m=1.2,
        prop_mass=0.0,
        i_deg=51.6,
        ω_deg=30.0,
        Ω_deg=25.0
    )]
    sc_mars_aerobrake = [make_spacecraft(
        mars;
        id=1,
        with_panel=false,
        ra_alt_m=4500e3,
        rp_alt_m=120e3,
        i_deg=93.0,
        ω_deg=80.0,
        Ω_deg=30.0,
        ν_deg=180.0
    )]
    earth_gram_point_density = _build_earth_gram_point_density()
    mars_gram_point_density = _build_mars_gram_point_density()
    earth_gram_surrogate_density = _build_earth_gram_surrogate_density()
    thermal_stress_density = SimulationModel.EnvironmentModels.PolynomialFitAtmosphereModel([-27.0])
    multi_scaling_effectors = (InverseSquaredGravityModel(), harmonics20)
    sc_thermal_stress = make_constellation(planet, 8; with_panel=true, panel_count=12)
    sc_thermal_aerobrake = [make_spacecraft(
        mars;
        id=1,
        with_panel=true,
        panel_count=16,
        ra_alt_m=4500e3,
        rp_alt_m=120e3,
        i_deg=93.0,
        ω_deg=80.0,
        Ω_deg=30.0,
        ν_deg=180.0,
        root_mass=340.0,
        root_area=5.0,
        panel_mass=6.0,
        panel_area=1.8,
        panel_offset_y=1.2
    )]
    sc_srp_heavy = [make_spacecraft(
        planet;
        id=1,
        with_panel=true,
        panel_count=8,
        ra_alt_m=700e3,
        rp_alt_m=680e3,
        i_deg=97.0,
        ω_deg=10.0,
        Ω_deg=45.0,
        ν_deg=30.0,
        root_mass=600.0,
        root_area=40.0,
        panel_mass=20.0,
        panel_area=35.0,
        panel_offset_y=2.5
    )]
    srp_heavy_effectors = (
        InverseSquaredGravityModel(),
        SolarRadiationPressureModel(1.2, 120.0),
        SolarRadiationPressureModel(1.8, 220.0)
    )
    sc_articulated_heavy = [make_spacecraft(
        planet;
        id=1,
        with_panel=true,
        panel_count=28,
        orientation_state=(q0, w0),
        ra_alt_m=520e3,
        rp_alt_m=500e3,
        ν_deg=160.0,
        root_mass=450.0,
        root_area=10.0,
        panel_mass=8.0,
        panel_area=3.0,
        panel_offset_y=2.0
    )]
    sc_multi_sat_control = make_constellation(planet, 8; with_panel=false)
    constellation_thruster = BaseThrusterModel(
        thrust=fill(0.18, 8),
        direction=fill(0.0, 8),
        Δv=fill(0.0, 8),
        start_burn_time=repeat([120.0, 540.0, 960.0, 1380.0], 2),
        stop_burn_time=repeat([180.0, 600.0, 1020.0, 1440.0], 2),
        Isp=fill(285.0, 8)
    )
    sc_long_constellation = make_constellation(planet, 12; with_panel=false)
    sc_effector_stress6 = make_constellation(planet, 6; with_panel=false)
    sc_effector_stress12 = make_constellation(planet, 12; with_panel=false)
    effector_stress_effectors = (
        InverseSquaredJ2GravityModel(),
        SolarRadiationPressureModel(1.2, 16.0),
        SolarRadiationPressureModel(1.6, 24.0)
    )
    sc_proximity_fullstack = [
        make_spacecraft(
            planet;
            id=1,
            with_panel=true,
            orientation_state=(q0, w0),
            ra_alt_m=520e3,
            rp_alt_m=515e3,
            ν_deg=168.0,
            root_mass=280.0,
            root_area=4.0,
            panel_mass=12.0,
            panel_area=1.8,
            panel_offset_y=1.0,
            prop_mass=18.0
        ),
        make_spacecraft(
            planet;
            id=2,
            with_panel=true,
            orientation_state=(q1, w1),
            ra_alt_m=520.15e3,
            rp_alt_m=515.1e3,
            ν_deg=168.18,
            root_mass=280.0,
            root_area=4.0,
            panel_mass=12.0,
            panel_area=1.8,
            panel_offset_y=1.0,
            prop_mass=18.0
        )
    ]
    proximity_thruster = BaseThrusterModel(
        thrust=fill(0.28, 2),
        direction=fill(0.0, 2),
        Δv=fill(2.5, 2),
        start_burn_time=fill(-1.0, 2),
        stop_burn_time=fill(-1.0, 2),
        Isp=fill(290.0, 2)
    )
    proximity_guidance = (
        AerobrakingCampaignPropulsiveManeuverGuidanceModel(
            maneuver_orbit_number=[1],
            maneuver_Δv=[2.5]
        ),
        AerobrakingCampaignPropulsiveManeuverGuidanceModel(
            maneuver_orbit_number=[1],
            maneuver_Δv=[2.5]
        )
    )

    cases = BenchmarkCase[
        BenchmarkCase(
            name="single_orientation_aero",
            category="orientation",
            description="1 spacecraft (5-body), orientation dynamics on, aerodynamic model active",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_orientation,
                mission_time_s=spec.mission_short_s,
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM())
            )
        ),
        BenchmarkCase(
            name="thermal_8sat_panel12_aero",
            category="thermal_stress",
            description="8 spacecraft (13 links each) with aerodynamic model and fixed polynomial atmosphere to stress thermal callback throughput",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_thermal_stress,
                mission_time_s=min(spec.mission_short_s, 3600.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(thermal_stress_density),
                dt_max_orbit=0.5
            )
        ),
        BenchmarkCase(
            name="thermal_aerobrake_mars_panel16",
            category="thermal_entry",
            description="1 articulated spacecraft (17 links) in Mars aerobraking regime (10 orbits, aero + MarsGRAM point-to-point) to stress thermal callback under entry-like heating",
            args_template=build_config(
                planet=mars,
                spacecraft=sc_thermal_aerobrake,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                mission_type=MissionOrbits,
                mission_keplerian=false,
                mission_orbits=10,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(mars_gram_point_density),
                dt_max_orbit=1.0
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="srp_heavy_high_area",
            category="srp_heavy",
            description="1 high-area spacecraft (9 links) with stacked SRP effectors to stress SRP-heavy workloads",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_srp_heavy,
                mission_time_s=spec.mission_short_s,
                orientation_sim=true,
                dynamic_effectors=srp_heavy_effectors,
                dt_max_orbit=1.0
            )
        ),
        BenchmarkCase(
            name="articulated_1sat_panel28_fullstack",
            category="articulated_multibody",
            description="1 heavily articulated spacecraft (29 links), orientation dynamics on, harmonics + aero active to stress multibody/link-level kernels",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_articulated_heavy,
                mission_time_s=min(spec.mission_short_s, 2400.0),
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20, AerodynamicCoefficientfM()),
                density_model=deepcopy(thermal_stress_density),
                dt_max_orbit=0.5
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_sat_control_8sat_thruster",
            category="multi_sat_control",
            description="8-spacecraft constellation with active thruster control callbacks to capture multi-satellite + active control behavior",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_multi_sat_control,
                mission_time_s=min(spec.mission_short_s, 2400.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(),),
                control_effectors=(constellation_thruster,),
                control_rates=[0.2],
                dt_max_orbit=0.5
            )
        ),
        BenchmarkCase(
            name="long_constellation_12sat",
            category="long_constellation",
            description="12-spacecraft long-duration constellation with L20 harmonics + SRP and GRAM surrogate density",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_long_constellation,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20, SolarRadiationPressureModel(1.2, 12.0)),
                density_model=deepcopy(earth_gram_surrogate_density),
                dt_max_orbit=2.0
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="effector_6sat_dual_srp_stack",
            category="effector_stress",
            description="6 spacecraft with J2 + dual SRP effectors (no atmosphere/control) to stress dynamic effector reduction with outer routing off (r2)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_effector_stress6,
                mission_time_s=min(spec.mission_short_s, 3600.0),
                orientation_sim=false,
                dynamic_effectors=effector_stress_effectors,
                density_model=NoAtmosphereModel(),
                dt_max_orbit=0.5
            )
        ),
        BenchmarkCase(
            name="effector_12sat_dual_srp_stack",
            category="effector_stress",
            description="12 spacecraft with J2 + dual SRP effectors (no atmosphere/control) for larger multi-satellite effector scaling",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_effector_stress12,
                mission_time_s=min(spec.mission_short_s, 3600.0),
                orientation_sim=false,
                dynamic_effectors=effector_stress_effectors,
                density_model=NoAtmosphereModel(),
                dt_max_orbit=0.5
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_entry_earth_shallow",
            category="entry",
            description="1 blunted-cone entry spacecraft, Earth shallow entry from gamma/v/h (target 1 entry interface downcrossing, drag + aero, Earth GRAM point-to-point)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_entry_shallow,
                mission_time_s=max(spec.mission_short_s, 2400.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(earth_gram_point_density),
                dt_max_orbit=0.5
            ),
            entry_target_count_override=1
        ),
        BenchmarkCase(
            name="single_entry_earth_nominal",
            category="entry",
            description="1 blunted-cone entry spacecraft, Earth nominal entry from gamma/v/h (target 1 entry interface downcrossing, drag + aero, Earth GRAM point-to-point)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_entry_nominal,
                mission_time_s=max(spec.mission_short_s, 1800.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(earth_gram_point_density),
                dt_max_orbit=0.5
            ),
            entry_target_count_override=1
        ),
        BenchmarkCase(
            name="single_entry_earth_steep",
            category="entry",
            description="1 blunted-cone entry spacecraft, Earth steep entry from gamma/v/h (target 1 entry interface downcrossing, drag + aero, Earth GRAM point-to-point)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_entry_steep,
                mission_time_s=max(spec.mission_short_s, 1200.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(earth_gram_point_density),
                dt_max_orbit=0.25
            ),
            entry_target_count_override=1
        ),
        BenchmarkCase(
            name="multi_4_gravity",
            category="satellite_scaling",
            description="4 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 4; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            )
        ),
        BenchmarkCase(
            name="multi_8_gravity",
            category="satellite_scaling",
            description="8 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 8; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_8_gravity_surrogate_cached",
            category="satellite_scaling",
            description="8 spacecraft, L20 harmonics with GRAM surrogate density and track-cache enabled",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 8; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false,
            env_overrides=Pair{String, String}[
                "SPACEAGORA_GRAM_TRACK_CACHE" => "on"
            ]
        ),
        BenchmarkCase(
            name="multi_16_gravity",
            category="satellite_scaling",
            description="16 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 16; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_32_gravity",
            category="satellite_scaling",
            description="32 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 32; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_64_gravity",
            category="satellite_scaling",
            description="64 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 64; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_j2",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square + J2 gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(),)
            )
        ),
        BenchmarkCase(
            name="single_nbody_sun_moon",
            category="dynamics_fidelity",
            description="1 spacecraft, J2 gravity + N-body Sun/Moon perturbations",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(), nbody_sun_moon),
                dt_max_orbit=10.0
            )
        ),
        BenchmarkCase(
            name="single_harmonics_l20",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=20",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20)
            )
        ),
        BenchmarkCase(
            name="single_harmonics_l50",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=50",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics50)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="proximity_2sat_orientation_fullstack_gnc_highrate",
            category="rpo_gnc",
            description="2-spacecraft close-proximity operations with orientation on, high-rate guidance, and BaseThrusterModel control callback",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_proximity_fullstack,
                mission_time_s=min(spec.mission_short_s, 2400.0),
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(),),
                guidance_effectors=proximity_guidance,
                guidance_rates=[0.1, 0.1],
                control_effectors=(proximity_thruster,),
                control_rates=[0.1],
                dt_max_orbit=0.2
            )
        ),
        BenchmarkCase(
            name="single_baseline_long_mission",
            category="mission_length",
            description="1 spacecraft Mars aerobraking-style long mission (10 orbits, atmospheric drag + aero, MarsGRAM point-to-point)",
            args_template=build_config(
                planet=mars,
                spacecraft=sc_mars_aerobrake,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                mission_type=MissionOrbits,
                mission_keplerian=false,
                mission_orbits=10,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(mars_gram_point_density),
                dt_max_orbit=1.0
            )
        )
    ]

    return cases
end

@inline function _split_variant_case_name(base::AbstractString, solver::AbstractString)::String
    return string(base, "__split_imex_", lowercase(String(solver)))
end

@inline function _split_base_scenario_name(name::AbstractString)::String
    token = "__split_imex_"
    idx = findfirst(token, String(name))
    if idx === nothing
        return String(name)
    end
    start_idx = first(idx)
    return String(name)[1:(start_idx - 1)]
end

@inline function _split_is_variant_case(name::AbstractString)::Bool
    return occursin("__split_imex_", String(name))
end

function _split_rollout_benchmark_cases(cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    if !_split_rollout_enabled()
        return cases
    end
    target_names = Set(_split_rollout_case_names())
    split_solvers = _split_rollout_solver_variants()
    expanded = copy(cases)
    for case in cases
        if !(case.name in target_names)
            continue
        end
        for split_solver in split_solvers
            push!(expanded, BenchmarkCase(
                name=_split_variant_case_name(case.name, split_solver),
                category=case.category,
                description=string(case.description, " [split_imex:", split_solver, "]"),
                args_template=case.args_template,
                run_in_quick=case.run_in_quick,
                solver_mode_override="split_imex",
                split_imex_solver_override=split_solver,
                entry_target_count_override=case.entry_target_count_override,
                env_overrides=copy(case.env_overrides)
            ))
        end
    end
    return expanded
end

@inline function _multirate_variant_case_name(base::AbstractString)::String
    return string(base, "__multirate")
end

function _multirate_rollout_benchmark_cases(cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    if !_multirate_rollout_enabled()
        return cases
    end
    target_names = Set(_multirate_rollout_case_names())
    expanded = copy(cases)
    for case in cases
        if !(case.name in target_names)
            continue
        end
        push!(expanded, BenchmarkCase(
            name=_multirate_variant_case_name(case.name),
            category=case.category,
            description=string(case.description, " [multirate]"),
            args_template=case.args_template,
            run_in_quick=case.run_in_quick,
            solver_mode_override="multirate",
            split_imex_solver_override=nothing,
            entry_target_count_override=case.entry_target_count_override,
            env_overrides=copy(case.env_overrides)
        ))
    end
    return expanded
end

