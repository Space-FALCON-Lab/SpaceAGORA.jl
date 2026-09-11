include(joinpath(PPC_REPO_ROOT, "src", "parallel", "routing", "parallel_profiles.jl"))
include(joinpath(PPC_REPO_ROOT, "src", "simulation", "runtime_services.jl"))
include(joinpath(PPC_REPO_ROOT, "src", "core", "simulation_model.jl"))
include(joinpath(PPC_REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))

using .SimulationModel

const PPC_SPICE_PATH = joinpath(PPC_REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
const PPC_EARTH_HARMONICS_FILE = joinpath(PPC_REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")

Base.@kwdef struct PPCCaseSpec
    name::String
    family::String
    description::String
    builder::Function
    montecarlo::Bool = false
    orientation::Bool = false
    default_samples::Int = 1
end

function ppc_spacecraft(
    planet;
    id::Int=1,
    ra_alt_m::Float64=550e3,
    rp_alt_m::Float64=500e3,
    i_deg::Float64=35.0,
    omega_deg::Float64=40.0,
    raan_deg::Float64=10.0,
    nu_deg::Float64=170.0,
    with_panel::Bool=false,
    panel_count::Int=1,
    root_mass::Float64=500.0,
    root_area::Float64=12.0,
    prop_mass::Float64=0.0,
    orientation_state=nothing
)
    root = Link{0}(root=true, m=root_mass, ref_area=root_area)
    links = Link[root]
    if with_panel
        for panel_idx in 1:panel_count
            theta = 2π * (panel_idx - 1) / max(1, panel_count)
            panel = Link{0}(
                root=false,
                m=8.0,
                ref_area=3.0,
                r=MVector{3, Float64}(1.8 * cos(theta), 1.8 * sin(theta), 0.0)
            )
            push!(links, panel)
        end
    end
    ra = planet.Rp_e + ra_alt_m
    rp = planet.Rp_e + rp_alt_m
    ic = if orientation_state === nothing
        InitialCondition(ra=ra, rp=rp, i=i_deg, ω=omega_deg, Ω=raan_deg, ν=nu_deg)
    else
        q0, w0 = orientation_state
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        InitialCondition(a, e, i_deg, omega_deg, raan_deg, nu_deg, q0, w0)
    end
    dry_mass = sum(link.m for link in links)
    return SpacecraftModel(Joint[], links, root, true, dry_mass, prop_mass, root.inertia, 0, 0, ic, id)
end

function ppc_constellation(planet, n::Int; with_panel::Bool=false, panel_count::Int=1)
    sats = SpacecraftModel[]
    for i in 1:n
        push!(sats, ppc_spacecraft(
            planet;
            id=i,
            ra_alt_m=540e3 + 2e3 * (i - 1),
            rp_alt_m=500e3 + 1e3 * (i - 1),
            nu_deg=120.0 + 240.0 * (i - 1) / max(1, n),
            with_panel=with_panel,
            panel_count=panel_count
        ))
    end
    return sats
end

function ppc_harmonics_model(planet, degree::Int)
    if isfile(PPC_EARTH_HARMONICS_FILE)
        return GravitationalHarmonicsModel(degree, degree, PPC_EARTH_HARMONICS_FILE, planet)
    end
    return InverseSquaredJ2GravityModel()
end

function ppc_build_config(;
    planet,
    spacecraft,
    mission_time_s::Float64,
    orientation_sim::Bool,
    dynamic_effectors::Tuple,
    density_model=NoAtmosphereModel(),
    guidance_effectors::Tuple=(),
    guidance_rates::Vector{Float64}=Float64[],
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    dt_max_orbit::Float64=10.0,
    reltol_orbit::Float64=1e-9,
    abstol_orbit::Float64=1e-9
)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=orientation_sim,
            num_steps_to_save=300
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

@inline function ppc_mission_time(
    profile::String;
    test::Float64=10.0,
    smoke::Float64=120.0,
    full::Float64=1800.0
)::Float64
    profile == "test" && return test
    profile == "smoke" && return smoke
    return full
end

function ppc_single_config(case_name::String, cfg::PPCConfig; seed::Int=cfg.seed, mc_index::Int=1)
    planet = Earth("", PPC_SPICE_PATH)
    mars = Mars("", PPC_SPICE_PATH)
    q0 = normalize(SVector{4, Float64}(0.15, -0.05, 0.2, 0.96))
    w0 = SVector{3, Float64}(0.01, -0.02, 0.015)
    rng = MersenneTwister(seed + 1009 * mc_index)

    if case_name == "single_inverse_square_vacuum"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet)],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(),)
        )
    elseif occursin(r"^gravity_[0-9]+sat_inverse_square_vacuum$", case_name)
        n = parse(Int, match(r"^gravity_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    elseif case_name == "single_harmonics_l20_vacuum"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet)],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=10.0
        )
    elseif occursin(r"^gravity_[0-9]+sat_l20_vacuum$", case_name)
        n = parse(Int, match(r"^gravity_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    elseif case_name == "single_harmonics_l50_vacuum" || case_name == "montecarlo_high_accuracy"
        ra_jitter = case_name == "montecarlo_high_accuracy" ? randn(rng) * 8e3 : 0.0
        rp_jitter = case_name == "montecarlo_high_accuracy" ? randn(rng) * 8e3 : 0.0
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet; ra_alt_m=550e3 + ra_jitter, rp_alt_m=500e3 + rp_jitter)],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=5.0,
            reltol_orbit=1e-10,
            abstol_orbit=1e-10
        )
    elseif case_name == "srp_heavy_high_area"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet; with_panel=true, panel_count=8, root_area=40.0, orientation_state=(q0, w0))],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=true,
            dynamic_effectors=(InverseSquaredGravityModel(), SolarRadiationPressureModel(1.2, 120.0), SolarRadiationPressureModel(1.8, 220.0)),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=2.0
        )
    elseif case_name == "articulated_1sat_fullstack"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet; with_panel=true, panel_count=28, orientation_state=(q0, w0), root_area=10.0)],
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=true,
            dynamic_effectors=(ppc_harmonics_model(planet, 20), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=1.0
        )
    elseif case_name == "multi_16_aero_surrogate_cached"
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 16),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=5.0
        )
    elseif case_name == "multi_64_high_fidelity" || case_name == "multi_128_high_fidelity"
        n = case_name == "multi_128_high_fidelity" ? 128 : 64
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20), SolarRadiationPressureModel(1.2, 12.0), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=10.0
        )
    elseif case_name == "callback_128_aero_thermal"
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 128),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=3.0
        )
    elseif case_name == "montecarlo_multi_sat"
        spacecraft = SpacecraftModel[]
        for i in 1:4
            push!(spacecraft, ppc_spacecraft(
                planet;
                id=i,
                ra_alt_m=540e3 + randn(rng) * 6e3,
                rp_alt_m=500e3 + randn(rng) * 6e3,
                nu_deg=130.0 + 45.0 * i + randn(rng) * 2.0
            ))
        end
        return ppc_build_config(
            planet=planet,
            spacecraft=spacecraft,
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=10.0
        )
    elseif case_name == "montecarlo_mars_aerobraking"
        return ppc_build_config(
            planet=mars,
            spacecraft=[ppc_spacecraft(
                mars;
                ra_alt_m=4500e3 + randn(rng) * 100e3,
                rp_alt_m=max(110e3, 135e3 + randn(rng) * 10e3),
                i_deg=93.0,
                omega_deg=80.0,
                raan_deg=30.0,
                nu_deg=180.0 + randn(rng) * 4.0
            )],
            mission_time_s=ppc_mission_time(cfg.profile; smoke=120.0, full=1800.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(mars),
            dt_max_orbit=1.0
        )
    end

    throw(ArgumentError("Unknown parallelization performance case '$case_name'."))
end

function ppc_case_catalog()::Dict{String, PPCCaseSpec}
    cases = Dict{String, PPCCaseSpec}()
    add!(name, family, description; montecarlo=false, orientation=false) = begin
        cases[name] = PPCCaseSpec(
            name=name,
            family=family,
            description=description,
            builder=ppc_single_config,
            montecarlo=montecarlo,
            orientation=orientation
        )
    end
    add!("single_inverse_square_vacuum", "gravity_only", "1 spacecraft, inverse-square gravity, no atmosphere")
    for n in (4, 16, 64, 256, 1024, 2048)
        add!("gravity_$(n)sat_inverse_square_vacuum", "gravity_only", "$(n) spacecraft, inverse-square gravity, no atmosphere")
        add!("gravity_$(n)sat_l20_vacuum", "gravity_only", "$(n) spacecraft, L20 harmonics, no atmosphere")
    end
    add!("single_harmonics_l20_vacuum", "gravity_only", "1 spacecraft, L20 harmonics, no atmosphere")
    add!("single_harmonics_l50_vacuum", "few_sat_high_fidelity", "1 spacecraft, L50 harmonics, no atmosphere")
    add!("srp_heavy_high_area", "few_sat_high_fidelity", "1 high-area articulated spacecraft with stacked SRP", orientation=true)
    add!("articulated_1sat_fullstack", "few_sat_high_fidelity", "1 articulated spacecraft with harmonics, aero, thermal, and attitude", orientation=true)
    add!("multi_16_aero_surrogate_cached", "many_sat_high_fidelity", "16 spacecraft with aero and analytic density")
    add!("multi_64_high_fidelity", "many_sat_high_fidelity", "64 spacecraft with harmonics, SRP, aero, and analytic density")
    add!("multi_128_high_fidelity", "many_sat_high_fidelity", "128 spacecraft with harmonics, SRP, aero, and analytic density")
    add!("callback_128_aero_thermal", "many_sat_high_fidelity", "128 spacecraft callback-heavy aero and thermal stress")
    add!("montecarlo_high_accuracy", "monte_carlo", "Monte Carlo high-accuracy gravity seeds", montecarlo=true)
    add!("montecarlo_multi_sat", "monte_carlo", "Monte Carlo 4-spacecraft gravity seeds", montecarlo=true)
    add!("montecarlo_mars_aerobraking", "monte_carlo", "Monte Carlo Mars aerobraking seeds", montecarlo=true)
    return cases
end

function ppc_resolve_cases(requested::Vector{String})::Vector{String}
    catalog = ppc_case_catalog()
    names = isempty(requested) ? sort(collect(keys(catalog))) : requested
    unknown = [name for name in names if !haskey(catalog, name)]
    isempty(unknown) || throw(ArgumentError("Unknown case(s): $(join(unknown, ", "))"))
    return names
end
