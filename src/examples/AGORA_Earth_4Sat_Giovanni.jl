if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using Dates
using SPICE
using StaticArrays

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    include("../simulation/execution/run_simulation.jl")
end
if !isdefined(@__MODULE__, :run_and_report)
    include("typed_example_utils.jl")
end

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
const MISSION_TIME_SEC = 24.0 * 3600.0

@inline function mean_to_true_anomaly_rad(M_rad::Float64, e::Float64; tol::Float64=1e-13, max_iter::Int=30)
    M_norm = mod(M_rad + pi, 2pi) - pi
    E = e < 0.8 ? M_norm : pi
    for _ in 1:max_iter
        f = E - e * sin(E) - M_norm
        fp = 1.0 - e * cos(E)
        Δ = f / fp
        E -= Δ
        if abs(Δ) < tol
            break
        end
    end
    ν = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0))
    return mod(ν, 2pi)
end

@inline mean_to_true_anomaly_deg(M_deg::Float64, e::Float64) = rad2deg(mean_to_true_anomaly_rad(deg2rad(M_deg), e))

@inline function eclipse_factor_cylindrical(r_sc::SVector{3, Float64}, r_sun::SVector{3, Float64}, rp::Float64)::Float64
    if dot(r_sc, r_sun) >= 0.0
        return 1.0
    end
    d = norm(cross(r_sc, r_sun)) / norm(r_sun)
    return d <= rp ? 0.0 : 1.0
end

function moon_pos_analytic_j2000(et::Float64)::SVector{3, Float64}
    d = et / 86_400.0
    Ω = deg2rad(mod(125.1228 - 0.0529538083 * d, 360.0))
    i = deg2rad(5.1454)
    ω = deg2rad(mod(318.0634 + 0.1643573223 * d, 360.0))
    a = 384_400_000.0
    e = 0.0549
    M = deg2rad(mod(115.3654 + 13.0649929509 * d, 360.0))

    E = M
    for _ in 1:20
        Δ = (E - e * sin(E) - M) / (1.0 - e * cos(E))
        E -= Δ
        abs(Δ) < 1e-13 && break
    end

    ν = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0))
    r = a * (1.0 - e * cos(E))
    u = ω + ν

    x_ecl = r * (cos(Ω) * cos(u) - sin(Ω) * sin(u) * cos(i))
    y_ecl = r * (sin(Ω) * cos(u) + cos(Ω) * sin(u) * cos(i))
    z_ecl = r * (sin(u) * sin(i))

    ϵ = deg2rad(23.439291)
    x_eq = x_ecl
    y_eq = y_ecl * cos(ϵ) - z_ecl * sin(ϵ)
    z_eq = y_ecl * sin(ϵ) + z_ecl * cos(ϵ)
    return SVector{3, Float64}(x_eq, y_eq, z_eq)
end

const _moon_fallback_warned = Ref(false)

function body_pos_wrt_earth_pci(body_name::String, et::Float64, planet::AbstractPlanet)::SVector{3, Float64}
    name = uppercase(body_name)
    if name == "MOON"
        try
            r_moon_j2000 = lock(SimulationModel.SPICE_LOCK) do
                SVector{3, Float64}(spkpos("MOON", et, "J2000", "NONE", "EARTH")[1]) * 1e3
            end
            return planet.J2000_to_pci * r_moon_j2000
        catch
            if !_moon_fallback_warned[]
                _moon_fallback_warned[] = true
                @warn "MOON SPICE ephemeris is unavailable; using analytic lunar orbit fallback for third-body acceleration."
            end
            return planet.J2000_to_pci * moon_pos_analytic_j2000(et)
        end
    end

    r_body_j2000 = lock(SimulationModel.SPICE_LOCK) do
        SVector{3, Float64}(spkpos(name, et, "J2000", "NONE", "EARTH")[1]) * 1e3
    end
    return planet.J2000_to_pci * r_body_j2000
end

struct SunMoonThirdBodyModel{P <: AbstractPlanet} <: AbstractForceTorqueModel
    planet::P
    body_names::NTuple{2, String}
    body_mus::NTuple{2, Float64}
end

SunMoonThirdBodyModel(planet::P) where {P <: AbstractPlanet} = SunMoonThirdBodyModel(
    planet,
    ("SUN", "MOON"),
    (1.3271244002331e20, 4.9028005821478e12)
)

function SimulationModel.calcForceTorque(
    model::SunMoonThirdBodyModel,
    x,
    p::ODEParams,
    i::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_sc = SVector{3, Float64}(x.pos)
    mass_sc = x.mass
    et = p.shared_buffers.et_start[] + p.shared_buffers.current_time[]
    a_tb = MVector{3, Float64}(0.0, 0.0, 0.0)

    @inbounds for k in eachindex(model.body_names)
        body_name = model.body_names[k]
        μ_k = model.body_mus[k]
        r_primary_k = body_pos_wrt_earth_pci(body_name, et, model.planet)
        r_sc_k = r_primary_k - pos_sc
        a_tb .+= μ_k * (r_sc_k / norm(r_sc_k)^3 - r_primary_k / norm(r_primary_k)^3)
    end

    return mass_sc * SVector{3, Float64}(a_tb), SVector{3, Float64}(0.0, 0.0, 0.0)
end

struct CannonballSRPModel{P <: AbstractPlanet} <: AbstractForceTorqueModel
    planet::P
    areas_m2::Vector{Float64}
    cr::Float64
    p_ref::Float64
    au_m::Float64
end

function SimulationModel.calcForceTorque(
    model::CannonballSRPModel,
    x,
    p::ODEParams,
    i::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_sc = SVector{3, Float64}(x.pos)
    et = p.shared_buffers.et_start[] + p.shared_buffers.current_time[]
    r_sun = body_pos_wrt_earth_pci("SUN", et, model.planet)
    r_sc_to_sun = r_sun - pos_sc
    d_sc_sun = norm(r_sc_to_sun)
    u_sun_to_sc = -r_sc_to_sun / d_sc_sun

    eclipse = eclipse_factor_cylindrical(pos_sc, r_sun, model.planet.Rp_e)
    p_srp = model.p_ref * (model.au_m / d_sc_sun)^2
    area = model.areas_m2[min(i, length(model.areas_m2))]
    force = p_srp * model.cr * area * eclipse * u_sun_to_sc
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function make_translational_spacecraft(
    id::Int64,
    a_m::Float64,
    e::Float64,
    i_deg::Float64,
    raan_deg::Float64,
    aop_deg::Float64,
    mean_anomaly_deg::Float64;
    mass_kg::Float64=220.0,
    area_m2::Float64=2.0
)::SpacecraftModel
    root = Link{0}(root=true, m=mass_kg, ref_area=area_m2, dims=MVector{3, Float64}(1.0, 1.0, 1.0))
    ν_deg = mean_to_true_anomaly_deg(mean_anomaly_deg, e)
    ic = InitialCondition(a_m, e, i_deg, aop_deg, raan_deg, ν_deg)
    return SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end

planet = Earth("", SPICE_PATH)
isfile(EARTH_HARMONICS_FILE) || error("Earth harmonics file not found: $EARTH_HARMONICS_FILE")

e = 1e-4
raan_deg = 10.0
aop_deg = 14.0

observer_a_m = 6_963.0e3
target_a_m = 6_964.0e3

observer_i_deg = 85.0
target_i_deg = 86.0

observer_mean_anomalies_deg = (288.0, 290.0, 292.0)
target_mean_anomaly_deg = observer_mean_anomalies_deg[1]

spacecraft = SpacecraftModel[
    make_translational_spacecraft(1, observer_a_m, e, observer_i_deg, raan_deg, aop_deg, observer_mean_anomalies_deg[1]),
    make_translational_spacecraft(2, observer_a_m, e, observer_i_deg, raan_deg, aop_deg, observer_mean_anomalies_deg[2]),
    make_translational_spacecraft(3, observer_a_m, e, observer_i_deg, raan_deg, aop_deg, observer_mean_anomalies_deg[3]),
    make_translational_spacecraft(4, target_a_m, e, target_i_deg, raan_deg, aop_deg, target_mean_anomaly_deg)
]

areas = [sc.root.ref_area for sc in spacecraft]
dynamic_effectors = (
    InverseSquaredGravityModel(),
    GravitationalHarmonicsModel(4, 0, EARTH_HARMONICS_FILE, planet), # up to J4 (zonal)
    SunMoonThirdBodyModel(planet),
    CannonballSRPModel(planet, areas, 1.3, 4.56e-6, 149_597_870_700.0)
)

now_utc = Dates.now(Dates.UTC)
args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true,
        verbose=true,
        generate_plots=false,
        results_directory=joinpath(REPO_ROOT, "output"),
        normalize=false
    ),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=MISSION_TIME_SEC,
        orientation_sim=false,
        num_steps_to_save=2_000
    ),
    environment_model=EnvironmentModel(
        planet=planet,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    ),
    dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=InitialTime(
        year=Dates.year(now_utc),
        month=Dates.month(now_utc),
        day=Dates.day(now_utc),
        hour=Dates.hour(now_utc),
        minute=Dates.minute(now_utc),
        second=Float32(Dates.second(now_utc))
    ),
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-9,
        abstol_orbit=1e-9,
        dt_max_orbit=20.0
    )
)

println("Mission start UTC: $(now_utc)")
println("Propagation time: $(MISSION_TIME_SEC / 3600.0) hours")
run_and_report(args)
