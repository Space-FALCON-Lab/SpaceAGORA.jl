if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using SPICE
using StaticArrays
using LinearAlgebra

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    include("../simulation/run_simulation.jl")
end
if !isdefined(@__MODULE__, :make_example_config)
    include("typed_example_utils.jl")
end

struct ConstantDensityModel <: AbstractDensityModel
    rho::Float64
    temp::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

struct TimedVelocityThrusterModel
    thrust::Float64
    start_time::Float64
    stop_time::Float64
end

function SimulationModel.calcControlForceTorque(
    model::TimedVelocityThrusterModel,
    u,
    p::ODEParams,
    i::Int64,
    t::Float64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if t < model.start_time || t > model.stop_time
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    vel = SVector{3, Float64}(u.vel)
    speed = norm(vel)
    if speed <= 1e-9
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    force = model.thrust * vel / speed
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlEffect!(
    model::TimedVelocityThrusterModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

function geo_getter(u, t, integrator)
    n = length(u.sc)
    planet = integrator.p.args.environment_model.planet
    out = Vector{NamedTuple{(:alt_m, :lat_rad, :lon_rad), Tuple{Float64, Float64, Float64}}}(undef, n)
    @inbounds for i in 1:n
        pos = SVector{3,Float64}(u.sc[i].pos)
        vel = SVector{3,Float64}(u.sc[i].vel)
        rp, _ = r_intor_p!(pos, vel, planet)
        alt, lat, lon = rtolatlong(rp, planet)
        out[i] = (alt_m=alt, lat_rad=lat, lon_rad=lon)
    end
    out
end

function oe_getter(u, t, integrator)
    n = length(u.sc)
    planet = integrator.p.args.environment_model.planet
    out = Vector{NamedTuple{(:a, :e, :i, :omega, :OMEGA, :nu), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}}(undef, n)
    @inbounds for i in 1:n
        pos = SVector{3,Float64}(u.sc[i].pos)
        vel = SVector{3,Float64}(u.sc[i].vel)
        out_temp = rvtoorbitalelement(pos, vel, planet)
        out[i] = (a=out_temp[1], e=out_temp[2], i=out_temp[3], omega=out_temp[4], OMEGA=out_temp[5], nu=out_temp[6])
    end
    out
end

function density_getter(u, t, integrator)
    return Vector{Float64}(integrator.p.shared_buffers.densities)
end


planet = Mars("", SPICE_PATH)
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"

ic = InitialCondition(
    ra=28_559.615e3,
    rp=planet.Rp_e + 95_000.0,
    i=93.522,
    ω=109.7454,
    Ω=28.1517,
    ν=180.1
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7),
    panel_dims=(0.01, 3.89 / 2.0, 1.7),
    bus_mass=391.0,
    panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 3.89 / 4.0,
    ic=ic,
    prop_mass=60.0,
    id=100
)

mars_harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "GMM3.csv")
dynamic_effectors = smoke_mode ? (InverseSquaredGravityModel(),) : (
    InverseSquaredGravityModel(),
    GravitationalHarmonicsModel(20, 20, mars_harmonics_file, planet),
    AerodynamicCoefficientfM(),
    NBodyGravityModel(body_names=("Sun",), primary_body_name=planet.name),
)
# density_model = smoke_mode ? NoAtmosphereModel() : ConstantDensityModel(2e-8, 180.0)
# density_model = NoAtmosphereModel()
density_model = GRAMAtmosphereModel(planet_name="mars")
base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    # mission_time=3_600.0*18.0*50.0, # 50 orbits
    orbits=100,
    initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0),
    dynamic_effectors=dynamic_effectors,
    density_model=density_model,
    orientation_sim=false,
    keplerian=smoke_mode,
    EI_km=160.0
)

# thruster = TimedVelocityThrusterModel(4.0, 900.0, 1_200.0)
thruster = BaseThrusterModel(thrust=[4.0], direction=[deg2rad(180.0)], Δv=[0.0], Isp=[300.0], start_burn_time=[0.0], stop_burn_time=[-0.1])
thruster_guidance = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
    maneuver_orbit_number=[7, 14, 26, 30, 35, 47, 54, 69, 72, 80, 87, 110, 128, 161, 179, 195, 211, 223, 239, 251, 263, 274, 287, 299, 311],
    maneuver_Δv=[0.15, -0.15, -0.1, -0.1, -0.2, -0.2, 0.3, -0.15, -0.15, -0.15, 0.15, 0.15, 1.0, 0.84, 0.6, 0.84, 0.6, 0.6, 1.2, 1.0, 1.2, 1.2, 1.0, 1.0, 1.2]
)
guidance = GuidanceModel(guidance_effectors=(thruster_guidance,), guidance_rates=[30.0])
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=base_args.mission_configuration,
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=guidance,
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(thruster,), control_rates=[30.0]),
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=30.0,
        reltol_atmosphere=1e-8,
        abstol_atmosphere=1e-8,
        dt_max_atmosphere=0.5
    )
)
save_fields = vcat(
    default_save_fields(args),   # keep existing fields
    [
        SaveField(:geo, geo_getter; per_satellite=true, column_prefix="geo"),
        SaveField(:orbital_elements, oe_getter; per_satellite=true, column_prefix="oe"),
        SaveField(:density, density_getter; per_satellite=true, column_prefix="density")
    ]
)

run_and_report(args; save_fields=save_fields)