if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using SPICE
using StaticArrays

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
    p,
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
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
    out = Vector{NamedTuple{(:a, :e, :i, :OMEGA, :omega, :nu), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}}(undef, n)
    @inbounds for i in 1:n
        pos = SVector{3,Float64}(u.sc[i].pos)
        vel = SVector{3,Float64}(u.sc[i].vel)
        out_temp = rvtoorbitalelement(pos, vel, planet)
        out[i] = (a=out_temp[1], e=out_temp[2], i=out_temp[3], OMEGA=out_temp[4], omega=out_temp[5], nu=out_temp[6])
    end
    out
end

function density_getter(u, t, integrator)
    return Vector{Float64}(integrator.p.shared_buffers.densities)
end

planet = Earth("", SPICE_PATH)

ic = InitialCondition(
    a=7000.0e3,
    e=0.00001,
    i=28.5,
    ω=0.0,
    Ω=45.0,
    ν=0.0
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 3.7, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    prop_mass=10.0,
    id=102,
)

# density_model = GRAMAtmosphereModel(planet_name="earth")
density_model = ConstantDensityModel(1e-20, 300.0) # for testing

dynamic_effectors = (
    InverseSquaredJ2GravityModel(),
    # GravitationalHarmonicsModel(2, 0, "Gravity_harmonics_data/EarthGGM05C_full.csv", planet),
    # AerodynamicCoefficientfM(),
    # NBodyGravityModel(body_ids=(10,), primary_body_id=planet.spice_id),
)
# guidance = GuidanceModel(guidance_effectors=(thruster_guidance,), guidance_rates=[30.0])
base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=1.0e6, # 1e6 seconds ~ 11.6 days
    data_rate=10.0,
    # orbits=50,
    initial_time=InitialTime(year=2026, month=1, day=1, hour=12, minute=0, second=0.0),
    dynamic_effectors=dynamic_effectors,
    density_model=density_model,
    orientation_sim=false,
    keplerian=false,
    EI_km=100.0
)
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=base_args.mission_configuration,
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=base_args.guidance_model,
    navigation_model=base_args.navigation_model,
    control_model=base_args.control_model,
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-12,
        abstol_orbit=1e-12,
        dt_max_orbit=10.0,
        reltol_atmosphere=1e-12,
        abstol_atmosphere=1e-12,
        dt_max_atmosphere=1.0
    )
)

# args = make_example_config(
#     planet=planet,
#     spacecraft=spacecraft,
#     orbits=50,
#     initial_time=InitialTime(
#         year=2014,
#         month=5,
#         day=19,
#         hour=14,
#         minute=7,
#         second=32.0,
#     ),
#     dynamic_effectors=(
#         InverseSquaredJ2GravityModel(),
#         AerodynamicCoefficientfM(),
#         NBodyGravityModel(body_names=("Sun",), primary_body_name=planet.name),
#     ),
#     density_model=GRAMAtmosphereModel(planet_name="venus"),
#     orientation_sim=false,
#     keplerian=false,
#     EI_km=250.0,
# )

save_fields = vcat(
    default_save_fields(args),   # keep existing fields
    [
        # SaveField(:geo, geo_getter; per_satellite=true, column_prefix="geo"),
        SaveField(:orbital_elements, oe_getter; per_satellite=true, column_prefix="oe"),
        # SaveField(:density, density_getter; per_satellite=true, column_prefix="density")
    ]
)

run_and_report(args; save_fields=save_fields)
