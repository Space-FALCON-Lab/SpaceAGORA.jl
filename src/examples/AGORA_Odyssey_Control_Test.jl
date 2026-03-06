if !isdefined(@__MODULE__, :SimulationModel)
    include("../core/simulation_model.jl")
end
using .SimulationModel
using SPICE
using StaticArrays
using LinearAlgebra

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include("../simulation/engine/simulation_engine.jl")
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end
if !isdefined(@__MODULE__, :make_example_config)
    include(joinpath(@__DIR__, "..", "core", "utils", "typed_example_utils.jl"))
end

struct BodyRateDamperModel
    gain::Float64
end

function SimulationModel.calcControlForceTorque(
    model::BodyRateDamperModel,
    u,
    p::ODEParams,
    i::Int64,
    t::Float64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    torque = -model.gain * SVector{3, Float64}(u.ω)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque
end

function SimulationModel.calcControlEffect!(
    model::BodyRateDamperModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

planet = Mars("", SPICE_PATH)

q0 = normalize(SVector{4, Float64}(0.0, -0.6321683, -0.07370895, 0.7713171))
w0 = SVector{3, Float64}(0.02, -0.01, 0.015)

ra = 28_559.615e3
rp = planet.Rp_e + 93_000.0
a = (ra + rp) / 2.0
e = (ra - rp) / (ra + rp)
ic = InitialCondition(a, e, 93.522, 109.7454, 28.1517, 175.0, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7),
    panel_dims=(0.01, 3.89 / 2.0, 1.7),
    bus_mass=391.0,
    panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 3.89 / 4.0,
    ic=ic,
    prop_mass=50.0,
    id=101
)

base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=1_800.0,
    initial_time=InitialTime(year=2001, month=11, day=6, hour=8, minute=30, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    density_model=NoAtmosphereModel(),
    orientation_sim=true,
    keplerian=true,
    EI_km=125.0
)

control_model = ControlModel(control_effectors=(BodyRateDamperModel(0.05),), control_rates=[1.0])
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=base_args.mission_configuration,
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=base_args.guidance_model,
    navigation_model=base_args.navigation_model,
    control_model=control_model,
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-9,
        abstol_orbit=1e-9,
        reltol_quaternion=1e-9,
        abstol_quaternion=1e-9,
        dt_max_orbit=5.0
    )
)

run_and_report(args)
