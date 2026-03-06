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

struct QuaternionPDControlModel
    kq::Float64
    kw::Float64
end

function SimulationModel.calcControlForceTorque(
    model::QuaternionPDControlModel,
    u,
    p::ODEParams,
    i::Int64,
    t::Float64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    q = SVector{4, Float64}(u.q)
    w = SVector{3, Float64}(u.ω)

    q_vec = q[1:3]
    if q[4] < 0.0
        q_vec = -q_vec
    end

    torque = -model.kq * q_vec - model.kw * w
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque
end

function SimulationModel.calcControlEffect!(
    model::QuaternionPDControlModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

planet = Earth("", SPICE_PATH)

q0 = normalize(SVector{4, Float64}(0.25, -0.2, 0.15, 0.92))
w0 = SVector{3, Float64}(0.015, -0.02, 0.01)

ra = planet.Rp_e + 1_000e3
rp = planet.Rp_e + 700e3
a = (ra + rp) / 2.0
e = (ra - rp) / (ra + rp)
ic = InitialCondition(a, e, 51.6, 0.0, 30.0, 170.0, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(1.8, 2.0, 1.5),
    panel_dims=(0.01, 1.2, 0.6),
    bus_mass=320.0,
    panel_mass_each=4.0,
    panel_offset_y=1.0,
    ic=ic,
    prop_mass=0.0,
    id=105
)

base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=2_000.0,
    initial_time=InitialTime(year=2023, month=4, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    density_model=NoAtmosphereModel(),
    orientation_sim=true,
    keplerian=true,
    EI_km=140.0
)

control_model = ControlModel(control_effectors=(QuaternionPDControlModel(0.15, 0.6),), control_rates=[0.5])
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
        dt_max_orbit=3.0
    )
)

run_and_report(args)
