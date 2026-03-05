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
    include("../simulation/execution/run_simulation.jl")
end
if !isdefined(@__MODULE__, :make_example_config)
    include("typed_example_utils.jl")
end

struct GravityGradientTorqueModel <: AbstractForceTorqueModel
    inertia::SMatrix{3, 3, Float64}
end

function SimulationModel.calcForceTorque(
    model::GravityGradientTorqueModel,
    x,
    p::ODEParams,
    i::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    r = SVector{3, Float64}(x.pos)
    rmag = norm(r)
    if rmag <= 1e-9
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    mu = p.args.environment_model.planet.μ
    rhat = r / rmag
    torque = 3.0 * mu / rmag^3 * cross(rhat, model.inertia * rhat)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque
end

planet = Earth("", SPICE_PATH)

q0 = normalize(SVector{4, Float64}(0.1, 0.2, -0.3, 0.9))
w0 = SVector{3, Float64}(0.001, -0.01, 0.03)

ra = planet.Rp_e + 15_000e3
rp = planet.Rp_e + 145e3
a = (ra + rp) / 2.0
e = (ra - rp) / (ra + rp)
ic = InitialCondition(a, e, 33.3, 347.8, 48.2, 85.3, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.0, 2.36643, 2.9664),
    panel_dims=(0.01, 1.0, 0.6),
    bus_mass=749.0,
    panel_mass_each=1.0,
    panel_offset_y=0.8,
    ic=ic,
    prop_mass=1.0,
    id=103
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=2_500.0,
    initial_time=InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(), GravityGradientTorqueModel(spacecraft.root.inertia)),
    density_model=NoAtmosphereModel(),
    orientation_sim=true,
    keplerian=true,
    EI_km=160.0
)

run_and_report(args)
