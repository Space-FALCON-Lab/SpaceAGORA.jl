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

planet = Earth("", SPICE_PATH)

q0 = normalize(SVector{4, Float64}(0.12, -0.22, 0.05, 0.96))
w0 = SVector{3, Float64}(0.01, -0.015, 0.02)

ra = 7_200e3
rp = 6_900e3
a = (ra + rp) / 2.0
e = (ra - rp) / (ra + rp)
ic = InitialCondition(a, e, 35.0, 40.0, 10.0, 170.0, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.6, 1.8, 1.5),
    panel_dims=(0.01, 1.5, 0.8),
    bus_mass=750.0,
    panel_mass_each=5.0,
    panel_offset_y=1.2,
    ic=ic,
    prop_mass=1.0
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=900.0,
    initial_time=InitialTime(year=2025, month=1, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    density_model=NoAtmosphereModel(),
    orientation_sim=true,
    keplerian=true,
    EI_km=120.0
)

run_and_report(args)
