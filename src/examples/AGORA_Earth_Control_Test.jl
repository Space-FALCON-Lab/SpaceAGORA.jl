include("../simulation_model/SimulationModel.jl")
using .SimulationModel
using SPICE
using StaticArrays
using LinearAlgebra

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
include("../simulation/run_simulation.jl")
include("typed_example_utils.jl")

planet = Earth("", SPICE_PATH)

q0 = normalize(SVector{4, Float64}(0.2, -0.1, 0.15, 0.95))
w0 = SVector{3, Float64}(0.003, -0.004, 0.005)

ra = 56_378.7978559e3
rp = 200_590.0 + planet.Rp_e
a = (ra + rp) / 2.0
e = (ra - rp) / (ra + rp)
ic = InitialCondition(a, e, 89.876, 75.505, 104.115, 175.0, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(3.7, 2.05, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=640.0,
    panel_mass_each=20.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    prop_mass=200.0
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=1200.0,
    initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(),),
    density_model=NoAtmosphereModel(),
    orientation_sim=true,
    keplerian=true,
    EI_km=300.0
)

run_and_report(args)
