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

planet = Titan("", SPICE_PATH)

q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
w0 = SVector{3, Float64}(0.005, -0.008, 0.01)

ra = 15_000e3 + planet.Rp_e
rp = 720_000.0 + planet.Rp_e
a = (ra + rp) / 2.0
e = (ra - rp) / (ra + rp)

ic = InitialCondition(a, e, 85.37, 90.0, 64.495, 175.0, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(3.7, 2.05, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=2_400.0,
    panel_mass_each=50.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    prop_mass=1_000.0
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=2.0 * 3600.0,
    initial_time=InitialTime(year=2031, month=10, day=15, hour=14, minute=21, second=0.0),
    dynamic_effectors=(InverseSquaredJ2GravityModel(),),
    density_model=NoAtmosphereModel(),
    orientation_sim=true,
    keplerian=true,
    EI_km=1500.0
)

run_and_report(args)
