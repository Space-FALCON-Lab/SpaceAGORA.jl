include(joinpath(@__DIR__, "common.jl"))
using SPICE
using StaticArrays
using LinearAlgebra


planet = Earth("", SPICE_PATH)

q0 = normalize(SVector{4, Float64}(0.0, 0.0, 0.3, 0.9539392))
w0 = SVector{3, Float64}(0.0, 0.0, 0.0)

ra = planet.Rp_e + 36_000e3
rp = planet.Rp_e + 600e3
a = (ra + rp) / 2.0
e = (ra - rp) / (ra + rp)
ic = InitialCondition(a, e, 45.0, 20.0, 10.0, 170.0, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(1.8, 2.2, 1.6),
    panel_dims=(0.01, 2.4, 0.8),
    bus_mass=420.0,
    panel_mass_each=8.0,
    panel_offset_y=1.5,
    ic=ic,
    prop_mass=0.0,
    id=104
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=3_000.0,
    initial_time=InitialTime(year=2022, month=1, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(
        InverseSquaredGravityModel(),
        SolarRadiationPressureModel(1.2, spacecraft.root.ref_area)
    ),
    density_model=NoAtmosphereModel(),
    orientation_sim=true,
    keplerian=true,
    EI_km=180.0
)

run_and_report(args)
