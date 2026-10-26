## A CubeSat's power budget: solar array, battery, constant load and a laser terminal.
#
# A 10 kg bus in a 500 km orbit holds its attitude, its 0.06 m² array facing the
# Sun at the start. A 5 Wh battery, half full, feeds a 5 W avionics load and a
# laser terminal transmitting 2 W optical at 20 % efficiency (10 W electrical).
# In sunlight the array's 24.5 W charges the battery until it is full; in
# eclipse the battery carries the 15 W until it is empty, the laser is shed,
# and when the Sun returns the laser comes back and the battery recharges.
#
# Checks (also run by the CI examples smoke): the battery stays within its
# capacity, fills in sunlight, empties in eclipse with the laser shed, and
# recovers with the laser back on.
include(joinpath(@__DIR__, "common.jl"))
using StaticArrays
using LinearAlgebra
using Printf

planet = Earth()
ephemerides = SimpleEphemeridesModel()
epoch = InitialTime(year=2024, month=3, day=20, hour=0, minute=0, second=0.0)
sun_at_start = SM.ephemerides_sun_ii(ephemerides, planet, SM.ephemerides_epoch_seconds(ephemerides, epoch))

radius = planet.Rp_e + 500e3
ic = CartesianInitialCondition(SVector(radius, 0.0, 0.0), SVector(0.0, sqrt(planet.μ / radius), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
bus = Link(root=true, m=10.0, dims=MVector(0.1, 0.2, 0.3))
spacecraft = SpacecraftModel(links=[bus], root=bus, inertia_tensor=bus.inertia, initial_condition=ic, id=1)
add_power!(spacecraft, bus, SolarArray(normal=sun_at_start, area_m2=0.06, efficiency=0.3))
add_power!(spacecraft, bus, Battery(capacity_wh=5.0, initial_soc=0.5))
add_power!(spacecraft, bus, PowerLoad(power_w=5.0, name="avionics"))
add_actuator!(spacecraft, bus, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=2.0, efficiency=0.2))

# The controller keeps the laser transmitting 2 W; the power system decides
# whether it may.
mutable struct TransmitContinuously <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::TransmitContinuously) = :mounted
SpaceAGORA.actuator_levels(c::TransmitContinuously, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::TransmitContinuously, u, p, t::Float64, i::Int64) = nothing

base = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=2 * 5677.0,            # two orbits, short enough for the smoke run
    initial_time=epoch,
    dynamic_effectors=(InverseSquaredGravityModel(),),
    ephemerides_model=ephemerides,
    orientation_sim=true,
    keplerian=false,
    EI_km=100.0,
    verbose=false,
    results=false,
)
args = SM.SimConfig._with_configuration(base;
    control_model=ControlModel(control_effectors=(TransmitContinuously([2.0]),), control_rates=[10.0]))
table = run_simulation(args; return_results=true).table

soc, laser = table.sc1_battery_soc, table.sc1_actuator_level_1
for k in 1:max(1, length(table.time) ÷ 12):length(table.time)
    @printf("t %6.0f s  state of charge %5.1f %%  generated %5.1f W  load %5.1f W  shed %4.1f W\n", table.time[k],
        100 * soc[k], table.sc1_power_generated_w[k], table.sc1_power_load_w[k], table.sc1_power_shed_w[k])
end

full = findall(>=(1.0 - 1e-9), soc)
empty = findall(<=(1e-9), soc)
@assert all(s -> -1e-9 <= s <= 1.0 + 1e-9, soc) "the battery left its capacity"
@assert !isempty(full) && all(i -> table.sc1_power_generated_w[i] > 15.0, full) "the battery did not fill in sunlight"
@assert !isempty(empty) && all(i -> laser[i] == 0.0, empty) "the laser was not shed with the battery empty"
@assert any(i -> i > maximum(empty) && soc[i] > 0.1 && laser[i] == 2.0, eachindex(soc)) "the laser did not come back after the eclipse"

println("Earth_CubeSat_Power_Budget ok")
