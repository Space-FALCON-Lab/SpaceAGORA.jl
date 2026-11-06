## Mounted actuators on a CubeSat: the wheels detumble, the torque rods unload them.
#
# The vehicle carries three reaction wheels, three torque rods and an optical
# crosslink terminal, each mounted on its bus with `add_actuator!`. One
# controller, in `:mounted` actuation, sets a level for every mounted actuator
# at each control tick, and the engine turns those levels into the wrench
# through each actuator's own record: each wheel's torque within its torque and
# momentum limits, each rod's dipole crossed with the run's magnetic field, the
# terminal's recoil. The run's results carry the level each actuator applied
# (`sc1_actuator_level_k`, in mounting order) and the power the actuators draw
# (`sc1_actuator_power`).
#
# Checks (also run by the CI examples smoke): the wheels damp the initial
# tumble, no wheel passes its momentum limit, no rod passes its dipole limit,
# and the power column is the terminal's draw. Over the full run the rods also
# take momentum out of the wheels.
include(joinpath(@__DIR__, "common.jl"))
using StaticArrays
using LinearAlgebra

const SMOKE = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
const MISSION_TIME_S = SMOKE ? 120.0 : 3000.0

include(joinpath(@__DIR__, "support", "mounted_actuator_case.jl"))
args = mounted_actuator_configuration(duration_s=MISSION_TIME_S)
result = run_simulation(args; return_results=true, return_solution=true)
states = [u.sc[1] for u in result.solution.u]
table = result.table

ω0, ω_end = norm(states[1].ω), norm(states[end].ω)
h_peak = maximum(s -> maximum(abs, s.h_wheels), states)
h_after_detumble = maximum(abs, states[findfirst(t -> t >= 100.0, result.solution.t)].h_wheels)
h_end = maximum(abs, states[end].h_wheels)
rod_peak = maximum(k -> maximum(abs, table[!, "sc1_actuator_level_$k"]), 4:6)
println("body rate:       ", round(ω0; sigdigits=3), " -> ", round(ω_end; sigdigits=3), " rad/s")
println("wheel momentum:  peak ", round(h_peak; sigdigits=3), " N m s (limit ", WHEEL_MOMENTUM_NMS,
    "), at 100 s ", round(h_after_detumble; sigdigits=3), ", at the end ", round(h_end; sigdigits=3))
println("rod dipole:      peak ", round(rod_peak; sigdigits=3), " A m^2 (limit ", ROD_DIPOLE_AM2, ")")
println("actuator power:  ", table.sc1_actuator_power[end], " W")

@assert ω_end < 0.05 * ω0 "the wheels did not damp the tumble"
@assert h_peak <= WHEEL_MOMENTUM_NMS * (1 + 1e-9) "a wheel passed its momentum limit"
@assert rod_peak <= ROD_DIPOLE_AM2 * (1 + 1e-12) "a rod passed its dipole limit"
@assert all(==(LASER_POWER_W / LASER_EFFICIENCY), table.sc1_actuator_power[2:end]) "the power column is not the terminal's draw"
SMOKE || @assert h_end < 0.5 * h_after_detumble "the rods did not unload the wheels"

println("Earth_CubeSat_Mounted_Actuators ok")
