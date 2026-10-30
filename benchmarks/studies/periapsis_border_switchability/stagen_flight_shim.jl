# Compatibility shim for the current main, loaded BEFORE the Stage F/H common files (which are frozen and not edited). Main renamed the guidance tick calcGuidanceEffect! to update_guidance! and removed the old
# name; the Stage F/H policies and the example's corridor policy still define methods under the old name. This file (1) loads SpaceAGORA through the examples environment, (2) recreates the old function in
# GuidanceHooks so those definitions load, and (3) bridges update_guidance! to it for any guidance model without a method of its own. Same arguments, same call times (the rename commit says only the names change).
# Verified by the bitwise reproduction of Stage I flights (stagen_flight_plan.md section 6, repro check).
include(joinpath(@__DIR__, "..", "..", "..", "examples", "common.jl"))
if !isdefined(SimulationModel.GuidanceHooks, :calcGuidanceEffect!)
    Core.eval(SimulationModel.GuidanceHooks, :(function calcGuidanceEffect! end))
    Core.eval(SimulationModel.GuidanceHooks, :(update_guidance!(m::$(SimulationModel.AbstractGuidanceModel), u, p, t::Float64, i::Int64) = calcGuidanceEffect!(m, u, p, t, i)))
end
# The same release (0.5.0) removed three ephemerides hook names used by the frozen study files; forward them to the renamed hooks (model first; same values).
if !isdefined(SimulationModel, :ephemerides_time_seconds)
    Core.eval(SimulationModel, :(ephemerides_time_seconds(initial_time, model) = ephemerides_epoch_seconds(model, initial_time)))
    Core.eval(SimulationModel, :(planet_frame_lpi(planet, et, model) = ephemerides_frame_lpi(model, planet, et)))
    Core.eval(SimulationModel, :(ephemerides_sun_direction_ii(planet, et, model) = ephemerides_sun_ii(model, planet, et)))
end
