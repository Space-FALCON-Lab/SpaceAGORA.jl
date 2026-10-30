# Sun position and Mars pole from the ENGINE's own ephemerides, tabulated once (no GRAM, no flight, one short process): used by geometry_rank.py for the Sun third-body channel.
#   julia --project=<repo> sun_table.jl --out=<dir>   ->  <dir>/sun_table.csv (daily) and <dir>/sun_check.csv (off-grid epochs), each with a sha256 in <dir>/manifest.sha256
# Sun: SimulationModel's own SPICE call (spice_position_j2000_m("sun", et, "mars"), the call behind ephemerides_sun_ii and the third-body model), kernels furnished by Mars("", SPICE_PATH).
# Pole: third row of planet_frame_lpi (the engine's J2000 -> Mars-fixed transform, which mars_rot uses for the study's element frame) = the body spin axis in J2000.
# Time: UTC days since 2000-01-01T12:00:00 (the geometry_rank.py convention); et from the engine's ephemerides_time_seconds.
using SHA, Dates, LinearAlgebra
include(joinpath(@__DIR__, "stagen_flight_shim.jl"))   # current-main compatibility; must precede common.jl
include(joinpath(@__DIR__, "common.jl"))
const OUT = abspath(let m = match(r"--out=(.*)", join(ARGS, " ")); m === nothing ? error("--out=<dir>") : m.captures[1]; end)
mkpath(OUT)
planet = Mars("", SPICE_PATH); eph = SpiceEphemeridesModel(); J2000_UTC = DateTime(2000, 1, 1, 12, 0, 0)
function row(d::DateTime)
    ep = InitialTime(year=Dates.year(d), month=Dates.month(d), day=Dates.day(d), hour=Dates.hour(d), minute=Dates.minute(d), second=Dates.second(d))
    et = SimulationModel.ephemerides_time_seconds(ep, eph)
    s = SimulationModel.EphemeridesModels.spice_position_j2000_m("sun", et, "mars"); l = SimulationModel.planet_frame_lpi(planet, et, eph)
    return (Dates.value(d - J2000_UTC) / 86_400_000, et, s[1], s[2], s[3], l[3, 1], l[3, 2], l[3, 3])
end
hdr = "days_j2000_utc,et_s,sun_x_m,sun_y_m,sun_z_m,pole_x,pole_y,pole_z"
write_tab(path, ds) = open(path, "w") do io; println(io, hdr); for d in ds; println(io, join(row(d), ",")); end; end
write_tab(joinpath(OUT, "sun_table.csv"), DateTime(2001, 10, 1):Day(1):DateTime(2003, 1, 1))
# off-grid check epochs (whole seconds; deterministic): 300 epochs spread over the same range
n0 = Dates.value(DateTime(2003, 1, 1) - DateTime(2001, 10, 1)) ÷ 1000
write_tab(joinpath(OUT, "sun_check.csv"), [DateTime(2001, 10, 1) + Second(mod(round(Int, 7919 * k * 3607.31), n0 - 86400)) for k in 1:300])
open(joinpath(OUT, "manifest.sha256"), "w") do io
    for f in ("sun_table.csv", "sun_check.csv"); println(io, bytes2hex(sha256(read(joinpath(OUT, f)))), "  ", f); end
end
println("wrote ", OUT)
