# Venus analogue of sun_table.jl: Sun position seen from Venus (J2000), Venus pole and prime-meridian axis from the ENGINE's own ephemerides, tabulated once
# (no GRAM, no flight, one short process).   julia --project=<repo> sun_table_venus.jl --out=<dir>
#   -> <dir>/sun_table.csv (daily, 2014-05-01 to 2015-07-01) and <dir>/sun_check.csv (300 off-grid epochs), sha256 in <dir>/manifest.sha256
# Sun: SimulationModel's own SPICE call (spice_position_j2000_m("sun", et, "venus")); kernels furnished by Venus("", SPICE_PATH).
# Pole and prime meridian: rows 3 and 1 of planet_frame_lpi (the engine's J2000 -> Venus-fixed transform); geometry_rank_venus.py derives the rotation angle from them.
# Time: UTC days since 2000-01-01T12:00:00 (the geometry_rank.py convention); et from the engine's ephemerides_time_seconds.
using SHA, Dates, LinearAlgebra
include(joinpath(@__DIR__, "stagen_flight_shim.jl"))
include(joinpath(@__DIR__, "common.jl"))
const OUT = abspath(let m = match(r"--out=(.*)", join(ARGS, " ")); m === nothing ? error("--out=<dir>") : m.captures[1]; end)
mkpath(OUT)
planet = Venus("", SPICE_PATH); eph = SpiceEphemeridesModel(); J2000_UTC = DateTime(2000, 1, 1, 12, 0, 0)
function row(d::DateTime)
    ep = InitialTime(year=Dates.year(d), month=Dates.month(d), day=Dates.day(d), hour=Dates.hour(d), minute=Dates.minute(d), second=Dates.second(d))
    et = SimulationModel.ephemerides_time_seconds(ep, eph)
    s = SimulationModel.EphemeridesModels.spice_position_j2000_m("sun", et, "venus"); l = SimulationModel.planet_frame_lpi(planet, et, eph)
    return (Dates.value(d - J2000_UTC) / 86_400_000, et, s[1], s[2], s[3], l[3, 1], l[3, 2], l[3, 3], l[1, 1], l[1, 2], l[1, 3])
end
hdr = "days_j2000_utc,et_s,sun_x_m,sun_y_m,sun_z_m,pole_x,pole_y,pole_z,pm_x,pm_y,pm_z"
write_tab(path, ds) = open(path, "w") do io; println(io, hdr); for d in ds; println(io, join(row(d), ",")); end; end
d0, d1 = DateTime(2014, 5, 1), DateTime(2015, 7, 1)
write_tab(joinpath(OUT, "sun_table.csv"), d0:Day(1):d1)
n0 = Dates.value(d1 - d0) ÷ 1000
write_tab(joinpath(OUT, "sun_check.csv"), [d0 + Second(mod(round(Int, 7919 * k * 3607.31), n0 - 86400)) for k in 1:300])
open(joinpath(OUT, "manifest.sha256"), "w") do io
    for f in ("sun_table.csv", "sun_check.csv"); println(io, bytes2hex(sha256(read(joinpath(OUT, f)))), "  ", f); end
end
println("wrote ", OUT)
