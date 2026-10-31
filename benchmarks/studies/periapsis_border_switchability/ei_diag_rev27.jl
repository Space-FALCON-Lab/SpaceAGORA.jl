# Revision 27 diagnostic (run AFTER rev27ei; point queries only, no flight): along the flown passes of ONE flight, compare the GRAM mean density with the planet-fit density (the model the simulation uses above the entry interface,
# EnvironmentModels.density_polyfit) at altitudes between EI_LOW (150 km) and 2000 km, and the drag delta-v of each pass above 150 km against the total of the pass.
# Orbit per pass: the osculating elements at the flown periapsis (passes.csv: a_m, e, i_rad, raan_rad, omega_rad, time_s) advanced as a Kepler ellipse (no J2 drift inside one orbit: an approximation, stated). Ground track from the same frame chain as the flight
# (elements in the body-equator inertial frame, rotated by mars_rot, then planet_frame_lpi at the point's time). Drag deceleration 0.5 rho v_rel^2 / beta (ballistic coefficient of the flight setup), integrated in time (trapezoid).
# Usage: systemd-run --user --scope -p MemoryMax=16G julia --project=. --threads=1 ei_diag_rev27.jl --flight=DIR --gid=ODY --tag=ODY_s13 --out=FILE.csv [--every=10] [--n=240] [--ei-km=150]
# Output columns per pass: pass, t_s, a_m, e, dv_above_gram_mps, dv_above_fit_mps, dv_below_gram_mps, frac_above_gram, frac_above_fit, med_ln_fit_over_gram_150_500, _500_1000, _1000_2000 (NaN when the band is empty). One Julia process, serial.
using TOML, CSV, DataFrames, Statistics, Printf
include(joinpath(@__DIR__, "planner_common.jl"))
const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const CLI = Dict(m.captures[1] => something(m.captures[2], "1") for a in ARGS for m in [match(r"^--([^=]+)(?:=(.*))?$", a)] if m !== nothing)
const FLIGHT = abspath(get(CLI, "flight", "")); const GID = get(CLI, "gid", "ODY"); const TAG = get(CLI, "tag", "ODY_s13"); const OUTF = abspath(get(CLI, "out", "ei_diag_rev27.csv"))
const EVERY = parse(Int, get(CLI, "every", "10")); const NPT = parse(Int, get(CLI, "n", "240")); const EI_LOW = 1000.0 * parse(Float64, get(CLI, "ei-km", "150"))
isfile(joinpath(FLIGHT, "passes.csv")) || error("--flight=DIR needs a passes.csv")
cfgp = joinpath(@__DIR__, "planner_config.toml"); meta = planner_meta(cfgp); meta, sc, geom, ep = stageg3_load(cfgp, GID); ms = stagef_setup(sc); epoch = parse_epoch(ep); planet = ms.planet; rot = mars_rot(planet, epoch)
et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); mu = Float64(ms.harm.gm_m3s2); polyfit = EM.PolynomialFitAtmosphereModel(planet); beta = Float64(ms.beta)
_corridor_gram_load!(); _corridor_gram_require_direct_sampling()
model = Base.invokelatest(SimulationModel.GRAMAtmosphereModel; planet_name=String(sc["atmosphere_gram"]["planet_name"]), initial_time=epoch, seed=stagef_seed("mars" * TAG))
P = CSV.read(joinpath(FLIGHT, "passes.csv"), DataFrame)

"""Kepler: eccentric anomaly of mean anomaly M (Newton)."""
function kepler_E(M, e); E = M + e * sin(M); for _ in 1:50; d = (E - e * sin(E) - M) / (1 - e * cos(E)); E -= d; abs(d) < 1e-13 && break; end; return E; end

function pass_profile(row)
    a, e, inc, raan, omega, tp = Float64(row.a_m), Float64(row.e), Float64(row.i_rad), Float64(row.raan_rad), Float64(row.omega_rad), Float64(row.time_s); n = sqrt(mu / a^3); p = a * (1 - e^2)
    r150 = planet.Rp_e + EI_LOW; cnu = (p / r150 - 1) / e; cnu > 1 && return nothing                              # the orbit never rises to EI_LOW (cannot happen for these flights)
    nu150 = acos(clamp(cnu, -1.0, 1.0)); E150 = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu150 / 2)); M150 = E150 - e * sin(E150)
    out = NamedTuple[]
    for M in vcat(collect(range(-M150, M150; length=21)), collect(range(M150, 2π - M150; length=NPT)))               # below the interface (dense) and above it
        E = kepler_E(M, e); nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)); t = tp + M / n
        rb, vb = elements_to_cartesian((a=a, e=e, i=inc, raan=raan, omega=omega, nu=nu); mu=mu)
        r = rot * SVector{3, Float64}(rb); v = rot * SVector{3, Float64}(vb); lpi = SimulationModel.planet_frame_lpi(planet, et0 + t, SpiceEphemeridesModel())
        rp = lpi * r; alt, lat = geodetic_altitude_latitude(rp, planet); lon = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rp), planet)[3]
        vrel = norm(lpi * v - cross(planet.ω, rp)); rg = Float64(EM._gram_point_density(model, alt, lat, lon, t, false)[1]); rf = Float64(EM._polyfit_density(polyfit, alt))
        push!(out, (; t, alt, vrel, rg, rf))
    end
    return out
end
trapz(ts, ys) = sum(0.5 * (ys[k] + ys[k + 1]) * (ts[k + 1] - ts[k]) for k in 1:(length(ts) - 1))
rows = NamedTuple[]
for k in 1:EVERY:size(P, 1)
    row = P[k, :]; isfinite(row.a_m) && isfinite(row.e) && isfinite(row.raan_rad) || continue; pr = pass_profile(row); pr === nothing && continue
    ts = [x.t for x in pr]; ag = [0.5 * x.rg * x.vrel^2 / beta for x in pr]; af = [0.5 * x.rf * x.vrel^2 / beta for x in pr]; hi = [x.alt >= EI_LOW for x in pr]
    seg(ys, sel) = (s = 0.0; for j in 1:(length(ts) - 1); (sel[j] && sel[j + 1]) && (s += 0.5 * (ys[j] + ys[j + 1]) * (ts[j + 1] - ts[j])); end; s)
    dv_ag = seg(ag, hi); dv_af = seg(af, hi); dv_bg = seg(ag, .!hi); bandmed(lo, hi_) = (v = [log(x.rf / x.rg) for x in pr if lo <= x.alt < hi_ && x.rg > 0 && x.rf > 0]; isempty(v) ? NaN : median(v))
    push!(rows, (; pass=Int(row.pass), t_s=Float64(row.time_s), a_m=Float64(row.a_m), e=Float64(row.e), dv_above_gram_mps=dv_ag, dv_above_fit_mps=dv_af, dv_below_gram_mps=dv_bg, frac_above_gram=dv_ag / (dv_ag + dv_bg), frac_above_fit=dv_af / (dv_af + dv_bg),
        med_ln_fit_over_gram_150_500=bandmed(EI_LOW, 5.0e5), med_ln_fit_over_gram_500_1000=bandmed(5.0e5, 1.0e6), med_ln_fit_over_gram_1000_2000=bandmed(1.0e6, 2.0e6)))
    println(@sprintf("pass %d: above-%.0f-km drag share GRAM %.3e FIT %.3e", Int(row.pass), EI_LOW / 1e3, rows[end].frac_above_gram, rows[end].frac_above_fit))
end
write_rows_m(OUTF, rows); println("wrote ", OUTF, " (", length(rows), " passes)")
