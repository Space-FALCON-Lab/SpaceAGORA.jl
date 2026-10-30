# Revision 21 (WP21-C): the Revision 17b dust-level rule at a generic geometry (copy of gram_dust_levels_rev17.jl, parametrized): find the constant MGCM dust optical depth whose mean ln(rho / rho_base) over the first NP periapsis
# points of the nominal flight is +TARGET (L2p) and -TARGET (L3p). Point queries only; the nominal flight is read from --flight=DIR (passes.csv) or, if that directory has none, flown first in this same process (arm OL, tag ODY_s13,
# --fly-days, needs --cal-root with the geometry's calibration). Serial, one process, one thread.
# Usage: systemd-run --user --scope -p MemoryMax=16G julia --project=. --threads=1 gram_dust_levels_rev21.jl --geom=G0975|G0063 --out=DIR [--flight=DIR] [--cal-root=DIR] [--fly-days=35] [--target=0.2231]
using TOML, CSV, DataFrames, Statistics, Printf, JSON
include(joinpath(@__DIR__, "planner_common.jl"))
const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const CLI = Dict(m.captures[1] => something(m.captures[2], "1") for a in ARGS for m in [match(r"^--([^=]+)(?:=(.*))?$", a)] if m !== nothing)
const OUT = abspath(get(CLI, "out", ".")); mkpath(OUT); const TARGET = parse(Float64, get(CLI, "target", "0.2231")); const NP = 40
const GEOM = get(CLI, "geom", "G0975"); const FLIGHT = abspath(get(CLI, "flight", joinpath(OUT, "nominal_flight")))
cfgp = joinpath(@__DIR__, "planner_config.toml"); meta = planner_meta(cfgp); gid = planner_gid(meta, GEOM == "G0975" ? "GC" : GEOM == "G0063" ? "GA" : GEOM)
meta, sc, geom, ep = stageg3_load(cfgp, gid); ms0 = stagef_setup(sc); epoch = parse_epoch(ep)
if !isfile(joinpath(FLIGHT, "passes.csv"))
    haskey(CLI, "cal-root") || error("no passes.csv in $FLIGHT and no --cal-root to fly the nominal flight")
    seed = (; config_path=cfgp, gid, grole=GEOM == "G0975" ? "GC" : "GA", role="flight", disp="", arm="OL", tag="ODY_s13", code=planner_arm_code("OL"), run_id=basename(FLIGHT), out=dirname(FLIGHT), cal_root=abspath(CLI["cal-root"]),
        horizon_days=parse(Float64, get(CLI, "fly-days", "35")), dt_atm_s=NaN, smoke_force=false, decision_speed=String(meta["stageg3"]["decision_speed"]), host="local")
    row = planner_run_task(seed); println("nominal flight: status=$(row.status) passes=$(row.passes) stop=$(row.stop_reason) $(row.error)"); row.status == "ok" || error("nominal flight failed")
end
P = CSV.read(joinpath(FLIGHT, "passes.csv"), DataFrame); size(P, 1) >= NP || error("nominal flight has only $(size(P, 1)) passes, need $NP")
pts = [(t=P.time_s[k], h=P.alt_m[k], lat=P.lat_rad[k], lon=P.lon_rad[k]) for k in 1:NP]
function dens(; kw...)
    _corridor_gram_load!(); _corridor_gram_require_direct_sampling()
    m = Base.invokelatest(SimulationModel.GRAMAtmosphereModel; planet_name=String(sc["atmosphere_gram"]["planet_name"]), initial_time=epoch, seed=stagef_seed("mars" * "ODY_s13"), kw...)
    return [Float64(EM._gram_point_density(m, p.h, p.lat, p.lon, p.t, false)[1]) for p in pts]
end
base = dens(); off(tau) = log.(dens(; mars_mgcm_dust_levels=(tau, 0.3, 1.0)) ./ base)
scan = Dict(string(t) => mean(off(t)) for t in (0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0)); println("scan tau => mean ln(rho / rho_base): ", sort(collect(scan); by=x -> parse(Float64, x[1])))
function bisect(lo, hi, tgt)
    (mean(off(lo)) - tgt) * (mean(off(hi)) - tgt) < 0 || return nothing                    # target not reachable inside [lo, hi]: the scan in the output shows the extreme offset
    for _ in 1:40; mid = (lo + hi) / 2; (mean(off(mid)) < tgt) == (tgt > 0) ? (tgt > 0 ? (lo = mid) : (hi = mid)) : (tgt > 0 ? (hi = mid) : (lo = mid)); end
    return round((lo + hi) / 2; digits=4)
end
res = Dict{String, Any}("scan_mean_ln" => scan, "geometry" => GEOM, "gid" => gid, "epoch_utc" => ep, "n_points" => NP, "flight" => basename(FLIGHT), "alt_m_first_last" => [pts[1].h, pts[end].h])
for (name, tgt, lo, hi) in (("L2p_up", TARGET, 0.30, 2.0), ("L3p_down", -TARGET, 0.02, 0.30))
    tau = bisect(lo, hi, tgt)
    if tau === nothing; res[name] = Dict("tau" => nothing, "target" => tgt, "note" => "not reachable with tau in [$lo, $hi]; see scan_mean_ln"); println("$name: target $tgt NOT reachable with tau in [$lo, $hi]"); continue; end
    d = off(tau); res[name] = Dict("tau" => tau, "mean_ln" => mean(d), "min_ln" => minimum(d), "max_ln" => maximum(d), "target" => tgt)
    println(@sprintf("%s: tau %.4f -> mean ln offset %+.4f (min %+.4f max %+.4f) over %d periapsis points; target %+.4f", name, tau, mean(d), minimum(d), maximum(d), NP, tgt))
    for dt in (-0.005, 0.005); d2 = off(tau + dt); println(@sprintf("   neighbor tau %.4f: mean %+.4f", tau + dt, mean(d2))); end
end
res["base_check_const_0.3"] = mean(off(0.3)); println("const 0.3 vs untouched: mean ln ", res["base_check_const_0.3"])
open(joinpath(OUT, "dust_levels_$(GEOM).json"), "w") do io; JSON.print(io, res, 2); end
