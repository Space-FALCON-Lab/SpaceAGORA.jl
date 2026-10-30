# Revision 17b: find the constant MGCM dust level that gives a persistent mean-density offset of about +/-TARGET ln at the periapsis of a nominal Odyssey flight (point queries only, no flights). Bisection on the mean over
# the first NP periapsis points; prints and writes the measured offset (mean, min, max over the points) for the chosen level and its neighbors. Serial, one process.
# Usage: systemd-run --user --scope -p MemoryMax=16G julia --project=. --threads=1 gram_dust_levels_rev17.jl --out=DIR [--target=0.223]
using TOML, CSV, DataFrames, Statistics, Printf, JSON
include(joinpath(@__DIR__, "planner_common.jl"))
const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const CLI = Dict(m.captures[1] => something(m.captures[2], "1") for a in ARGS for m in [match(r"^--([^=]+)(?:=(.*))?$", a)] if m !== nothing)
const OUT = abspath(get(CLI, "out", ".")); mkpath(OUT); const TARGET = parse(Float64, get(CLI, "target", "0.223")); const NP = 40
const FLIGHT = joinpath(get(ENV, "PLANNER_DATA_ROOT", "."), "planner/runs/ODY_PL_L3_pert_ODY_s13")
cfgp = joinpath(@__DIR__, "planner_config.toml"); meta = planner_meta(cfgp); gid = planner_gid(meta, "ODY")
meta, sc, geom, ep = stageg3_load(cfgp, gid); ms0 = stagef_setup(sc); epoch = parse_epoch(ep)
P = CSV.read(joinpath(FLIGHT, "passes.csv"), DataFrame); pts = [(t=P.time_s[k], h=P.alt_m[k], lat=P.lat_rad[k], lon=P.lon_rad[k]) for k in 1:NP]
function dens(; kw...)
    _corridor_gram_load!(); _corridor_gram_require_direct_sampling()
    m = Base.invokelatest(SimulationModel.GRAMAtmosphereModel; planet_name=String(sc["atmosphere_gram"]["planet_name"]), initial_time=epoch, seed=stagef_seed("mars" * "ODY_s13"), kw...)
    return [Float64(EM._gram_point_density(m, p.h, p.lat, p.lon, p.t, false)[1]) for p in pts]
end
base = dens(); off(tau) = log.(dens(; mars_mgcm_dust_levels=(tau, 0.3, 1.0)) ./ base)
function bisect(lo, hi, tgt)
    for _ in 1:40; mid = (lo + hi) / 2; (mean(off(mid)) < tgt) == (tgt > 0) ? (tgt > 0 ? (lo = mid) : (hi = mid)) : (tgt > 0 ? (hi = mid) : (lo = mid)); end
    return round((lo + hi) / 2; digits=4)
end
res = Dict{String, Any}()
for (name, tgt, lo, hi) in (("L2_up", TARGET, 0.30, 0.45), ("L3_down", -TARGET, 0.15, 0.30))
    tau = bisect(lo, hi, tgt); d = off(tau); res[name] = Dict("tau" => tau, "mean_ln" => mean(d), "min_ln" => minimum(d), "max_ln" => maximum(d), "target" => tgt)
    println(@sprintf("%s: tau %.4f -> mean ln offset %+.4f (min %+.4f max %+.4f) over %d periapsis points; target %+.3f", name, tau, mean(d), minimum(d), maximum(d), NP, tgt))
    for dt in (-0.005, 0.005); d2 = off(tau + dt); println(@sprintf("   neighbor tau %.4f: mean %+.4f", tau + dt, mean(d2))); end
end
res["base_check_const_0.3"] = mean(off(0.3)); println("const 0.3 vs untouched: mean ln ", res["base_check_const_0.3"])
open(joinpath(OUT, "dust_levels.json"), "w") do io; JSON.print(io, res, 2); end
