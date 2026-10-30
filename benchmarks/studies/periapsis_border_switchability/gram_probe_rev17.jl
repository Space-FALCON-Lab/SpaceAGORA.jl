# Revision 17 GRAM probe (no flights): point densities at the periapsis of a nominal Odyssey flight for each exposed persistent knob, and the spread of the walk ratio for two perturbation scales.
# Usage (one Julia process, serial; native GRAM is not thread-safe so there is no thread route and no second process):
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. benchmarks/studies/periapsis_border_switchability/gram_probe_rev17.jl --out=DIR [--nseeds=60]
# Points: the periapsis events (time_s, alt_m, lat, lon) of passes 1..NP of the nominal held flight ODY_PL_L3_pert_ODY_s13 (path below, SOURCED local data). Mean density = what the flight RHS uses (EM._gram_point_density).
using TOML, CSV, DataFrames, Statistics, Printf
include(joinpath(@__DIR__, "planner_common.jl"))
const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const CLI = Dict(m.captures[1] => something(m.captures[2], "1") for a in ARGS for m in [match(r"^--([^=]+)(?:=(.*))?$", a)] if m !== nothing)
const OUT = abspath(get(CLI, "out", ".")); mkpath(OUT)
const NP = parse(Int, get(CLI, "np", "40")); const NSEEDS = parse(Int, get(CLI, "nseeds", "60"))
const FLIGHT = get(CLI, "flight", joinpath(get(ENV, "PLANNER_DATA_ROOT", "."), "planner/runs/ODY_PL_L3_pert_ODY_s13"))
cfgp = joinpath(@__DIR__, "planner_config.toml"); meta = planner_meta(cfgp)
gid = planner_gid(meta, "ODY")
meta, sc, geom, ep = stageg3_load(cfgp, gid); ms0 = stagef_setup(sc); epoch = parse_epoch(ep)
P = CSV.read(joinpath(FLIGHT, "passes.csv"), DataFrame); M = CSV.read(joinpath(FLIGHT, "measured.csv"), DataFrame)
pts = [(t=P.time_s[k], h=P.alt_m[k], lat=P.lat_rad[k], lon=P.lon_rad[k]) for k in 1:NP]
seed0 = stagef_seed("mars" * "ODY_s13")
function build(; seed=seed0, kw...)
    _corridor_gram_load!(); _corridor_gram_require_direct_sampling()
    return Base.invokelatest(SimulationModel.GRAMAtmosphereModel; planet_name=String(sc["atmosphere_gram"]["planet_name"]), initial_time=epoch, seed=seed, kw...)
end
mean_at(model, p) = Float64(EM._gram_point_density(model, p.h, p.lat, p.lon, p.t, false)[1])
variants = Pair{String, NamedTuple}["base" => (;),
    "f107=100" => (; mars_f107=100.0), "f107=130" => (; mars_f107=130.0), "f107=200" => (; mars_f107=200.0), "f107=68(explicit default)" => (; mars_f107=68.0),
    "dust const 0.3" => (; mars_mgcm_dust_levels=(0.3, 0.3, 1.0)), "dust const 1.0" => (; mars_mgcm_dust_levels=(1.0, 0.3, 1.0)), "dust const 2.0" => (; mars_mgcm_dust_levels=(2.0, 0.3, 1.0)),
    "f107=50" => (; mars_f107=50.0), "f107=30" => (; mars_f107=30.0),
    "dust const 0.1" => (; mars_mgcm_dust_levels=(0.1, 0.3, 1.0)), "dust const 0.2" => (; mars_mgcm_dust_levels=(0.2, 0.3, 1.0)), "dust const 0.35" => (; mars_mgcm_dust_levels=(0.35, 0.3, 1.0)), "dust const 0.4" => (; mars_mgcm_dust_levels=(0.4, 0.3, 1.0)),
    "dust const 0.5" => (; mars_mgcm_dust_levels=(0.5, 0.3, 1.0)), "dust const 0.6" => (; mars_mgcm_dust_levels=(0.6, 0.3, 1.0)), "dust const 0.8" => (; mars_mgcm_dust_levels=(0.8, 0.3, 1.0)),
    "dust seasonal, min=max=0.3" => (; mars_mgcm_dust_levels=(0.0, 0.3, 0.3)), "dust seasonal, min=max=1.0" => (; mars_mgcm_dust_levels=(0.0, 1.0, 1.0)), "dust seasonal, 0.3..1.0 (explicit default)" => (; mars_mgcm_dust_levels=(0.0, 0.3, 1.0)), "map_year=0" => (; mars_map_year=0),
    "map_year=1" => (; mars_map_year=1), "map_year=2" => (; mars_map_year=2), "min_max=0" => (; mars_min_max=0), "scales=(2,1,1,1)" => (; gram_perturbation_scales=(2.0, 1.0, 1.0, 1.0))]
rows = NamedTuple[]; const BASE = Ref(Float64[])
for (name, kw) in variants
    t0 = time(); rho = try
        m = build(; kw...); [mean_at(m, p) for p in pts]
    catch e
        println("variant $name failed: ", sprint(showerror, e)); fill(NaN, NP)
    end
    name == "base" && (BASE[] = rho)
    d = log.(rho ./ BASE[]); fd = filter(isfinite, d); fd = isempty(fd) ? [NaN] : fd
    println(@sprintf("%-28s mean ln(rho/rho_base) over %d periapsis points = %+.4f  (min %+.4f  max %+.4f)  [%.1f s]", name, NP, mean(fd), minimum(fd), maximum(fd), time() - t0))
    for k in 1:NP; push!(rows, (variant=name, pass=k, h_m=pts[k].h, rho=rho[k], ln_ratio=d[k])); end
end
CSV.write(joinpath(OUT, "probe_means.csv"), DataFrame(rows))
# walk ratio spread at the periapsis point of pass 1: fresh model per seed, one query, read perturbedDensity / density (the same ratio the pass mode applies)
srows = NamedTuple[]
for sc_ in (1.0, 1.25, 1.5, 1.75, 2.0, 2.5)
    lr = Float64[]
    for k in 1:NSEEDS
        m = build(; seed=1000 + k, gram_perturbation_scales=(sc_, 1.0, 1.0, 1.0)); core = m.core
        mean_at(m, pts[1]); ds = Base.invokelatest(getproperty(core.gram, :get_density_state), core.gram_atmosphere)
        push!(lr, log(Float64(ds.perturbedDensity) / Float64(ds.density)))
    end
    push!(srows, (scale=sc_, n=NSEEDS, sd_ln=std(lr), mean_ln=mean(lr), frac_r_gt_2_6=mean(lr .> log(2.6)), frac_r_gt_1_4=mean(lr .> log(1.4)), min_ln=minimum(lr), max_ln=maximum(lr)))
    println(@sprintf("  fractions r>2.6: %.3f, r>1.4: %.3f, min ln %.3f max ln %.3f", mean(lr .> log(2.6)), mean(lr .> log(1.4)), minimum(lr), maximum(lr)))
    println(@sprintf("density scale %.1f: ln(perturbed/mean) at the pass-1 periapsis point, %d fresh seeds: sd %.4f mean %+.4f", sc_, NSEEDS, std(lr), mean(lr)))
end
CSV.write(joinpath(OUT, "probe_noise.csv"), DataFrame(srows))
# flown record: ln r (measured.csv r = perturbed/mean at the pass) over the passes of this nominal flight
println(@sprintf("flown nominal (scale 1) ln r over the passes of %s: sd %.4f mean %+.4f (n=%d)", basename(FLIGHT), std(log.(M.r)), mean(log.(M.r)), nrow(M)))
