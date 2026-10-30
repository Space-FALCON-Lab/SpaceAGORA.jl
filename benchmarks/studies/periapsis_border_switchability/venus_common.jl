# Planet generality at Venus: the Venus setup for the Mars truth machinery (mars_common.jl, unchanged) and the frozen estimators.
# mars_common.jl is config-driven except for the planet; here the planet is the SPICE-backed Venus() with the atmosphere constants replaced by the
# DERIVED Venus Express single-layer exponential (see venus_generality_config.toml [atmosphere]), so that the engine's ExponentialAtmosphereModel(planet)
# and the drag fit agree exactly (as at Mars). No src/ or examples/ change. Included by run_venus_generality.jl and by process workers.
using Statistics
include(joinpath(@__DIR__, "mars_common.jl"))

const VENUS_CONFIG_PATH = joinpath(@__DIR__, "venus_generality_config.toml")
load_venus(path::AbstractString=VENUS_CONFIG_PATH) = TOML.parsefile(path)

"""DERIVED exponential constants `(rho_ref, h_ref, H)` from the repo's Venus Express per-pass table (rule in the config)."""
function venus_atmosphere(vc, planet, mu)
    a = vc["atmosphere"]; tbl = CSV.read(joinpath(REPO_ROOT, String(a["pass_table"])), DataFrame)
    Hs = unique(tbl.scale_height_km); length(Hs) == 1 || throw(ErrorException("VEX table has several scale heights: $Hs"))
    H = 1000.0 * Hs[1]; h = 1000.0 .* tbl.h_peri_km; h_ref = median(h); ra = Float64(a["vex_apoapsis_radius_m"])
    rp = planet.Rp_e .+ h; sma = 0.5 .* (ra .+ rp); vp2 = mu .* (2 ./ rp .- 1 ./ sma)
    rho = 2.0 .* tbl.pdyn_Nm2 ./ vp2
    return (rho_ref=median(rho .* exp.((h .- h_ref) ./ H)), h_ref=h_ref, H=H, n_passes=length(h))
end

"""Venus bundle (a `MarsSetup`; the struct is planet-agnostic): SPICE Venus with the derived atmosphere, MGNP180U harmonics, exponential fit, beta."""
function venus_setup(vc=load_venus(); degree::Int=Int(vc["gravity"]["degree"]))
    p0 = Venus("", SPICE_PATH)
    harm0 = GravitationalHarmonicsModel(degree, Int(vc["gravity"]["order"]), joinpath(REPO_ROOT, String(vc["gravity"]["coefficient_file"])), p0)
    at = venus_atmosphere(vc, p0, harm0.gm_m3s2)
    planet = Venus(; (f => getfield(p0, f) for f in fieldnames(typeof(p0)))..., ρ_ref=at.rho_ref, h_ref=at.h_ref, H=at.H)
    harm = GravitationalHarmonicsModel(degree, Int(vc["gravity"]["order"]), joinpath(REPO_ROOT, String(vc["gravity"]["coefficient_file"])), planet)
    top = Float64(vc["atmosphere"]["fit_top_altitude_m"])
    fit = CorridorAtmosphereFit([at.h_ref, top], [at.rho_ref, at.rho_ref * exp(-(top - at.h_ref) / at.H)], [at.H], :extend, :extend, Float64(vc["atmosphere"]["temperature_k"]),
        "Venus Express table (vex_final/step2pp.csv): single-layer exponential, rho_ref/h_ref/H derived", Float64[], Float64[])
    v = vc["vehicle"]
    area = Float64(v["bus_dims_m"][1]) * Float64(v["bus_dims_m"][3]) + 2 * Float64(v["panel_dims_m"][2]) * Float64(v["panel_dims_m"][3])
    dry = Float64(v["bus_mass_kg"]) + 2 * Float64(v["panel_mass_each_kg"]); wet = dry + Float64(v["prop_mass_kg"])
    beta_der = wet / (Float64(v["cd"]) * area); bo = Float64(v["beta_override_kg_m2"])
    return MarsSetup(vc, planet, harm, fit, bo > 0 ? bo : beta_der, area, wet, dry, Float64(vc["orbit"]["apoapsis_radius_m"]) - planet.Rp_e)
end

"""One Venus cell: calibrate, then fly the horizon (the Mars cell task with the Venus setup). Plain-data row; never throws."""
function venus_cell_task(seed)
    t0 = time()
    try
        vc = load_venus(seed.config_path); ms = venus_setup(vc)
        epoch = parse_epoch(vc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
        dir = joinpath(seed.out, seed.tag); mkpath(dir)
        seed_alt, recs = mars_calibrate(ms, rot, epoch, dir; inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega)
        write_rows_m(joinpath(dir, "calibration.csv"), recs)
        cal_wall = time() - t0
        r = mars_run(ms, rot, epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir; seed_alt_m=seed_alt, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega)
        write_rows_m(joinpath(dir, "passes_3b.csv"), r.rows)
        return (; tag=seed.tag, variant="full", inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx, status="ok",
            seed_alt_m=seed_alt, cal_wall_s=cal_wall, run_wall_s=r.wall_s, retcode=r.retcode, passes=length(r.rows),
            t_end_days=isempty(r.rows) ? NaN : r.rows[end].time_s / SECONDS_PER_DAY, error="")
    catch err
        return (; tag=seed.tag, variant="full", inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, epoch_idx=seed.epoch_idx, status="run_error",
            seed_alt_m=NaN, cal_wall_s=NaN, run_wall_s=NaN, retcode="", passes=0, t_end_days=NaN, error=replace(first(sprint(showerror, err), 300), "\n" => " "))
    end
end
