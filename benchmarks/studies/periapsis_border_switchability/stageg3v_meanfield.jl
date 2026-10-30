# Venus gate V-G(V), prerequisite of R8.3: GRAM mean density along the predicted periapsis track of each gate geometry (the Venus counterpart of meanfield_sample.jl, T10; that file is not edited).
# No dynamics and no flight: native GRAM point calls only (mean field, perturbations off). For every pass of every geometry:
#   mf_<geometry>.csv  three altitudes h_ref + alt_offsets_m  (cand,k,alt_idx,alt_m,lat_rad,lon_rad,t_s,rho_kg_m3): the format geometry_rank.load_meanfield reads (H per pass from the secant, ln rho at h_ref)
#   pf_<geometry>.csv  the profile sweep profile_alt_min_m .. profile_alt_max_m every profile_alt_step_m at the periapsis point (cand,k,alt_m,lat_rad,lon_rad,t_s,rho_kg_m3): H(h) and the profile shape
# Tracks: python3 -I stageg3v_tracks.py <out>/tracks.csv.   julia --project=. stageg3v_meanfield.jl --out=<dir> [--workers=N] [--dry-run]
# Rules: a Julia job (uptime and ps first; systemd-run --user --scope -p MemoryMax=24G). Native GRAM is not thread-safe: every GRAM call goes through run_gram_campaign (run_monte_carlo, route :process asserted, one process per
# worker). Each worker builds one GRAM model per epoch and checks the model's core epoch against the requested one (the GRAM epoch defect). Resumable: a geometry file on disk is skipped. Wall times are contended.
using TOML, SHA, Dates
const DRY = "--dry-run" in ARGS
const HERE = @__DIR__
const CFGV = joinpath(HERE, "stageg3v_config.toml")
const C3 = TOML.parsefile(CFGV)["stageg3v"]; const SM = C3["sampling"]
const CLI = let kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    kv
end
DRY || include(joinpath(HERE, "stageg3v_common.jl"))   # SpaceAGORA, GRAMSuite, gram_rho, run_gram_campaign (the infrastructure the gate flights use)

core_epoch_string(y, mo, d, h, mi, s) = string(y, "-", mo, "-", d, " ", h, ":", mi, ":", s)
profile_alts() = collect(Float64(SM["profile_alt_min_m"]):Float64(SM["profile_alt_step_m"]):Float64(SM["profile_alt_max_m"]))

const MF_MODELS = Dict{String, Any}()   # per worker process: one GRAM model per epoch
function mf_env(epoch_utc::String)
    get!(MF_MODELS, epoch_utc) do
        meta = TOML.parsefile(CFGV); sc = load_stagef("venus"; path=CFGV); stagef_env!(false, meta)   # perturbations off: the mean field
        epoch = parse_epoch(epoch_utc); planet = stagef_planet("venus")
        model = stagef_gram_model(sc, epoch, stagef_seed(String(SM["gram_seed_string"])))
        pargs = gram_pargs(planet, Float64(sc["atmosphere_gram"]["entry_interface_km"])); t = model.core.initial_time
        (; model, pargs, core=core_epoch_string(t.year, t.month, t.day, t.hour, t.minute, t.second))
    end
end

"""One case = one geometry. seed: (; run_id, epoch_utc, out, points=[(; cand, k, t_s, lat, lon, h_center)]). Plain-data row; never throws."""
function mf_task(seed)
    t0 = time()
    try
        env = mf_env(String(seed.epoch_utc)); off = Float64.(SM["alt_offsets_m"]); alts = profile_alts(); bad = 0; mf = String[]; pf = String[]; ms = Float64[]
        d = DateTime(String(seed.epoch_utc)); expect = core_epoch_string(year(d), month(d), day(d), hour(d), minute(d), float(second(d)))
        for p in seed.points
            for (ia, o) in enumerate(off)
                alt = p.h_center + o; tc = time_ns(); r = gram_rho(env.model, env.pargs, alt, p.lat, p.lon, p.t_s); push!(ms, (time_ns() - tc) / 1e6)
                (isfinite(r) && r > 0) || (bad += 1); push!(mf, string(p.cand, ",", p.k, ",", ia, ",", alt, ",", p.lat, ",", p.lon, ",", p.t_s, ",", r))
            end
            for alt in alts
                r = gram_rho(env.model, env.pargs, alt, p.lat, p.lon, p.t_s); (isfinite(r) && r > 0) || (bad += 1); push!(pf, string(p.cand, ",", p.k, ",", alt, ",", p.lat, ",", p.lon, ",", p.t_s, ",", r))
            end
        end
        mkpath(seed.out); cand = seed.points[1].cand
        for (nm, hdr, rows) in (("mf", "cand,k,alt_idx,alt_m,lat_rad,lon_rad,t_s,rho_kg_m3", mf), ("pf", "cand,k,alt_m,lat_rad,lon_rad,t_s,rho_kg_m3", pf))
            tmp = joinpath(seed.out, "$(nm)_$(cand).part"); open(tmp, "w") do io; println(io, hdr); foreach(r -> println(io, r), rows); end
            bad == 0 && mv(tmp, joinpath(seed.out, "$(nm)_$(cand).csv"); force=true)
        end
        return (; status = bad == 0 ? "ok" : "bad_rho", run_id=seed.run_id, n_calls=length(mf) + length(pf), n_bad=bad, ms_median_mf=isempty(ms) ? NaN : sort(ms)[(length(ms) + 1) ÷ 2], core_epoch=env.core, epoch_ok=(env.core == expect), wall_s=time() - t0, error="")
    catch err
        return (; status="error", run_id=seed.run_id, n_calls=0, n_bad=0, ms_median_mf=NaN, core_epoch="", epoch_ok=false, wall_s=time() - t0, error=replace(first(sprint(showerror, err), 500), "\n" => " "))
    end
end

write_summary(path, rows) = open(path, "w") do io
    ks = collect(keys(rows[1])); println(io, join(ks, ",")); for r in rows; println(io, join((replace(string(r[k]), "," => ";") for k in ks), ",")); end
end

campaign(seeds, workers) = run_gram_campaign(mf_task, seeds, workers; mission_time_s=Float64(SM["route_hint_mission_s"]), worker_files=[abspath(@__FILE__)])

function dry_run(out)
    tr = joinpath(out, "tracks.csv"); n = isfile(tr) ? count(_ -> true, eachline(tr)) - 1 : 0; na = length(SM["alt_offsets_m"]); np = length(profile_alts())
    println("tracks file $(isfile(tr) ? tr : "MISSING") rows=$n; mf altitudes $(SM["alt_offsets_m"]) around h_ref; profile altitudes $(np) ($(SM["profile_alt_min_m"]) to $(SM["profile_alt_max_m"]) every $(SM["profile_alt_step_m"]) m)")
    println("calls: $(n * (na + np)) native GRAM point calls (mean field, perturbations off); 3 geometries = 3 cases; no GRAM, no SpaceAGORA loaded, nothing flown")
end

function main()
    out = abspath(get(CLI, "out", "")); isempty(get(CLI, "out", "")) && error("--out=<dir> required")
    DRY && return dry_run(out)
    tracks = joinpath(out, "tracks.csv"); isfile(tracks) || error("$tracks missing (python3 -I stageg3v_tracks.py $tracks)")
    T = CSV.read(tracks, DataFrame; types=Dict(:epoch_utc => String, :cand => String)); seeds = NamedTuple[]; skipped = 0
    for g in groupby(T, :cand)
        cand = String(g.cand[1]); isfile(joinpath(out, "mf_$(cand).csv")) && isfile(joinpath(out, "pf_$(cand).csv")) && (skipped += 1; continue)
        push!(seeds, (; run_id="mf_$(cand)", epoch_utc=String(g.epoch_utc[1]), out=out, points=[(; cand, k=Int(r.k), t_s=Float64(r.t_s), lat=Float64(r.lat_gd_rad), lon=Float64(r.lon_rad), h_center=Float64(r.h_center_m)) for r in eachrow(g)]))
    end
    println("sample: $(length(seeds)) geometries to run ($skipped on disk)"); length(seeds) >= 2 || error("a one-case campaign plans route none")
    res = campaign(seeds, parse(Int, get(CLI, "workers", string(SM["workers"]))))
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    res.route === :process || error("route was $(res.route), not :process")
    rows = [s.value for s in res.successful]; isempty(rows) || write_summary(joinpath(out, "summary_sample.csv"), rows)
    open(joinpath(out, "mc_sample.toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"")
    end
    bad = count(r -> r.status != "ok" || !r.epoch_ok, rows); println("rows=$(length(rows)) not ok or epoch mismatch=$bad")
    foreach(r -> println("  $(r.run_id): status=$(r.status) calls=$(r.n_calls) ms_median_mf=$(round(r.ms_median_mf; digits=3)) core_epoch=$(r.core_epoch) epoch_ok=$(r.epoch_ok) $(r.error)"), rows)
    (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
    files = sort(filter(f -> endswith(f, ".csv"), readdir(out)))
    open(joinpath(out, "manifest.sha256"), "w") do io; for f in files; println(io, bytes2hex(sha256(read(joinpath(out, f)))), "  ", f); end; end
    println("manifest.sha256 written for $(length(files)) files")
end
(abspath(PROGRAM_FILE) == @__FILE__) && main()
