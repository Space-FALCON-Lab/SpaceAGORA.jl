# WP-G1: GRAM mean-density sampler along each candidate geometry's predicted periapsis track (plan R2.3 "the one Julia job"; field term T10 of the geometry ranking).
# No dynamics and no flight: native GRAM point calls only (mean field, perturbations off), at three altitudes around h_ref_m for every pass of every candidate.
#
#   python3 -I geometry_rank.py --tracks <dir>/tracks.csv          # candidate tracks (cand, epoch, k, t_s, lat, lon, in_cut256/384)
#   julia --project=. meanfield_sample.jl --phase=probe  --out=<dir>                       # ~1,000 calls first: per-call cost, decides the grid cut
#   julia --project=. meanfield_sample.jl --phase=sample --out=<dir> --tracks=<dir>/tracks.csv [--cut=none|256|384]
#   julia --project=. meanfield_sample.jl --dry-run [--tracks=<file>]      # counts and plan only: loads no GRAM, no SpaceAGORA, flies nothing
# Rules (plan R2.3, project rules): a Julia job, so the host must be quiet (uptime; ps -eo pcpu,nlwp,args | grep julia) and it waits for the coordinator's RELEASE; memory capped
#   systemd-run --user --scope -p MemoryMax=24G julia ...
# Native GRAM is not thread-safe: every GRAM call goes through run_gram_campaign (run_monte_carlo, route :process asserted, one process per worker, --threads=1). The probe is two
# cases (a one-case campaign plans route none). Each worker builds one GRAM model per candidate epoch and checks the model's core epoch against the requested one (GRAM epoch defect).
# Output <out>/mf_<epochtag>_<chunk>.csv: cand,k,alt_idx,alt_m,lat_rad,lon_rad,t_s,rho_kg_m3 (resumable: a chunk file on disk is skipped), summary_<phase>.csv, mc_<phase>.toml,
# probe.toml (probe), manifest.sha256 (sample). Wall times are contended: never performance numbers.
using TOML, SHA, Dates
const DRY = "--dry-run" in ARGS
const HERE = @__DIR__
const CFGR = TOML.parsefile(joinpath(HERE, "geometry_rank_config.toml")); const SM = CFGR["sampling"]
const CFGO = joinpath(HERE, "stageo_config.toml")
const CLI = let kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    kv
end
DRY || include(joinpath(HERE, "stageo_common.jl"))   # SpaceAGORA, GRAMSuite, the compatibility shim, gram_rho, run_gram_campaign (the infrastructure the Stage O runs use)

# ---- pure helpers (no GRAM) -----------------------------------------------------------------------------------------------------------------------

"""GRAM core epoch string in the form the Stage runs record: year-month-day hour:minute:second, no zero padding."""
core_epoch_string(y, mo, d, h, mi, s) = string(y, "-", mo, "-", d, " ", h, ":", mi, ":", s)
epoch_tag(epoch_utc::AbstractString) = replace(epoch_utc, r"[^0-9]" => "")

"""Deterministic probe points (no tracks needed): `n` calls spread over latitude, longitude and the first 70 days, altitudes cycling through the three offsets."""
function probe_points(n::Int)
    off = Float64.(SM["alt_offsets_m"]); h0 = Float64(CFGR["orbit"]["h_ref_m"])
    return [(; lat=deg2rad(-80 + 160 * mod(0.6180339887 * j, 1)), lon=2pi * mod(0.7548776662 * j, 1), t_s=86400.0 * 70 * mod(0.5698402910 * j, 1), alt_m=h0 + off[mod1(j, length(off))]) for j in 1:n]
end

median_of(v) = (s = sort(v); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1]))

"""Plan counts from the config, optionally from a tracks file (header cols cand, in_cut256, in_cut384; one row per pass of a candidate)."""
function plan_counts(tracks::Union{Nothing, String})
    g = CFGR["grid"]; na = length(SM["alt_offsets_m"]); np = Int(SM["passes_nominal"])
    ncand = length(0:g["omega0_step_deg"]:359.999) * length(0:g["raan0_step_deg"]:359.999) * length(g["epoch_offsets_days"]) * length(g["inclinations_deg"])
    c = Dict{String, Any}("grid_candidates" => ncand, "calls_nominal" => ncand * np * na, "alt" => na, "passes_nominal" => np)
    if tracks !== nothing && isfile(tracks)
        hdr = Dict{String, Int}(); full = 0; c256 = 0; c384 = 0; anchor = 0; cands = Set{String}()
        for (i, l) in enumerate(eachline(tracks))
            f = split(l, ","); if i == 1; for (j, h) in enumerate(f); hdr[h] = j; end; continue; end
            id = f[hdr["cand"]]; push!(cands, id)
            id == "ODY" && (anchor += 1; continue)
            full += 1; f[hdr["in_cut256"]] == "1" && (c256 += 1); f[hdr["in_cut384"]] == "1" && (c384 += 1)
        end
        merge!(c, Dict("tracks_grid_rows" => full, "tracks_anchor_rows" => anchor, "tracks_calls_full" => full * na, "tracks_calls_cut256" => c256 * na, "tracks_calls_cut384" => c384 * na,
                       "tracks_candidates" => length(cands)))
    end
    return c
end

function dry_run()
    c = plan_counts(get(CLI, "tracks", nothing)); lo, hi = SM["est_ms_per_call_low"], SM["est_ms_per_call_high"]
    println("candidates (R2.3 grid): $(c["grid_candidates"]) (plan: 1,536) -> ", c["grid_candidates"] == 1536 ? "MATCH" : "MISMATCH")
    println("nominal calls: $(c["grid_candidates"]) x $(c["passes_nominal"]) passes x $(c["alt"]) altitudes = $(c["calls_nominal"]) (plan: about 0.92 million) -> ", round(c["calls_nominal"] / 1e6, digits=2) == 0.92 ? "MATCH" : "MISMATCH")
    haskey(c, "tracks_calls_full") ? println("tracks file: $(c["tracks_candidates"]) candidates incl. the Odyssey anchor; grid rows $(c["tracks_grid_rows"]) -> calls full $(c["tracks_calls_full"]), cut256 $(c["tracks_calls_cut256"]), cut384 $(c["tracks_calls_cut384"]); anchor rows $(c["tracks_anchor_rows"]) -> $(c["tracks_anchor_rows"] * c["alt"]) calls (always sampled)") :
        println("no tracks file given (--tracks=...): real per-candidate pass counts come from geometry_rank.py --tracks")
    println("probe: $(SM["probe_calls"]) calls = $(SM["probe_workers"]) cases of $(SM["probe_calls"] ÷ SM["probe_workers"]) (a one-case campaign plans route none); first $(SM["probe_warmup_calls"]) calls per case excluded; limit $(SM["probe_limit_ms"]) ms/call, above it the full grid is refused (--cut=256 or 384)")
    for (nm, calls) in (("full", c["calls_nominal"]),)
        println("estimated job ($nm, ASSUMED $lo to $hi ms/call): $(round(calls * lo / 6e4, digits=1)) to $(round(calls * hi / 6e4, digits=1)) min single process = $(round(calls * lo / 3.6e6, digits=2)) to $(round(calls * hi / 3.6e6, digits=2)) worker-h; $(SM["workers"]) workers, $(SM["cands_per_chunk"]) candidates per chunk")
    end
    println("dry run: no GRAM, no SpaceAGORA loaded, nothing flown")
end

# ---- worker task -----------------------------------------------------------------------------------------------------------------------------------

const MF_MODELS = Dict{String, Any}()   # per worker process: one GRAM model per epoch
function mf_env(epoch_utc::String)
    get!(MF_MODELS, epoch_utc) do
        meta = TOML.parsefile(CFGO); sc = load_stagef("mars"; path=CFGO); stagef_env!(false, meta)   # perturbations off: the mean field
        epoch = parse_epoch(epoch_utc); planet = stagef_planet("mars")
        model = stagef_gram_model(sc, epoch, stagef_seed(String(SM["gram_seed_string"])))
        pargs = gram_pargs(planet, Float64(sc["atmosphere_gram"]["entry_interface_km"])); t = model.core.initial_time
        (; model, pargs, core=core_epoch_string(t.year, t.month, t.day, t.hour, t.minute, t.second))
    end
end

"""One case. seed: (; kind="probe"|"chunk", run_id, epoch_utc, out, points=[(; cand, k, t_s, lat, lon, h_center)] for chunk | probe points for probe). Plain-data row; never throws."""
function mf_task(seed)
    t0 = time()
    try
        env = mf_env(String(seed.epoch_utc)); off = Float64.(SM["alt_offsets_m"]); warm = Int(SM["probe_warmup_calls"]); ms = Float64[]; first_ms = NaN; bad = 0; rows = String[]
        d = DateTime(String(seed.epoch_utc)); expect = core_epoch_string(year(d), month(d), day(d), hour(d), minute(d), float(second(d)))
        if seed.kind == "probe"
            for (j, p) in enumerate(seed.points)
                tc = time_ns(); r = gram_rho(env.model, env.pargs, p.alt_m, p.lat, p.lon, p.t_s); dt = (time_ns() - tc) / 1e6
                j == 1 && (first_ms = dt); j > warm && push!(ms, dt); (isfinite(r) && r > 0) || (bad += 1)
            end
        else
            for p in seed.points, (ia, o) in enumerate(off)
                alt = p.h_center + o; tc = time_ns(); r = gram_rho(env.model, env.pargs, alt, p.lat, p.lon, p.t_s); dt = (time_ns() - tc) / 1e6
                isempty(ms) && (first_ms = dt); push!(ms, dt); (isfinite(r) && r > 0) || (bad += 1)
                push!(rows, string(p.cand, ",", p.k, ",", ia, ",", alt, ",", p.lat, ",", p.lon, ",", p.t_s, ",", r))
            end
            mkpath(seed.out); tmp = joinpath(seed.out, seed.run_id * ".part")
            open(tmp, "w") do io; println(io, "cand,k,alt_idx,alt_m,lat_rad,lon_rad,t_s,rho_kg_m3"); foreach(r -> println(io, r), rows); end
            bad == 0 && mv(tmp, joinpath(seed.out, seed.run_id * ".csv"); force=true)
        end
        return (; status = bad == 0 ? "ok" : "bad_rho", run_id=seed.run_id, kind=seed.kind, n_calls=length(seed.points) * (seed.kind == "probe" ? 1 : length(off)), n_bad=bad,
                 ms_median=isempty(ms) ? NaN : median_of(ms), ms_mean=isempty(ms) ? NaN : sum(ms) / length(ms), ms_max=isempty(ms) ? NaN : maximum(ms), first_call_ms=first_ms,
                 core_epoch=env.core, epoch_ok=(env.core == expect), wall_s=time() - t0, error="")
    catch err
        return (; status="error", run_id=seed.run_id, kind=seed.kind, n_calls=0, n_bad=0, ms_median=NaN, ms_mean=NaN, ms_max=NaN, first_call_ms=NaN, core_epoch="", epoch_ok=false,
                 wall_s=time() - t0, error=replace(first(sprint(showerror, err), 500), "\n" => " "))
    end
end

# ---- driver ----------------------------------------------------------------------------------------------------------------------------------------

write_summary(path, rows) = open(path, "w") do io
    ks = collect(keys(rows[1])); println(io, join(ks, ",")); for r in rows; println(io, join((replace(string(r[k]), "," => ";") for k in ks), ",")); end
end

campaign(seeds, workers) = run_gram_campaign(mf_task, seeds, workers; mission_time_s=Float64(SM["route_hint_mission_s"]), worker_files=[abspath(@__FILE__)])

function finish(res, phase, out)
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    res.route === :process || error("route was $(res.route), not :process")
    rows = [s.value for s in res.successful]; mkpath(out); isempty(rows) || write_summary(joinpath(out, "summary_$(phase).csv"), rows)
    open(joinpath(out, "mc_$(phase).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"")
    end
    bad = count(r -> r.status != "ok" || !r.epoch_ok, rows); println("rows=$(length(rows)) not ok or epoch mismatch=$bad")
    foreach(r -> (r.status != "ok" || !r.epoch_ok) && println("  $(r.run_id): status=$(r.status) core_epoch=$(r.core_epoch) epoch_ok=$(r.epoch_ok) $(r.error)"), rows)
    (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
    return rows
end

function phase_probe(out)
    n = Int(SM["probe_calls"]); w = Int(SM["probe_workers"]); per = n ÷ w; epoch = String(CFGR["anchor"]["epoch_utc"])
    seeds = [(; kind="probe", run_id="probe_$(i)", epoch_utc=epoch, out=out, points=probe_points(per)) for i in 1:w]
    rows = finish(campaign(seeds, w), "probe", out)
    med = median_of([r.ms_median for r in rows]); worst = maximum(r.ms_median for r in rows); lim = Float64(SM["probe_limit_ms"])
    println("probe: median ms/call per case $(round.([r.ms_median for r in rows]; digits=3)), first call $(round.([r.first_call_ms for r in rows]; digits=1)) ms; limit $lim ms -> ", worst > lim ? "CUT THE GRID (--cut=256)" : "full grid allowed")
    open(joinpath(out, "probe.toml"), "w") do io
        println(io, "calls = $n\ncases = $w\nmedian_ms_per_call = $med\nworst_case_median_ms_per_call = $worst\nlimit_ms = $lim\ncut_required = $(worst > lim)\nroute = \"process\"\nnote = \"contended (2 concurrent workers), not a performance number\"")
    end
end

function phase_sample(out)
    probe = joinpath(out, "probe.toml"); isfile(probe) || error("run --phase=probe first ($probe missing)")
    cut = get(CLI, "cut", "none"); req = TOML.parsefile(probe)["cut_required"]
    req && cut == "none" && error("the probe exceeded $(SM["probe_limit_ms"]) ms/call: rerun with --cut=256 (or 384)")
    cut in ("none", "256", "384") || error("--cut must be none, 256 or 384")
    tracks = get(CLI, "tracks", ""); isfile(tracks) || error("--tracks=<tracks.csv> required (geometry_rank.py --tracks)")
    T = CSV.read(tracks, DataFrame; types=Dict(:epoch_utc => String, :cand => String)); cut == "none" || (T = T[(T.cand .== "ODY") .| (T[!, "in_cut" * cut] .== 1), :])
    seeds = NamedTuple[]; cpc = Int(SM["cands_per_chunk"]); h = Float64(CFGR["orbit"]["h_ref_m"]); skipped = 0
    for (ep, g) in pairs(groupby(T, :epoch_utc))
        ids = sort(unique(g.cand)); tag = epoch_tag(String(ep.epoch_utc))
        for (ci, chunk) in enumerate(Iterators.partition(ids, cpc))
            rid = "mf_$(tag)_$(lpad(ci, 3, '0'))"; isfile(joinpath(out, rid * ".csv")) && (skipped += 1; continue)
            sub = g[in.(g.cand, Ref(Set(chunk))), :]
            push!(seeds, (; kind="chunk", run_id=rid, epoch_utc=String(ep.epoch_utc), out=out,
                           points=[(; cand=String(r.cand), k=Int(r.k), t_s=Float64(r.t_s), lat=Float64(r.lat_gd_rad), lon=Float64(r.lon_rad), h_center=h) for r in eachrow(sub)]))
        end
    end
    println("sample: $(length(seeds)) chunks to run ($skipped already on disk), cut=$cut, $(sum(length(s.points) for s in seeds) * length(SM["alt_offsets_m"])) calls")
    isempty(seeds) && error("no chunks (all on disk)")
    length(seeds) >= 2 || error("a one-case campaign plans route none")
    finish(campaign(seeds, Int(SM["workers"])), "sample", out)
    files = sort(filter(f -> endswith(f, ".csv"), readdir(out)))
    open(joinpath(out, "manifest.sha256"), "w") do io; for f in files; println(io, bytes2hex(sha256(read(joinpath(out, f)))), "  ", f); end; end
    println("manifest.sha256 written for $(length(files)) files")
end

function main()
    DRY && return dry_run()
    phase = get(CLI, "phase", ""); out = abspath(get(CLI, "out", "")); isempty(get(CLI, "out", "")) && error("--out=<dir> required")
    phase == "probe" && return phase_probe(out); phase == "sample" && return phase_sample(out)
    error("--phase must be probe or sample (or --dry-run)")
end
(abspath(PROGRAM_FILE) == @__FILE__) && main()
