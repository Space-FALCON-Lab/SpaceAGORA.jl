# Stage O driver: the periapsis-corridor keeper against the flown Mars Odyssey main phase (stageo_config.toml). One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. stageo_run.jl --phase=selftest|cal|smoke|gate|dev1|dev2|held|heldp|repro|xhost --results-root=DIR [--workers=N] --host=NAME [--dry-run]
# --dry-run flies nothing (it loads the common file, so the arm codes come from the real code) and prints the case counts of the phase (or of every phase when --phase is omitted) against the plan's section 6 counts.
# Trees below <root>/stageo/mars: cal/<band>/{guess,nodes,calibration}.csv (phase cal), runs/<code>_pert_<tag>/ (dev1, dev2, held, heldp), smoke/, gate/, repro/, xhost_<host>/, calibrated_arms.csv
# (columns arm,g_a,t_mult for arms "N2-1" and "N2-I"; written by the analysis calibrate step and read by dev2/held/repro/xhost). The Kalman constants are READ from <root>/stagen_flight/mars/kalman_params.csv (Stage N-flight,
# unchanged). Runs already on disk (passes.csv present) are skipped. Rows: <root>/stageo/summary_<phase>_mars_<host>.csv (the Stage N2 columns plus stop_reason; a summary is never overwritten: use a DISTINCT --host stamp per call).
# Arms (plan 4.4): F-O, N2-T, N2-1, N2-I fly the Odyssey band; F-P, N2-P fly the paper band.
#   cal    2 cases, one fit-only calibration per band (pass 1 at the band's target Q)       smoke  F-O and N2-T, seed ODY_s1, 3 days (2 cases: one case plans route none)
#   gate   G1: F-O, ODY_s1, full horizon at dt_max_atmosphere 1.0 s and 0.2 s                dev1   12 dev seeds x {F-O, N2-T}        dev2  12 dev seeds x {N2-1, N2-I}
#   held   24 held seeds x {F-O, N2-1, N2-I}                                                 heldp  24 held seeds x {F-P, N2-P}
#   repro  the first two held seeds of N2-1 into repro/ (bitwise check against the held runs) xhost  N2-1 on the cross-host seed, plus one 3-day F-O padding flight (a one-case campaign plans route none)
using TOML
const CLI = let kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    kv
end
const DRY = haskey(CLI, "dry-run")
const CFGO = joinpath(@__DIR__, "stageo_config.toml")
const META = TOML.parsefile(CFGO); const SO = META["stageo"]
include(joinpath(@__DIR__, "stageo_common.jl"))
const LENIENT = DRY || get(CLI, "phase", "") == "selftest"   # no calibrated_arms.csv / kalman_params.csv needed to count cases
const ROOT = abspath(get(CLI, "results-root", ".")); const SFS = joinpath(ROOT, String(META["stagef"]["results_subdir"])); const SFO = joinpath(SFS, "mars"); const SFK = joinpath(ROOT, "stagen_flight")
const HOST = get(CLI, "host", "local")
const PHASES = ["cal", "smoke", "gate", "dev1", "dev2", "held", "heldp", "repro", "xhost"]
const ENS_PHASES = ["ens", "ensn", "xh0"]   # step-size ensemble, outside the 176-flight plan
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(META["stagef"]["workers"])
function seedrange(r)
    rng = Int(r[1]):Int(r[2])
    haskey(CLI, "seeds") || return rng
    lo, hi = parse.(Int, split(CLI["seeds"], ":")); return max(first(rng), lo):min(last(rng), hi)   # --seeds=a:b splits a phase by seed across hosts
end
horizon_days(key) = (A = SO["arms"][key]; haskey(A, "horizon_days") ? Float64(A["horizon_days"]) : Float64(SO["horizon_factor"]) * Float64(A["predicted_passes"]) * Float64(SO["flown"]["days"]) / Float64(SO["flown"]["passes"]))   # DERIVED
dtlabel(dt) = replace(string(dt), "." => "p")

"""Kalman constants (unchanged Stage N-flight values); NaN placeholders only for --dry-run."""
function kconst()
    f = joinpath(SFK, "mars", "kalman_params.csv")
    if !isfile(f)
        LENIENT && return (; kr=NaN, kq=NaN, km0=NaN, kv0=NaN)
        error("missing $f (Stage N-flight kalman_params.csv)")
    end
    P = CSV.read(f, DataFrame)[1, :]; return (; kr=Float64(P.r), kq=Float64(P.q), km0=Float64(P.m0), kv0=Float64(P.v0))
end
"""(g_a, t_mult) of a calibrated arm from calibrated_arms.csv; --dry-run falls back to the N2-T values (counts only)."""
function calibrated_pair(name)
    f = joinpath(SFO, "calibrated_arms.csv")
    if !isfile(f)
        LENIENT && return (Float64(SO["arms"]["N2T"]["g_a"]), Float64(SO["arms"]["N2T"]["t_mult"]))
        error("missing $f (the analysis calibrate step writes it: columns arm,g_a,t_mult)")
    end
    D = CSV.read(f, DataFrame); i = findfirst(==(name), String.(D.arm)); i === nothing && error("arm $name not in $f"); return (Float64(D.g_a[i]), Float64(D.t_mult[i]))
end

"""One flight seed for arm `key` and seed index k (or tag override). `dir` is the output directory; `role` is "flight" or "pad" (route padding, not counted as a plan flight)."""
function mkseed(key, k, dir; horizon=nothing, dt=NaN, suffix="", role="flight", tag=nothing)
    A = SO["arms"][key]; fixed = A["ctrl"] == "fixed"; K = kconst()
    g, t = fixed ? (NaN, NaN) : get(A, "calibrated", false) ? calibrated_pair(A["name"]) : (Float64(A["g_a"]), Float64(A["t_mult"]))
    eps = fixed ? NaN : Float64(A["eps"]); w = fixed ? NaN : Float64(SO["n2"]["w_cap_m"]); margin = Float64(SO["n2"]["fixed_margin_m"])
    code = stageo_arm_code(key, A["ctrl"], margin, eps, t, g, w); tg = tag === nothing ? String(SO["seed_prefix"]) * string(k) : tag
    return (; config_path=CFGO, band=String(A["band"]), tag=tg, inc=Float64(SO["inclination_deg"]), raan=Float64(SO["raan_deg"]), omega=Float64(SO["omega_deg"]), arm=key, ctrl=String(A["ctrl"]), code,
        eps, w_m=w, aim=g, t_mult=t, kr=K.kr, kq=K.kq, km0=K.km0, kv0=K.kv0, run_id="$(code)_pert_$(tg)$(suffix)", out=dir, cal_root=joinpath(SFO, "cal"),
        horizon_days=horizon === nothing ? horizon_days(key) : Float64(horizon), dt_atm_s=Float64(dt), role=role, name=String(A["name"]))
end

"""Cases of a phase (flights; the cal phase returns calibration seeds). Skipping of runs already on disk is done by the caller."""
function cases(phase)
    runs = joinpath(SFO, "runs"); sd = Float64(SO["smoke_days"]); dev = seedrange(SO["dev_seeds"]); held = seedrange(SO["held_seeds"]); out = NamedTuple[]
    if phase == "cal"
        for band in ("odyssey", "paper")
            push!(out, (; config_path=CFGO, band, gram_tag=String(SO["cal_gram_tag"]), inc=Float64(SO["inclination_deg"]), raan=Float64(SO["raan_deg"]), omega=Float64(SO["omega_deg"]), out=joinpath(SFO, "cal"), run_id="cal_$band", role="cal"))
        end
    elseif phase == "smoke"
        for key in ("FO", "N2T"); push!(out, mkseed(key, 1, joinpath(SFO, "smoke"); horizon=sd)); end
    elseif phase == "gate"
        for dt in SO["gate_dt_s"]; push!(out, mkseed("FO", Int(SO["gate_seed"]), joinpath(SFO, "gate"); dt=Float64(dt), suffix="_dt" * dtlabel(Float64(dt)))); end
    elseif phase == "dev1"
        for k in dev, key in ("FO", "N2T"); push!(out, mkseed(key, k, runs)); end
    elseif phase == "dev2"
        for k in dev, key in ("N21", "N2I"); push!(out, mkseed(key, k, runs)); end
    elseif phase == "held"
        for k in held, key in ("FO", "N21", "N2I"); push!(out, mkseed(key, k, runs)); end
    elseif phase == "heldp"
        for k in held, key in ("FP", "N2P"); push!(out, mkseed(key, k, runs)); end
    elseif phase == "repro"
        for k in first(held, Int(SO["repro_n"])); push!(out, mkseed("N21", k, joinpath(SFO, "repro"))); end
    elseif phase == "xhost"
        d = joinpath(SFO, "xhost_" * HOST)
        push!(out, mkseed("N21", Int(SO["xhost_seed"]), d)); push!(out, mkseed("FO", 1, d; horizon=sd, suffix="_pad", role="pad"))
    elseif phase == "xh0"   # cross-host determinism check before the production phases: the full-length F-O ODY_s1 flight (reference: the local gate/ensemble flight) plus a 3-day padding flight (a one-case campaign plans route none)
        d = joinpath(SFO, "xhost0_" * HOST)
        push!(out, mkseed("FO", 1, d; horizon=Float64(SO["ensemble"]["horizon_days"]))); push!(out, mkseed("FO", 1, d; horizon=sd, suffix="_pad", role="pad"))
    elseif phase == "ens" || phase == "ensn"
        key = phase == "ens" ? "FO" : "N2T"; E = SO["ensemble"]; d = joinpath(SFO, "ensemble")
        for dt in sort(Float64.(E["dts"])), k in seedrange(E["seeds"]); push!(out, mkseed(key, k, d; horizon=Float64(E["horizon_days"]), dt=Float64(dt), suffix="_dt" * dtlabel(Float64(dt)))); end
    else
        error("unknown phase $phase")
    end
    return out
end
flights(cs) = count(c -> c.role != "pad" && c.role != "cal", cs)
ondisk(c) = c.role == "cal" ? isfile(joinpath(c.out, c.band, "guess.csv")) : isfile(joinpath(c.out, c.run_id, "passes.csv"))

function dry_run(which)
    pc = SO["plan_counts"]; total = 0; ok = true
    for ph in which
        cs = cases(ph); n = flights(cs); pads = count(c -> c.role == "pad", cs); done = count(ondisk, cs)
        if ph == "cal"
            println(rpad("cal", 6), "cases=$(length(cs)) (plan: $(pc["cal_cases"]) calibration cases, one per band)", done > 0 ? "  [$done on disk]" : ""); ok &= length(cs) == pc["cal_cases"]; continue
        end
        arms = join(sort(unique(["$(c.name)@$(c.band)" for c in cs if c.role != "pad"])), ", "); hz = join(sort(unique(round.(getfield.(filter(c -> c.role != "pad", cs), :horizon_days); digits=1))), "/")
        # xhost: the plan counts 2 flights over BOTH hosts (1 per host); the padding flight is extra and not counted
        expect = ph == "xhost" ? pc["xhost"] ÷ 2 : pc[ph]
        println(rpad(ph, 6), "flights=$n (plan: $(ph == "xhost" ? "$(pc["xhost"]) over both hosts, $expect per host" : string(pc[ph])))", pads > 0 ? " + $pads padding flight(s) not in the plan count" : "", "  horizon_days=$hz  arms=[$arms]", done > 0 ? "  [$done on disk]" : "")
        ok &= n == expect; total += ph == "xhost" ? 2expect : n
    end
    if length(which) == length(PHASES)
        println("total flights over all phases (xhost counted for both hosts) = $total (plan: $(pc["total"]))"); ok &= total == pc["total"]
    end
    println(ok ? "dry run: counts match plan section 6; nothing flown" : "dry run: COUNT MISMATCH against plan section 6; nothing flown"); ok || exit(1)
end

"""Exact checks without any flight: config parse and band tables, arm codes, case counts against plan section 6, pass bookkeeping and the stop rule, vis-viva periapsis speed, Odyssey start geometry."""
function selftest()
    # config
    bo = SO["bands"]["odyssey"]; bp = SO["bands"]["paper"]; cb = META["common"]["band"]
    @assert abs(sqrt(bo["q_low_w_cm2"] * bo["q_high_w_cm2"]) / bo["q_target_w_cm2"] - 1) < 0.01 "Odyssey target is not the geometric mean"
    @assert all(bp[k] == cb[k] for k in keys(cb)) "paper band differs from [common.band]"
    for band in ("odyssey", "paper")
        meta, sc = stageo_load(CFGO, band; dt_atm_s=0.2); scc = stageo_cal_view(sc)
        @assert sc["orbit"]["apoapsis_radius_m"] == 28533000.0 && sc["integration"]["dt_max_atmosphere_s"] == 0.2 && sc["band"]["q_high_w_cm2"] == SO["bands"][band]["q_high_w_cm2"]
        @assert abs(sqrt(scc["band"]["q_low_w_cm2"] * scc["band"]["q_high_w_cm2"]) / scc["band"]["q_target_w_cm2"] - 1) < 1e-12 && sc["band"]["q_low_w_cm2"] == SO["bands"][band]["q_low_w_cm2"]
        ms = stagef_setup(sc); @assert abs(ms.apo_alt_m + ms.planet.Rp_e - SO["ra0_m"]) < 1e-6
    end
    @assert stageo_load(CFGO, "odyssey")[2]["integration"]["dt_max_atmosphere_s"] == 1.0
    # arm codes and case tables
    @assert stageo_arm_code("FO", "fixed", 1000.0, NaN, NaN, NaN, NaN) == "FO_fixed1000" && stageo_arm_code("N2T", "cc", 1000.0, 0.01, 1.52, 0.10, 4000.0) == "N2T_cc01_a010_t152_w04000"
    for ph in PHASES[2:end]
        cs = cases(ph); ids = [(c.out, c.run_id) for c in cs]; @assert length(unique(ids)) == length(ids) "duplicate run ids in $ph"
    end
    pc = SO["plan_counts"]; tot = sum(flights(cases(ph)) for ph in PHASES[2:end] if ph != "xhost") + pc["xhost"]
    @assert tot == pc["total"] && flights(cases("xhost")) == pc["xhost"] ÷ 2 && length(cases("cal")) == pc["cal_cases"]
    @assert all(c.horizon_days == SO["smoke_days"] for c in cases("smoke")) && all(c.band == "odyssey" for c in cases("held")) && all(c.band == "paper" for c in cases("heldp"))
    @assert horizon_days("FO") == 130.0 && horizon_days("N2P") == 195.0
    # pass bookkeeping: two passes, the second closes at a radius below the target
    pk = PassBook(7.408e6)
    book_sample!(pk, 1e-9); book_periapsis!(pk, 3e-9); book_sample!(pk, 2e-9)
    @assert !book_apoapsis!(pk, 2.8e7) && pk.recs == [(; pass=1, rho_peak_kg_m3=3e-9, ra_m=2.8e7)] && pk.cur_max == 0.0
    book_periapsis!(pk, 5e-9); @assert book_apoapsis!(pk, 7.4e6) && pk.stopped && pk.recs[end] == (; pass=2, rho_peak_kg_m3=5e-9, ra_m=7.4e6)
    pk2 = PassBook(7.408e6); @assert !book_apoapsis!(pk2, 1.0e6) && isempty(pk2.recs)   # no pass flown yet: never stops, logs nothing
    pk3 = PassBook(7.408e6); book_periapsis!(pk3, 4e-9); book_flush!(pk3); @assert length(pk3.recs) == 1 && isnan(pk3.recs[1].ra_m) && pk3.recs[1].rho_peak_kg_m3 == 4e-9
    # vis-viva periapsis speed: start at the apoapsis of an orbit with the Odyssey radii, v tangential
    mu = 42828.375e9; ra = 28533e3; rp = 3.49e6; a = (ra + rp) / 2; va = sqrt(mu * (2 / ra - 1 / a)); vp = sqrt(mu * (2 / rp - 1 / a))
    @assert abs(periapsis_speed([ra, 0.0, 0.0], [0.0, va, 0.0], mu) / vp - 1) < 1e-12 && abs(periapsis_speed([rp, 0.0, 0.0], [0.0, vp, 0.0], mu) / vp - 1) < 1e-12
    a2 = (7.408e6 + rp) / 2; vp2 = sqrt(mu * (2 / rp - 1 / a2))
    println("stageo selftest ok: counts match plan section 6 (total $tot); periapsis speed $(round(vp)) -> $(round(vp2)) m/s over the main phase (ratio $(round(vp2 / vp; digits=3)), cube $(round((vp2 / vp)^3; digits=3)); the track_vrel issue)")
end

"""Smoke-phase stop rule (plan section 9): the first-pass Q of every smoke flight must lie in [0.05, 0.54] W/cm2."""
function smoke_check(res)
    bad = false
    for s in res.successful
        f = joinpath(SFO, "smoke", s.value.run_id, "passes.csv"); isfile(f) || continue
        r = CSV.read(f, DataFrame)[1, :]; q = Float64(r.q_eff_w_cm2); println("smoke $(s.value.run_id): first-pass q_eff = $(round(q; digits=4)) W/cm2 (mean-density Q $(round(Float64(r.q_w_cm2); digits=4)))")
        0.05 <= q <= 0.54 || (bad = true; println("STOP RULE: first-pass Q outside [0.05, 0.54]"))
    end
    bad && error("smoke first-pass Q outside [0.05, 0.54]")
end

function finish(res, phase)
    tag = phase * "_mars_" * HOST
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1)) (contended wall time, not a performance number)")
    foreach(f -> println("FAILED $(f.seed.run_id): $(first(sprint(showerror, f.error), 600))"), res.failed)
    rows = [merge(s.value, (; host=HOST)) for s in res.successful]
    write_rows_m(joinpath(SFS, "summary_$(tag).csv"), rows)
    open(joinpath(SFS, "mc_$(tag).toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_density_family = \"gram\"\nhost = \"$HOST\"")
    end
    res.route === :process || error("route was $(res.route), not :process")
    bad = count(r -> r.status != "ok", rows); println("rows=$(length(rows)) not ok=$bad")
    for r in rows; r.status != "ok" && println("  $(get(r, :run_id, r.tag)): $(r.error)"); end
    phase == "cal" ? foreach(r -> r.status == "ok" && println("  cal $(r.tag): seed_alt_m=$(round(r.seed_alt_m; digits=1)) q_first=$(round(r.q_first_w_cm2; digits=4)) target=$(r.q_target_w_cm2) vrel_first=$(round(r.vrel_first_mps; digits=1))"), rows) :
    for r in rows; r.status == "ok" && println("  $(r.run_id): passes=$(r.passes) t_end_days=$(round(r.t_end_days; digits=2)) stop_reason=$(r.stop_reason) dv=$(round(r.total_dv_mps; digits=4)) burns=$(r.burns)"); end
    phase == "smoke" && smoke_check(res)
    (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
end
campaign(seeds, mission_s) = run_gram_campaign(seeds[1] isa NamedTuple && haskey(seeds[1], :gram_tag) ? stageo_cal_task : stageo_run_task, seeds, workers(); mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "stageo_common.jl")])

function main()
    phase = get(CLI, "phase", "")
    if DRY
        phase in ENS_PHASES && (cs = cases(phase); println("$phase flights=$(length(cs)) (ensemble, outside the plan count)  horizon_days=$(cs[1].horizon_days)  runs=$(join([c.run_id for c in cs], ", "))"); return)
        phase in PHASES && return dry_run([phase]); phase == "" && return dry_run(PHASES); error("--phase must be one of $(PHASES) with --dry-run (or omitted for all)")
    end
    phase == "selftest" && return selftest()
    mkpath(SFO)
    phase in PHASES || phase in ENS_PHASES || error("--phase must be selftest|" * join(vcat(PHASES, ENS_PHASES), "|"))
    summary = joinpath(SFS, "summary_$(phase)_mars_$(HOST).csv"); isfile(summary) && error("$summary exists: use a distinct --host stamp per call")
    seeds = filter(c -> !ondisk(c), cases(phase)); println("$phase: $(length(seeds)) cases to fly ($(length(cases(phase)) - length(seeds)) already on disk)")
    isempty(seeds) && error("no cases (all runs already on disk)")
    finish(campaign(seeds, phase == "cal" ? 86400.0 : maximum(c.horizon_days for c in seeds) * 86400.0), phase)
end
main()
