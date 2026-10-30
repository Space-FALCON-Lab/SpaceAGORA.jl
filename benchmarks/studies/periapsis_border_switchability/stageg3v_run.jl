# Venus gate V-G(V) driver (plan R8.3, R8.4; config stageg3v_config.toml; the Mars gate stageg3_run.jl, never edited, with the Venus config view). One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=24G julia --project=. stageg3v_run.jl --phase=selftest|cal|smoke|gate --results-root=DIR [--workers=N] --host=NAME [--dry-run]
# --dry-run flies nothing (it loads the common file, so the arm codes and run ids come from the real code) and prints the case counts of the phase (or of every phase when --phase is omitted) against the plan's R8.3 counts, the flight list, and the cost estimate from the Stage L Venus flight timings ([stageg3v.cost]).
# NO GATE FLIGHT IS STARTED BY THIS FILE UNLESS --phase is cal or gate without --dry-run; the coordinator releases those phases (Freeze A(V) first). The smoke phase may be run after the build.
# Trees below <root>/stageg3v: cal/<gid>/{guess,nodes,calibration}.csv (phase cal), runs/<gid>_<code>_pert_<tag>/ (gate), smoke/, summary_<phase>_venus_<host>.csv (the Stage O columns plus geometry, role, keeper, epoch_utc,
# elements, decision_speed and the safety columns; a summary is never overwritten: use a DISTINCT --host stamp per call). Runs already on disk (passes.csv present) are skipped.
#   cal    3 fit-only calibration cases, one per geometry (pass 1 at the Odyssey band's target Q at that geometry and epoch)
#   smoke  forced safety-stop test (F4): a zero-burn flight at VG1 with the Q safety threshold lowered to smoke_forced_safety_q_w_cm2 (it must stop with safety_q), 3 d, plus one keeper padding flight (a one-case campaign plans route none)
#   gate   R8.3: per geometry 2 keeper (F-O) flights (seeds 101, 102); at the zero_burn_roles (VG1) also 2 zero-burn flights, paired GRAM seeds; 8 flights
using TOML
const CLI = let kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    kv
end
const DRY = haskey(CLI, "dry-run")
const CFG3 = joinpath(@__DIR__, "stageg3v_config.toml")
include(joinpath(@__DIR__, "stageg3v_common.jl"))
const META = stageg3v_meta(CFG3); const G3 = META["stageg3v"]   # geometry list from the hash-verified table
const ROOT = abspath(get(CLI, "results-root", ".")); const SFS = joinpath(ROOT, String(META["stagef"]["results_subdir"]))
const HOST = get(CLI, "host", "local")
const PHASES = ["cal", "smoke", "gate"]
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(G3["workers"])
geoms() = [String(g["id"]) for g in G3["geometry"]]
role_of(gid) = String(stageg3v_geometry(META, gid)["role"])
seed_tag(k) = String(G3["seed_tag_prefix"]) * string(k)

"""One flight seed. keeper "fixed" (F-O) or "zero_burn"; k the GRAM seed index; `role` is "flight" or "pad" (route padding, not counted)."""
function mkseed(gid, keeper, k, dir; horizon=nothing, suffix="", role="flight", smoke_force=false)
    code = stageg3_arm_code(keeper, Float64(G3["keeper_margin_m"])); tg = seed_tag(k); name = keeper == "fixed" ? "F-O" : "ZB"
    return (; config_path=CFG3, gid, grole=role_of(gid), role, keeper, tag=tg, arm=keeper == "fixed" ? "FO" : "ZB", code, run_id="$(gid)_$(code)_pert_$(tg)$(suffix)", out=dir, cal_root=joinpath(SFS, "cal"),
        horizon_days=horizon === nothing ? Float64(G3["horizon_days"]) : Float64(horizon), dt_atm_s=NaN, smoke_force, decision_speed=String(G3["decision_speed"]), name)
end
# stageg3v_run_task reads seed.role as the geometry role; the case table's own flight/pad flag is seed.role in stageo's convention, so the two are kept under distinct names and mapped here
task_seed(c) = merge(c, (; role=c.grole))

"""Cases of a phase (flights; the cal phase returns calibration seeds)."""
function cases(phase)
    out = NamedTuple[]; seeds = Int.(G3["seeds"]); sd = Float64(G3["smoke_days"])
    if phase == "cal"
        for gid in geoms(); push!(out, (; config_path=CFG3, gid, gram_tag=String(G3["cal_gram_tag"]), out=joinpath(SFS, "cal"), run_id="cal_$gid", role="cal", grole=role_of(gid))); end
    elseif phase == "smoke"
        g1 = first(filter(g -> String(g["role"]) == "VG1", G3["geometry"]))["id"]
        push!(out, mkseed(String(g1), "zero_burn", seeds[1], joinpath(SFS, "smoke"); horizon=sd, suffix="_smoke", smoke_force=true))
        push!(out, mkseed(String(g1), "fixed", seeds[1], joinpath(SFS, "smoke"); horizon=sd, suffix="_pad", role="pad"))
    elseif phase == "gate"
        dir = joinpath(SFS, "runs")
        for gid in geoms()
            for k in first(seeds, Int(G3["keeper_seeds_per_geometry"])); push!(out, mkseed(gid, "fixed", k, dir)); end
            if role_of(gid) in String.(G3["zero_burn_roles"])
                for k in first(seeds, Int(G3["zero_burn_seeds"])); push!(out, mkseed(gid, "zero_burn", k, dir)); end
            end
        end
    else
        error("unknown phase $phase")
    end
    return out
end
expect(key) = Int(G3["plan_counts"][key])
flights(cs) = count(c -> c.role != "pad" && c.role != "cal", cs)
ondisk(c) = c.role == "cal" ? isfile(joinpath(c.out, c.gid, "guess.csv")) : isfile(joinpath(c.out, c.run_id, "passes.csv"))

"""Epoch line of a geometry for the dry run, with the GRAM epoch assertion: the epoch the flights will use (stageg3v_load) must equal the table epoch, and the GRAM core epoch the run guard will demand is printed."""
function gram_epoch_line(gid)
    _, sc, g, ep = stageg3v_load(CFG3, gid); used = String(sc["grid"]["epochs_utc"][1]); @assert used == ep == String(g["epoch_utc"]) == stageg3_epoch_string(String(G3["base_epoch_utc"]), Float64(g["epoch_off_d"]))
    return "$(ep) (+$(g["epoch_off_d"]) d) GRAM_core_epoch_required=\"$(stageg3_gram_prefix(used))\""
end

"""Cost estimate of the dry run from the Stage L Venus flight timings ([stageg3v.cost]); DERIVED: flights x model passes x median s/pass, plus calibration and a per-phase start-up. Contended, not a performance number."""
function cost_lines()
    c = G3["cost"]; np = maximum(Int(g["n_passes_model"]) for g in JSON.parsefile(joinpath(stageg3_task_dir(), String(G3["table"]["path"])))["geometries"]); sp = Float64(c["s_per_pass_median"]); spx = Float64(c["s_per_pass_max"])
    nf = Int(G3["plan_counts"]["gate"]); fl = nf * np * sp; flx = nf * np * spx * 2.0; su = Float64(c["phase_startup_s"]); nc = Int(G3["plan_counts"]["cal_cases"]); cal = nc * Float64(c["cal_s_median"])
    sm = 3 * sp * 2
    println("cost estimate (DERIVED; source: Stage L Venus flights, $(c["s_per_pass_median"]) s per pass median, $(c["s_per_pass_max"]) max, sf3, dt_max_atmosphere 1.0 s; $(np) model passes per flight; contended, not a performance number):")
    println("   gate: $(nf) flights x $(np) passes x $(sp) s = $(round(fl; digits=0)) s of flight work = $(round(fl / 3600; digits=3)) worker-h; upper bound (max s/pass, 2x passes) $(round(flx / 3600; digits=3)) worker-h; plus about $(round(su; digits=0)) s of worker start-up per phase")
    println("   cal: $(nc) cases x $(c["cal_s_median"]) s (median; Stage I Venus calibrations, max 32.5 s) + start-up; smoke: 2 flights of 3 d (3 passes) about $(round(sm; digits=0)) s + start-up")
    println("   total: about $(round((fl + cal + sm + 3 * su) / 3600; digits=3)) worker-h (range to $(round((flx + 3 * su + nc * 32.5) / 3600; digits=2)) worker-h), a few minutes of wall at $(workers()) workers")
end

function dry_run(which)
    pc = G3["plan_counts"]; ok = true
    for ph in which
        cs = cases(ph); n = flights(cs); pads = count(c -> c.role == "pad", cs); done = count(ondisk, cs)
        if ph == "cal"
            println(rpad("cal", 6), "cases=$(length(cs)) (plan: $(pc["cal_cases"]) calibration cases, one per geometry)  geometries=[$(join(geoms(), ", "))]", done > 0 ? "  [$done on disk]" : ""); ok &= length(cs) == expect("cal_cases"); continue
        end
        if ph == "smoke"
            println(rpad("smoke", 6), "flights=$n (+ $pads padding flight(s)); plan smoke cases=$(pc["smoke"]) (outside the gate's $(expect("gate")))  runs=[$(join([c.run_id for c in cs], ", "))]"); ok &= length(cs) == pc["smoke"]; continue
        end
        nk = count(c -> c.keeper == "fixed", cs); nz = count(c -> c.keeper == "zero_burn", cs)
        println(rpad(ph, 6), "flights=$n (plan: $(expect("gate"))): keeper F-O $nk (plan $(expect("keeper"))), zero-burn $nz (plan $(pc["zero_burn"])); horizon_days=$(join(sort(unique(getfield.(cs, :horizon_days))), "/"))")
        for gid in geoms()
            g = filter(c -> c.gid == gid, cs); println("   $(rpad(gid, 6)) role=$(rpad(role_of(gid), 7)) i=$(stageg3v_geometry(META, gid)["inc_deg"]) w0=$(stageg3v_geometry(META, gid)["omega_deg"]) O0=$(stageg3v_geometry(META, gid)["raan_deg"]) epoch=$(gram_epoch_line(gid)) keeper=$(count(c -> c.keeper == "fixed", g)) zero_burn=$(count(c -> c.keeper == "zero_burn", g))")
            foreach(c -> println("      $(c.run_id)  keeper=$(c.keeper) gram_seed_string=venus$(c.tag) horizon_d=$(c.horizon_days)"), g)
        end
        ok &= n == expect("gate") && nk == expect("keeper") && nz == pc["zero_burn"]
    end
    which == PHASES && cost_lines()
    println(ok ? "dry run: counts match plan R8.3 (gate $(expect("gate")) flights = $(expect("keeper")) keeper + $(pc["zero_burn"]) zero-burn; $(expect("cal_cases")) calibration cases; $(pc["smoke"])-case smoke), GRAM epochs follow the geometry epochs, geometry table hash verified ($(G3["table_hash"])), burn_response_ratio $(G3["burn_response_ratio"]); nothing flown" : "dry run: COUNT MISMATCH against plan R8.3; nothing flown"); ok || exit(1)
end

"""Exact checks without any flight: config and geometry table, epochs, arm codes and case tables, the U15 identity (rotation correction = omega r_p cos i along track) and V_rel against the closed form at the Venus rotation, the safety state machine."""
function selftest()
    pc = G3["plan_counts"]; gs = G3["geometry"]
    @assert [String(g["id"]) for g in gs] == ["VG1", "VG2", "VGref"] && [String(g["role"]) for g in gs] == ["VG1", "VG2", "VGref"] && [Float64(g["inc_deg"]) for g in gs] == [90.0, 63.4, 89.876]
    @assert all(String(g["epoch_utc"]) == "2014-05-19T04:00:00" && Float64(g["epoch_off_d"]) == 0.0 for g in gs) && G3["table_hash"] == G3["table"]["sha256"]
    @assert stageg3_gram_prefix("2014-05-19T04:00:00") == "2014-5-19 4:0:0" && stageg3_gram_guard((year=2014, month=5, day=19, hour=4, minute=0, second=0.0), "2014-05-19T04:00:00") == "2014-5-19 4:0:0.0"
    @assert (try; stageg3_gram_guard((year=2000, month=1, day=1, hour=0, minute=0, second=0.0), "2014-05-19T04:00:00"); false; catch; true; end)   # a wrong GRAM epoch (the 2000-01-01 defect) is an error
    bad = mktempdir() do d; f = joinpath(d, "bad.toml"); write(f, replace(read(CFG3, String), String(G3["table"]["sha256"]) => "0"^64)); try; stageg3v_meta(f); false; catch e; occursin("mismatch", sprint(showerror, e)); end; end
    @assert bad "a wrong table hash must be an error"
    @assert G3["burn_response_ratio"] == 1.0 && G3["decision_speed"] == "relative" && G3["band_name"] == "venus_carried" && G3["safety"]["q_stop_w_cm2"] == 0.65 && G3["safety"]["n_hloc"] == 3.0
    @assert stageg3_epoch_string("2014-05-19T04:00:00", 0.0) == "2014-05-19T04:00:00" && stageg3_epoch_string("2014-05-19T04:00:00", 90.0) == "2014-08-17T04:00:00"
    for g in gs
        meta, sc, geom, ep = stageg3v_load(CFG3, String(g["id"]); dt_atm_s=1.0)
        @assert sc["orbit"]["apoapsis_radius_m"] == 72651770.0 && sc["integration"]["dt_max_atmosphere_s"] == 1.0 && sc["band"]["q_high_w_cm2"] == 0.32 && sc["grid"]["epochs_utc"] == [ep] && sc["atmosphere_gram"]["planet_name"] == "venus" && sc["planet"]["name"] == "Venus"
        @assert ep == stageg3_epoch_string(String(G3["base_epoch_utc"]), Float64(g["epoch_off_d"])) && meta["stageg3v"]["table_start"]["ra_target_m"] == G3["start"]["ra_target_m"]
    end
    @assert stageg3_arm_code("fixed", 1000.0) == "FO_fixed1000" && stageg3_arm_code("zero_burn", 1000.0) == "ZB_zeroburn"
    for ph in PHASES[2:end]; cs = cases(ph); ids = [(c.out, c.run_id) for c in cs]; @assert length(unique(ids)) == length(ids) "duplicate run ids in $ph"; end
    gate = cases("gate"); @assert flights(gate) == expect("gate") && count(c -> c.keeper == "fixed", gate) == expect("keeper") && count(c -> c.keeper == "zero_burn", gate) == pc["zero_burn"]
    @assert all(c.keeper == "fixed" || c.grole == "VG1" for c in gate) && all(c.role == "flight" && !c.smoke_force for c in gate) && length(cases("cal")) == expect("cal_cases")
    sm = cases("smoke"); @assert length(sm) == pc["smoke"] && sm[1].smoke_force && sm[1].keeper == "zero_burn" && sm[1].gid == "VG1" && sm[2].role == "pad" && flights(sm) == 1
    # keeper and zero-burn flights of one geometry share the GRAM seed strings (paired perturbations)
    for gid in geoms(); @assert Set(c.tag for c in gate if c.gid == gid && c.keeper == "zero_burn") <= Set(c.tag for c in gate if c.gid == gid && c.keeper == "fixed"); end
    # U15 identities on synthetic osculating states at Venus: omega x r_p . t_hat = omega r_p cos i, and V_rel against the closed form including the cos(latitude) term
    mu = Float64(G3["rotation"]["mu_m3s2"]); w = Float64(G3["rotation"]["w_rot_rad_s"]); omega = [0.0, 0.0, w]; ra = Float64(G3["start"]["ra0_m"]); rp0 = 6051800.0 + 135150.0
    worst = 0.0; worst_v = 0.0; cnt = 0
    for inc in (39.2, 63.4, 89.876, 90.0, 116.6, 150.0), raan in (0.0, 90.0, 270.0), argp in (0.0, 45.0, 90.0, 200.0), nu in (0.0, 180.0)
        a = (ra + rp0) / 2; e = (ra - rp0) / (ra + rp0); i = deg2rad(inc); O = deg2rad(raan); om = deg2rad(argp); f = deg2rad(nu); p = a * (1 - e^2); rr = p / (1 + e * cos(f))
        rq = [rr * cos(f), rr * sin(f), 0.0]; vq = sqrt(mu / p) .* [-sin(f), e + cos(f), 0.0]
        R = [cos(O) -sin(O) 0; sin(O) cos(O) 0; 0 0 1] * [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)] * [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1]; r = R * rq; v = R * vq
        pv = periapsis_vectors(r, v, mu); @assert abs(pv.rp - rp0) < 1e-3 && abs(pv.vp - sqrt(mu * (2 / rp0 - 1 / a))) < 1e-6
        alt = along_track_rotation(r, v, mu, omega); worst = max(worst, abs(alt - w * rp0 * cos(i)))
        sphi = sin(i) * sin(om); closed = sqrt(pv.vp^2 - 2 * w * rp0 * pv.vp * cos(i) + (w * rp0)^2 * (1 - sphi^2)); worst_v = max(worst_v, abs(periapsis_vrel(r, v, mu, omega) - closed)); cnt += 1
    end
    @assert worst < 1e-6 "rotation correction differs from omega r_p cos i by $worst m/s"; @assert worst_v < 1e-6 "V_rel differs from the closed form by $worst_v m/s"
    # Venus rotation is retrograde and slow: V_rel / V_p = 1 - w r_p cos(i) / V_p is above 1 for prograde (63.4) and below 1 for retrograde (116.6) orbits, by under 3e-4 (omega r_p about 1.9 m/s against V_p about 10 km/s)
    vp0 = sqrt(mu * (2 / rp0 - 2 / (ra + rp0))); r63 = (vp0 - w * rp0 * cos(deg2rad(63.4))) / vp0; r117 = (vp0 - w * rp0 * cos(deg2rad(116.6))) / vp0
    @assert r63 > 1 && r117 < 1 && abs(r63 - 1) < 3e-4 && abs(r117 - 1) < 3e-4 "V_rel/V_p at Venus: $r63, $r117"
    # safety state machine: Q on V_rel above the threshold at the closing apoapsis; altitude stop below alt0 - n_hloc H; disabled never fires
    ss = SafetyStop(true, 0.65, 3.0, 6000.0); @assert !safety_periapsis!(ss, 100000.0, 4700.0) && ss.alt0_m == 100000.0 && !safety_periapsis!(ss, 82500.0, 4700.0) && safety_periapsis!(ss, 81999.0, 4700.0) && ss.reason == "safety_alt"
    ss = SafetyStop(true, 0.65, 3.0, 6000.0); safety_periapsis!(ss, 100000.0, 4700.0); rho_stop = 2 * 0.65e4 / 4700.0^3
    @assert !safety_apoapsis!(ss, 0.99 * rho_stop) && safety_apoapsis!(ss, 1.01 * rho_stop) && ss.reason == "safety_q"
    sd = SafetyStop(false, 0.65, 3.0, 6000.0); safety_periapsis!(sd, 100000.0, 4700.0); @assert !safety_apoapsis!(sd, 1.0) && !safety_periapsis!(sd, 1.0, 4700.0)
    @assert abs(qs_relative_w_cm2(rho_stop, 4700.0) - 0.65) < 1e-12
    println("stageg3v selftest ok: $(length(gs)) geometries, gate $(flights(gate)) flights ($(pc["keeper"]) keeper + $(pc["zero_burn"]) zero-burn), epochs, arm codes, paired seeds; rotation correction = omega r_p cos i along track (max error $(worst) m/s over $(cnt) states), V_rel closed form (max error $(worst_v) m/s), Venus V_rel/V_p signs, safety stop")
end

"""Smoke-phase check (F4): the forced zero-burn flight must have stopped on safety_q after at least one pass, with no burn, and the first-pass Q on V_rel must lie in [0.05, 0.54] W/cm2."""
function smoke_check(res)
    bad = false
    for s in res.successful
        v = s.value; f = joinpath(SFS, "smoke", v.run_id, "passes.csv"); isfile(f) || continue
        r = CSV.read(f, DataFrame)[1, :]; q = Float64(r.q_eff_w_cm2); println("smoke $(v.run_id): first-pass q_eff = $(round(q; digits=4)) W/cm2 (mean-density Q $(round(Float64(r.q_w_cm2); digits=4)), vrel $(round(Float64(r.vrel_mps); digits=1)) m/s)")
        0.05 <= q <= 0.54 || (bad = true; println("STOP RULE: first-pass Q outside [0.05, 0.54]"))
        if v.keeper == "zero_burn"
            v.stop_reason == "safety_q" || (bad = true; println("SMOKE FAIL: forced zero-burn flight stopped with $(v.stop_reason), not safety_q"))
            v.burns == 0 || (bad = true; println("SMOKE FAIL: zero-burn flight has $(v.burns) burns"))
        end
    end
    bad && error("smoke check failed")
end

function finish(res, phase)
    tag = phase * "_venus_" * HOST
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
    phase == "cal" ? foreach(r -> r.status == "ok" && println("  cal $(r.tag): seed_alt_m=$(round(r.seed_alt_m; digits=1)) q_first=$(round(r.q_first_w_cm2; digits=4)) target=$(r.q_target_w_cm2) vrel_first=$(round(r.vrel_first_mps; digits=1)) hloc=$(round(r.hloc_at_hg_m; digits=1)) core_epoch=$(r.gram_core_epoch)"), rows) :
    for r in rows; r.status == "ok" && println("  $(r.run_id): passes=$(r.passes) t_end_days=$(round(r.t_end_days; digits=2)) stop_reason=$(r.stop_reason) dv=$(round(r.total_dv_mps; digits=4)) burns=$(r.burns) core_epoch=$(r.gram_core_epoch)"); end
    phase == "smoke" && smoke_check(res)
    (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
end
campaign(seeds, mission_s) = run_gram_campaign(seeds[1] isa NamedTuple && haskey(seeds[1], :gram_tag) ? stageg3v_cal_task : stageg3v_run_task, [haskey(s, :gram_tag) ? s : task_seed(s) for s in seeds], workers();
    mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "stageg3v_common.jl")])

function main()
    phase = get(CLI, "phase", "")
    if DRY
        phase in PHASES && return dry_run([phase]); phase == "" && return dry_run(PHASES); error("--phase must be one of $(PHASES) with --dry-run (or omitted for all)")
    end
    phase == "selftest" && return selftest()
    mkpath(SFS)
    phase in PHASES || error("--phase must be selftest|" * join(PHASES, "|"))
    summary = joinpath(SFS, "summary_$(phase)_venus_$(HOST).csv"); isfile(summary) && error("$summary exists: use a distinct --host stamp per call")
    seeds = filter(c -> !ondisk(c), cases(phase)); println("$phase: $(length(seeds)) cases to fly ($(length(cases(phase)) - length(seeds)) already on disk)")
    isempty(seeds) && error("no cases (all runs already on disk)")
    finish(campaign(seeds, phase == "cal" ? 86400.0 : maximum(c.horizon_days for c in seeds) * 86400.0), phase)
end
main()
