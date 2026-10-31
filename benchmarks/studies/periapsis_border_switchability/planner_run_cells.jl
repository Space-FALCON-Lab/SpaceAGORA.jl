# planner_run_cells.jl: a COPY of planner_run.jl (regenerated after Revision 16) with two added options: --cells=ARM@GEOM[,ARM@GEOM] and --seedlist=k[,k] restricting the flights of a phase; used to split blocks across hosts and to screen
# remedies on chosen dev seeds. Everything else, including the config, is identical.
# Planner driver (WP-P2; config planner_config.toml). One run_monte_carlo campaign per call (process route asserted); memory capped:
#   systemd-run --user --scope -p MemoryMax=24G julia --project=. planner_run.jl --phase=selftest|cal|smoke|smoke_ol|repro|xhost|dev|held --results-root=DIR [--workers=N] --host=NAME [--hostrole=ws|sf3] [--blocks=a,b] [--dry-run]
# --dry-run flies nothing (it loads the common file, so the arm codes and run ids come from the real code) and prints the case counts of the phase (or of every phase when --phase is omitted) against the plan's 274.
# NO FLIGHT IS STARTED BY THIS FILE UNLESS --phase is cal, smoke, smoke_ol, repro, xhost, dev or held without --dry-run; the coordinator releases those phases (Freeze B first). NOT RUN: written without executing Julia;
# the first command to run is --phase=selftest (no GRAM, no flight), then --dry-run.
# The Stage O baselines F-O and N2-I are flown by stageo_run.jl (frozen) and the 10 gate flights by stageg3_run.jl (frozen); this driver lists them in --dry-run (so the total is 274) but flies only the planner arms.
# Trees below <root>/planner: cal/<gid>/{guess,nodes,calibration}.csv (phase cal; gid ODY and the G1 id), runs/ (held), dev/, smoke/, repro/, xhost/, summary_<phase>_mars_<host>.csv (a summary is never overwritten: use a
# DISTINCT --host stamp per call). Runs already on disk (passes.csv present) are skipped.
using TOML
const CLI = let kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)(?:=(.*))?$", a); m === nothing || (kv[m.captures[1]] = something(m.captures[2], "1")); end
    kv
end
const DRY = haskey(CLI, "dry-run")
const CFGP = joinpath(@__DIR__, "planner_config.toml")
include(joinpath(@__DIR__, "planner_common.jl"))
const META = planner_meta(CFGP); const PL = META["planner"]; const G3 = META["stageg3"]
const ROOT = abspath(get(CLI, "results-root", ".")); const SFS = joinpath(ROOT, String(META["stagef"]["results_subdir"]))
const HOST = get(CLI, "host", "local"); const HOSTROLE = get(CLI, "hostrole", "")
const PHASES = ["cal", "smoke", "smoke_ol", "repro", "xhost", "dev", "held", "sens", "devT", "heldT", "devR", "calG", "gen1", "genG", "rev17", "calx0", "rev20olt", "rev20val", "rev21c", "rev21v", "calR23", "rev23geo", "calR25", "rev25low", "mf25", "rev27ei", "rev27low", "rev27b", "rev26val", "rev26val2", "rev28t", "calR28", "rev29aim", "rev30seed", "rev31pos", "rev31node", "rev31g1230", "calR31", "rev31cls", "rev31ext", "rev32dev", "rev31ext2"]
const PHASE_BLOCKS = Dict("smoke" => ["smoke"], "smoke_ol" => ["smoke_ol"], "repro" => ["repro"], "xhost" => ["xhost"], "dev" => ["dev"], "sens" => ["sens"], "devT" => ["devT"], "heldT" => ["heldT"], "devR" => ["devR"], "rev17" => ["rev17"], "gen1" => ["gen1"], "genG" => ["genG"], "rev20olt" => ["rev20olt"], "rev20val" => ["rev20val"], "rev21c" => ["rev21c"], "rev21v" => ["rev21v"], "rev23geo" => ["rev23geo"], "rev25low" => ["rev25low"], "rev27ei" => ["rev27ei"], "rev27b" => ["rev27b"], "rev26val" => ["rev26val"], "rev26val2" => ["rev26val2"], "rev28t" => ["rev28t"], "rev29aim" => ["rev29aim"], "rev30seed" => ["rev30seed"], "rev31pos" => ["rev31pos"], "rev31node" => ["rev31node"], "rev31g1230" => ["rev31g1230"], "rev31cls" => ["rev31cls"], "rev31ext" => ["rev31ext"], "rev32dev" => ["rev32dev"], "rev31ext2" => ["rev31ext2"], "rev27low" => ["rev27low"],
    "held" => ["held_baselines", "held_geometry", "held_ladder", "held_ol", "held_olr"])
const PLANNER_ARMS = [k for (k, a) in PL["arms"] if String(a["runner"]) == "planner"]
workers() = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(PL["workers"])
seed_tag(k) = String(PL["seed_tag_prefix"]) * string(k)
seedset(name) = Int.(PL["seedsets"][name])

"""Config rows (flight rows, plus the extra rows for smoke_ol) of the blocks of a phase, optionally filtered by --hostrole and --blocks."""
function plan_rows(phase)
    blocks = PHASE_BLOCKS[phase]; haskey(CLI, "blocks") && (blocks = filter(b -> b in split(CLI["blocks"], ","), blocks))
    rows = phase == "smoke_ol" ? PL["extra"] : phase == "sens" ? PL["sensitivity"] : phase in ("devT", "heldT", "devR", "gen1", "genG", "rev17", "rev20olt", "rev20val", "rev21c", "rev21v", "rev23geo", "rev25low", "rev27ei", "rev27low", "rev27b", "rev26val", "rev26val2", "rev28t", "rev29aim", "rev30seed", "rev31pos", "rev31node", "rev31g1230", "rev31cls", "rev31ext", "rev32dev", "rev31ext2") ? PL["terminal_flight"] : PL["flight"]
    return [r for r in rows if String(r["block"]) in blocks]
end

"""One flight seed of a cell (ARM@GEOM). Output dir: runs (held), else the phase name; run id <gid>_<code>_pert_<tag><suffix>."""
function mkseed(cell, k, phase, host)
    parts = String.(split(cell, "@")); arm, tok = parts[1], parts[2]; disp = length(parts) >= 3 ? parts[3] : ""; gid = planner_gid(META, tok); code = planner_arm_code(arm); tg = seed_tag(k)
    suffix = !isempty(disp) ? "_$(disp)" : phase in ("held", "dev", "devT", "heldT", "devR", "gen1", "genG") ? "" : phase == "xhost" ? "_xhost_$(host)" : "_$(phase)"
    sm = phase in ("smoke", "smoke_ol")
    return (; config_path=CFGP, gid, grole=tok, role="flight", disp, arm, tag=tg, code, run_id="$(gid)_$(code)_pert_$(tg)$(suffix)", out=joinpath(SFS, phase == "held" ? "runs" : phase),
        cal_root=joinpath(SFS, "cal"), horizon_days=sm ? Float64(PL["smoke_days"]) : planner_horizon_days(META, gid), dt_atm_s=NaN, smoke_force=false, decision_speed=String(G3["decision_speed"]), host)
end

"""All flights of a phase: planner-arm cells only (baseline cells are listed by the dry run but flown by stageo_run.jl)."""
function cases(phase)
    out = NamedTuple[]
    if phase in ("cal", "calG", "calR23", "calR25", "calR28", "calR31")
        for tok in (phase == "cal" ? ("ODY", "G1") : phase == "calR23" ? String.(PL["rev23"]["cal_tokens"]) : phase == "calR25" ? String.(PL["rev25"]["cal_tokens"]) : phase == "calR28" ? String.(PL["rev28"]["cal_tokens"]) : phase == "calR31" ? String.(PL["rev31"]["cal_tokens"]) : ("GA", "GB", "GC", "GD")); gid = planner_gid(META, tok)
            push!(out, (; config_path=CFGP, gid, gram_tag=String(META["stageo"]["cal_gram_tag"]), out=joinpath(SFS, "cal"), run_id="cal_$gid", role="cal", grole=tok)); end
        return out
    end
    if phase == "calx0"                                           # Revision 20: S-free calibrations (first pass calibrated to exp(x0) of the pinned IC plan)
        for tok in String.(PL["ic"]["cal_tokens"]); gid = planner_gid(META, tok); ic = planner_ic_plan(META, planner_ic_token(gid, tok), "S-free")
            push!(out, (; config_path=CFGP, gid, gram_tag=String(META["stageo"]["cal_gram_tag"]), out=joinpath(SFS, "cal"), run_id="calx0_$gid", role="cal", grole=tok, cal_name=planner_ic_cal_name(gid, ic.x0), x0=ic.x0)); end
        return out
    end
    for r in plan_rows(phase), host in String.(r["hosts"])
        (HOSTROLE == "" || host == HOSTROLE) || continue
        for cell in String.(r["cells"]), k in (haskey(r, "cell_seeds") ? Int.(r["cell_seeds"][cell]) : seedset(String(r["seeds"])))      # Revision 20b: cell_seeds gives a cell its own seed list
            (!haskey(CLI, "cells") || cell in split(CLI["cells"], ",")) || continue                  # --cells=OLRT1@ODY,... selects cells (new in this driver)
            (!haskey(CLI, "seedlist") || k in parse.(Int, split(CLI["seedlist"], ","))) || continue         # --seedlist=1,2,3 selects seeds (new in this driver)
            arm = String(split(cell, "@")[1]); arm in PLANNER_ARMS && push!(out, mkseed(cell, k, phase, host))
        end
    end
    return out
end
flights(cs) = count(c -> c.role != "cal", cs)
ondisk(c) = c.role == "cal" ? isfile(joinpath(c.out, get(c, :cal_name, c.gid), "guess.csv")) : isfile(joinpath(c.out, c.run_id, "passes.csv"))

"""All 274 flights of the plan by block (config rows), planner or not."""
plan_total() = sum(Int(r["count"]) for r in PL["flight"])

function dry_run(which)
    pc = PL["plan_counts"]; ok = true
    for ph in which
        cs = try; cases(ph); catch err; ok = false; println(rpad(ph, 9), "ERROR: ", first(sprint(showerror, err), 300)); continue; end; n = flights(cs); done = count(ondisk, cs)
        if ph in ("cal", "calG", "calx0", "calR23", "calR25", "calR28", "calR31")
            println(rpad(ph, 9), "cases=$(length(cs)) (fit-only calibrations)  gids=[$(join([c.gid for c in cs], ", "))]", ph == "calx0" ? "  cal dirs=[$(join([c.cal_name for c in cs], ", "))]" : "", done > 0 ? "  [$done on disk]" : ""); ok &= length(cs) == (ph == "cal" ? 2 : ph == "calG" ? 4 : ph == "calR23" ? length(PL["rev23"]["cal_tokens"]) : ph == "calR25" ? length(PL["rev25"]["cal_tokens"]) : ph == "calR28" ? length(PL["rev28"]["cal_tokens"]) : ph == "calR31" ? length(PL["rev31"]["cal_tokens"]) : length(PL["ic"]["cal_tokens"])); continue
        end
        arms = sort(unique(c.arm for c in cs)); println(rpad(ph, 9), "planner flights=$n  arms=$(arms)  hosts=$(sort(unique(c.host for c in cs)))  horizon_days=$(join(sort(unique(getfield.(cs, :horizon_days))), "/"))", done > 0 ? "  [$done on disk]" : "")
        ids = [(c.out, c.run_id) for c in cs]; ok &= length(unique(ids)) == length(ids)
    end
    allrows = [(String(r["block"]), String(c), String(h), k) for r in PL["flight"] for c in r["cells"] for h in r["hosts"] for k in seedset(String(r["seeds"]))]
    nbase = count(x -> x[1] != "gate" && (startswith(x[2], "FO@") || startswith(x[2], "N2I@")), allrows); ngate = count(x -> x[1] == "gate", allrows)
    println("plan: $(length(allrows)) flights in all (plan $(pc["total"])): gate $(ngate) (stageg3_run.jl), Stage O baselines $(nbase) (stageo_run.jl), planner runner $(length(allrows) - ngate - nbase)")
    ok &= length(allrows) == Int(pc["total"]) && ngate == Int(pc["gate"])
    HOSTROLE == "" && (ok &= sum(flights(cases(p)) for p in ("smoke", "repro", "xhost", "dev", "held")) == length(allrows) - ngate - nbase)
    println("sensitivity (penalty 0): $(length(cases("sens"))) planner flights; total with sensitivity $(Int(pc["total_with_sens"]))")
    ok &= length(cases("sens")) == Int(pc["sens"])
    println(ok ? "dry run: counts match plan R4.5 amended by R7.3 / R8.2 ($(pc["total"]) flights; $(pc["total_with_sens"]) with the sensitivity arms); nothing flown" : "dry run: COUNT MISMATCH against the plan; not flying")
    return ok
end

"""Exact checks without any flight (no GRAM): the lever, the DP against the Python values, the plan inputs, the OL schedule, the B-ODY rule, the table loader, the geometry injection, the safety stop with n_hloc = Inf, the flight list."""
function selftest()
    mu = 42828.375e9; approx(a, b; tol=1e-9) = abs(a - b) <= tol * max(1.0, abs(b))
    # lever and weight
    a, e = 1.6e7, 0.78; @assert approx(planner_lever_m_per_mps(mu, a, e), 4 * a^2 * sqrt(mu * (1 - e) / ((1 + e) * a)) / mu); @assert approx(planner_weight(6450.0, 1.0, mu, a, e), 6450.0 / planner_lever_m_per_mps(mu, a, e))
    @assert approx(planner_weight(6450.0, 0.5, mu, a, e), 2 * planner_weight(6450.0, 1.0, mu, a, e))
    # the DP: values from planner_lib.cross_language_vector (Python)
    level = [0.00, 0.30, 0.55, 0.70, 0.40, 0.10, -0.20, -0.35]; w = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4]
    pl = planner_plan(level, fill(-0.5, 8), fill(0.5, 8), w; eps=1.0e-6); @assert pl !== nothing
    @assert all(abs.(pl.nu .- [-0.2, -0.2, -0.2, -0.2, -0.15, -0.15, -0.15, -0.15]) .< 1e-9) "DP offsets differ from planner_lib.py: $(pl.nu)"
    @assert approx(pl.dv_cost, 0.54; tol=1e-9)
    z = planner_plan(zeros(5), fill(-1.0, 5), fill(1.0, 5), ones(5)); @assert z.dv_cost == 0.0 && all(z.nu .== 0.0)
    z = planner_plan([0.0, 1.5, 1.5], fill(-1.0, 3), fill(1.0, 3), ones(3)); @assert approx(z.dv_cost, 0.5) && z.nu[1] == 0.0
    z = planner_plan([0.0, 0.0, 1.5], fill(-1.0, 3), fill(1.0, 3), [1.0, 2.0, 4.0]); @assert z.du[1] != 0.0 && approx(z.dv_cost, 0.5)      # rising weights: prepay
    lvl = [0.0, 0.6, 1.2, 1.8]; a0 = planner_plan(lvl, fill(-1.0, 4), fill(0.5, 4), ones(4)); a1 = planner_plan(lvl, fill(-1.0, 4), fill(0.5, 4), ones(4); pen=0.3)      # burn penalty: same dV, one burn instead of three
    @assert count(x -> abs(x) > 1e-12, a0.du) == 3 && count(x -> abs(x) > 1e-12, a1.du) == 1 && approx(a1.dv_cost, 1.3) && approx(a0.dv_cost, 1.3)
    @assert planner_dp([1.0], [0.5], [1.0]) === nothing                                                                                       # A > B is infeasible
    # constants and plan inputs (same numbers as planner_lib.selftest)
    P, B = planner_constants(META, 1.0e4, 6450.0); hi0 = log(P.q_cap_mean_w_cm2)
    n = 10; tab = (; n, q_req=0.1623, W=0.72, d_slow=fill(0.05, 9), fast=zeros(n), w=[1.0 + 0.1 * k for k in 0:(n - 1)], path=[0.05 * k for k in 0:(n - 1)], slow=[0.05 * k for k in 0:(n - 1)], H=fill(6450.0, n),
        a=fill(1.6e7, n), e=fill(0.7, n), t=[1.0e5 * (k + 1) for k in 0:(n - 1)])
    swm(k) = planner_switches(META, k)
    i0 = planner_inputs(swm("L0"), tab, P, 2, -1.8, 2.0); @assert i0.level == [-1.8] && approx(i0.hi[1], hi0) && approx(i0.lo[1], hi0 - P.l0_width_ln)
    i1 = planner_inputs(swm("L1"), tab, P, 2, -1.8, 2.0); @assert approx(i1.lo[1], hi0 - 0.72)
    i3 = planner_inputs((; slack=true, ff=true, early=true, horizon=4), tab, P, 2, -1.8, 2.0)
    @assert length(i3.level) == 4 && approx(i3.level[3], -1.8 + 0.1) && approx(i3.hi[2], hi0 - P.fast_reserve_ln / 2) && i3.w[2:4] == [tab.w[3], tab.w[4], tab.w[5]]
    i2 = planner_inputs((; slack=true, ff=true, early=false, horizon=4), tab, P, 2, -1.8, 2.0); @assert i2.w == fill(2.0, 4)
    ib = planner_inputs((; slack=true, ff=true, early=true, horizon=4), tab, P, 2, -1.8, 2.0; bias=0.01); @assert approx(ib.level[2], -1.8 + 0.05 + 0.01)
    function dv_of(sw, j, z0, w0; bias=0.0)
        inp = planner_inputs(sw, tab, P, j, z0, w0; bias=bias); pl = planner_plan(inp.level, inp.lo, inp.hi, inp.w; eps=P.dp_tie_eps, pen=P.dp_burn_penalty_mps); dv = -pl.nu[1] * w0; return abs(dv) < P.dp_min_burn_mps ? 0.0 : dv
    end
    @assert approx(dv_of(swm("L1"), 3, hi0 + 0.1, 2.0), 0.2) && approx(dv_of(swm("L1"), 3, hi0 - 0.72 - 0.1, 2.0), -0.2) && dv_of(swm("L1"), 3, hi0 - 0.3, 2.0) == 0.0     # reset to the nearest edge
    # ramp table: OL schedule and L3 / L2 decisions (values from planner_lib.py)
    nr = 40; ramp = (; n=nr, q_req=0.1623, W=0.72, d_slow=fill(0.04, nr - 1), fast=zeros(nr), w=[1.0 + 0.05 * k for k in 0:(nr - 1)], path=[0.04 * k for k in 0:(nr - 1)], slow=[0.04 * k for k in 0:(nr - 1)],
        H=fill(6450.0, nr), a=fill(1.6e7, nr), e=fill(0.7, nr), t=[1.0e5 * (k + 1) for k in 0:(nr - 1)])
    sched = planner_ol_schedule(ramp, P, hi0 - 0.3); @assert length(sched) == 10 && approx(sched[2], 0.348075; tol=1e-6) && approx(sched[5], 0.144; tol=1e-6) && approx(sum(abs, values(sched)), 2.236875; tol=1e-6) "OL schedule differs from planner_lib.py (burn penalty 0.01)"
    inp = planner_inputs((; slack=true, ff=true, early=true, horizon=30), ramp, P, 5, hi0 + 0.05, 2.0); plr = planner_plan(inp.level, inp.lo, inp.hi, inp.w; eps=P.dp_tie_eps); @assert approx(-plr.nu[1] * 2.0, 0.1; tol=1e-6)
    # B-ODY rule
    zin = log(0.25 / B.qs_over_qmean); @assert planner_body_dv(B, P, 10, 5.0e5, -Inf, zin, 2.0) == 0.0
    zhi = log(0.40 / B.qs_over_qmean); @assert approx(planner_body_dv(B, P, 10, 5.0e5, -Inf, zhi, 2.0), 2.0 * (zhi - log(0.24 / B.qs_over_qmean))) && planner_body_dv(B, P, 10, 5.0e5, 5.0e5 - 3600.0, zhi, 2.0) == 0.0
    @assert approx(planner_body_dv(B, P, 150, 5.0e5, -Inf, zhi, 2.0), 2.0 * (zhi - log(0.20 / B.qs_over_qmean)))
    # table loader (synthetic file, hash not checked when `file` is given) and geometry injection
    mktempdir() do d
        f = joinpath(d, "t.json"); write(f, JSON.json(Dict("geometries" => Dict("X" => Dict("n" => 3, "q_req" => 0.16, "W" => 0.7, "path" => [0.0, 0.1, 0.2], "slow" => [0.0, 0.1, 0.2], "d_slow" => [0.1, 0.1], "fast" => [0.0, 0.0, 0.0],
            "w" => [0.0, 1.0, 2.0], "H" => [6450.0, 6450.0, 6450.0], "a" => [1.6e7, 1.6e7, 1.6e7], "e" => [0.7, 0.7, 0.7], "t" => [1.0, 2.0, 3.0])))))
        t = planner_table(META, "X"; file=f); @assert t.n == 3 && t.W == 0.7 && t.d_slow == [0.1, 0.1] && planner_tget(t.w, 5) == 2.0 && planner_dslow(t, 5) == 0.0
        @assert (try; planner_table(META, "Y"; file=f); false; catch; true; end)
    end
    @assert planner_gid(META, "ODY") == "ODY" && planner_gid(META, "G1") == "G5863" && planner_gid(META, "G2") == "G2348" && planner_gid(META, "Gworst") == "G3503"
    meta, sc, g, ep = stageg3_load(CFGP, "ODY"); @assert ep == "2001-11-06T19:00:32" && g["inc_deg"] == 93.522 && sc["orbit"]["apoapsis_radius_m"] == 28533000.0 && sc["grid"]["epochs_utc"] == [ep]
    @assert META["planner"]["burn_response_ratio"] == G3["burn_response_ratio"] == 1.0 && String(META["stagef"]["results_subdir"]) == "planner"
    # safety stop with n_hloc = Inf never fires on altitude, still fires on Q
    ss = SafetyStop(true, 0.65, Inf, 6000.0); @assert !safety_periapsis!(ss, 100000.0, 4700.0) && !safety_periapsis!(ss, 1000.0, 4700.0) && safety_apoapsis!(ss, 1.01 * 2 * 0.65e4 / 4700.0^3)
    # Revision 14 terminal aim (same numbers as planner_lib.py)
    Pt, _ = planner_constants(META, 1.0e4, 6450.0, "OLRT"); ta = planner_terminal_aim(50.0e3, [2.0e5, 2.2e5, 1.8e5, 2.0e5], Pt)
    @assert ta !== nothing && ta.N == 2 && approx(ta.dstar, 31250.0) && approx(ta.lnq, log(31250.0 / 2.0e5)) && planner_terminal_aim(200.0e3, fill(2.0e5, 4), Pt) === nothing && planner_terminal_aim(50.0e3, [2.0e5, 2.0e5], Pt) === nothing
    @assert approx(planner_terminal_aim(10.0e3, fill(8.0e5, 4), Pt).lnq, log(Pt.q_cap_mean_w_cm2) + log(Pt.term_min_q_factor)) && planner_terminal_aim(100.0e3, fill(1.0e3, 4), Pt).lnq == log(Pt.q_cap_mean_w_cm2)
    @assert planner_switches(META, "OLRT").terminal && !planner_switches(META, "OLR").terminal && planner_switches(META, "OLRT").update && planner_median([3.0, 1.0, 2.0, 10.0]) == 2.5 && length(cases("devT")) == 16 && length(cases("heldT")) == 32 && length(cases("devR")) == 96 && length(cases("gen1")) == 20 && length(cases("genG")) == 2 && length(cases("calG")) == 4
    # Revision 15 remedies: switches default off everywhere except their own arm; lower_margin shifts only the lower tube edge
    @assert all(planner_switches(META, k).rev_deadband == 0 && planner_switches(META, k).dev_thr == 0 && planner_switches(META, k).smooth_beta == 0 && planner_switches(META, k).lower_margin == 0 && planner_switches(META, k).min_burn_arm == 0 for k in ("L3", "OLR", "OLRT"))
    @assert planner_switches(META, "OLRT1").rev_deadband == 0.25 && planner_switches(META, "OLRT1").rev_window == 4 && planner_switches(META, "OLRT2").min_burn_arm == 0.05 && planner_switches(META, "OLRT3").dev_thr == 0.15 &&
        planner_switches(META, "OLRT4").smooth_beta == 0.5 && planner_switches(META, "OLRT5").lower_margin == 0.30 && planner_constants(META, 1.0e4, 6450.0, "OLRT2")[1].dp_burn_penalty_mps == 0.05 && planner_switches(META, "OLRT1").terminal
    ia = planner_inputs(swm("L1"), tab, P, 2, -1.8, 2.0); ib = planner_inputs(swm("L1"), tab, P, 2, -1.8, 2.0; lower_margin=0.3); @assert approx(ib.lo[1] - ia.lo[1], 0.3) && ib.hi == ia.hi
    # arms and the flight list
    @assert Set(PLANNER_ARMS) == Set(["BODY", "L0", "L1", "L2", "L3", "OL", "OLR", "L3Z", "OLRZ", "OLRT", "OLRT1", "OLRT2", "OLRT3", "OLRT4", "OLRT5", "OLRT35", "BODYT", "N2IT", "N2I0", "OLT", "OLRT3-W", "OLRT3-WF", "OLRT3-W3", "OLRT3-WF3", "OLRT3-W3E", "OLRT3-W4", "OLRT3-W4L", "OLRT3-W3L", "OLRT3-W4LN", "OLRT3-W3T", "OLRT3-W3E5", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5", "OLRT3-W3E7"]) && [(planner_switches(META, k).slack, planner_switches(META, k).ff, planner_switches(META, k).early) for k in ("L0", "L1", "L2", "L3")] ==
        [(false, false, false), (true, false, false), (true, true, false), (true, true, true)]
    @assert plan_total() == Int(PL["plan_counts"]["total"]) == 274
    @assert sum(Int(r["count"]) for r in PL["sensitivity"]) == Int(PL["plan_counts"]["sens"]) == 32 && Int(PL["plan_counts"]["total_with_sens"]) == 306 && length(cases("sens")) == 32
    P0, _ = planner_constants(META, 1.0e4, 6450.0, "L3Z"); P1, _ = planner_constants(META, 1.0e4, 6450.0, "L3"); P2, _ = planner_constants(META, 1.0e4, 6450.0)
    @assert P0.dp_burn_penalty_mps == 0.0 && P1.dp_burn_penalty_mps == 0.01 && P2.dp_burn_penalty_mps == 0.01 && planner_switches(META, "L3Z") == planner_switches(META, "L3") && planner_switches(META, "OLRZ") == planner_switches(META, "OLR")
    # Revision 17b: per-level GRAM keywords and the 60-flight block
    @assert planner_disp_kwargs(META, "") == (;) && planner_disp_kwargs(META, "L0") == (;) && planner_disp_kwargs(META, "L1") == (; gram_perturbation_scales=(1.5, 1.0, 1.0, 1.0)) && planner_disp_kwargs(META, "L2") == (; mars_mgcm_dust_levels=(0.3385, 0.3, 1.0)) &&
        planner_disp_kwargs(META, "L3") == (; mars_mgcm_dust_levels=(0.2298, 0.3, 1.0)) && planner_disp_kwargs(META, "L4") == (; gram_perturbation_scales=(1.5, 1.0, 1.0, 1.0), mars_mgcm_dust_levels=(0.3385, 0.3, 1.0))
    c17 = cases("rev17"); @assert length(c17) == 60 && length(unique(c.run_id for c in c17)) == 60 && all(c.gid == "ODY" && c.disp in ("L0", "L1", "L2", "L3", "L4") && endswith(c.run_id, "_" * c.disp) for c in c17) && Set(c.arm for c in c17) == Set(["OL", "OLR", "OLRT3"]) && sort(unique(parse(Int, split(c.tag, "s")[end]) for c in c17)) == [13, 14, 15, 16]
    @assert all(!hasproperty(c, :disp) || c.disp == "" for c in cases("gen1")) && cases("gen1")[1].run_id == "G0975_PL_BODYT_pert_ODY_s13"
    # Revision 20b / 20 / 19 (WP5): the OLT arm, the IC and FB switches; every pre-existing arm keeps the old switches and decisions
    @assert all(!planner_switches(META, k).fb_term && planner_switches(META, k).ic_scenario == "ref" for k in keys(PL["arms"]) if String(PL["arms"][k]["runner"]) == "planner" && !(k in ("OLRT3-W", "OLRT3-WF", "OLRT3-W3", "OLRT3-WF3", "OLRT3-W3E", "OLRT3-W4", "OLRT3-W4L", "OLRT3-W3L", "OLRT3-W4LN", "OLRT3-W3T", "OLRT3-W3E5", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5", "OLRT3-W3E7"))) && "OLT" in PLANNER_ARMS
    sO = planner_switches(META, "OL"); sT = planner_switches(META, "OLT"); @assert sT.terminal && !sO.terminal && (sT.slack, sT.ff, sT.early, sT.update, sT.horizon) == (sO.slack, sO.ff, sO.early, sO.update, sO.horizon) && String(PL["arms"]["OLT"]["mode"]) == "ol" && String(PL["arms"]["OLT"]["safety"]) == String(PL["arms"]["OL"]["safety"])
    @assert planner_constants(META, 1.0e4, 6450.0, "OLT")[1] == planner_constants(META, 1.0e4, 6450.0, "OL")[1]
    c20 = cases("rev20olt"); @assert length(c20) == 13 && length(unique(c.run_id for c in c20)) == 13 && all(c.gid == "ODY" && c.arm == "OLT" && endswith(c.run_id, "_" * c.disp) && startswith(c.run_id, "ODY_PL_OLT_pert_") for c in c20)
    sd(d) = sort([parse(Int, split(c.tag, "s")[end]) for c in c20 if c.disp == d]); @assert sd("L0") == [13, 14, 15, 16] && sd("L1") == [14, 15, 16] && sd("L2") == [14, 15, 16] && sd("L3") == [13, 15, 16]
    # linear cost in the DP: nothing = zeros = -0.0 are bit-identical, a nonzero entry moves the plan (values from planner_lib.cross_language_vector)
    lvx = [0.00, 0.30, 0.55, 0.70, 0.40, 0.10, -0.20, -0.35]; wx = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4]; Ax = fill(-0.5, 8) .- lvx; Bx = fill(0.5, 8) .- lvx; Lx = planner_tie_weights(wx, 1.0e-6)
    d_none = planner_dp(Ax, Bx, Lx); @assert planner_dp(Ax, Bx, Lx; lin=zeros(8)) == d_none && planner_dp(Ax, Bx, Lx; lin=fill(-0.0, 8)) == d_none && planner_dp(Ax, Bx, Lx; pen=0.3, lin=zeros(8)) == planner_dp(Ax, Bx, Lx; pen=0.3)
    p_lin = planner_plan(lvx, fill(-0.5, 8), fill(0.5, 8), wx; lin=[0, 0, 0, 0, 0, 1.5, 0, 0]); @assert all(abs.(p_lin.nu .- [-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.15, -0.15]) .< 1e-9) && approx(p_lin.dv_cost, 0.56; tol=1e-9)
    # first-burn slope identity (c_T / H = 1 / lever = dFB/drp at fixed desired periapsis) and the plan-input hook
    for rp_ in (3.480e6, 3.484e6, 3.490e6)
        fd = (planner_first_burn_dv(mu, 7.408e6, rp_ - 250.0, rp_ + 150e3) - planner_first_burn_dv(mu, 7.408e6, rp_ + 250.0, rp_ + 150e3)) / 500.0; @assert approx(planner_first_burn_slope(mu, 7.408e6, rp_, 6450.0) / 6450.0, fd; tol=1e-6)
    end
    tabx = (; tab..., a=fill(1.0e7, n), e=fill(0.65, n), H=fill(6450.0, n)); fbx = (; weight=2.0); lvz = planner_fb_lin(tabx, fbx, mu, 7.408e6, 2, 5, 9)
    @assert planner_fb_lin(tabx, fbx, mu, 7.408e6, 2, 1, 5) === nothing && planner_fb_lin(tabx, fbx, mu, 7.408e6, 2, 8, 5) === nothing && count(!iszero, lvz) == 1 && lvz[4] > 0 && planner_fb_lin(tabx, (; weight=0.0), mu, 7.408e6, 2, 5, 9) == zeros(9)
    ioff = planner_inputs((; slack=true, ff=true, early=true, horizon=4), tabx, P, 2, -1.8, 2.0); ion = planner_inputs((; slack=true, ff=true, early=true, horizon=4), tabx, P, 2, -1.8, 2.0; fb=fbx, fb_N=4, mu=mu)
    @assert !hasproperty(ioff, :lin) && ion.level == ioff.level && ion.lo == ioff.lo && ion.hi == ioff.hi && ion.w == ioff.w && length(ion.lin) == 4 && count(!iszero, ion.lin) == 1 && ion.lin[3] > 0
    # the Odyssey and G1 ground tables: the term off (and weight 0) reproduces the existing schedule exactly; the term on gives N 216 -> 231 and 106 -> 140 burns at Odyssey (planner_lib.py / rev20_data/ic/headline.json crosscheck), G1 unchanged
    P0x = (; P..., dp_burn_penalty_mps=0.0); muA = Float64(META["planner"]["analysis"]["mu_mars_m3s2"]); kapx = META["planner"]["terminal_periapsis"]["kappa"]; q1x = log(0.212)
    for (gid, kap, N0, N1, b0, b1) in (("ODY", kapx["ODY"], 216, 231, 106, 140), ("G5863", kapx["G5863"], 265, 265, 4, 4))
        gt = planner_table(META, gid); fbp = (; mu=muA, kappa=Float64(kap), weight=1.0)
        s_old = planner_ol_schedule(gt, P0x, q1x); o_none = planner_ol_plan(gt, P0x, q1x); o_zero = planner_ol_plan(gt, P0x, q1x; fb=(; fbp..., weight=0.0)); o_on = planner_ol_plan(gt, P0x, q1x; fb=fbp)
        @assert o_none.sched == s_old && o_zero.sched == s_old && o_zero.nu == o_none.nu && planner_ol_schedule(gt, P, q1x) == planner_ol_plan(gt, P, q1x).sched
        @assert o_zero.N == N0 && o_on.N == N1 && length(s_old) == b0 && length(o_on.sched) == b1 && !o_on.dropped && o_on.converged "FB term: N $(o_zero.N) -> $(o_on.N), burns $(length(s_old)) -> $(length(o_on.sched))"
    end
    # Revision 20 initial-condition plans from the pinned dump: ref -> nothing; S-free / S-fixed read x0, nu_pre, w_pre (no number retyped); the signed pre-burn raises the periapsis when nu_pre < 0; the file is hash-pinned
    @assert planner_ic_plan(META, "ODY", "ref") === nothing
    icF = planner_ic_plan(META, "ODY", "S-free"); icX = planner_ic_plan(META, "ODY", "S-fixed"); icG = planner_ic_plan(META, "G1", "S-fixed")
    @assert icF.x0 == icX.x0 && icF.nu_pre == icX.nu_pre && icF.dv_pre_charged == 0.0 && approx(icX.dv_pre_charged, abs(icX.nu_pre) * icX.w_pre) && icX.dv_pre_signed == -icX.nu_pre * icX.w_pre && (icX.nu_pre < 0) == (icX.dv_pre_signed > 0) && icG.nu_pre > 0 && icG.dv_pre_signed < 0
    @assert abs(icF.x0 - (-1.7707)) < 1e-3 && abs(icX.dv_pre_charged - 0.1043) < 1e-3 && abs(icG.x0 - (-1.5199)) < 1e-3 && planner_ic_cal_name("ODY", icF.x0) == "ODY__x0_$(round(icF.x0; digits=5))" && planner_ic_token("G5863", "G1") == "G1" && planner_ic_token("G0975", "GC") == "G0975"
    @assert (try; planner_ic_plan(META, "ODY", "S-free"; file=joinpath(mktempdir(), "none.json")); false; catch; true; end)
    # Revision 20c (WP5 Part B): guard v2 and the validation block; default off everywhere except the validation arm
    @assert all(!planner_switches(META, k).ttg_guard for k in keys(PL["arms"]) if String(PL["arms"][k]["runner"]) == "planner" && !(k in ("OLRT3-W", "OLRT3-WF", "OLRT3-W3", "OLRT3-WF3", "OLRT3-W3E", "OLRT3-W4", "OLRT3-W4L", "OLRT3-W3L", "OLRT3-W4LN", "OLRT3-W3T", "OLRT3-W3E5", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5", "OLRT3-W3E7"))) && planner_switches(META, "OLRT3-W").ttg_guard && planner_switches(META, "OLRT3-W").ic_scenario == "S-fixed" && planner_switches(META, "OLRT3-W").terminal
    s3 = planner_switches(META, "OLRT3"); sW = planner_switches(META, "OLRT3-W"); @assert Base.structdiff(s3, (; ttg_guard=false, ic_scenario="")) == Base.structdiff(sW, (; ttg_guard=false, ic_scenario="")) && planner_constants(META, 1.0e4, 6500.0, "OLRT3")[1] == planner_constants(META, 1.0e4, 6500.0, "OLRT3-W")[1]
    cv = cases("rev20val"); @assert length(cv) == 29 && length(unique(c.run_id for c in cv)) == 29 && count(c -> c.arm == "OLRT3" && c.gid == "ODY" && c.disp == "L4", cv) == 1 && count(c -> c.arm == "OLRT3-W" && c.gid == "ODY", cv) == 20 && count(c -> c.arm == "OLRT3-W" && c.gid == "G5863", cv) == 4 && count(c -> c.arm == "OLRT3" && c.gid == "G5863", cv) == 4 &&
        all(c.disp == "L0" for c in cv if c.gid == "G5863") && sort(unique(c.disp for c in cv if c.gid == "ODY")) == ["L0", "L1", "L2", "L3", "L4"] && sort(unique(parse(Int, split(c.tag, "s")[end]) for c in cv)) == [13, 14, 15, 16] && all(startswith(c.run_id, c.gid * "_PL_" * c.arm * "_pert_ODY_s") for c in cv)
    # Revision 21 (WP21-C): arm OLRT3-WF, blocks rev21c (36) and rev21v (28), generic-geometry dispersion level L2p, guard-set selector
    sF = planner_switches(META, "OLRT3-WF"); @assert sF.ttg_guard && sF.ic_scenario == "S-free" && Base.structdiff(sF, (; ic_scenario="")) == Base.structdiff(planner_switches(META, "OLRT3-W"), (; ic_scenario="")) && String(PL["arms"]["OLRT3-W"]["guard_set"]) == "v2" == String(PL["arms"]["OLRT3-WF"]["guard_set"])
    c21 = cases("rev21c"); @assert length(c21) == 36 && length(unique(c.run_id for c in c21)) == 36 && count(c -> c.gid == "G0975", c21) == 16 && count(c -> c.gid == "ODY", c21) == 8 && count(c -> c.gid == "G0063", c21) == 12 && all(c.arm in ("OLRT3", "BODYT", "N2IT") for c in c21) &&
        sort(unique(parse(Int, split(c.tag, "s")[end]) for c in c21 if c.gid == "ODY")) == [17, 18, 19, 20] && sort(unique(c.disp for c in c21 if c.gid == "G0975")) == ["L0", "L2p"] && all(c.disp == "L0" for c in c21 if c.gid == "G0063") && count(c -> c.arm == "OLRT3" && c.gid == "G0975" && c.disp == "L2p", c21) == 4
    v21 = cases("rev21v"); @assert length(v21) == 52 && length(unique(c.run_id for c in v21)) == 52 && all(c.arm in ("OLRT3-W3", "OLRT3-WF3") for c in v21) && count(c -> c.arm == "OLRT3-W3" && c.gid == "ODY", v21) == 28 && count(c -> c.arm == "OLRT3-WF3", v21) == 8 && count(c -> c.gid == "G0063", v21) == 4 && count(c -> c.gid == "G0975", v21) == 8 &&
        count(c -> c.arm == "OLRT3-W3" && c.gid == "G5863", v21) == 4 && count(c -> c.arm == "OLRT3-WF3" && c.gid == "G5863", v21) == 4 && count(c -> c.gid == "ODY" && c.disp in ("L0", "L2"), v21) == 20 && sort(unique(c.disp for c in v21 if c.gid == "ODY")) == ["L0", "L1", "L2", "L3", "L4"]
    # Revision 23b: geometry entries, mean-field tables and IC plans (pinned), block rev23geo (12), guard t* of the three geometries
    c23 = cases("rev23geo"); @assert length(c23) == 12 && length(unique(c.run_id for c in c23)) == 12 && all(c.arm == "OLRT3-W3" && c.disp == "L0" for c in c23) && sort(unique(c.gid for c in c23)) == ["G0063M", "G0528", "G0890"] && all(sort([parse(Int, split(c.tag, "s")[end]) for c in c23 if c.gid == g]) == [13, 14, 15, 16] for g in ("G0528", "G0890", "G0063M"))
    cc23 = cases("calR23"); @assert sort([c.gid for c in cc23]) == ["G0063M", "G0528", "G0890"]
    @assert planner_gid(META, "GE") == "G0528" && planner_gid(META, "GF") == "G0890" && planner_gid(META, "GAM") == "G0063M" && planner_gid(META, "GA") == "G0063" && planner_gid(META, "GC") == "G0975"
    muT = Float64(META["planner"]["analysis"]["mu_mars_m3s2"]); PT = planner_constants(META, 1.0e4, 6500.0, "OLRT3")[1]; pc = JSON.parsefile(joinpath(stageg3_task_dir(), "rev23_data", "prep_check.json"))["tstar"]
    for (gid, ref) in (("G0528", pc["G0528"]["guard_tstar_d"]), ("G0890", pc["G0890"]["guard_tstar_d"]), ("G0063M", 65.48))
        tb = planner_table(META, gid); icx = planner_ic_plan(META, gid, "S-fixed"); tg = planner_ttg_setup(META, gid, tb, PT, muT, log(0.211971) + icx.nu_pre)          # x0 as the Python prep check: planner_rev20_replay.tstar_for with the q_first of G0975 (0.211971) registered for the new tokens
        @assert tb.n > 100 && icx.nu_pre == 0.0 && abs(tg.tstar_nominal_d - ref) < (gid == "G0063M" ? 0.02 : 1.0e-6) && tg.kappa == 1.0 && tg.tstar_s == tg.tstar_nominal_d * 86400.0 "$gid t* $(tg.tstar_nominal_d) vs $ref"
    end
    @assert planner_table(META, "G0975").q_req != planner_table(META, "G0975M").q_req                                                              # the preview and the mean-field G0975 tables are different, both load
    # Revision 25 (WP25-J): the apoapsis target as a parameter, one-sided tolerance (above the target), event, node check, IC variables, geometry/target entries, block rev25low; all default-off
    gtq = planner_table(META, "ODY"); T0 = planner_target_table(META); P0r, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3"); @assert T0["ra_target_m"] == Float64(META["stageo"]["ra_target_m"]) "default target radius"
    @assert Float64(META["stageo"]["ra_target_m"]) == P0r.ra_target_m && P0r.ra_target_m == P0r.ra_stop_m "default stop radius"
    @assert P0r.tail_passes == Int(PL["tail_passes"]) && T0["tail_passes"] == P0r.tail_passes "default tail"
    @assert !P0r.apo_one_sided && P0r.term_aim_below_m == 12500.0 "default aim below"
    @assert isequal(planner_constants(META, 1.0e4, 6500.0, "OLRT3"; target="default")[1], P0r) "default target by name"
    @assert String(T0["event"]) == "apoapsis" "event"
    Pt2, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3"; target="T2000"); tj = JSON.parsefile(joinpath(stageg3_task_dir(), "target_orbits.json")); trad = Dict(String(t["target_id"]) => 1000.0 * Float64(t["circle_r_km"]) for t in tj["targets"])
    @assert Pt2.ra_target_m == 5396190.0 && Pt2.ra_stop_m == Pt2.ra_target_m && Pt2.tail_passes == 50 && Pt2.term_step_m == P0r.term_step_m && planner_horizon_days(META, "G0486") == Float64(PL["horizon_days"]) && planner_horizon_days(META, "ODY") == Float64(PL["horizon_days"]) &&
        all(planner_target_table(META, k)["ra_target_m"] == trad[String(planner_target_table(META, k)["target_id"])] for k in ("T2000", "T1000", "T500")) && [planner_geometry_target(META, g) for g in ("G0486", "G1230", "G1038", "ODY", "G5863")] == ["T2000", "T1000", "T500", nothing, nothing]
    @assert planner_gid(META, "GG") == "G0486" && planner_gid(META, "GH") == "G1230" && planner_gid(META, "GI") == "G1038"
    for (g, k) in (("G0486", "T2000"), ("G1230", "T1000"), ("G1038", "T500")); e = only(o for (_, o) in PL["geometry"] if o["id"] == g); @assert e["target"] == k; end
    m2 = deepcopy(META); m2["planner"]["target"]["event"] = "periapsis"; @assert (try; planner_constants(m2, 1.0e4, 6500.0, "OLRT3"); false; catch err; occursin("Not implemented", sprint(showerror, err)); end)
    m3 = deepcopy(META); m3["planner"]["target"]["apo_one_sided"] = true; P1s, _ = planner_constants(m3, 1.0e4, 6500.0, "OLRT3"); @assert P1s.apo_one_sided && P1s.ra_stop_m == P1s.ra_target_m + 25000.0 && P1s.term_aim_below_m == -12500.0 && P1s.ra_target_m == P0r.ra_target_m
    ta_def = planner_terminal_aim(60.0e3, fill(2.0e5, 4), P0r); ta_one = planner_terminal_aim(60.0e3, fill(2.0e5, 4), P1s); @assert ta_def.N == 2 && approx(ta_def.dstar, (60.0e3 + 12500.0) / 2) && approx(ta_one.dstar, (60.0e3 - 12500.0) / 2) && planner_terminal_aim(10.0e3, fill(2.0e5, 4), P1s) === nothing
    @assert planner_apo_slope(0.056, 40000.0) == 0.056 * 40.0 && planner_term_lin(gtq, (; weight=0.0), P0r, 42828.375e9, P0r.ra_target_m, 2, 100, 150) == zeros(150) && sum(planner_term_lin(gtq, (; weight=0.0), P1s, 42828.375e9, P1s.ra_target_m, 2, 100, 150)) > 0
    @assert planner_ol_plan(gtq, P0r, log(0.212)).sched == planner_ol_schedule(gtq, P0r, log(0.212)) && planner_ol_plan(gtq, P0r, log(0.212); fb=(; mu=42828.375e9, kappa=0.99, weight=0.0)).sched == planner_ol_schedule(gtq, P0r, log(0.212))
    # node check: pure functions, default off, no tolerance: an error when enabled without one
    @assert planner_wrap_deg(190.0) == -170.0 && planner_wrap_deg(-190.0) == 170.0 && planner_wrap_deg(180.0) == 180.0 && planner_wrap_deg(-180.0) == 180.0 && planner_wrap_deg(720.0) == 0.0
    rt = planner_j2_node_rate(42828.375e9, 0.001958616126805748, 3394200.0, 1.0e7, 0.65, deg2rad(93.5)); @assert rt > 0 && planner_j2_node_rate(42828.375e9, 0.001958616126805748, 3394200.0, 1.0e7, 0.65, deg2rad(25.2)) < 0
    @assert approx(planner_node_error_deg(10.0, 1.0e-5, 0.0, 1.0e6, 20.0), 0.0) && approx(planner_node_aim(10.0, 1.0e-5, 0.0, 20.0, 5.0e5), 1.0e6) && approx(planner_node_aim(10.0, 1.0e-5, 0.0, 20.0, 3.0e7), 1.0e6 + 360.0 / 1.0e-5 * round((3.0e7 - 1.0e6) / (360.0 / 1.0e-5))) && planner_node_aim(10.0, 0.0, 0.0, 20.0, 5.0) == 5.0
    @assert !planner_guard_constants(META).node_check && !planner_guard_constants(META; set="v3").node_check && isnan(planner_guard_constants(META).node_tol_deg)
    m4 = deepcopy(META); m4["planner"]["guard"]["node_check"] = true; @assert (try; planner_guard_constants(m4); false; catch err; occursin("no tolerance", sprint(showerror, err)); end)
    m4["planner"]["guard"]["node_tol_deg"] = 1.0; Gn = planner_guard_constants(m4); @assert Gn.node_check && Gn.node_tol_deg == 1.0
    @assert (try; planner_ttg_decide(Gn, gtq, P0r, 42828.375e9, 1.0, 40, fill(-1.0, 30), fill(-0.5, 30), fill(0.5, 30), fill(1.0, 30), zeros(30), 1.2e7, 3.49e6, 2.0e6, 5.0e6, nothing, 0.0, 6500.0, nothing, NaN); false; catch err; occursin("no node state", sprint(showerror, err)); end)
    # IC variables (default off): no insertion node / t* in the Revision 23 plans; a t* selector finds nothing in them
    icq = planner_ic_plan(META, "ODY", "S-fixed"); @assert icq.insertion_node_deg === nothing && icq.t_star_d === nothing && !Bool(get(PL["ic"], "use_insertion_node", false)) && !Bool(get(PL["ic"], "use_plan_t_star", false))
    m5 = deepcopy(META); m5["planner"]["ic"]["t_star_d"] = 70.0; @assert (try; planner_ic_plan(m5, "ODY", "S-fixed"); false; catch err; occursin("0 entries", sprint(showerror, err)); end)
    t_nom = planner_ttg_setup(META, "G0528", planner_table(META, "G0528"), P0r, 42828.375e9, log(0.212)); t_ovr = planner_ttg_setup(META, "G0528", planner_table(META, "G0528"), P0r, 42828.375e9, log(0.212); tstar_d=70.0); @assert t_ovr.tstar_nominal_d == 70.0 && t_ovr.tstar_s == 70.0 * 86400.0 && t_nom.tstar_nominal_d != 70.0
    c25 = cases("rev25low"); @assert length(c25) == 13 && length(unique(c.run_id for c in c25)) == 13 && count(c -> c.arm == "OLRT3-W3" && c.disp == "L0", c25) == 12 && sort(unique(c.gid for c in c25 if c.arm == "OLRT3-W3")) == ["G0486", "G1038", "G1230"] && count(c -> c.arm == "OLRT3" && c.gid == "ODY" && c.disp == "L4" && endswith(c.tag, "s14"), c25) == 1
    # Revision 25 pins: tables (with ra_target_m equal to the target radii), IC plans, guard t* of WP25-P (rev25_data/guard_tstar_rev25mf.json)
    gtr = JSON.parsefile(joinpath(stageg3_task_dir(), String(PL["ground_table_rev25"]["path"])))["geometries"]; gts = JSON.parsefile(joinpath(stageg3_task_dir(), "rev25_data", "guard_tstar_rev25mf.json"))
    for (gid, tk) in (("G0486", "T2000"), ("G1230", "T1000"), ("G1038", "T500"))
        @assert Float64(gtr[gid]["ra_target_m"]) == planner_target_table(META, tk)["ra_target_m"] && String(gtr[gid]["target_id"]) == String(planner_target_table(META, tk)["target_id"]) && planner_table(META, gid).n == Int(gtr[gid]["n"])
        icr = planner_ic_plan(META, gid, "S-fixed"); Pg, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3"; target=tk); tg = planner_ttg_setup(META, gid, planner_table(META, gid), Pg, 42828.375e9, log(0.211971) + icr.nu_pre)
        @assert Pg.ra_target_m == Float64(gtr[gid]["ra_target_m"]) && abs(tg.tstar_nominal_d - gts[gid]["guard_tstar_d"]) < 1.0e-6 && tg.kappa == 1.0 && tg.tstar_s == tg.tstar_nominal_d * 86400.0 "$gid t* $(tg.tstar_nominal_d) vs $(gts[gid]["guard_tstar_d"])"
    end
    @assert planner_ic_plan(META, "G0528", "S-fixed") !== nothing && planner_ic_plan(META, "ODY", "S-fixed").x0 == planner_ic_plan(META, "ODY", "S-fixed").x0
    # Revision 27: entry interface as a study parameter (default untouched), arm OLRT3-W3E, blocks rev27ei (4) and rev27low (13)
    @assert all(!haskey(PL["arms"][k], "ei_km") for k in keys(PL["arms"]) if !(k in ("OLRT3-W3E", "OLRT3-W3T", "OLRT3-W3E5", "OLRT3-W3E7"))) && PL["arms"]["OLRT3-W3E"]["ei_km"] == 150.0 && planner_switches(META, "OLRT3-W3E") == planner_switches(META, "OLRT3-W3") && String(PL["arms"]["OLRT3-W3E"]["guard_set"]) == "v3"
    @assert planner_ei_km(META, PL["arms"]["OLRT3-W3"], nothing) === nothing && planner_ei_km(META, PL["arms"]["OLRT3-W3E"], "T2000") == 150.0 && planner_ei_km(META, PL["arms"]["OLRT3"], "T500") === nothing
    m6 = deepcopy(META); m6["planner"]["targets"]["T500"]["ei_km"] = 300.0; @assert planner_ei_km(m6, PL["arms"]["OLRT3-W3"], "T500") == 300.0 && planner_ei_km(m6, PL["arms"]["OLRT3-W3E"], "T500") == 150.0 && planner_ei_km(m6, PL["arms"]["OLRT3-W3"], "T2000") === nothing
    c27 = cases("rev27ei"); @assert length(c27) == 4 && all(c.arm == "OLRT3-W3E" && c.gid == "ODY" && c.disp == "L0" for c in c27) && sort([parse(Int, split(c.tag, "s")[end]) for c in c27]) == [13, 14, 15, 16]
    c27l = cases("rev27low"); @assert length(c27l) == 13 && length(unique(c.run_id for c in c27l)) == 13 && count(c -> c.arm == "OLRT3-W3E" && c.disp == "L0", c27l) == 12 && sort(unique(c.gid for c in c27l if c.arm == "OLRT3-W3E")) == ["G0486", "G1038", "G1230"] && count(c -> c.arm == "OLRT3" && c.gid == "ODY" && c.disp == "L4" && endswith(c.tag, "s14"), c27l) == 1
    cb = cases("rev27b"); @assert length(cb) == 2 && Set(c.arm for c in cb) == Set(["OLRT3-W3", "OLRT3-W3E"]) && all(c.gid == "ODY" && c.disp == "mean" && endswith(c.run_id, "_s13_mean") for c in cb)
    @assert planner_disp_perturbed(META, "") && planner_disp_perturbed(META, "L0") && planner_disp_perturbed(META, "L2") && planner_disp_perturbed(META, "L2p") && !planner_disp_perturbed(META, "mean") && planner_disp_kwargs(META, "mean"; gid="ODY") == (;)
    # Revision 27 / user 2026-10-30: EI 150 is the default for NEW flights only (arms with ei_default = "new"); every existing block resolves as before (2000, or its explicit ei_km)
    @assert PL["ei_km_new_default"] == 150.0 && all(String(get(a, "ei_default", "")) == "" for (k, a) in PL["arms"] if !(k in ("OLRT3-W4", "OLRT3-W4L", "OLRT3-W3L", "OLRT3-W4LN", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5")))
    for rows in (PL["flight"], PL["terminal_flight"], PL["extra"], PL["sensitivity"]), r in rows, cell in String.(r["cells"])
        pp = String.(split(cell, "@")); haskey(PL["arms"], pp[1]) || continue; eik = planner_ei_km(META, PL["arms"][pp[1]], planner_geometry_target(META, planner_gid(META, pp[2])))
        @assert (eik === nothing) == !(pp[1] in ("OLRT3-W3E", "OLRT3-W3T", "OLRT3-W3E5", "OLRT3-W3E7", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5", "OLRT3-W4", "OLRT3-W4L", "OLRT3-W3L", "OLRT3-W4LN")) && (eik === nothing || eik == 150.0) "block $(r["block"]) cell $cell resolves to EI $eik"
    end
    m7 = deepcopy(META); m7["planner"]["arms"]["OLRT3-W3"]["ei_default"] = "new"; @assert planner_ei_km(m7, m7["planner"]["arms"]["OLRT3-W3"], nothing) == 150.0 && planner_ei_km(m7, m7["planner"]["arms"]["OLRT3"], nothing) === nothing
    m7["planner"]["arms"]["OLRT3-W3"]["ei_km"] = 300.0; @assert planner_ei_km(m7, m7["planner"]["arms"]["OLRT3-W3"], nothing) == 300.0          # an explicit ei_km wins
    # Revision 26 (WP26-B): guard v4 (position stage), L1, default t*, node check ON (1.0 deg) in the new arms only; all default off
    NEW26 = ("OLRT3-W4", "OLRT3-W4L", "OLRT3-W3L", "OLRT3-W4LN", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5")
    @assert all(String(get(PL["arms"][k], "ei_default", "")) == (k in NEW26 ? "new" : "") for k in keys(PL["arms"])) && all(planner_ei_km(META, PL["arms"][k], nothing) == 150.0 for k in NEW26)
    sw26 = Dict(k => planner_switches(META, k) for k in NEW26); @assert all(sw26[k].ttg_guard && sw26[k].ic_scenario == "S-fixed" && sw26[k].terminal && !sw26[k].default_tstar for k in NEW26) && sw26["OLRT3-W4L"].l1 && sw26["OLRT3-W3L"].l1 && sw26["OLRT3-W4LN"].l1 && !sw26["OLRT3-W4"].l1
    @assert all(!planner_switches(META, k).l1 && !planner_switches(META, k).default_tstar for k in keys(PL["arms"]) if String(PL["arms"][k]["runner"]) == "planner" && !(k in NEW26) && k != "OLRT3-W3T")
    G4 = planner_guard_constants(META; set="v4"); G3n = planner_guard_constants(META; set="v3n"); G3p = planner_guard_constants(META; set="v3")
    @assert G4.pos_stage && !G4.node_check && G4.node_tol_deg == 1.0 && !G3n.node_check && G3n.node_tol_deg == 1.0 && !G3n.pos_stage && !G3p.pos_stage && !G3p.node_check && !planner_guard_constants(META).node_check && G3p.pos === nothing &&
        G4.pos.pos_band_deg == 0.3 && G4.pos.kC_deg_per_d == 0.03 && G4.pos.M_frac == 0.5 && G4.pos.min_steps == 12 && G4.pos.shift_scale == 0.3 && G4.pos.fd_step == 0.3 && G4.pos.hw == 0.15 && G4.pos.nit == 3 && G4.pos.max_solves == 14 &&
        G4.s_frac == G3p.s_frac && G4.term_model == G3p.term_model && G4.t_offset_d == G3p.t_offset_d && G4.pen_g == G3p.pen_g
    c26b = cases("rev26val2"); @assert length(c26b) == 17 && [c.run_id for c in c26b] == [c.run_id for c in cases("rev26val")] && all(endswith(c.out, "rev26val2") for c in c26b)
    c26 = cases("rev26val"); @assert length(c26) == 17 && length(unique(c.run_id for c in c26)) == 17 && count(c -> c.arm == "OLRT3-W4L" && c.gid == "ODY", c26) == 8 && count(c -> c.arm == "OLRT3-W3L" && c.gid == "G0528" && c.disp == "L0", c26) == 4 && count(c -> c.arm == "OLRT3-W4L" && c.gid == "G5863", c26) == 4 &&
        count(c -> c.arm == "OLRT3" && c.gid == "ODY" && c.disp == "L4" && endswith(c.tag, "s14"), c26) == 1 && sort(unique(c.disp for c in c26 if c.arm == "OLRT3-W4L" && c.gid == "ODY")) == ["L0", "L2"]
    # pure functions: geometry, J2 track, L1, default t*
    u1 = planner_apo_hat(0.7, 1.1, 1.6); u2 = planner_apo_hat(0.8, 1.15, 1.6); @assert approx(norm(u1), 1.0) && approx(norm(planner_logmap(u1, u2)), planner_ang(u1, u2)) && abs(dot(planner_tproj(u1, u2), u1)) < 1e-12 && norm(planner_logmap(u1, u1)) == 0.0
    muJ = 42828.375e9; J2k = 0.001958616126805748; Rk = 3394200.0; tt_ = [1.0e5 * k for k in 1:20]; trk = planner_j2_track(tt_, fill(1.6e7, 20), fill(0.7, 20), fill(1.0, 20), 0.5, 1.0, 1.0, muJ, J2k, Rk)
    @assert trk.n == 20 && approx(planner_j2_at(trk, 0.0)[1], 0.5) && planner_j2_at(trk, 1.0e5 * 20)[1] > planner_j2_at(trk, 1.0e5 * 10)[1] && approx(norm(planner_j2_apo(trk, 5.0e5)), 1.0) && norm(planner_j2_rate_vec(trk, 5.0e5)) > 0
    E_ = SVector(0.01, 0.02, 0.0); s_ = SVector(1.0, 0.0, 0.0); l1a = planner_l1(E_, s_, 1.0, 0.0, 0.5); @assert approx(l1a.along_deg, rad2deg(0.01)) && 0.0 <= l1a.tau_d <= 0.5 && l1a.err_deg <= rad2deg(norm(E_)) + 1e-12 && planner_l1(E_, s_, 1.0, 5.0, 0.5).tau_d == 0.0 && planner_l1(E_, s_, 1.0, 5.0, 0.5).err_deg == rad2deg(norm(E_))
    gO = planner_table(META, "ODY"); PO, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3-W4L"); kO = META["planner"]["terminal_periapsis"]["kappa"]["ODY"]; rg = planner_arrival_range(gO, PO, 42828.375e9, kO, log(0.212)); dg = planner_default_tstar(gO, PO, 42828.375e9, kO, log(0.212), 6.0)
    @assert rg.fast_d < planner_ttg_setup(META, "ODY", gO, PO, 42828.375e9, log(0.212)).tstar_nominal_d < rg.slow_d && dg.tstar_d == clamp(rg.fast_d + 6.0, rg.fast_d, rg.slow_d) && META["planner"]["tstar_margin_d"] == 6.0
    nomO = planner_nominal_track(gO, PO, 42828.375e9, kO, log(0.212), deg2rad(93.522), deg2rad(109.7454), deg2rad(28.1517), J2k, Rk); @assert nomO.days == planner_arrival(gO, 42828.375e9, PO.ra_stop_m, vcat(log(0.212), planner_ol_plan(gO, PO, log(0.212)).level .+ planner_ol_plan(gO, PO, log(0.212)).nu), log(0.212), kO).days && nomO.N > 100
    # Revision 26c: node check v2 (J2-track forecast, aim shift clamped to the window), default off; arm OLRT3-W4LN
    G4n = planner_guard_constants(META; set="v4n"); @assert G4n.node_check_v2 && !G4n.node_check && G4n.node_tol_deg == 1.0 && G4n.pos_stage && !G4.node_check_v2 && !G3n.node_check_v2 && !planner_guard_constants(META).node_check_v2 && !planner_guard_constants(META; set="v3").node_check_v2
    m8 = deepcopy(META); m8["planner"]["guard_sets"]["v4n"]["node_check"] = true; @assert (try; planner_guard_constants(m8; set="v4n"); false; catch err; occursin("mutually exclusive", sprint(showerror, err)); end)
    m8b = deepcopy(META); delete!(m8b["planner"]["guard_sets"]["v4n"], "node_tol_deg"); @assert (try; planner_guard_constants(m8b; set="v4n"); false; catch err; occursin("no tolerance", sprint(showerror, err)); end)
    @assert planner_node_v2_aim(0.5, 1e-6, 1e6, 1e6, 3e4, 1.0) == (; aim=1e6, ok=true) && planner_node_v2_aim(5.0, 1e-6, 1e6, 1e6, 3e4, 1.0) == (; aim=1e6 - 3e4, ok=false) && planner_node_v2_aim(-5.0, 1e-6, 1e6, 1e6, 3e4, 1.0) == (; aim=1e6 + 3e4, ok=false) && planner_node_v2_aim(2.0, 0.0, 1e6, 1e6, 3e4, 1.0) == (; aim=1e6, ok=false)
    # growing-rate synthetic orbit: the J2 node rate grows as the orbit shrinks, so the linear extrapolation of the CURRENT rate (the rev26val defect) under-predicts the node motion that the forecast sees
    tg = [8.0e4 * k for k in 1:80]; ag = [1.6e7 * (1 - 0.0075 * k) for k in 1:80]; trg = planner_j2_track(tg, ag, fill(0.7, 80), fill(1.2, 80), 0.3, 0.2, 1.2, muJ, J2k, Rk)
    err_f, rate_f = planner_node_forecast(trg, tg[end], 0.0); lin_rate = rad2deg(planner_j2_node_rate(muJ, J2k, Rk, ag[10], 0.7, 1.2)); Om_f = rad2deg(planner_j2_at(trg, tg[end])[2]); Om_lin = rad2deg(planner_j2_at(trg, tg[10])[2]) + lin_rate * (tg[end] - tg[10])
    @assert abs(rate_f) > 2.0 * abs(lin_rate) && abs(planner_wrap_deg(Om_f - Om_lin)) > 0.05 && approx(err_f, planner_wrap_deg(Om_f)) "growing-rate forecast: rate_f $(rate_f * 86400) deg/d vs current $(lin_rate * 86400)"
    # Revision 28A: arm OLRT3-W3T (W3E + default t*), geometry G0975M (role GCM), blocks rev28t (37) and calR28; the default t* each geometry will get
    sT = planner_switches(META, "OLRT3-W3T"); sE = planner_switches(META, "OLRT3-W3E"); @assert sT.default_tstar && !sE.default_tstar && Base.structdiff(sT, (; default_tstar=false)) == Base.structdiff(sE, (; default_tstar=false)) && PL["arms"]["OLRT3-W3T"]["ei_km"] == 150.0 && String(PL["arms"]["OLRT3-W3T"]["guard_set"]) == "v3" && !haskey(PL["arms"]["OLRT3-W3T"], "tstar_d") && PL["tstar_margin_d"] == 6.0
    c28 = cases("rev28t"); @assert length(c28) == 37 && length(unique(c.run_id for c in c28)) == 37 && count(c -> c.arm == "OLRT3-W3T", c28) == 24 && count(c -> c.arm == "OLRT3-W3E", c28) == 12 && count(c -> c.arm == "OLRT3", c28) == 1 && count(c -> c.gid == "G0975M", c28) == 8 && count(c -> c.gid == "G1038", c28) == 4 && count(c -> c.gid == "G1230", c28) == 12 && count(c -> c.gid == "ODY", c28) == 13
    @assert planner_gid(META, "GCM") == "G0975M" && Set(cases("calR28")[k].gid for k in 1:1) == Set(["G0975M"]) && planner_table(META, "G0975M").n > 200 && planner_ic_plan(META, "G0975M", "S-fixed") !== nothing && planner_fb_params(META, "G0975M", 42828.375e9).kappa == planner_fb_params(META, "G0975", 42828.375e9).kappa
    println("default t* (model feasible arrival range, kappa, chosen = fast end + $(PL["tstar_margin_d"]) d), S-fixed anchor ln 0.211971 + nu_pre:")
    for (gid, tk) in (("ODY", nothing), ("G0975M", nothing), ("G1038", "T500"), ("G1230", "T1000"))
        tbl = planner_table(META, gid); Pq, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3T"; target=tk); tok = gid == "ODY" ? "ODY" : gid; icq = planner_ic_plan(META, tok, "S-fixed"); kq = planner_fb_params(META, gid, 42828.375e9).kappa; x0q = log(0.211971) + icq.nu_pre
        dq = planner_default_tstar(tbl, Pq, 42828.375e9, kq, x0q, PL["tstar_margin_d"]); nq = planner_ttg_setup(META, gid, tbl, Pq, 42828.375e9, x0q).tstar_nominal_d
        println("  ", rpad(gid, 7), " fast ", round(dq.fast_d; digits=3), " d  slow ", round(dq.slow_d; digits=3), " d  chosen ", round(dq.tstar_d; digits=3), " d  (nominal ", round(nq; digits=3), " d, kappa ", round(kq; digits=4), ", x0 ", round(x0q; digits=4), ")")
        @assert dq.fast_d < dq.tstar_d <= dq.slow_d && dq.tstar_d == clamp(dq.fast_d + PL["tstar_margin_d"], dq.fast_d, dq.slow_d) && dq.fast_d < nq
    end
    dO = planner_default_tstar(planner_table(META, "ODY"), planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3T")[1], 42828.375e9, planner_fb_params(META, "ODY", 42828.375e9).kappa, log(0.211971), 6.0); @assert abs(dO.tstar_d - 53.9) < 0.1
    # Revision 28a: term_aim_below_m 5,000 m for arms flagged aim_default = "new" (none yet); every existing cell keeps 12,500 m
    @assert PL["terminal"]["term_aim_below_m_new"] == 5000.0 && all(!haskey(a, "aim_default") for (k, a) in PL["arms"] if !(k in ("OLRT3-W3E5", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5", "OLRT3-W3E7")))
    for rows in (PL["flight"], PL["terminal_flight"], PL["extra"], PL["sensitivity"]), r in rows, cell in String.(r["cells"])
        pp = String.(split(cell, "@")); haskey(PL["arms"], pp[1]) || continue; tkq = planner_geometry_target(META, planner_gid(META, pp[2]))
        @assert planner_constants(META, 1.0e4, 6500.0, pp[1]; target=tkq)[1].term_aim_below_m == (pp[1] == "OLRT3-W3E7" ? 7500.0 : pp[1] in ("OLRT3-W3E5", "OLRT3-W3L5", "OLRT3-W4L5", "OLRT3-W4LN5") ? 5000.0 : 12500.0) "block $(r["block"]) cell $cell aims $(planner_constants(META, 1.0e4, 6500.0, pp[1]; target=tkq)[1].term_aim_below_m)"
    end
    m9 = deepcopy(META); m9["planner"]["arms"]["OLRT3-W3"]["aim_default"] = "new"; @assert planner_constants(m9, 1.0e4, 6500.0, "OLRT3-W3")[1].term_aim_below_m == 5000.0 && planner_constants(m9, 1.0e4, 6500.0, "OLRT3")[1].term_aim_below_m == 12500.0 && planner_constants(m9, 1.0e4, 6500.0, "OLRT3-W3"; target="T2000")[1].term_aim_below_m == 5000.0
    m9["planner"]["target"]["apo_one_sided"] = true; @assert planner_constants(m9, 1.0e4, 6500.0, "OLRT3-W3")[1].term_aim_below_m == -12500.0          # the one-sided switch keeps its own aim
    # Revision 29: arm OLRT3-W3E5 = OLRT3-W3E with the 5 km terminal aim; block rev29aim (21)
    @assert planner_switches(META, "OLRT3-W3E5") == planner_switches(META, "OLRT3-W3E")
    aE = Dict{String, Any}(PL["arms"]["OLRT3-W3E"]); aE5 = Dict{String, Any}(PL["arms"]["OLRT3-W3E5"]); delete!(aE, "note"); delete!(aE5, "note"); delete!(aE5, "aim_default"); @assert aE == aE5 && PL["arms"]["OLRT3-W3E5"]["aim_default"] == "new" && !haskey(PL["arms"]["OLRT3-W3E"], "aim_default")
    @assert planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3E5")[1].term_aim_below_m == 5000.0 && planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3E")[1].term_aim_below_m == 12500.0 && planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3E5"; target="T500")[1].term_aim_below_m == 5000.0 && planner_ei_km(META, PL["arms"]["OLRT3-W3E5"], nothing) == 150.0
    @assert Base.structdiff(planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3E5")[1], (; term_aim_below_m=0.0)) == Base.structdiff(planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3E")[1], (; term_aim_below_m=0.0))
    c29 = cases("rev29aim"); @assert length(c29) == 21 && length(unique(c.run_id for c in c29)) == 21 && count(c -> c.arm == "OLRT3-W3E5" && c.gid == "ODY", c29) == 8 && count(c -> c.arm == "OLRT3-W3E5" && c.gid == "G1038", c29) == 8 && count(c -> c.arm == "OLRT3-W3E" && c.gid == "G1038", c29) == 4 && count(c -> c.arm == "OLRT3" && c.gid == "ODY" && endswith(c.tag, "s14") && c.disp == "L4", c29) == 1
    @assert all(planner_constants(META, 1.0e4, 6500.0, c.arm; target=planner_geometry_target(META, c.gid))[1].term_aim_below_m == (c.arm == "OLRT3-W3E5" ? 5000.0 : 12500.0) for c in c29)
    # Revision 30.3: arms OLRT3-W3L5 / OLRT3-W4L5 (parents W3L / W4L plus the 5 km aim), block rev30seed (25)
    for (par, ch) in (("OLRT3-W3L", "OLRT3-W3L5"), ("OLRT3-W4L", "OLRT3-W4L5"), ("OLRT3-W3E", "OLRT3-W3E5"))
        ap = Dict{String, Any}(PL["arms"][par]); ac = Dict{String, Any}(PL["arms"][ch]); delete!(ap, "note"); delete!(ac, "note"); @assert ac["aim_default"] == "new" && !haskey(ap, "aim_default"); delete!(ac, "aim_default"); @assert ap == ac "arm $ch differs from $par beyond the aim"
        @assert planner_switches(META, ch) == planner_switches(META, par) && planner_constants(META, 1.0e4, 6500.0, ch)[1].term_aim_below_m == 5000.0 && planner_constants(META, 1.0e4, 6500.0, par)[1].term_aim_below_m == 12500.0 &&
            Base.structdiff(planner_constants(META, 1.0e4, 6500.0, ch)[1], (; term_aim_below_m=0.0)) == Base.structdiff(planner_constants(META, 1.0e4, 6500.0, par)[1], (; term_aim_below_m=0.0)) && planner_ei_km(META, PL["arms"][ch], nothing) == planner_ei_km(META, PL["arms"][par], nothing)
    end
    c30 = cases("rev30seed"); @assert length(c30) == 25 && length(unique(c.run_id for c in c30)) == 25 && count(c -> c.arm == "OLRT3-W3L5" && c.gid == "G0528", c30) == 8 && count(c -> c.arm == "OLRT3-W3E5" && c.gid == "G1230", c30) == 8 && count(c -> c.arm == "OLRT3-W4L5" && c.gid == "ODY", c30) == 8 && count(c -> c.arm == "OLRT3" && c.disp == "L4" && endswith(c.tag, "s14"), c30) == 1
    @assert all(planner_constants(META, 1.0e4, 6500.0, c.arm; target=planner_geometry_target(META, c.gid))[1].term_aim_below_m == (c.arm == "OLRT3" ? 12500.0 : 5000.0) for c in c30)
    # Revision 31: arms OLRT3-W4LN5 (W4LN + 5 km aim) and OLRT3-W3E7 (W3E5 + the arm's own 7.5 km aim, ASSUMED), blocks rev31pos (45), rev31node (25), rev31g1230 (9)
    for (par, ch, extra) in (("OLRT3-W4LN", "OLRT3-W4LN5", "aim_default"), ("OLRT3-W3E5", "OLRT3-W3E7", "term_aim_below_m"))
        ap = Dict{String, Any}(PL["arms"][par]); ac = Dict{String, Any}(PL["arms"][ch]); delete!(ap, "note"); delete!(ac, "note"); @assert haskey(ac, extra) && (extra == "term_aim_below_m" || !haskey(ap, extra)); delete!(ac, extra); @assert ap == ac "arm $ch differs from $par beyond $extra"
        @assert planner_switches(META, ch) == planner_switches(META, par) && planner_ei_km(META, PL["arms"][ch], nothing) == 150.0 && Base.structdiff(planner_constants(META, 1.0e4, 6500.0, ch)[1], (; term_aim_below_m=0.0)) == Base.structdiff(planner_constants(META, 1.0e4, 6500.0, par)[1], (; term_aim_below_m=0.0))
    end
    @assert planner_constants(META, 1.0e4, 6500.0, "OLRT3-W4LN5")[1].term_aim_below_m == 5000.0 && planner_constants(META, 1.0e4, 6500.0, "OLRT3-W4LN")[1].term_aim_below_m == 12500.0 && planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3E7")[1].term_aim_below_m == 7500.0 && planner_constants(META, 1.0e4, 6500.0, "OLRT3-W3E7"; target="T1000")[1].term_aim_below_m == 7500.0
    m10 = deepcopy(META); m10["planner"]["targets"]["T1000"]["term_aim_below_m"] = 9000.0; @assert planner_constants(m10, 1.0e4, 6500.0, "OLRT3-W3E5"; target="T1000")[1].term_aim_below_m == 9000.0 && planner_constants(m10, 1.0e4, 6500.0, "OLRT3-W3E7"; target="T1000")[1].term_aim_below_m == 7500.0          # arm key > target > aim_default
    G4nq = planner_guard_constants(META; set=String(PL["arms"]["OLRT3-W4LN5"]["guard_set"])); @assert G4nq.node_check_v2 && G4nq.node_tol_deg == 1.0 && !G4nq.node_check && PL["arms"]["OLRT3-W4LN5"]["l1"] && PL["arms"]["OLRT3-W4LN5"]["ei_default"] == "new"
    cp = cases("rev31pos"); cn = cases("rev31node"); cg = cases("rev31g1230"); @assert length(cp) == 45 && length(cn) == 25 && length(cg) == 9 && length(unique(c.run_id for c in cp)) == 45 && length(unique(c.run_id for c in cn)) == 25
    @assert count(c -> c.arm == "OLRT3-W4L5", cp) == 24 && count(c -> c.arm == "OLRT3-W3E5", cp) == 20 && sort(unique(c.gid for c in cp if c.arm == "OLRT3-W4L5")) == sort(["G5863", "G0975M", "G0063M", "G0890", "G0486", "G1038"]) && sort(unique(c.gid for c in cp if c.arm == "OLRT3-W3E5")) == sort(["G5863", "G0975M", "G0063M", "G0890", "G0486"])
    @assert count(c -> c.arm == "OLRT3-W4LN5", cn) == 16 && count(c -> c.arm == "OLRT3-W4L5", cn) == 8 && all(c.gid in ("G5863", "G0975M") for c in cn if c.arm != "OLRT3") && count(c -> c.arm == "OLRT3-W3E7" && c.gid == "G1230", cg) == 8 && count(c -> c.arm == "OLRT3", cp) == 1 && count(c -> c.arm == "OLRT3", cn) == 1 && count(c -> c.arm == "OLRT3", cg) == 1
    @assert all(planner_constants(META, 1.0e4, 6500.0, c.arm; target=planner_geometry_target(META, c.gid))[1].term_aim_below_m == (c.arm == "OLRT3" ? 12500.0 : c.arm == "OLRT3-W3E7" ? 7500.0 : 5000.0) for c in vcat(cp, cn, cg))
    for gid in ("G5863", "G0975M", "G0063M", "G0890", "G0486", "G1038", "G1230"); @assert isfile(joinpath(stageg3_task_dir(), "planner_ground_tables_rev23.json")) && planner_table(META, gid).n > 100 && planner_ic_plan(META, gid == "G5863" ? "G1" : gid, "S-fixed") !== nothing; end
    cc25 = cases("calR25"); @assert sort([c.gid for c in cc25]) == ["G0486", "G1038", "G1230"]
    # Revision 31.3: fine-grid geometries GJ/GK/GL/GM, phase calR31, block rev31cls (17); the table/IC pins are checked once [planner.ground_table_rev31] exists
    cls = cases("rev31cls"); ccl = cases("calR31"); @assert length(cls) == 17 && length(unique(c.run_id for c in cls)) == 17 && count(c -> c.arm == "OLRT3-W4L5" && c.disp == "L0", cls) == 16 && count(c -> c.arm == "OLRT3", cls) == 1
    @assert sort([c.gid for c in ccl]) == ["G0042_A2000", "G0140_AODY", "G3161_A2000", "G3257_AODY"] && sort(unique(c.gid for c in cls if c.arm == "OLRT3-W4L5")) == ["G0042_A2000", "G0140_AODY", "G3161_A2000", "G3257_AODY"]
    @assert [planner_geometry_target(META, g) === nothing for g in ("G0140_AODY", "G3257_AODY")] == [true, true] && [planner_constants(META, 1.0e4, 6500.0, "OLRT3-W4L5"; target=planner_geometry_target(META, g))[1].ra_target_m for g in ("G0042_A2000", "G3161_A2000")] == [5396190.0, 5396190.0]
    cx = cases("rev31ext"); @assert length(cx) == 13 && length(unique(c.run_id for c in cx)) == 13 && count(c -> c.arm == "OLRT3-W4L5" && c.disp == "L0", cx) == 8 && count(c -> c.arm == "OLRT3-W3E5" && c.gid == "G0486", cx) == 4 && count(c -> c.arm == "OLRT3", cx) == 1 && sort(unique(c.gid for c in cx if c.arm == "OLRT3-W4L5")) == ["G0486", "G1038"]
    cv = cases("rev32dev"); @assert length(cv) == 13 && length(unique(c.run_id for c in cv)) == 13 && count(c -> c.arm == "OLRT3-W3E5" && c.disp == "L0", cv) == 12 && count(c -> c.arm == "OLRT3", cv) == 1 && sort(unique(c.gid for c in cv if c.arm == "OLRT3-W3E5")) == ["G0063M", "G0890", "G0975M"] && sort(unique(c.tag for c in cv if c.arm == "OLRT3-W3E5")) == sort(unique(seed_tag.(1:4)))
    c2 = cases("rev31ext2"); @assert length(c2) == 9 && length(unique(c.run_id for c in c2)) == 9 && count(c -> c.arm == "OLRT3-W4L5" && c.disp == "L0", c2) == 8 && count(c -> c.arm == "OLRT3", c2) == 1 && sort(unique(c.gid for c in c2 if c.arm == "OLRT3-W4L5")) == ["G0140_AODY", "G3257_AODY"]
    haskey(PL, "ground_table_rev31") && foreach(g -> (@assert planner_table(META, g).n > 100 && planner_ic_plan(META, g, "S-fixed") !== nothing), ("G0140_AODY", "G0042_A2000", "G3257_AODY", "G3161_A2000"))
    # guard v3: separate constant set (s_frac 0.25, term_model), v2 untouched; the terminal-aim emulation only acts through the set
    G3c = planner_guard_constants(META; set="v3"); G2c = planner_guard_constants(META); @assert G3c.s_frac == 0.25 && G3c.term_model && !G2c.term_model && G2c.s_frac == 0.5 && (G3c.pen_g, G3c.rel_slack, G3c.anchor_k, G3c.bias_steps, G3c.k_window, G3c.ema_beta, G3c.dev_thr_g, G3c.t_offset_d) == (G2c.pen_g, G2c.rel_slack, G2c.anchor_k, G2c.bias_steps, G2c.k_window, G2c.ema_beta, G2c.dev_thr_g, G2c.t_offset_d)
    @assert String(PL["arms"]["OLRT3-W3"]["guard_set"]) == "v3" == String(PL["arms"]["OLRT3-WF3"]["guard_set"]) && planner_switches(META, "OLRT3-W3") == planner_switches(META, "OLRT3-W") && planner_switches(META, "OLRT3-WF3") == planner_switches(META, "OLRT3-WF")
    @assert planner_disp_kwargs(META, "L2p"; gid="G0975") == (; mars_mgcm_dust_levels=(1.0042, 0.3, 1.0)) && planner_disp_kwargs(META, "L2"; gid="ODY") == planner_disp_kwargs(META, "L2") && planner_disp_kwargs(META, "L0"; gid="G0975") == (;) && (try; planner_disp_kwargs(META, "L2p"; gid="ODY"); false; catch; true; end)
    @assert isequal(planner_guard_constants(META; set="v2"), planner_guard_constants(META)) && (try; planner_guard_constants(META; set="v9"); false; catch; true; end)
    # G0975 kappa from the calibration file (MEASURED, WP21-B) and the guard's t* with it (planner_rev21_position.nominal("G0975").days, S-fixed IC, q_first 0.211971)
    @assert planner_fb_params(META, "G0975", 42828.375e9).kappa == JSON.parsefile(joinpath(stageg3_task_dir(), String(META["planner"]["terminal_periapsis"]["kappa_file"])))["G0975"]["kappa_med"] && planner_fb_params(META, "G0063", 42828.375e9).kappa == 1.0 && planner_fb_params(META, "ODY", 42828.375e9).kappa == META["planner"]["terminal_periapsis"]["kappa"]["ODY"]
    tG = planner_ttg_setup(META, "G0975", planner_table(META, "G0975"), planner_constants(META, 1.0e4, 6500.0, "OLRT3")[1], 42828.375e9, planner_anchor_ln(0.211971, planner_ic_plan(META, "G0975", "S-fixed"))); @assert abs(tG.tstar_nominal_d - 87.09) < 0.01 "G0975 t* $(tG.tstar_nominal_d)"
    # Revision 21: the repinned IC file keeps every ODY / G1 / G0975 entry of the Revision 20 file value for value (so the flown rev20val selections are unchanged) and adds G0063
    let old = JSON.parsefile(joinpath(stageg3_task_dir(), "rev20_data", "ic", "headline.json")), new = JSON.parsefile(joinpath(stageg3_task_dir(), String(META["planner"]["ic"]["path"])))
        ko = p -> (p["geom"], p["scenario"], p["fb_term"], p["horizon"], p["delta_ln"], p["h_des_m"], p["t_star_d"]); dn = Dict(ko(p) => p for p in new["plans"])
        @assert all(haskey(dn, ko(p)) && dn[ko(p)] == p for p in old["plans"]) && any(p["geom"] == "G0063" for p in new["plans"])
        for tok in ("ODY", "G1", "G0975"), sc in ("S-fixed", "S-free"); ic = planner_ic_plan(META, tok, sc); po = only(p for p in old["plans"] if p["geom"] == tok && p["scenario"] == sc && !p["fb_term"] && p["horizon"] == "arrival" && p["delta_ln"] == 0.0 && p["h_des_m"] == 250000.0 && p["t_star_d"] === nothing); @assert ic.x0 == po["x0"] && ic.nu_pre == po["nu_pre"] && ic.w_pre == po["w_pre"]; end
        @assert planner_ic_plan(META, "G0063", "S-fixed") !== nothing
    end
    Gc = planner_guard_constants(META); @assert Gc.pen_g == 0.1 && Gc.rel_slack == 0.15 && Gc.anchor_k == 15 && Gc.bias_steps == 60 && Gc.k_window == 20 && Gc.ema_beta == 0.5 && Gc.dev_thr_g == 0.05 && Gc.max_solves == 12 && Gc.t_offset_d["ODY"] == 0.15 && Gc.t_offset_d["G5863"] == 0.0
    # guard core on synthetic inputs: an in-window prediction leaves the planner's plan unchanged; a late / early prediction tightens the tube on the needed side only and never raises the upper edge; at most max_solves solves
    gt2 = planner_table(META, "ODY"); P3, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3"); mu2 = 42828.375e9; kap2 = META["planner"]["terminal_periapsis"]["kappa"]["ODY"]
    hi2 = log(P3.q_cap_mean_w_cm2); n2 = 60; lvl2 = fill(hi2 - 0.5, n2); lo2 = fill(hi2 - gt2.W + 0.0, n2); hi2v = fill(hi2 - 0.0, n2); A2 = lo2 .- lvl2; B2 = hi2v .- lvl2; L2 = planner_tie_weights(fill(0.5, n2), P3.dp_tie_eps)
    nu02 = zeros(n2); ra2 = 1.2e7; rp2 = 3.49e6; t2 = 20.0 * 86400.0
    tpred = planner_ttg_decide(Gc, gt2, P3, mu2, kap2, 40, lvl2, A2, B2, L2, nu02, ra2, rp2, t2, 1.0e12, nothing, 0.0, 6500.0, nothing, NaN).ta
    r_in = planner_ttg_decide(Gc, gt2, P3, mu2, kap2, 40, lvl2, A2, B2, L2, nu02, ra2, rp2, t2, tpred, nothing, 0.0, 6500.0, nothing, NaN); @assert !r_in.engaged && r_in.nu == nu02 && !isnan(r_in.ema) && r_in.ta == tpred
    r_late = planner_ttg_decide(Gc, gt2, P3, mu2, kap2, 40, lvl2, A2, B2, L2, nu02, ra2, rp2, t2, tpred - 40.0 * 86400.0, nothing, 0.0, 6500.0, nothing, NaN)
    @assert r_late.engaged && r_late.late && all(r_late.nu .<= B2 .+ 1e-12) && all(r_late.nu .>= A2 .- 1e-12) && r_late.solves <= Gc.max_solves && any(r_late.nu .!= 0.0)
    r_early = planner_ttg_decide(Gc, gt2, P3, mu2, kap2, 40, lvl2, A2, B2, L2, nu02, ra2, rp2, t2, tpred + 40.0 * 86400.0, nothing, 0.0, 6500.0, nothing, NaN); @assert r_early.engaged && !r_early.late && all(r_early.nu .<= 1e-12) && r_early.solves <= Gc.max_solves
    @assert planner_ttg_kappa_hat(zeros(3), zeros(3), gt2, 20) === nothing && planner_ttg_lvl_alt(Gc, gt2, 5, lvl2, zeros(4), zeros(4), 0.0, 0.0) === nothing
    # guard policy-level: switch off = no ra_hist, no guard columns; on = the planner's own burn while the prediction is inside the window
    # Revision 20b OLT before activation equals OL: same decisions and logged fields (a policy built the way planner_run_task builds it, synthetic measurements, no flight)
    olt_identity_selftest(); ttg_policy_selftest(); v4_policy_selftest(); node_v2_policy_selftest()
    println("planner selftest ok: lever, DP offsets equal planner_lib.py, plan inputs, OL schedule, B-ODY rule, table loader, geometry injection (ODY), safety stop, OLT = OL before activation (policy-level), linear-cost DP (lin = zeros bit-identical), FB term (ODY N 216 -> 231, G1 unchanged), IC plans, rev20olt 13 flights, guard v2 core + rev20val 28 flights, $(plan_total()) flights in the plan")
end

"""Minimal PlannerPolicy of arm `arm` over synthetic measurements (no flight, no GRAM): the structure of planner_run_task with the cal files replaced by constants."""
function selftest_policy(arm::AbstractString, apriori)
    meta, sc, g, ep = stageg3_load(CFGP, "ODY"); ms = stagef_setup(sc); G3x = meta["stageg3"]; w2m = Float64(sc["band"]["w_cm2_to_w_m2"])
    cfg, hq, hl = stagef_corridor_config(ms, 4700.0, Float64(G3x["keeper_margin_m"])); c = cfg.values["comparison"]
    inner = CorridorGuidancePolicy(cfg, :favorable, Symbol(meta["stagef"]["policy"]), Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]), Float64(c["total_delta_v_budget_mps"]),
        Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]), deg2rad(Float64(c["planning_lead_deg"])), 0, false, nothing, NamedTuple[], NamedTuple[], String[])
    hh = collect(0.0:5000.0:400000.0); pwf = pw_fit(hh, 2.0e-8 .* exp.(-(hh .- 100000.0) ./ 6500.0), Float64(sc["atmosphere"]["temperature_k"]), "piecewise")
    meas = NamedTuple[]; vlog = NamedTuple[]
    fpol = StageFPolicy(inner, nothing, :fixed, pwf, 4700.0, 0.14 * w2m, 0.32 * w2m, 1000.0, meas, "", 0.0, 0.0, 100.0, NamedTuple[])
    P, B = planner_constants(meta, w2m, 6500.0, arm); sw = planner_switches(meta, arm); tab = planner_table(meta, "ODY")
    return PlannerPolicy(fpol, Float64(ms.harm.gm_m3s2), t -> zeros(3), 4700.0, NaN, Symbol(String(meta["planner"]["arms"][arm]["mode"])), sw, P, B, tab, vlog, 6500.0, 0.212, 0.0, [0.0], apriori, -Inf, NaN, Float64[], 0.0, NamedTuple[], NaN, 0.0, 0, nothing, 0, nothing, nothing, Float64[], NaN, 0, nothing, nothing, nothing), P
end

function olt_identity_selftest()
    gt = planner_table(META, "ODY"); Pb, _ = planner_constants(META, 1.0e4, 6500.0, "OL"); sched = planner_ol_schedule(gt, Pb, log(0.212)); @assert !isempty(sched)
    wO, P = selftest_policy("OL", sched); wT, _ = selftest_policy("OLT", copy(sched)); @assert wO.mode === :ol && wT.mode === :ol && wT.sw.terminal && !wO.sw.terminal
    mu = wO.mu; ra_target = P.ra_target_m; rp = 3.49e6; ras = [2.5e7, 2.0e7, 1.2e7, 9.0e6, ra_target + 4.0e5, ra_target + 3.0e5, ra_target + 2.0e5, ra_target + 1.3e5]      # far from the target: inactive
    for (i, ra) in enumerate(ras)
        j = i + 1; a = (ra + rp) / 2; e = (ra - rp) / (ra + rp)
        push!(wO.f.meas, (; rho_mean=1.0e-8, alt_m=1.1e5, lon_rad=0.1 * i)); push!(wO.vlog, (; pass=i, vrel_p_mps=4700.0, vp_inertial_mps=4500.0))
        push!(wT.f.meas, (; rho_mean=1.0e-8, alt_m=1.1e5, lon_rad=0.1 * i)); push!(wT.vlog, (; pass=i, vrel_p_mps=4700.0, vp_inertial_mps=4500.0))
        dO = planner_decide(wO, j, 1.0e5 * j, nothing, 4700.0, a, e); dT = planner_decide(wT, j, 1.0e5 * j, nothing, 4700.0, a, e)
        @assert isequal(dO, dT) "OLT differs from OL before activation at pass $j: $dO vs $dT"
        @assert dT.reason in ("ol_hold", "ol_scheduled")
    end
    @assert length(wT.c_hist) == length(ras) - 1 && isempty(wO.c_hist) && wT.last_ra == ras[end] && wO.nu_cum == wT.nu_cum == 0.0 && wO.apriori == wT.apriori
    # activation: 40 km above the target with >= 3 samples; the aim replaces the schedule from here on (predictor and density fit are synthetic)
    ra = ra_target + 4.0e4; a = (ra + rp) / 2; e = (ra - rp) / (ra + rp); pred = dv -> 1.1e5
    push!(wT.f.meas, (; rho_mean=1.0e-8, alt_m=1.1e5, lon_rad=1.0)); push!(wT.vlog, (; pass=9, vrel_p_mps=4700.0, vp_inertial_mps=4500.0))
    dT = planner_decide(wT, 10, 1.0e6, pred, 4700.0, a, e); @assert dT.reason in ("terminal_aim", "terminal_hold") && isfinite(dT.z0) && dT.n_plan >= 1
end

function ttg_policy_selftest()
    gt = planner_table(META, "ODY"); Pb, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3"); sched = planner_ol_schedule(gt, Pb, log(0.212)); mu = Float64(META["planner"]["analysis"]["mu_mars_m3s2"])
    w3, P = selftest_policy("OLRT3", copy(sched)); wW, _ = selftest_policy("OLRT3-W", copy(sched)); @assert wW.sw.ttg_guard && !w3.sw.ttg_guard
    t0 = planner_ttg_setup(META, "ODY", gt, P, mu, log(0.212)); @assert t0.tstar_nominal_d > 40 && t0.tstar_nominal_d < 90 && t0.tstar_s == (t0.tstar_nominal_d + 0.15) * 86400.0
    wW.ttg = (; t0..., G=(; t0.G..., s_frac=1.0e6))                                          # an unreachably wide window: the guard never engages and the decisions must equal OLRT3's
    rp = 3.49e6; pred = dv -> 1.1e5 + 300.0 * dv; nb = 0
    for i in 1:30
        ra = 2.4e7 - 4.5e5 * i; a = (ra + rp) / 2; e = (ra - rp) / (ra + rp); rho = 1.0e-8 * (1 + 0.3 * sin(0.7 * i))
        for w in (w3, wW); push!(w.f.meas, (; t_s=1.5e5 * i, rho_mean=rho, alt_m=1.1e5 + 50 * i, lon_rad=0.1 * i)); push!(w.vlog, (; pass=i, vrel_p_mps=4700.0, vp_inertial_mps=4500.0)); end
        d3 = planner_decide(w3, i + 1, 1.5e5 * i, pred, 4700.0, a, e); dW = planner_decide(wW, i + 1, 1.5e5 * i, pred, 4700.0, a, e)
        @assert d3.dv == dW.dv && d3.z0 == dW.z0 && d3.reason == dW.reason && dW.ttg_engaged == false && hasproperty(dW, :ttg_dv_ref) && !hasproperty(d3, :ttg_engaged) "guard-off/idle mismatch at pass $(i + 1)"
        nb += d3.dv != 0.0
        for (w, d) in ((w3, d3), (wW, dW)); d.dv != 0.0 && (w.nu_cum += -d.dv / d.w0); push!(w.nu_hist, w.nu_cum); end
    end
    @assert length(wW.ra_hist) == 30 && isempty(w3.ra_hist) && isfinite(wW.ta_ema) && wW.ttg_N > 0
    # an engaged guard: aim the guard at a time 30 d before its own prediction at the next decision
    wW.ttg = (; t0..., tstar_s=wW.ta_ema - 30.0 * 86400.0); ra = 2.4e7 - 4.5e5 * 31; a = (ra + rp) / 2; e = (ra - rp) / (ra + rp)
    push!(wW.f.meas, (; t_s=1.5e5 * 31, rho_mean=1.0e-8, alt_m=1.1e5, lon_rad=0.1)); push!(wW.vlog, (; pass=31, vrel_p_mps=4700.0, vp_inertial_mps=4500.0))
    dE = planner_decide(wW, 32, 1.5e5 * 31, pred, 4700.0, a, e); @assert dE.ttg_engaged && dE.ttg_solves >= 2 && dE.ttg_solves <= 12
end

function v4_policy_selftest()
    gt = planner_table(META, "ODY"); Pb, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3-W4"); sched = planner_ol_schedule(gt, Pb, log(0.212)); mu = Float64(META["planner"]["analysis"]["mu_mars_m3s2"]); nc = META["planner"]["node_check"]
    w, P = selftest_policy("OLRT3-W4", copy(sched)); @assert w.sw.ttg_guard && !w.sw.l1
    kO = META["planner"]["terminal_periapsis"]["kappa"]["ODY"]; metaN = deepcopy(META); metaN["planner"]["guard_sets"]["v4"]["node_check"] = true; tt0 = planner_ttg_setup(metaN, "ODY", gt, P, mu, log(0.212); set="v4"); @assert tt0.G.pos_stage && tt0.G.node_check          # node check switched ON on a copy: exercises the node path
    gc = (; om0=deg2rad(109.7454), Om0=deg2rad(28.1517), i0=deg2rad(93.522), J2=Float64(nc["J2"]), R=Float64(nc["R_m"]), rot=nothing)
    nom = planner_nominal_track(gt, P, mu, kO, log(0.212), gc.i0, gc.om0, gc.Om0, gc.J2, gc.R)
    w.ttg = tt0; w.pos = PlannerPosCtx(nom, gc, nom.days * 86400.0, Float64[], Float64[], Float64[], Float64[], gc.i0, (; fast_d=NaN, slow_d=NaN, tstar_d=nom.days, source="nominal"), true, NamedTuple[])
    nctx = planner_node_ctx(metaN, nothing, nothing, w.sw, metaN["planner"]["arms"]["OLRT3-W4"]; pos=w.pos); @assert nctx.omega_target_deg == mod(rad2deg(planner_j2_at(nom.track, nom.days * 86400.0)[2]), 360.0) && nctx.J2 == gc.J2
    rp = 3.49e6; pred = dv -> 1.1e5 + 300.0 * dv; nst = 0
    for i in 1:60
        ra = 2.4e7 - 3.4e5 * i; a = (ra + rp) / 2; e = (ra - rp) / (ra + rp); rho = 1.0e-8 * (1 + 0.3 * sin(0.7 * i))
        push!(w.f.meas, (; t_s=1.2e5 * i, rho_mean=rho, alt_m=1.1e5 + 50 * i, lon_rad=0.1 * i)); push!(w.vlog, (; pass=i, vrel_p_mps=4700.0, vp_inertial_mps=4500.0))
        w.node_now = (; omega_now_deg=28.0 + 0.01 * i, rate_deg_s=-1.0e-6, omega_target_deg=nctx.omega_target_deg)
        d = planner_decide(w, i + 1, 1.2e5 * i, pred, 4700.0, a, e); @assert hasproperty(d, :ttg_engaged) && isfinite(d.dv)
        d.dv != 0.0 && (w.nu_cum += -d.dv / d.w0); push!(w.nu_hist, w.nu_cum)
    end
    ps = w.pos; @assert length(ps.ka) == 60 && length(ps.de) == 60 && length(ps.inc) == 60 && length(ps.rp_hist) == 60 && all(isfinite, ps.ka) && all(isfinite, ps.de) && length(w.ra_hist) == 60 && all(0.8 .< ps.ka .< 1.2)
    # result columns of a finished flight (synthetic last row) and the all-NaN default
    last = (; a_m=4.45e6, e=0.215, i_rad=gc.i0, raan_rad=gc.Om0 + 0.3, omega_rad=gc.om0 + 0.2, time_s=nom.days * 86400.0 - 9.0e3, ra_m=5.4e6, rp_m=3.5e6)
    pc = planner_pos_final(ps, mu, last); @assert isfinite(pc.pos_err_deg) && isfinite(pc.l1_err_deg) && pc.l1_err_deg <= pc.pos_err_deg + 1e-9 && isfinite(pc.pos_t_arr_d) && pc.pos_tstar_d == nom.days && pc.tstar_source == "nominal" && planner_pos_final(ps, mu, nothing).pos_err_deg |> isnan
    @assert all(isnan(getfield(PLANNER_POS_NAN, k)) for k in keys(PLANNER_POS_NAN) if getfield(PLANNER_POS_NAN, k) isa Float64)
end

function node_v2_policy_selftest()
    gt = planner_table(META, "ODY"); P, _ = planner_constants(META, 1.0e4, 6500.0, "OLRT3-W4LN"); sched = planner_ol_schedule(gt, P, log(0.212)); mu = Float64(META["planner"]["analysis"]["mu_mars_m3s2"]); nc = META["planner"]["node_check"]
    w, _ = selftest_policy("OLRT3-W4LN", copy(sched)); kO = META["planner"]["terminal_periapsis"]["kappa"]["ODY"]; tt0 = planner_ttg_setup(META, "ODY", gt, P, mu, log(0.212); set="v4n"); @assert tt0.G.node_check_v2 && !tt0.G.node_check
    gc = (; om0=deg2rad(109.7454), Om0=deg2rad(28.1517), i0=deg2rad(93.522), J2=Float64(nc["J2"]), R=Float64(nc["R_m"]), rot=nothing); nom = planner_nominal_track(gt, P, mu, kO, log(0.212), gc.i0, gc.om0, gc.Om0, gc.J2, gc.R)
    w.ttg = tt0; w.pos = PlannerPosCtx(nom, gc, nom.days * 86400.0, Float64[], Float64[], Float64[], Float64[], gc.i0, (; fast_d=NaN, slow_d=NaN, tstar_d=nom.days, source="nominal"), true, NamedTuple[])
    W3 = tt0.G.w_max_periods * planner_kepler_period(mu, P.ra_stop_m, tt0.G.w_rp_m); rp = 3.49e6; pred = dv -> 1.1e5 + 300.0 * dv
    for i in 1:60
        ra = 2.4e7 - 3.4e5 * i; a = (ra + rp) / 2; e = (ra - rp) / (ra + rp); rho = 1.0e-8 * (1 + 0.3 * sin(0.7 * i))
        push!(w.f.meas, (; t_s=1.2e5 * i, rho_mean=rho, alt_m=1.1e5 + 50 * i, lon_rad=0.1 * i)); push!(w.vlog, (; pass=i, vrel_p_mps=4700.0, vp_inertial_mps=4500.0))
        d = planner_decide(w, i + 1, 1.2e5 * i, pred, 4700.0, a, e); @assert hasproperty(d, :ttg_engaged) && isfinite(d.dv) && abs(d.ttg_tt_s - tt0.tstar_s) <= W3 + 1e-6 "aim $(d.ttg_tt_s) outside the window around $(tt0.tstar_s)"
        d.dv != 0.0 && (w.nu_cum += -d.dv / d.w0); push!(w.nu_hist, w.nu_cum)
    end
    @assert length(w.pos.ka) == 60
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
    phase in ("cal", "calG", "calx0", "calR23", "calR25", "calR28", "calR31") ? foreach(r -> r.status == "ok" && println("  cal $(r.tag): seed_alt_m=$(round(r.seed_alt_m; digits=1)) q_first=$(round(r.q_first_w_cm2; digits=4)) target=$(r.q_target_w_cm2) vrel_first=$(round(r.vrel_first_mps; digits=1)) hloc=$(round(r.hloc_at_hg_m; digits=1))"), rows) :
        foreach(r -> r.status == "ok" && println("  $(r.run_id): passes=$(r.passes) t_end_days=$(round(r.t_end_days; digits=2)) stop_reason=$(r.stop_reason) dv=$(round(r.total_dv_mps; digits=4)) burns=$(r.burns) core_epoch=$(r.gram_core_epoch)"), rows)
    (isempty(res.failed) && bad == 0) || error("failed or not-ok cases present")
end
campaign(seeds, mission_s) = run_gram_campaign(haskey(seeds[1], :x0) ? planner_cal_x0_task : haskey(seeds[1], :gram_tag) ? planner_cal_task : planner_run_task, seeds, workers(); mission_time_s=mission_s, worker_files=[joinpath(@__DIR__, "planner_common.jl")])

"""Mean-field GRAM sampling of the Revision 25 low-altitude geometries: runs meanfield_sample.jl (unchanged; track-driven, so the target-aware track length comes from the tracks file) as the sampling_request.json of WP25-P says:
--request=<sampling_request.json> (the tracks csv beside it, sha256 checked against the request) and --mf-out=<dir>; the command is `julia --project=. meanfield_sample.jl --phase=sample --out=DIR --tracks=CSV --cut=none`. With --dry-run it only prints the command."""
function mf25()
    rq = get(CLI, "request", ""); out = get(CLI, "mf-out", ""); (isempty(rq) || isempty(out)) && error("--phase=mf25 needs --request=<sampling_request.json> and --mf-out=<dir>")
    J = JSON.parsefile(rq); tr = joinpath(dirname(abspath(rq)), String(J["tracks_csv"])); bytes2hex(sha256(read(tr))) == String(J["tracks_sha256"]) || error("tracks csv sha256 differs from the request: $tr")
    cmd = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --threads=1 $(joinpath(@__DIR__, "meanfield_sample.jl")) --phase=sample --out=$out --tracks=$tr --cut=none`
    DRY ? println(cmd) : run(cmd)
end

function main()
    phase = get(CLI, "phase", "")
    phase == "mf25" && return mf25()
    if DRY
        phase in PHASES && return dry_run([phase]); phase == "" && return dry_run(PHASES); error("--phase must be one of $(PHASES) with --dry-run (or omitted for all)")
    end
    phase == "selftest" && return selftest()
    mkpath(SFS)
    phase in PHASES || error("--phase must be selftest|" * join(PHASES, "|"))
    summary = joinpath(SFS, "summary_$(phase)_mars_$(HOST).csv"); isfile(summary) && error("$summary exists: use a distinct --host stamp per call")
    seeds = filter(c -> !ondisk(c), cases(phase)); println("$phase: $(length(seeds)) cases to fly ($(length(cases(phase)) - length(seeds)) already on disk)")
    isempty(seeds) && error("no cases (all runs already on disk)")
    finish(campaign(seeds, phase in ("cal", "calG", "calx0", "calR23", "calR25", "calR28", "calR31") ? 86400.0 : maximum(c.horizon_days for c in seeds) * 86400.0), phase)
end
main()
