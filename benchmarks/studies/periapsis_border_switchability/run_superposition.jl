# Part 3 (task B1): superposition check. Replays maintenance schedules (computed by the A1 python library from a passive or a
# closed-loop truth history) as burns in the full truth model (Sun+Moon via SPICE, J2/J3, US76 drag) on the process route.
#   --phase=null     zero-kick replay (checks the callback and determinism against the stored passive run) + finite-difference
#                    lever probes at three apoapses of the passive run
#   --phase=replay   IMPULSIVE replay of schedules/<tag>__<sched>__it<N>.csv (time-event kicks at apoapsis; --iter=N --sched=NAME)
#   --phase=finite   FINITE-BURN replay (the engine's thruster path: 4 N, Isp 300 s, burn centered on apoapsis) of the same file
#   --pilot          one cell, serial, timed (phase null)
#   julia --project=. run_superposition.jl --phase=... [--results-root=DIR] [--workers=N] [--only=tag,tag] [--iter=N] [--sched=lazy_s0]
# Reads <root>/superposition/sup_cells.csv (analyze_superposition.py select) and the passive truth under <root>/border{,2}/grid.
# Kick machinery: a ContinuousCallback written here (no src/ change). Time-event kicks, not the restart-from-state machinery of
# run_switch.jl: one continuous run per schedule, so the kicks and the drift act in one solve.
include(joinpath(@__DIR__, "run_static.jl"))
import DiffEqBase

const SUP_CONFIG_PATH = joinpath(@__DIR__, "superposition_config.toml")
sup_old_dir(tag) = tag[7:end]     # I0500_R000_W000_E1 -> R000_W000_E1

passive_dir(root, tag, inc) = inc == 50.0 ? joinpath(root, "border", "grid", sup_old_dir(tag)) : joinpath(root, "border2", "grid", tag)

# ---- impulsive tangential kick at apoapsis ------------------------------------------------------------------------------------------

"""
Continuous callback: at the j-th apoapsis crossing (radial-velocity downcrossing of r.v = OFFSET, a hair before the true apoapsis so
the root is not shared with the engine's apsis callback at r.v = 0) apply the signed tangential impulse `schedule[j + 1]` (m/s,
along the inertial velocity; schedule[k] is the burn before pass k; the run starts at apoapsis, so pass 1's burn must be 0). The condition
is held at -1 for `guard_s` (5000 s, ASSUMED: far shorter than the orbit period, long enough that r.v is hugely negative on release) after a kick, so the kicked root cannot be found twice.
`log` receives one row per firing.
"""
function apoapsis_kick_callback(schedule::Vector{Float64}, log::Vector{NamedTuple}; offset::Float64=1.0, guard_s::Float64=5000.0)
    last_t = Ref(-Inf); nfired = Ref(0)
    function condition(u, t, integrator)
        t - last_t[] < guard_s && return -1.0
        return dot(u.sc[1].pos, u.sc[1].vel) - offset
    end
    function affect!(integrator)
        sc = integrator.u.sc[1]
        r = SVector{3, Float64}(sc[1], sc[2], sc[3]); v = SVector{3, Float64}(sc[4], sc[5], sc[6])
        last_t[] = Float64(integrator.t)
        nfired[] += 1; k = nfired[] + 1
        dv = k <= length(schedule) ? schedule[k] : 0.0
        mu = integrator.p.args.environment_model.planet.μ
        el0 = _corridor_element_summary(r, v, mu)
        if dv != 0.0
            vn = v + dv * normalize(v)
            for i in 1:3; sc[3 + i] = vn[i]; end
            buffers = integrator.p.shared_buffers
            hasproperty(buffers, :density_sample_t) && (buffers.density_sample_t[1] = NaN)
            hasproperty(buffers, :in_atmosphere_sample_t) && (buffers.in_atmosphere_sample_t[1] = NaN)
            applicable(DiffEqBase.derivative_discontinuity!, integrator, true) && DiffEqBase.derivative_discontinuity!(integrator, true)
        else
            vn = v
        end
        el1 = _corridor_element_summary(r, vn, mu)
        push!(log, (; firing=nfired[], pass=k, t_s=Float64(integrator.t), dv_mps=dv, a_before_m=el0.a, e_before=el0.e, a_after_m=el1.a, e_after=el1.e,
            rp_before_m=el0.a * (1 - el0.e), rp_after_m=el1.a * (1 - el1.e), nu_before_rad=el0.nu))
        return nothing
    end
    initialize(cb, u, t, integrator) = (last_t[] = -Inf; nfired[] = 0; empty!(log); nothing)
    return ContinuousCallback(condition, nothing, affect!; initialize=initialize)
end

"""Impulsive replay of `schedule`; returns `(rows, kicks, wall_s, retcode, sol, args)`; rows = the periapsis table with band columns."""
function impulsive_replay(cfg, study, epoch, horizon_s, dir, schedule)
    args = build_args(cfg, study, epoch, Float64(horizon_s), dir; third_bodies=true)
    log = NamedTuple[]
    cb = apoapsis_kick_callback(schedule, log)
    t0 = time()
    sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(cb,))
    wall = time() - t0
    rows = with_band(periapsis_table(sol, args, cfg, study), cfg, study)
    return (; rows, kicks=log, wall_s=wall, retcode=string(sol.retcode), sol, args)
end

# ---- finite-burn replay -------------------------------------------------------------------------------------------------------------

"""Pre-apoapsis callback (the corridor example's, 12 deg lead) that commands `schedule[j + 1]` at the j-th pre-apoapsis crossing and logs the burn plan."""
function finite_burn_callback(ctrl::CorridorFixedBurnGuidance, schedule::Vector{Float64}, plans::Vector{NamedTuple}, lead_rad::Float64)
    count = Ref(0); target_nu = π - lead_rad
    function true_anomaly(u, planet)
        return _corridor_element_summary(SVector{3, Float64}(u.sc[1].pos), SVector{3, Float64}(u.sc[1].vel), planet.μ).nu
    end
    condition(u, t, integrator) = sin(true_anomaly(u, integrator.p.args.environment_model.planet) - target_nu)
    function affect!(integrator)
        dot(integrator.u.sc[1].pos, integrator.u.sc[1].vel) > 0.0 || return nothing
        count[] += 1; j = count[]; k = j + 1
        dv = k <= length(schedule) ? schedule[k] : 0.0
        ctrl.delta_v_mps = dv; ctrl.source_orbit = j; ctrl.armed = true
        SimulationModel.SimulationCallbacks.schedule_event_driven_thruster_controls!(integrator, 1)
        ctrl.armed = false
        if dv != 0.0
            plan = integrator.p.shared_buffers.maneuver_burn_plans[1]
            plan.valid && plan.source_orbit == j && plan.stop_burn_s > plan.start_burn_s ||
                throw(ErrorException("finite-burn plan rejected at t=$(integrator.t) for firing $j (dv=$dv)"))
            push!(plans, (; firing=j, pass=k, request_t_s=Float64(integrator.t), start_s=plan.start_burn_s, stop_s=plan.stop_burn_s,
                dv_cmd_mps=dv, direction_rad=plan.direction_rad, propellant_kg=plan.propellant_required_kg))
        end
        return nothing
    end
    return ContinuousCallback(condition, affect!, nothing; save_positions=(false, false))
end

function finite_burn_replay(cfg, study, epoch, horizon_s, dir, schedule)
    args0 = build_args(cfg, study, epoch, Float64(horizon_s), dir; third_bodies=true)
    v = cfg.values["vehicle"]
    ctrl = CorridorFixedBurnGuidance(0.0, 1, false)
    thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
    args = SimulationConfiguration(
        file_paths=args0.file_paths, simulation_settings=args0.simulation_settings, mission_configuration=args0.mission_configuration,
        environment_model=args0.environment_model, dynamics_model=args0.dynamics_model,
        guidance_model=GuidanceModel((ctrl,), [1.0e12]), navigation_model=args0.navigation_model, control_model=ControlModel((thruster,), [1.0e12]),
        initial_time=args0.initial_time, integration_tolerances=args0.integration_tolerances, solver_config=args0.solver_config)
    plans = NamedTuple[]
    cb = finite_burn_callback(ctrl, schedule, plans, deg2rad(Float64(cfg.values["comparison"]["planning_lead_deg"])))
    t0 = time()
    sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(cb,))
    wall = time() - t0
    rows = with_band(periapsis_table(sol, args, cfg, study), cfg, study)
    return (; rows, plans, wall_s=wall, retcode=string(sol.retcode), sol, args)
end

# ---- finite-difference lever probe -----------------------------------------------------------------------------------------------

"""Central-difference lever at the apoapsis before pass `k` of a solved passive run: m of periapsis altitude per m/s, tangential kick of +-dv."""
function fd_lever_at(sol, args, cfg, study, epoch, dir, k_apo::Int, dv::Float64)
    geo = _geodetic_of(args); ev = physical_apsis_events(sol, geo)
    apo = [e for e in ev if e.kind == "apoapsis"]
    k_apo <= length(apo) || return (; ok=false)
    ta = apo[k_apo].time_s
    nxt = [p.time_s for p in ev if p.kind == "periapsis" && p.time_s > ta]
    isempty(nxt) && return (; ok=false)
    tp = first(nxt)
    s = state_at(sol, ta); ep = _shift_initial_time(epoch, ta); horizon = tp - ta + 600.0
    ref = run_truth(cfg, study, ep, horizon, dir; third_bodies=true, state=(s.r, s.v), mass_kg=s.m)
    up = run_truth(cfg, study, ep, horizon, dir; third_bodies=true, state=(s.r, s.v + dv * normalize(s.v)), mass_kg=s.m)
    dn = run_truth(cfg, study, ep, horizon, dir; third_bodies=true, state=(s.r, s.v - dv * normalize(s.v)), mass_kg=s.m)
    (isempty(ref.rows) || isempty(up.rows) || isempty(dn.rows)) && return (; ok=false)
    el = _corridor_element_summary(s.r, s.v, ref.args.environment_model.planet.μ)
    return (; ok=true, apo_index=k_apo, t_apo_days=ta / SECONDS_PER_DAY, a_m=el.a, e=el.e, dv_mps=dv, ref_alt_m=ref.rows[1].alt_m,
        dalt_per_mps=(up.rows[1].alt_m - dn.rows[1].alt_m) / (2dv), asym_m=(up.rows[1].alt_m + dn.rows[1].alt_m) / 2 - ref.rows[1].alt_m)
end

# ---- tasks ----------------------------------------------------------------------------------------------------------------------------

function sup_setup(seed)
    study = load_study(seed.study_path); nbody_planet()
    pdir = passive_dir(seed.root, seed.tag, seed.inc)
    cfg = cell_config(study, seed.inc, seed.raan, seed.omega; seed_alt=seed_alt_of(pdir))
    epoch = parse_epoch(study["grid"]["epochs_utc"][seed.epoch_idx])
    return study, cfg, epoch, pdir
end

read_schedule(path) = (d = CSV.read(path, DataFrame); v = zeros(maximum(d.pass)); v[d.pass] .= d.dv_mps; v)

function null_task(seed)
    study, cfg, epoch, pdir = sup_setup(seed)
    horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
    dir = joinpath(seed.out, "runs", "$(seed.tag)__null")
    r = impulsive_replay(cfg, study, epoch, horizon, dir, zeros(0))
    write_rows(joinpath(dir, "passes_cl.csv"), r.rows)
    stored = CSV.read(joinpath(pdir, "passes_3b.csv"), DataFrame)
    n = min(nrow(stored), length(r.rows))
    dalt = maximum(abs.([r.rows[k].alt_m - stored.alt_m[k] for k in 1:n]))
    fds = NamedTuple[]
    npass = length(r.rows)
    for k in unique(clamp.([2, round(Int, npass / 3), round(Int, 2 * npass / 3)], 2, npass - 2))
        f = try fd_lever_at(r.sol, r.args, cfg, study, epoch, dir, k, Float64(seed.fd_dv_mps)) catch err; (; ok=false, error=first(sprint(showerror, err), 200)) end
        push!(fds, merge((; tag=seed.tag, k_apoapsis=k), f))
    end
    write_rows(joinpath(dir, "fd_lever.csv"), [merge(f, (; tag=seed.tag)) for f in fds])
    return (; tag=seed.tag, wall_s=r.wall_s, retcode=r.retcode, passes=length(r.rows), stored_passes=nrow(stored), max_alt_diff_to_stored_m=dalt,
        fd_ok=count(f -> f.ok, fds))
end

function replay_task(seed)
    study, cfg, epoch, pdir = sup_setup(seed)
    horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
    name = "$(seed.tag)__$(seed.sched)__it$(seed.iter)"
    sched = read_schedule(joinpath(seed.out, "schedules", name * ".csv"))
    dir = joinpath(seed.out, "runs", name)
    r = impulsive_replay(cfg, study, epoch, horizon, dir, sched)
    write_rows(joinpath(dir, "passes_cl.csv"), r.rows); write_rows(joinpath(dir, "kicks.csv"), r.kicks)
    return (; name, tag=seed.tag, sched=seed.sched, iter=seed.iter, wall_s=r.wall_s, retcode=r.retcode, passes=length(r.rows), n_sched=count(!=(0.0), sched),
        n_fired=count(k -> k.dv_mps != 0.0, r.kicks), dv_sched_mps=sum(abs, sched), dv_applied_mps=sum(abs(k.dv_mps) for k in r.kicks; init=0.0))
end

function finite_task(seed)
    study, cfg, epoch, pdir = sup_setup(seed)
    horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
    name = "$(seed.tag)__$(seed.sched)__it$(seed.iter)"
    sched = read_schedule(joinpath(seed.out, "schedules", name * ".csv"))
    dir = joinpath(seed.out, "runs", name * "__finite")
    r = finite_burn_replay(cfg, study, epoch, horizon, dir, sched)
    write_rows(joinpath(dir, "passes_cl.csv"), r.rows); write_rows(joinpath(dir, "burn_plans.csv"), r.plans)
    geo = _geodetic_of(r.args); ev = physical_apsis_events(r.sol, geo)
    write_rows(joinpath(dir, "apsis_events.csv"), [(; kind=e.kind, time_s=e.time_s, alt_m=e.physical_altitude_m, mass_kg=e.mass_kg) for e in ev])
    return (; name=name * "__finite", tag=seed.tag, wall_s=r.wall_s, retcode=r.retcode, passes=length(r.rows), n_plans=length(r.plans),
        dv_sched_mps=sum(abs, sched), dv_cmd_mps=sum(abs(p.dv_cmd_mps) for p in r.plans; init=0.0), propellant_kg=sum(p.propellant_kg for p in r.plans; init=0.0))
end

# ---- driver ---------------------------------------------------------------------------------------------------------------------------

function sup_cli(args)
    o = Dict{String, Any}("root" => nothing, "phase" => "null", "pilot" => false, "workers" => nothing, "only" => nothing, "iter" => 0, "sched" => nothing)
    for a in args
        a == "--pilot" && (o["pilot"] = true)
        startswith(a, "--results-root=") && (o["root"] = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--phase=") && (o["phase"] = split(a, "=", limit=2)[2])
        startswith(a, "--workers=") && (o["workers"] = parse(Int, split(a, "=", limit=2)[2]))
        startswith(a, "--iter=") && (o["iter"] = parse(Int, split(a, "=", limit=2)[2]))
        startswith(a, "--sched=") && (o["sched"] = split(a, "=", limit=2)[2])
        startswith(a, "--only=") && (o["only"] = split(split(a, "=", limit=2)[2], ","))
    end
    return o
end

function sup_main(args=ARGS)
    o = sup_cli(args); study = load_study(); sc = TOML.parsefile(SUP_CONFIG_PATH)
    root = something(o["root"], results_root(study)); out = joinpath(root, "superposition"); mkpath(out)
    cells = CSV.read(joinpath(out, "sup_cells.csv"), DataFrame)
    o["only"] === nothing || (cells = cells[in.(cells.tag, Ref(o["only"])), :])
    hd = Float64(sc["run"]["horizon_days"]); workers = something(o["workers"], Int(sc["run"]["workers"]))
    base(r; kw...) = (; study_path=STUDY_CONFIG_PATH, out, root, tag=r.tag, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg, epoch_idx=Int(r.epoch),
        horizon_days=hd, fd_dv_mps=Float64(sc["fd_lever"]["probe_dv_mps"]), iter=o["iter"], sched=o["sched"], kw...)
    files = [joinpath(@__DIR__, "run_superposition.jl")]
    phase = o["phase"]; mt = hd * SECONDS_PER_DAY
    if o["pilot"]
        t0 = time(); r = null_task(base(first(eachrow(cells))))
        println("PILOT wall_s=$(round(time() - t0; digits=1)) $(r)"); return
    end
    if phase == "null"
        res = run_process_campaign(null_task, [base(r) for r in eachrow(cells)], workers; mission_time_s=mt, worker_files=files)
        report(res); write_rows(joinpath(out, "null_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "null")
    elseif phase == "replay"
        scheds = o["sched"] === nothing ? String.(sc["schedules"]["names"]) : [o["sched"]]
        seeds = [base(r; sched=s) for s in scheds for r in eachrow(cells)]
        res = run_process_campaign(replay_task, seeds, workers; mission_time_s=mt, worker_files=files)
        report(res); write_rows(joinpath(out, "replay_it$(o["iter"])_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "replay_it$(o["iter"])")
    elseif phase == "finite"
        o["sched"] === nothing && throw(ArgumentError("--sched is required for --phase=finite"))
        res = run_process_campaign(finite_task, [base(r) for r in eachrow(cells)], workers; mission_time_s=mt, worker_files=files)
        report(res); write_rows(joinpath(out, "finite_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "finite")
    else
        throw(ArgumentError("--phase must be null, replay or finite"))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    sup_main()
end
