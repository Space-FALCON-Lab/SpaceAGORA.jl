# Half-lever test I1: single-burn vacuum flight at the Stage O Odyssey-band start (tasks/2026-10-15_periapsis-region-border-targeting/halflever_results.md, "What would settle it").
# Atmosphere OFF (no drag effector; variant nosun_nodrag = Mars50c 20x20 only, nodrag = + Sun), ONE prograde burn through the SAME scheduler the Stage O keepers use
# (a BaseThrusterModel with the engine's apoapsis-centered window, armed by the pre-apoapsis crossing at the Stage O planning lead; the callback is the one of run_superposition.jl
# finite_burn_callback / _preapoapsis_guidance_callback, copied because including run_superposition.jl would re-include run_static.jl). No src/ or examples/ change, no new burn model.
# Start: Stage O (stageo_config.toml [stageo], odyssey band): epoch, apoapsis radius 28,533 km, i/RAAN/omega, calibrated first-pass periapsis altitude from the Stage O odyssey cal guess.csv (SOURCED file).
# The start is at apoapsis (t = 0); the burn is the one centered on the NEXT apoapsis (about one period later), as in the flights.
# Writes <out>/<run>/{events.csv, grid.csv, burn.csv} and <out>/meta.csv. Frame of all states: Mars body-equator inertial (rot' * J2000).
#   systemd-run --user --scope -p MemoryMax=24G julia --project=. halflever_vacuum.jl --out=DIR [--cal=guess.csv] [--only=nom,none]
using OrdinaryDiffEq
include(joinpath(@__DIR__, "stageo_common.jl"))
import DiffEqBase

const HLV_CLI = Dict(String(m.captures[1]) => String(something(m.captures[2], "1")) for a in ARGS for m in [match(r"^--([^=]+)(?:=(.*))?$", a)] if m !== nothing)
const HLV_OUT = abspath(get(HLV_CLI, "out", joinpath(@__DIR__, "..", "..", "..", "output", "halflever_vacuum")))
const HLV_CAL = get(HLV_CLI, "cal", joinpath(get(ENV, "PLANNER_TASK_ROOT", "."), "stageo_data/phase_cal_smoke_gate/cal/odyssey/guess.csv"))
# (name, variant, dv_mps, burn-window shift s). shift = ASSUMED +-60 s per the task (timing sensitivity); applied to the scheduler's own window after it is set.
const HLV_RUNS = [("nom", "nosun_nodrag", 0.25, 0.0), ("none", "nosun_nodrag", 0.0, 0.0), ("dv014", "nosun_nodrag", 0.14, 0.0),
    ("tp60", "nosun_nodrag", 0.25, 60.0), ("tm60", "nosun_nodrag", 0.25, -60.0), ("nom_sun", "nodrag", 0.25, 0.0), ("none_sun", "nodrag", 0.0, 0.0)]

function hlv_flight(name, variant, dv, shift, ms, sc, epoch, rot, seed_alt, SO, lead_rad, horizon_s, fine_center)
    dir = joinpath(HLV_OUT, name); mkpath(dir)
    state0 = mars_initial_state(ms, rot, seed_alt, Float64(SO["inclination_deg"]), Float64(SO["raan_deg"]), Float64(SO["omega_deg"]))
    ctrl = CorridorFixedBurnGuidance(0.0, 1, false)
    v = sc["vehicle"]; thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
    args = mars_args(ms, epoch, horizon_s, dir; state=state0, variant=variant, controller=ctrl, thruster=thruster)
    Rt = transpose(rot); events = NamedTuple[]; burn = NamedTuple[]; count = Ref(0); target_nu = π - lead_rad; mu = args.environment_model.planet.μ
    # apsis events (radial velocity zero crossings), body-equator frame state
    function rec!(integrator, kind)
        integrator.t > 1e-3 || return nothing
        s = integrator.u.sc[1]; r = Rt * SVector{3, Float64}(s.pos); w = Rt * SVector{3, Float64}(s.vel)
        push!(events, (; kind, t_s=Float64(integrator.t), rx=r[1], ry=r[2], rz=r[3], vx=w[1], vy=w[2], vz=w[3], mass_kg=Float64(s.mass))); return nothing
    end
    apsis = ContinuousCallback((u, t, i) -> dot(u.sc[1].pos, u.sc[1].vel), i -> rec!(i, "periapsis"), i -> rec!(i, "apoapsis"); save_positions=(false, false))
    nu_of(u) = _corridor_element_summary(SVector{3, Float64}(u.sc[1].pos), SVector{3, Float64}(u.sc[1].vel), mu).nu
    function affect!(integrator)
        dot(integrator.u.sc[1].pos, integrator.u.sc[1].vel) > 0.0 || return nothing
        count[] += 1; j = count[]
        d = j == 1 ? dv : 0.0; ctrl.delta_v_mps = d; ctrl.source_orbit = j; ctrl.armed = true
        SimulationModel.SimulationCallbacks.schedule_event_driven_thruster_controls!(integrator, 1)
        ctrl.armed = false
        if d != 0.0
            plan = integrator.p.shared_buffers.maneuver_burn_plans[1]
            plan.valid && plan.source_orbit == j && plan.stop_burn_s > plan.start_burn_s || throw(ErrorException("burn plan rejected at t=$(integrator.t) dv=$d"))
            th = integrator.p.args.control_model.control_effectors[1]
            if shift != 0.0   # deliberate off-apoapsis burn: shift the scheduler's window; the burn-edge callback re-registers the stops (model window has precedence over the plan)
                th.start_burn_time[1] += shift; th.stop_burn_time[1] += shift
            end
            push!(burn, (; request_t_s=Float64(integrator.t), nu_request_rad=nu_of(integrator.u), plan_start_s=plan.start_burn_s, plan_stop_s=plan.stop_burn_s, win_start_s=th.start_burn_time[1], win_stop_s=th.stop_burn_time[1],
                dv_cmd_mps=d, direction_rad=plan.direction_rad, thrust_n=plan.thrust_n, isp_s=plan.isp_s, mass_kg=Float64(integrator.u.sc[1].mass), shift_s=shift))
        end
        return nothing
    end
    pre = ContinuousCallback((u, t, i) -> sin(nu_of(u) - target_nu), affect!, nothing; save_positions=(false, false))
    t0 = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(apsis, pre)); wall = time() - t0
    # dense grid (same times in every run so burn and no-burn states compare at equal t): 5 s within 700 s of the Keplerian period (the second apoapsis), 120 s elsewhere
    tend = Float64(sol.t[end])
    ts = Float64[]; t = 1.0
    while t <= tend
        push!(ts, t); t += abs(t - fine_center) <= 700.0 ? 5.0 : 120.0
    end
    rows = NamedTuple[]
    for t in ts
        s = sol(t).sc[1]; r = Rt * SVector{3, Float64}(s.pos); w = Rt * SVector{3, Float64}(s.vel)
        push!(rows, (; t_s=t, rx=r[1], ry=r[2], rz=r[3], vx=w[1], vy=w[2], vz=w[3], mass_kg=Float64(s.mass)))
    end
    write_rows_m(joinpath(dir, "events.csv"), events); write_rows_m(joinpath(dir, "grid.csv"), rows)
    isempty(burn) || write_rows_m(joinpath(dir, "burn.csv"), burn)
    return (; name, variant, dv_mps=dv, shift_s=shift, retcode=string(sol.retcode), t_end_s=tend, n_events=length(events), n_burns=length(burn), wall_s=wall, mu_engine=mu)
end

function hlv_main()
    meta, sc = stageo_load(joinpath(@__DIR__, "stageo_config.toml"), "odyssey"); SO = meta["stageo"]; ms = stagef_setup(sc)
    sc["environment"] = Dict{String, Any}("EI_km" => Float64(sc["atmosphere_gram"]["entry_interface_km"]))   # key read by mars_args; the Stage O value (gram_args)
    epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); rot = mars_rot(ms.planet, epoch)
    seed_alt = Float64(CSV.read(HLV_CAL, DataFrame)[1, :].seed_alt_m); vrel = Float64(CSV.read(HLV_CAL, DataFrame)[1, :].vrel_first_mps)
    cfg, _, _ = stagef_corridor_config(ms, vrel, Float64(meta["stagef"]["fixed_margin_m"])); lead = deg2rad(Float64(cfg.values["comparison"]["planning_lead_deg"]))
    mu = ms.harm.gm_m3s2; ra = Float64(SO["ra0_m"]); rp = ms.planet.Rp_e + seed_alt; a = (ra + rp) / 2; T = 2π * sqrt(a^3 / mu); horizon = 1.7 * T
    println("start: seed_alt_m=$seed_alt rp_m=$rp ra_m=$ra a_m=$a T_s=$T lead_deg=$(rad2deg(lead)) horizon_s=$horizon mu_harm=$mu planet_mu=$(ms.planet.μ) gravity=$(sc["gravity"]["degree"])x$(sc["gravity"]["order"])")
    only = haskey(HLV_CLI, "only") ? Set(split(HLV_CLI["only"], ",")) : nothing
    mkpath(HLV_OUT); res = NamedTuple[]
    for (name, variant, dv, shift) in HLV_RUNS
        only === nothing || name in only || continue
        r = hlv_flight(name, variant, dv, shift, ms, sc, epoch, rot, seed_alt, SO, lead, horizon, T); push!(res, r); println(r)
    end
    write_rows_m(joinpath(HLV_OUT, "summary.csv"), res)
    write_rows_m(joinpath(HLV_OUT, "meta.csv"), [(; seed_alt_m=seed_alt, rp0_m=rp, ra0_m=ra, a0_m=a, period_s=T, lead_deg=rad2deg(lead), horizon_s=horizon, mu_harm=mu, mu_planet=ms.planet.μ, rp_eq_m=ms.planet.Rp_e,
        gravity_degree=Int(sc["gravity"]["degree"]), gravity_order=Int(sc["gravity"]["order"]), epoch=string(sc["grid"]["epochs_utc"][1]))])
end
hlv_main()
