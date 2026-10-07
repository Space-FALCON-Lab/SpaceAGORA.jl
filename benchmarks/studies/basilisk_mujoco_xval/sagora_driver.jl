# SpaceAGORA side of the Basilisk MJScene cross-validation (see README.md).
#   julia --project=<env with SpaceAGORA and SpaceAGORAMuJoCo> sagora_driver.jl <case.toml> <outdir> [duration_s]
# Runs the case once per dt in `sagora_dt_s` and writes <outdir>/sagora_dt<dt>.csv.
using TOML, StaticArrays, LinearAlgebra
import SpaceAGORA
using SpaceAGORAMuJoCo
using SpaceAGORAMuJoCo: Binding
const B = Binding
const SM = SpaceAGORA.SimulationModel

const OBJ_JOINT = 3
const OBJ_ACTUATOR = 19

ctrl_value(sched, t) = (i = findlast(<=(t + 1e-9), sched["times"]); i === nothing ? 0.0 : Float64(sched["values"][i]))

# Instability guard. MuJoCo answers a NaN/huge QACC with a warning in MUJOCO_LOG.TXT (working directory) and an
# mj_resetData, and carries on, so the run looks valid. The Stage 1 binding exposes no warning counters, so the
# harness checks (a) that the log file did not grow, (b) that mjData.time still equals the scene time (a reset sets it
# to zero), (c) that every recorded value is finite.
const MJLOG = joinpath(pwd(), "MUJOCO_LOG.TXT")
logsize() = isfile(MJLOG) ? filesize(MJLOG) : 0

function run_one(cfg, dt, duration, casedir)
    mu = cfg["mu"]; a = cfg["orbit_radius_m"]; inc = deg2rad(cfg["inclination_deg"])
    planet = SM.make_no_gram_planet(:earth)
    isapprox(planet.μ, mu; rtol = 1e-14) || error("case mu $mu differs from SpaceAGORA Earth mu $(planet.μ)")
    r0 = SVector(a, 0.0, 0.0)
    v0 = sqrt(mu / a) * SVector(0.0, cos(inc), sin(inc))
    states = [SceneBodyState(b["name"], r0 + SVector{3}(b["dr"]...), v0 + SVector{3}(b["dv"]...); ω = Tuple(b["omega"]))
              for b in cfg["body"]]
    scene = ProximityScene(; mjcf_path = joinpath(casedir, cfg["mjcf"]), dt, planet,
        gravity_effectors = (SM.InverseSquaredGravityModel(),), initial_states = states,
        integrator = Symbol(get(cfg, "sagora_integrator", "implicitfast")))
    m = scene.model
    qadr = B.jnt_qposadr(m); dadr = B.jnt_dofadr(m); jtype = B.jnt_type(m)
    jnames = [B.id2name(m, OBJ_JOINT, j) for j in 0:B.njnt(m)-1]
    scalar = [j for j in 1:B.njnt(m) if jtype[j] != B.JNT_FREE]
    # Initial joint state. The chief and the free-body states are untouched: every free body keeps its
    # state, only the (jointed) child bodies move, so the absolute free-body states stay as given.
    for jc in get(cfg, "joint", [])
        j = findfirst(==(jc["name"]), jnames); j === nothing && error("no joint $(jc["name"])")
        scene.qpos[qadr[j] + 1] = jc["qpos"]; scene.qvel[dadr[j] + 1] = jc["qvel"]
    end
    scene.fresh = false
    sched = [(B.name2id(m, OBJ_ACTUATOR, c["actuator"]) + 1, c) for c in get(cfg, "ctrl", [])]
    jadr = B.body_jntadr(m)
    bodies = scene_body_names(scene)
    cols = String["t"]
    for b in bodies; append!(cols, ["$(b)_r$k" for k in "xyz"], ["$(b)_v$k" for k in "xyz"]); end
    for b in scene.names[scene.roots]
        append!(cols, ["$(b)_q$k" for k in "wxyz"], ["$(b)_w$k" for k in "xyz"])
    end
    for j in scalar; append!(cols, ["q_$(jnames[j])", "qd_$(jnames[j])"]); end
    nsamp = round(Int, duration / cfg["sample_dt_s"]); per = round(Int, cfg["sample_dt_s"] / dt)
    abs(per * dt - cfg["sample_dt_s"]) < 1e-12 || error("sample_dt_s must be a multiple of dt=$dt")
    rows = Vector{Vector{Float64}}()
    function record()
        row = Float64[scene_time(scene)]
        for i in 1:scene.nb
            s = scene_body_state(scene, i); append!(row, s.r, s.v)
        end
        for id in scene.roots
            j = jadr[id + 1] + 1
            append!(row, scene_body_state(scene, id).q, [scene.qvel[dadr[j] + 3 + k] for k in 1:3])
        end
        for j in scalar; push!(row, scene.qpos[qadr[j] + 1], scene.qvel[dadr[j] + 1]); end
        push!(rows, row)
    end
    record()
    for s in 0:nsamp-1, k in 0:per-1
        t = (s * per + k) * dt
        for (idx, c) in sched; scene_ctrl(scene)[idx] = ctrl_value(c, t); end
        scene_step!(scene)
        if k == per - 1
            abs(B.data_time(scene.data) - scene_time(scene)) <= 1e-9 * max(1.0, scene_time(scene)) ||
                error("unstable: mjData.time=$(B.data_time(scene.data)) != scene time $(scene_time(scene)) (MuJoCo reset the state), dt=$dt")
            record()
        end
    end
    all(all(isfinite, r) for r in rows) || error("unstable: non-finite value recorded, dt=$dt")
    return cols, rows
end

function main(args)
    cfgfile, outdir = args[1], args[2]
    cfg = TOML.parsefile(cfgfile)
    duration = length(args) >= 3 ? parse(Float64, args[3]) : cfg["duration_s"]
    mkpath(outdir)
    for dt in cfg["sagora_dt_s"]
        size0 = logsize()
        cols, rows = run_one(cfg, dt, duration, dirname(abspath(cfgfile)))
        logsize() == size0 || error("MuJoCo wrote a warning to $MJLOG during dt=$dt (instability); level aborted, no output written")
        open(joinpath(outdir, "sagora_dt$(dt).csv"), "w") do io
            println(io, join(cols, ","))
            for r in rows; println(io, join(r, ",")); end
        end
        println("sagora dt=$dt: $(length(rows)) samples")
    end
end
main(ARGS)
