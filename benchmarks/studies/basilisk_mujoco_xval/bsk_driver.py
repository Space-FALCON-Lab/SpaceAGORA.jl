#!/usr/bin/env python3
"""Basilisk 2.12 MJScene side of the MuJoCo cross-validation (see README.md).

  python bsk_driver.py <case.toml> <outdir> [duration_s]

Runs the case once per step in `bsk_dt_s` (Basilisk's fixed-step RK4 on the MJScene dynamics, task step = integration
step) and writes <outdir>/bsk_dt<dt>.csv with the same columns as sagora_driver.jl. Requires a Basilisk build with
BUILD_MUJOCO=ON.
"""
import os, sys, tomllib
import xml.etree.ElementTree as ET
import numpy as np

from Basilisk.simulation import mujoco, svIntegrators
from Basilisk.utilities import SimulationBaseClass, macros, simIncludeGravBody
from Basilisk.architecture import messaging


def mjcf_layout(path):
    """Bodies (document order) and scalar/free joints (document order) of an MJCF file."""
    root = ET.parse(path).getroot()
    bodies, joints = [], []  # joints: (name or None, kind, body)

    def walk(el, body):
        for ch in el:
            if ch.tag == "body":
                bodies.append(ch.get("name")); walk(ch, ch.get("name"))
            elif ch.tag == "freejoint" or (ch.tag == "joint" and ch.get("type") == "free"):
                joints.append((ch.get("name"), "free", body))
            elif ch.tag == "joint":
                joints.append((ch.get("name"), "hinge", body))
            else:
                walk(ch, body)
    walk(root.find("worldbody"), None)
    return bodies, joints


MJLOG = os.path.join(os.getcwd(), "MUJOCO_LOG.TXT")


def logsize():
    return os.path.getsize(MJLOG) if os.path.exists(MJLOG) else 0


def ctrl_value(sched, t):
    v = 0.0
    for tt, vv in zip(sched["times"], sched["values"]):
        if tt <= t + 1e-9:
            v = vv
    return float(v)


def run_one(cfg, dt, duration, casedir):
    path = os.path.join(casedir, cfg["mjcf"])
    bodies, joints = mjcf_layout(path)
    mu, a, inc = cfg["mu"], cfg["orbit_radius_m"], np.deg2rad(cfg["inclination_deg"])
    r0 = np.array([a, 0.0, 0.0]); v0 = np.sqrt(mu / a) * np.array([0.0, np.cos(inc), np.sin(inc)])

    sim = SimulationBaseClass.SimBaseClass()
    proc = sim.CreateNewProcess("p")
    proc.addTask(sim.CreateNewTask("t", macros.sec2nano(dt)))
    scene = mujoco.MJScene.fromFile(path)
    scene.extraEoMCall = True
    sim.AddModelToTask("t", scene)
    scene.setIntegrator(svIntegrators.svIntegratorRK4(scene))

    fac = simIncludeGravBody.gravBodyFactory()
    earth = fac.createEarth()
    earth.isCentralBody = True
    earth.mu = mu                      # SpaceAGORA's Earth mu, so both tools share one point-mass field
    fac.addBodiesTo(scene)             # NBodyGravity: m*g(r_com) at every body COM

    msgs = []
    for c in cfg.get("ctrl", []):
        m = messaging.SingleActuatorMsg()
        scene.getSingleActuator(c["actuator"]).actuatorInMsg.subscribeTo(m)
        m.write(messaging.SingleActuatorMsgPayload(input=ctrl_value(c, 0.0)))
        msgs.append((m, c))

    sim.InitializeSimulation()
    for b in cfg["body"]:
        body = scene.getBody(b["name"])
        body.setPosition(r0 + np.array(b["dr"])); body.setVelocity(v0 + np.array(b["dv"]))
        if any(b["omega"]):
            body.setAttitudeRate(np.array(b["omega"]))
    jobj = {}
    for name, kind, body in joints:
        if kind == "hinge":
            jobj[name] = scene.getBody(body).getScalarJoint(name)
    for jc in cfg.get("joint", []):
        jobj[jc["name"]].setPosition(jc["qpos"]); jobj[jc["name"]].setVelocity(jc["qvel"])

    free = [(n, b) for n, k, b in joints if k == "free"]
    scalar = [(n, b) for n, k, b in joints if k == "hinge"]
    coms = {b: scene.getBody(b).getCenterOfMass().stateOutMsg for b in bodies}
    origins = {b: scene.getBody(b).getOrigin().stateOutMsg for b in bodies}   # site frame = body frame

    cols = ["t"]
    for b in bodies:
        cols += [f"{b}_r{k}" for k in "xyz"] + [f"{b}_v{k}" for k in "xyz"]
    for n, b in free:
        cols += [f"{b}_q{k}" for k in "wxyz"] + [f"{b}_w{k}" for k in "xyz"]
    for n, b in scalar:
        cols += [f"q_{n}", f"qd_{n}"]

    def mrp_to_q_mj(sig):
        """MJSite writes sigma_BN = MRP of the MuJoCo site matrix site_xmat (body-to-world) itself (MJSite.cpp,
        `mrpd = rot`), so the MuJoCo quaternion (w, x, y, z) is the plain MRP -> Euler-parameter map, no conjugate."""
        s2 = float(np.dot(sig, sig))
        return [(1 - s2) / (1 + s2)] + list(2 * np.asarray(sig) / (1 + s2))

    def record(t):
        row = [t]
        for b in bodies:
            p = coms[b].read()
            row += list(np.array(p.r_BN_N).ravel()) + list(np.array(p.v_BN_N).ravel())
        for n, b in free:
            o = origins[b].read()
            row += mrp_to_q_mj(np.array(o.sigma_BN).ravel()) + list(np.array(o.omega_BN_B).ravel())
        for n, b in scalar:
            row += [jobj[n].stateOutMsg.read().state, jobj[n].stateDotOutMsg.read().state]
        return row

    nsamp = int(round(duration / cfg["sample_dt_s"]))
    rows = []
    for s in range(nsamp + 1):
        t = s * cfg["sample_dt_s"]
        if s == 0:
            sim.ConfigureStopTime(0); sim.ExecuteSimulation()
        else:
            sim.ConfigureStopTime(macros.sec2nano(t)); sim.ExecuteSimulation()
        rows.append(record(t))
        for m, c in msgs:   # value held over the next sample interval (switch times are multiples of sample_dt_s)
            m.write(messaging.SingleActuatorMsgPayload(input=ctrl_value(c, t)))
    return cols, rows


def main():
    cfgfile, outdir = sys.argv[1], sys.argv[2]
    cfg = tomllib.load(open(cfgfile, "rb"))
    duration = float(sys.argv[3]) if len(sys.argv) > 3 else cfg["duration_s"]
    os.makedirs(outdir, exist_ok=True)
    for dt in cfg["bsk_dt_s"]:
        size0 = logsize()
        cols, rows = run_one(cfg, dt, duration, os.path.dirname(os.path.abspath(cfgfile)))
        # MJScene throws on NaN accelerations itself; any MuJoCo warning in the log or a non-finite value aborts the level.
        if logsize() != size0 or not np.all(np.isfinite(np.array(rows))):
            sys.exit(f"unstable Basilisk run at dt={dt}: MuJoCo warning in {MJLOG} or non-finite output; level aborted")
        with open(os.path.join(outdir, f"bsk_dt{dt}.csv"), "w") as f:
            f.write(",".join(cols) + "\n")
            for r in rows:
                f.write(",".join(repr(float(x)) for x in r) + "\n")
        print(f"bsk dt={dt}: {len(rows)} samples", flush=True)


if __name__ == "__main__":
    main()
