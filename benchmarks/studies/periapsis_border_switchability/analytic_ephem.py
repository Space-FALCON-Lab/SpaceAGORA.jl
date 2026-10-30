"""Minimal pure-Python ephemeris access for the analytic third-body term (WP-A3, T8): SPK type 2 reader, UTC -> ET, IAU pole, body-equator frame.

Everything is read from the repo's SPICE bundle (SPICE_DIR below; SOURCED: de440s.bsp, naif0012.tls leap seconds, pck00011.tpc poles, de_403_masses.tpc GM).
Frame convention (SOURCED, mars_common.jl mars_rot): the "inertial" frame of the study pass tables is the body-equator frame, z = body pole at the epoch,
x = node of the J2000 equator on the body equator, y = z x x.  The pole is the pck00011 IAU value (RA, Dec with the linear rate terms); for Earth the engine may
use a higher-precision frame (precession-nutation), which differs by < ~1 deg in the x axis direction (ASSUMED negligible for the tide direction).
"""
import math, os, pathlib, re, struct
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent


def spice_dir():
    env = os.environ.get("SPACEAGORA_SPICE_DIR")
    if env: return pathlib.Path(env)
    for q in [HERE, *HERE.parents]:
        p = q / "data" / "GRAMSuite.jl" / "GRAM Suite 2.0" / "SPICE"
        if (p / "spk" / "planets" / "de440s.bsp").exists(): return p
    raise FileNotFoundError("SPICE dir not found; set SPACEAGORA_SPICE_DIR")


class SPK:
    """SPK type 2 (Chebyshev, position only) reader: segments[(target, center)] -> list of (et0, et1, init, intlen, rsize, n, data offset)."""
    def __init__(self, path):
        self.f = open(path, "rb"); self.mm = np.memmap(path, dtype="<f8", mode="r")
        h = self.f.read(1024); nd, ni = struct.unpack("<ii", h[8:16]); fward = struct.unpack("<i", h[76:80])[0]
        self.seg = {}; rec = fward
        while rec:
            self.f.seek((rec - 1) * 1024); r = self.f.read(1024); nxt, _, ns = struct.unpack("<ddd", r[:24])
            for k in range(int(ns)):
                s = r[24 + k * 40: 24 + (k + 1) * 40]; et0, et1 = struct.unpack("<dd", s[:16]); tgt, cen, frm, typ, ini, fin = struct.unpack("<6i", s[16:40])
                if typ != 2: continue
                tail = self.mm[fin - 4:fin]; init, intlen, rsize, n = tail[0], tail[1], int(tail[2]), int(tail[3])
                self.seg.setdefault((tgt, cen), []).append((et0, et1, init, intlen, rsize, n, ini - 1))
            rec = int(nxt)

    def pos(self, tgt, cen, et):
        for et0, et1, init, intlen, rsize, n, off in self.seg[(tgt, cen)]:
            if et0 <= et <= et1:
                idx = min(int((et - init) // intlen), n - 1); rec = self.mm[off + idx * rsize: off + (idx + 1) * rsize]
                mid, rad = rec[0], rec[1]; x = (et - mid) / rad; nc = (rsize - 2) // 3
                T = np.empty(nc); T[0] = 1.0; T[1] = x
                for k in range(2, nc): T[k] = 2 * x * T[k - 1] - T[k - 2]
                return np.array([rec[2 + j * nc: 2 + (j + 1) * nc] @ T for j in range(3)]) * 1e3
        raise ValueError(f"no segment for {tgt}/{cen} at et={et}")


_CHAIN = {10: [(10, 0)], 301: [(301, 3), (3, 0)], 399: [(399, 3), (3, 0)], 299: [(299, 2), (2, 0)], 4: [(4, 0)], 499: [(4, 0)]}
_spk = None


def body_pos_ssb(body, et):
    global _spk
    if _spk is None: _spk = SPK(spice_dir() / "spk" / "planets" / "de440s.bsp")
    return sum((_spk.pos(t, c, et) for t, c in _CHAIN[body]), np.zeros(3))


def rel_pos_j2000(body, center, et):
    """Position of `body` relative to `center` in J2000 [m] (NAIF ids; Mars 499 = Mars barycenter, ~cm from the center)."""
    return body_pos_ssb(body, et) - body_pos_ssb(center, et)


def utc_to_et(y, mo, d, h, mi, s):
    """UTC -> ET. TAI-UTC from the leap-second kernel (naif0012.tls: 32 s for 1999-2005, 35 s from 2012-07, 37 s from 2017); ET - TAI = 32.184 s (periodic term < 2 ms ignored)."""
    jd = lambda Y, M, D: 367 * Y - (7 * (Y + (M + 9) // 12)) // 4 + (275 * M) // 9 + D + 1721013.5
    days = jd(y, mo, d) - 2451545.0
    ymd = (y, mo)   # TAI-UTC steps (SOURCED naif0012.tls): 2006-01 33, 2009-01 34, 2012-07 35, 2015-07 36, 2017-01 37
    tai_utc = 37.0 if ymd >= (2017, 1) else 36.0 if ymd >= (2015, 7) else 35.0 if ymd >= (2012, 7) else 34.0 if ymd >= (2009, 1) else 33.0 if ymd >= (2006, 1) else 32.0
    return days * 86400 + h * 3600 + mi * 60 + s + tai_utc + 32.184


_POLES = {  # NAIF id: (RA0, RA rate/century, Dec0, Dec rate/century) deg; SOURCED pck00010.tpc / pck00008.tpc BODY499 (the Mars values; pck00011 has 317.269/54.433 and W0 176.0499).
    # The Mars choice is DERIVED: the F20 pass tables' lon_rad matches W0 = 176.630 (pck00010/00008) to 2e-5 rad, not pck00011's 176.0499. Earth, Venus: pck00011 = pck00010
    399: (0.0, -0.641, 90.0, -0.557), 499: (317.68143, -0.1061, 52.8865, -0.0609), 299: (272.76, 0.0, 67.16, 0.0)}
_PCK_CENTER = {"earth": 399, "mars": 499, "venus": 299}
_SPK_CENTER = {"earth": 399, "mars": 4, "venus": 299}
GM = {10: 132712440023.310e9, 301: 4902.799e9}   # SOURCED de_403_masses.tpc BODY10_GM, BODY301_GM (km^3/s^2 -> m^3/s^2)


def body_frame(planet, et):
    """3x3 matrix whose columns are the body-equator frame axes in J2000 (v_J2000 = M @ v_body)."""
    ra0, rar, d0, ddr = _POLES[_PCK_CENTER[planet]]; T = et / (86400 * 36525); ra = math.radians(ra0 + rar * T); dec = math.radians(d0 + ddr * T)
    z = np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    node = np.array([-z[1], z[0], 0.0]); x = node / np.linalg.norm(node); y = np.cross(z, x)
    return np.column_stack([x, y, z])


def third_body_in_frame(planet, body, et):
    """Position of the perturber (10 Sun, 301 Moon) relative to the planet center in the body-equator frame [m]."""
    return body_frame(planet, et).T @ rel_pos_j2000(body, _SPK_CENTER[planet], et)
