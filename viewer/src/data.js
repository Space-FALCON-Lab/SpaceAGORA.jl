// Decoding of the payload the Julia bundler embeds, plus frame interpolation
// and the small geometry helpers the viewer shares.
//
// Payload contract (see docs/architecture/interactive_visualization_plan.md):
//   frames.t_s      Float[N]           elapsed seconds per kept row (Float64 when t_dtype == "f64")
//   frames.pos_km   Float[N*S*3]       inertial position, km, row-major (frame, sat, xyz);
//                                      Float64 when pos_dtype == "f64" (small runs), else Float32
//   frames.vel_kms  Float32[N*S*3]     inertial velocity, km/s (optional)
//   frames.q        Float32[N*S*4]     attitude, scalar-last, inertial->body passive (optional)
//   frames.mass_kg  Float32[N*S]       total mass (optional)
//   frames.link_pose {stride, counts[S], offsets[S], total, data: Float32[N*total]} (optional)
//   frames.sun_dir Float32[N*3]       unit vector planet center -> Sun, inertial (optional)
//   frames.channels [{name, label, unit, digits, log, data: Float32[N*S]}] (optional)
// Every block is base64 of little-endian floats.

export function decodeBytes(b64) {
  if (b64 == null || b64.length === 0) return new Uint8Array(0);
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return bytes;
}

// Blocks arrive as base64 strings from the Julia bundler, or as typed arrays
// (or plain arrays) when a page builds its payload in the browser
// (viewer/src/standalone.js); both forms are accepted everywhere.
export function decodeFloat32(b64) {
  if (b64 instanceof Float32Array) return b64;
  if (ArrayBuffer.isView(b64) || Array.isArray(b64)) return Float32Array.from(b64);
  const bytes = decodeBytes(b64);
  return new Float32Array(bytes.buffer, 0, bytes.length >> 2);
}

export function decodeFloat64(b64) {
  if (b64 instanceof Float64Array) return b64;
  if (ArrayBuffer.isView(b64) || Array.isArray(b64)) return Float64Array.from(b64);
  const bytes = decodeBytes(b64);
  return new Float64Array(bytes.buffer, 0, bytes.length >> 3);
}

export function decodeFloatBlock(b64, dtype) {
  return dtype === 'f64' ? decodeFloat64(b64) : decodeFloat32(b64);
}

// Cubic Hermite position between two samples `a` and `b` (flat xyz offsets
// into P) using the saved velocities V (same offsets, per second) over the
// sample spacing dt, at blend f in [0, 1]. Linear interpolation cuts a chord
// inside the arc: at a 25 s cadence and 4.7 km/s near a periapsis the chord
// sags a few hundred meters mid-segment and the follow camera, which rides
// the spacecraft, shows every other object wobbling by that much. The Hermite
// curve keeps the sag well under a meter at the same cadence. Falls back to
// the chord when velocities are missing, non-finite, or so long against the
// chord that the samples cannot resolve them (a gap or an impulsive burn).
export function hermitePosition(P, a, b, V, dt, f, out, o = 0) {
  const dx = P[b] - P[a], dy = P[b + 1] - P[a + 1], dz = P[b + 2] - P[a + 2];
  if (V && dt > 0) {
    const ux = dt * V[a], uy = dt * V[a + 1], uz = dt * V[a + 2];
    const wx = dt * V[b], wy = dt * V[b + 1], wz = dt * V[b + 2];
    const chord2 = dx * dx + dy * dy + dz * dz;
    const u2 = ux * ux + uy * uy + uz * uz, w2 = wx * wx + wy * wy + wz * wz;
    if (Number.isFinite(u2 + w2) && Math.max(u2, w2) <= 16 * chord2) {
      const f2 = f * f, f3 = f2 * f;
      const h00 = 2 * f3 - 3 * f2 + 1, h10 = f3 - 2 * f2 + f, h01 = -2 * f3 + 3 * f2, h11 = f3 - f2;
      out[o] = h00 * P[a] + h10 * ux + h01 * P[b] + h11 * wx;
      out[o + 1] = h00 * P[a + 1] + h10 * uy + h01 * P[b + 1] + h11 * wy;
      out[o + 2] = h00 * P[a + 2] + h10 * uz + h01 * P[b + 2] + h11 * wz;
      return out;
    }
  }
  out[o] = P[a] + f * dx; out[o + 1] = P[a + 1] + f * dy; out[o + 2] = P[a + 2] + f * dz;
  return out;
}

export class FrameData {
  constructor(frames) {
    this.count = frames.count;
    this.sats = frames.sats;
    this.t = decodeFloatBlock(frames.t_s, frames.t_dtype || 'f32');
    this.pos = decodeFloatBlock(frames.pos_km, frames.pos_dtype || 'f32');
    this.posDtype = frames.pos_dtype || 'f32';
    this.vel = frames.vel_kms ? decodeFloat32(frames.vel_kms) : null;
    this.q = frames.q ? decodeFloat32(frames.q) : null;
    this.mass = frames.mass_kg ? decodeFloat32(frames.mass_kg) : null;
    this.density = frames.density_kg_m3 ? decodeFloat32(frames.density_kg_m3) : null;
    this.heatRate = frames.heat_rate_w_m2 ? decodeFloat32(frames.heat_rate_w_m2) : null;
    this.drag = frames.drag_n ? decodeFloat32(frames.drag_n) : null;
    this.wind = frames.wind_ms ? decodeFloat32(frames.wind_ms) : null;
    // Plume-surface interaction (viewer/dust.js): one Float32 array per quantity,
    // frame-major then spacecraft. Absent unless the run carried a plume effector.
    this.plume = null;
    if (frames.plume) {
      this.plume = {};
      for (const key of Object.keys(frames.plume)) this.plume[key] = decodeFloat32(frames.plume[key]);
    }
    const specs = frames.channels ?? [];
    if (!Array.isArray(specs)) throw new Error('channels must be an array');
    const names = new Set();
    this.channels = specs.map((c) => {
      if (!c || typeof c.name !== 'string' || !c.name || names.has(c.name)) throw new Error('invalid or duplicate channel name');
      names.add(c.name);
      const label = c.label ?? c.name.replaceAll('_', ' '), unit = c.unit ?? '';
      const digits = c.digits ?? 3, log = c.log ?? false;
      if (typeof label !== 'string' || typeof unit !== 'string') throw new Error('channel label and unit must be strings');
      if (!Number.isInteger(digits) || digits < 0 || digits > 100) throw new Error('channel digits must be an integer from 0 through 100');
      if (log !== true && log !== false && log !== 'auto') throw new Error('invalid channel log option');
      let y;
      if (typeof c.data === 'string') {
        const bytes = decodeBytes(c.data);
        if (bytes.length !== 4 * this.count * this.sats) throw new Error('channel data length does not match frames');
        y = new Float32Array(bytes.buffer);
      } else {
        if ((!Array.isArray(c.data) && !ArrayBuffer.isView(c.data)) ||
            !Array.from(c.data).every((v) => typeof v === 'number')) throw new Error('channel data must contain numbers');
        y = decodeFloat32(c.data);
      }
      if (y.length !== this.count * this.sats) throw new Error('channel data length does not match frames');
      if (y.some((v) => !Number.isFinite(v) && !Number.isNaN(v))) throw new Error('channel samples must not be infinite or overflow Float32');
      return { name: c.name, label, unit, digits, log, y };
    });
    // Unit vector from the planet center to the Sun, inertial, one per frame
    // (not per spacecraft); null when the run could not resolve the Sun.
    this.sunDir = frames.sun_dir ? decodeFloat32(frames.sun_dir) : null;
    this.armPose = null;
    if (frames.arm_pose && frames.arm_pose.total > 0) {
      this.armPose = {
        stride: frames.arm_pose.stride,
        counts: frames.arm_pose.counts,
        offsets: frames.arm_pose.offsets,
        total: frames.arm_pose.total,
        data: decodeFloat32(frames.arm_pose.data),
      };
    }
    this.linkPose = null;
    if (frames.link_pose && frames.link_pose.total > 0) {
      this.linkPose = {
        stride: frames.link_pose.stride,
        counts: frames.link_pose.counts,
        offsets: frames.link_pose.offsets,
        total: frames.link_pose.total,
        data: decodeFloat32(frames.link_pose.data),
      };
    }
    this.strideRows = frames.stride_rows;
    this.sourceRows = frames.source_rows;
    this.tStart = this.count > 0 ? this.t[0] : 0;
    this.tEnd = this.count > 0 ? this.t[this.count - 1] : 0;
  }

  // Index i such that t[i] <= time < t[i+1], clamped to [0, count-2]; and the blend factor.
  locate(time) {
    const n = this.count;
    if (n < 2) return { i: 0, f: 0 };
    if (time <= this.t[0]) return { i: 0, f: 0 };
    if (time >= this.t[n - 1]) return { i: n - 2, f: 1 };
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (this.t[mid] <= time) lo = mid; else hi = mid;
    }
    const dt = this.t[hi] - this.t[lo];
    return { i: lo, f: dt > 0 ? (time - this.t[lo]) / dt : 0 };
  }

  // Every spacecraft position at `time` into `out` (S*3 floats): cubic Hermite
  // through the saved velocities, the chord where the run saved none.
  positionsAt(time, out) {
    const { i, f } = this.locate(time);
    const S = this.sats, a = i * S * 3, b = (i + 1) * S * 3;
    if (this.count < 2) { for (let k = 0; k < S * 3; k++) out[k] = this.pos[k]; return out; }
    const dt = this.t[i + 1] - this.t[i];
    for (let s = 0; s < S; s++) hermitePosition(this.pos, a + 3 * s, b + 3 * s, this.vel, dt, f, out, 3 * s);
    return out;
  }

  positionAt(time, sat, out) {
    const { i, f } = this.locate(time);
    const S = this.sats, a = (i * S + sat) * 3, b = ((i + 1) * S + sat) * 3;
    if (this.count < 2) { out[0] = this.pos[a]; out[1] = this.pos[a + 1]; out[2] = this.pos[a + 2]; return out; }
    return hermitePosition(this.pos, a, b, this.vel, this.t[i + 1] - this.t[i], f, out);
  }

  // Velocity from the saved block, or a finite difference of positions when the run had none.
  velocityAt(time, sat, out) {
    const { i, f } = this.locate(time);
    const S = this.sats, a = (i * S + sat) * 3, b = ((i + 1) * S + sat) * 3;
    if (this.vel) {
      if (this.count < 2) { out[0] = this.vel[a]; out[1] = this.vel[a + 1]; out[2] = this.vel[a + 2]; return out; }
      for (let k = 0; k < 3; k++) out[k] = this.vel[a + k] + f * (this.vel[b + k] - this.vel[a + k]);
      return out;
    }
    if (this.count < 2) { out[0] = 0; out[1] = 0; out[2] = 0; return out; }
    const dt = this.t[i + 1] - this.t[i] || 1;
    for (let k = 0; k < 3; k++) out[k] = (this.pos[b + k] - this.pos[a + k]) / dt;
    return out;
  }

  massAt(time, sat) {
    if (!this.mass) return null;
    const { i, f } = this.locate(time);
    const S = this.sats, a = i * S + sat, b = (i + 1) * S + sat;
    if (this.count < 2) return this.mass[a];
    return this.mass[a] + f * (this.mass[b] - this.mass[a]);
  }

  // Scalar diagnostics per (frame index k, spacecraft s). `kind` is one of
  // 'heat_rate' (W/m^2), 'density' (kg/m^3), 'dynamic_pressure' (Pa, from
  // density and speed), 'altitude' (km above the equatorial radius, needs Re),
  // 'speed' (km/s), 'drag' (N), 'wind' (m/s). NaN when the block is absent.
  scalarAtIndex(kind, k, s, Re = 0) {
    const S = this.sats, i = k * S + s, b = i * 3;
    switch (kind) {
      case 'heat_rate': return this.heatRate ? this.heatRate[i] : NaN;
      case 'density': return this.density ? this.density[i] : NaN;
      case 'drag': return this.drag ? this.drag[i] : NaN;
      case 'wind': return this.wind ? Math.hypot(this.wind[b], this.wind[b + 1], this.wind[b + 2]) : NaN;
      case 'speed': return this.vel ? Math.hypot(this.vel[b], this.vel[b + 1], this.vel[b + 2]) : NaN;
      case 'altitude': return Math.hypot(this.pos[b], this.pos[b + 1], this.pos[b + 2]) - Re;
      case 'dynamic_pressure': {
        if (!this.density) return NaN;
        const v = this.vel ? Math.hypot(this.vel[b], this.vel[b + 1], this.vel[b + 2]) * 1000 : NaN;
        return 0.5 * this.density[i] * v * v;
      }
      default: return NaN;
    }
  }

  // Interpolated scalar at time t (linear between the bracketing frames).
  scalarAt(kind, t, s, Re = 0) {
    const { i, f } = this.locate(t);
    if (this.count < 2) return this.scalarAtIndex(kind, 0, s, Re);
    const a = this.scalarAtIndex(kind, i, s, Re), b = this.scalarAtIndex(kind, i + 1, s, Re);
    return a + f * (b - a);
  }

  hasScalar(kind) {
    switch (kind) {
      case 'heat_rate': return !!this.heatRate;
      case 'density': return !!this.density;
      case 'dynamic_pressure': return !!this.density && !!this.vel;
      case 'drag': return !!this.drag;
      case 'wind': return !!this.wind;
      case 'speed': return !!this.vel;
      case 'altitude': return true;
      default: return false;
    }
  }

  // Values are linear in time even when the plot uses a logarithmic axis.
  channelAt(k, time, sat) {
    const c = Number.isInteger(k) ? this.channels[k] : null;
    if (!c || this.count === 0 || !Number.isFinite(time) || !Number.isInteger(sat) || sat < 0 || sat >= this.sats) return NaN;
    if (this.count === 1) return c.y[sat];
    const { i, f } = this.locate(time);
    const a = c.y[i * this.sats + sat], b = c.y[(i + 1) * this.sats + sat];
    // A neighbouring gap must not hide a saved endpoint.
    return f === 0 ? a : f === 1 ? b : a + f * (b - a);
  }

  hasPlume() { return !!this.plume; }

  // One plume quantity ('height_m', 'shear_pa', 'pressure_pa', 'erosion_kg_s',
  // 'eroded_kg', 'ejecta_mps', 'ground_effect_n') for spacecraft `sat`, linearly
  // interpolated at `time`; NaN when the block is absent.
  plumeAt(name, time, sat) {
    const a = this.plume ? this.plume[name] : null;
    if (!a) return NaN;
    const { i, f } = this.locate(time);
    const S = this.sats;
    if (this.count < 2) return a[sat];
    const p = a[i * S + sat], q = a[(i + 1) * S + sat];
    return p + f * (q - p);
  }

  // Attitude quaternion [x, y, z, w] via slerp; null when the run had no orientation state.
  quaternionAt(time, sat, out) {
    if (!this.q) return null;
    const { i, f } = this.locate(time);
    const S = this.sats, a = (i * S + sat) * 4, b = ((i + 1) * S + sat) * 4;
    if (this.count < 2) { for (let k = 0; k < 4; k++) out[k] = this.q[a + k]; return out; }
    return slerp(this.q, a, this.q, b, f, out);
  }

  linkCount(sat) { return this.linkPose ? this.linkPose.counts[sat] : 0; }
  armLinkCount(sat) { return this.armPose ? this.armPose.counts[sat] : 0; }

  // Arm link poses (COM relative to the spacecraft, meters, inertial; quaternion inertial).
  armPosesAt(time, sat, out) { return this._posesAt(this.armPose, time, sat, out); }

  // Link poses of one spacecraft at `time`: [7 * n_links] floats into `out`, or null.
  linkPosesAt(time, sat, out) { return this._posesAt(this.linkPose, time, sat, out); }

  _posesAt(lp, time, sat, out) {
    if (!lp || lp.counts[sat] === 0) return null;
    const { i, f } = this.locate(time);
    const n = lp.counts[sat] * lp.stride;
    const a = i * lp.total + lp.offsets[sat], b = (i + 1) * lp.total + lp.offsets[sat];
    if (this.count < 2) { for (let k = 0; k < n; k++) out[k] = lp.data[a + k]; return out; }
    for (let k = 0; k < n; k += lp.stride) {
      for (let c = 0; c < 3; c++) out[k + c] = lp.data[a + k + c] + f * (lp.data[b + k + c] - lp.data[a + k + c]);
      slerp(lp.data, a + k + 3, lp.data, b + k + 3, f, out.subarray(k + 3, k + 7));
    }
    return out;
  }
}

// Spherical linear interpolation between scalar-last quaternions stored in flat arrays.
export function slerp(A, ia, B, ib, f, out) {
  let ax = A[ia], ay = A[ia + 1], az = A[ia + 2], aw = A[ia + 3];
  let bx = B[ib], by = B[ib + 1], bz = B[ib + 2], bw = B[ib + 3];
  let dot = ax * bx + ay * by + az * bz + aw * bw;
  if (dot < 0) { bx = -bx; by = -by; bz = -bz; bw = -bw; dot = -dot; }
  let s0, s1;
  if (dot > 0.9995) {
    s0 = 1 - f; s1 = f;
  } else {
    const th = Math.acos(Math.min(1, dot)), sn = Math.sin(th);
    s0 = Math.sin((1 - f) * th) / sn; s1 = Math.sin(f * th) / sn;
  }
  let x = s0 * ax + s1 * bx, y = s0 * ay + s1 * by, z = s0 * az + s1 * bz, w = s0 * aw + s1 * bw;
  const n = Math.hypot(x, y, z, w) || 1;
  out[0] = x / n; out[1] = y / n; out[2] = z / n; out[3] = w / n;
  return out;
}

// Attitude for runs without orientation state: body +x along the velocity,
// body +z along nadir with its velocity component removed. Returns the
// scalar-last body-to-inertial quaternion, the same convention as the saved
// `sc{i}_q` columns (see SceneVisualization.velocity_aligned_quaternion).
export function velocityAlignedQuaternion(r, v, out) {
  const speed = Math.hypot(v[0], v[1], v[2]);
  if (speed === 0) { out[0] = 0; out[1] = 0; out[2] = 0; out[3] = 1; return out; }
  const xb = [v[0] / speed, v[1] / speed, v[2] / speed];
  let nz = [-r[0], -r[1], -r[2]];
  let d = nz[0] * xb[0] + nz[1] * xb[1] + nz[2] * xb[2];
  let zb = [nz[0] - d * xb[0], nz[1] - d * xb[1], nz[2] - d * xb[2]];
  let zn = Math.hypot(zb[0], zb[1], zb[2]);
  if (zn <= 1e-12 * Math.max(1, Math.hypot(r[0], r[1], r[2]))) {
    const h = Math.abs(xb[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    d = h[0] * xb[0] + h[1] * xb[1] + h[2] * xb[2];
    zb = [h[0] - d * xb[0], h[1] - d * xb[1], h[2] - d * xb[2]];
    zn = Math.hypot(zb[0], zb[1], zb[2]);
  }
  zb = [zb[0] / zn, zb[1] / zn, zb[2] / zn];
  const yb = [zb[1] * xb[2] - zb[2] * xb[1], zb[2] * xb[0] - zb[0] * xb[2], zb[0] * xb[1] - zb[1] * xb[0]];
  // Columns of the body-to-inertial matrix are the body axes; convert to a quaternion.
  return matrixToQuaternion(
    xb[0], yb[0], zb[0],
    xb[1], yb[1], zb[1],
    xb[2], yb[2], zb[2],
    out,
  );
}

// Row-major 3x3 rotation matrix to scalar-last quaternion.
export function matrixToQuaternion(m11, m12, m13, m21, m22, m23, m31, m32, m33, out) {
  const tr = m11 + m22 + m33;
  let x, y, z, w;
  if (tr > 0) {
    const s = 0.5 / Math.sqrt(tr + 1);
    w = 0.25 / s; x = (m32 - m23) * s; y = (m13 - m31) * s; z = (m21 - m12) * s;
  } else if (m11 > m22 && m11 > m33) {
    const s = 2 * Math.sqrt(1 + m11 - m22 - m33);
    w = (m32 - m23) / s; x = 0.25 * s; y = (m12 + m21) / s; z = (m13 + m31) / s;
  } else if (m22 > m33) {
    const s = 2 * Math.sqrt(1 + m22 - m11 - m33);
    w = (m13 - m31) / s; x = (m12 + m21) / s; y = 0.25 * s; z = (m23 + m32) / s;
  } else {
    const s = 2 * Math.sqrt(1 + m33 - m11 - m22);
    w = (m21 - m12) / s; x = (m13 + m31) / s; y = (m23 + m32) / s; z = 0.25 * s;
  }
  const n = Math.hypot(x, y, z, w) || 1;
  out[0] = x / n; out[1] = y / n; out[2] = z / n; out[3] = w / n;
  return out;
}

// Rotate a vector by the conjugate of a scalar-last quaternion (inertial -> body).
export function rotateByConjugate(q, v, out) {
  // Float32 storage can move a unit quaternion off the unit sphere. Normalize
  // before q* v q so a rotation cannot scale the radius and terrain clearance.
  const n = Math.hypot(q[0], q[1], q[2], q[3]) || 1;
  const x = -q[0] / n, y = -q[1] / n, z = -q[2] / n, w = q[3] / n;
  const ix = w * v[0] + y * v[2] - z * v[1];
  const iy = w * v[1] + z * v[0] - x * v[2];
  const iz = w * v[2] + x * v[1] - y * v[0];
  const iw = -x * v[0] - y * v[1] - z * v[2];
  out[0] = ix * w + iw * -x + iy * -z - iz * -y;
  out[1] = iy * w + iw * -y + iz * -x - ix * -z;
  out[2] = iz * w + iw * -z + ix * -y - iy * -x;
  return out;
}

// Geodetic latitude, longitude (degrees) and altitude (same unit as the
// inputs) for a body-fixed position on an oblate spheroid.
export function geodetic(rBody, Re, Rp) {
  const x = rBody[0], y = rBody[1], z = rBody[2];
  const lon = Math.atan2(y, x);
  const p = Math.hypot(x, y);
  const e2 = 1 - (Rp * Rp) / (Re * Re);
  let lat = Math.atan2(z, p * (1 - e2));
  let N = Re;
  for (let k = 0; k < 5; k++) {
    const s = Math.sin(lat);
    N = Re / Math.sqrt(1 - e2 * s * s);
    lat = Math.atan2(z + e2 * N * s, p);
  }
  const c = Math.cos(lat);
  const alt = Math.abs(c) > 1e-9 ? p / c - N : Math.abs(z) - Rp;
  return { latDeg: lat * 180 / Math.PI, lonDeg: lon * 180 / Math.PI, alt };
}

// Piecewise-slerp lookup in the planet rotation table, extrapolated with the spin rate
// about the body z axis past the last sample.
export class RotationTable {
  constructor(rotation, spinRadS) {
    this.t = Float64Array.from(rotation.t_s);
    this.q = new Float32Array(this.t.length * 4);
    rotation.q_pi.forEach((q, k) => { this.q.set(q, 4 * k); });
    this.spin = spinRadS[2];
  }

  at(time, out) {
    const n = this.t.length;
    if (n === 0) { out[0] = 0; out[1] = 0; out[2] = 0; out[3] = 1; return out; }
    if (time <= this.t[0]) { out.set(this.q.subarray(0, 4)); return out; }
    if (time >= this.t[n - 1]) {
      const last = this.q.subarray(4 * (n - 1), 4 * n);
      const half = 0.5 * this.spin * (time - this.t[n - 1]);
      // q_last (x) q_z(half): extra spin about the body's own pole
      const s = Math.sin(half), c = Math.cos(half);
      const x = last[0], y = last[1], z = last[2], w = last[3];
      out[0] = x * c + y * s;
      out[1] = y * c - x * s;
      out[2] = z * c + w * s;
      out[3] = w * c - z * s;
      return out;
    }
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (this.t[mid] <= time) lo = mid; else hi = mid; }
    const f = (time - this.t[lo]) / (this.t[hi] - this.t[lo]);
    return slerp(this.q, 4 * lo, this.q, 4 * hi, f, out);
  }
}
