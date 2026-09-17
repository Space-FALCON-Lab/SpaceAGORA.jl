// Playback clock: elapsed mission time, play/pause, speed, scrubbing.
export class Timeline {
  constructor(tStart, tEnd, options = {}) {
    this.tStart = tStart;
    this.tEnd = tEnd;
    this.t = tStart;
    this.playing = options.autoplay ?? true;
    this.speed = options.speed ?? 60; // simulated seconds per wall second
    this.loop = options.loop ?? false;   // playback stops at the end; Play from the end starts over
    this.listeners = new Set();
  }

  onChange(fn) { this.listeners.add(fn); return () => this.listeners.delete(fn); }

  emit() { for (const fn of this.listeners) fn(this); }

  tick(dtWall) {
    if (!this.playing) return;
    let t = this.t + this.speed * dtWall;
    if (t > this.tEnd) {
      if (this.loop) t = this.tStart + ((t - this.tStart) % Math.max(1e-9, this.tEnd - this.tStart));
      else { t = this.tEnd; this.playing = false; }
    }
    this.t = t;
    this.emit();
  }

  seek(t) { this.t = Math.min(this.tEnd, Math.max(this.tStart, t)); this.emit(); }
  seekFraction(f) { this.seek(this.tStart + f * (this.tEnd - this.tStart)); }
  fraction() { const d = this.tEnd - this.tStart; return d > 0 ? (this.t - this.tStart) / d : 0; }
  togglePlay() {
    if (!this.playing && this.t >= this.tEnd) this.t = this.tStart;
    this.playing = !this.playing;
    this.emit();
  }
  setSpeed(s) { this.speed = s; this.emit(); }
}

// Wall-clock string for an elapsed time given the sidecar epoch (ISO 8601 UTC).
export function utcString(epochUtc, elapsedSeconds) {
  const base = Date.parse(epochUtc);
  if (Number.isNaN(base)) return '';
  return new Date(base + 1000 * elapsedSeconds).toISOString().replace('.000Z', 'Z');
}

export function elapsedString(seconds) {
  const s = Math.max(0, seconds);
  const d = Math.floor(s / 86400), h = Math.floor((s % 86400) / 3600), m = Math.floor((s % 3600) / 60), sec = s % 60;
  const parts = [];
  if (d > 0) parts.push(`${d}d`);
  parts.push(`${String(h).padStart(2, '0')}:${String(m).padStart(2, '0')}:${sec.toFixed(1).padStart(4, '0')}`);
  return parts.join(' ');
}
