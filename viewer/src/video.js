// Video export: render the animation frame by frame at a fixed simulated
// cadence, encode with the browser's VideoEncoder (H.264) and mux into an
// MP4 in memory (mp4-muxer), then hand the file to the user. Frames are
// rendered deterministically through `viewer.renderAt(t)`, so the output is
// independent of the machine's frame rate and of the playback speed the
// page was showing. Where WebCodecs H.264 is not available the fallback
// records the live canvas with MediaRecorder into WebM.
//
// Downloads are what the page starts itself; a host that blocks page-started
// downloads (the claude.ai artifact sandbox does) will refuse them, in which
// case the page says so and the offline/standalone page is the place to
// record.
import { Muxer, ArrayBufferTarget } from 'mp4-muxer';
import { elapsedString } from 'viewer/timeline.js';

const PRESETS = [
  { label: 'canvas size', width: 0, height: 0 },
  { label: '1280 × 720', width: 1280, height: 720 },
  { label: '1920 × 1080', width: 1920, height: 1080 },
  { label: '2560 × 1440', width: 2560, height: 1440 },
];

function even(n) { return Math.max(2, 2 * Math.floor(n / 2)); }

async function h264Codec(width, height, fps) {
  if (typeof VideoEncoder === 'undefined' || typeof VideoFrame === 'undefined') return null;
  // Baseline 3.1 up to 720p30, High 4.0 above; the first supported wins.
  const candidates = width * height <= 1280 * 720 ? ['avc1.42001f', 'avc1.4d001f', 'avc1.640028'] : ['avc1.640028', 'avc1.640032', 'avc1.42001f'];
  for (const codec of candidates) {
    try {
      const { supported } = await VideoEncoder.isConfigSupported({ codec, width, height, framerate: fps, bitrate: 8e6 });
      if (supported) return codec;
    } catch (e) { /* try the next */ }
  }
  return null;
}

function saveBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(url); a.remove(); }, 2000);
}

// Record the animation between two elapsed times into an MP4.
//   viewer: { renderer, timeline, frames, renderAt(t), setRecording(v), resizeTo(w, h) }
//   options: { tStart, tEnd, fps, speed (simulated seconds per video second), width, height, bitrate, filename, onProgress }
export async function recordMp4(viewer, options) {
  const { renderer, timeline } = viewer;
  const fps = Math.max(1, Math.round(options.fps || 30));
  const speed = options.speed > 0 ? options.speed : 60;
  const t0 = Math.max(viewer.frames.tStart, options.tStart ?? viewer.frames.tStart);
  const t1 = Math.min(viewer.frames.tEnd, options.tEnd ?? viewer.frames.tEnd);
  if (!(t1 > t0)) throw new Error('the end time must be after the start time');
  const canvas = renderer.domElement;
  const width = even(options.width || canvas.width), height = even(options.height || canvas.height);
  const codec = await h264Codec(width, height, fps);
  if (!codec) throw new Error('this browser has no H.264 video encoder (WebCodecs); use the WebM recording instead');
  const nFrames = Math.floor((t1 - t0) / (speed / fps)) + 1;
  const dt = speed / fps;
  const wasPlaying = timeline.playing;
  const tBefore = timeline.t;
  const sizeBefore = { w: canvas.clientWidth || canvas.width, h: canvas.clientHeight || canvas.height };
  const pixelRatioBefore = renderer.getPixelRatio();
  viewer.setRecording(true);
  timeline.playing = false;
  renderer.setPixelRatio(1);
  viewer.resizeTo(width, height);
  const muxer = new Muxer({ target: new ArrayBufferTarget(), video: { codec: 'avc', width, height }, fastStart: 'in-memory', firstTimestampBehavior: 'offset' });
  let failure = null;
  const encoder = new VideoEncoder({ output: (chunk, meta) => muxer.addVideoChunk(chunk, meta), error: (e) => { failure = e; } });
  encoder.configure({ codec, width, height, framerate: fps, bitrate: options.bitrate || Math.round(0.12 * width * height * fps / 8) * 8 });
  try {
    for (let k = 0; k < nFrames; k++) {
      if (failure) throw failure;
      if (options.cancelled && options.cancelled()) throw new Error('cancelled');
      const t = Math.min(t1, t0 + k * dt);
      viewer.renderAt(t);
      const frame = new VideoFrame(canvas, { timestamp: Math.round(k * 1e6 / fps), duration: Math.round(1e6 / fps) });
      encoder.encode(frame, { keyFrame: k % (2 * fps) === 0 });
      frame.close();
      while (encoder.encodeQueueSize > 6) await new Promise((r) => setTimeout(r, 4));
      if (k % 5 === 0) {
        options.onProgress && options.onProgress(k + 1, nFrames);
        await new Promise((r) => setTimeout(r, 0));
      }
    }
    await encoder.flush();
    muxer.finalize();
    const blob = new Blob([muxer.target.buffer], { type: 'video/mp4' });
    saveBlob(blob, options.filename || 'spaceagora.mp4');
    options.onProgress && options.onProgress(nFrames, nFrames);
    return { frames: nFrames, bytes: blob.size, width, height, codec };
  } finally {
    try { encoder.state !== 'closed' && encoder.close(); } catch (e) { /* ignore */ }
    renderer.setPixelRatio(pixelRatioBefore);
    viewer.resizeTo(sizeBefore.w, sizeBefore.h);
    timeline.seek(tBefore);
    timeline.playing = wasPlaying;
    viewer.setRecording(false);
  }
}

// Fallback: record the live playback for `seconds` of wall time into WebM.
export function recordWebm(viewer, { seconds = 10, fps = 30, filename = 'spaceagora.webm', onProgress } = {}) {
  return new Promise((resolve, reject) => {
    if (typeof MediaRecorder === 'undefined') { reject(new Error('this browser has no MediaRecorder')); return; }
    const stream = viewer.renderer.domElement.captureStream(fps);
    const mime = ['video/webm;codecs=vp9', 'video/webm;codecs=vp8', 'video/webm'].find((m) => MediaRecorder.isTypeSupported(m));
    if (!mime) { reject(new Error('this browser cannot record WebM')); return; }
    const rec = new MediaRecorder(stream, { mimeType: mime, videoBitsPerSecond: 8e6 });
    const chunks = [];
    rec.ondataavailable = (e) => { if (e.data.size) chunks.push(e.data); };
    rec.onerror = (e) => reject(e.error || new Error('recording failed'));
    rec.onstop = () => {
      const blob = new Blob(chunks, { type: 'video/webm' });
      saveBlob(blob, filename);
      resolve({ bytes: blob.size, mime });
    };
    const wasPlaying = viewer.timeline.playing;
    viewer.timeline.playing = true;
    rec.start(500);
    const started = performance.now();
    const tick = () => {
      const s = (performance.now() - started) / 1000;
      onProgress && onProgress(Math.min(s, seconds), seconds);
      if (s >= seconds) { rec.stop(); viewer.timeline.playing = wasPlaying; } else setTimeout(tick, 250);
    };
    tick();
  });
}

// Dialog: times, cadence, size, then run the export with a progress line.
export function createVideoDialog(viewer, container) {
  const box = document.createElement('div');
  box.className = 'sa-panel sa-video';
  box.hidden = true;
  box.style.cssText = 'left: 50%; top: 50%; transform: translate(-50%, -50%); max-width: 46ch; z-index: 6;';
  const span = viewer.frames.tEnd - viewer.frames.tStart;
  box.innerHTML = `
    <h1>Save a video</h1>
    <dl style="grid-template-columns: max-content 1fr; gap: 4px 10px">
      <dt>from</dt><dd><input type="number" data-role="t0" value="${viewer.frames.tStart.toFixed(0)}" step="any" style="width:12ch"> s</dd>
      <dt>to</dt><dd><input type="number" data-role="t1" value="${viewer.frames.tEnd.toFixed(0)}" step="any" style="width:12ch"> s</dd>
      <dt>simulated s per video s</dt><dd><input type="number" data-role="speed" value="${Math.max(1, Math.round(span / 30))}" step="any" style="width:10ch"></dd>
      <dt>frames per second</dt><dd><select data-role="fps"><option>24</option><option selected>30</option><option>60</option></select></dd>
      <dt>size</dt><dd><select data-role="size">${PRESETS.map((p, i) => `<option value="${i}">${p.label}</option>`).join('')}</select></dd>
      <dt>file</dt><dd><input type="text" data-role="name" value="spaceagora.mp4" style="width:18ch"></dd>
    </dl>
    <p class="sa-hint" data-role="estimate"></p>
    <div class="sa-row" style="gap:8px; margin-top:8px">
      <button data-role="go">Save MP4</button>
      <button data-role="webm" title="Record the live playback for 10 s of wall time (no H.264 needed)">Record WebM (10 s)</button>
      <button data-role="cancel">Cancel</button>
    </div>
    <p class="sa-hint" data-role="status" style="min-height:1.4em"></p>`;
  container.appendChild(box);
  const q = (r) => box.querySelector(`[data-role="${r}"]`);
  const status = q('status'), estimate = q('estimate');
  let cancelled = false, running = false;
  function refreshEstimate() {
    const t0 = Number(q('t0').value), t1 = Number(q('t1').value), speed = Number(q('speed').value), fps = Number(q('fps').value);
    if (!(t1 > t0) || !(speed > 0)) { estimate.textContent = 'end after start, speed positive'; return; }
    const seconds = (t1 - t0) / speed;
    estimate.textContent = `${elapsedString(t1 - t0)} of simulation → ${seconds.toFixed(1)} s of video, ${Math.round(seconds * fps)} frames`;
  }
  for (const r of ['t0', 't1', 'speed', 'fps']) q(r).addEventListener('input', refreshEstimate);
  refreshEstimate();
  q('cancel').addEventListener('click', () => { if (running) cancelled = true; else box.hidden = true; });
  q('go').addEventListener('click', async () => {
    if (running) return;
    running = true; cancelled = false;
    const preset = PRESETS[Number(q('size').value)];
    const name = q('name').value.trim() || 'spaceagora.mp4';
    status.textContent = 'starting…';
    try {
      const result = await recordMp4(viewer, {
        tStart: Number(q('t0').value), tEnd: Number(q('t1').value), speed: Number(q('speed').value), fps: Number(q('fps').value),
        width: preset.width, height: preset.height, filename: name.endsWith('.mp4') ? name : `${name}.mp4`,
        cancelled: () => cancelled,
        onProgress: (k, n) => { status.textContent = `frame ${k} of ${n}`; },
      });
      status.textContent = `saved ${result.frames} frames, ${(result.bytes / 1e6).toFixed(1)} MB (${result.width}×${result.height}, ${result.codec}). If nothing downloaded, this host blocks page-started downloads; use the offline page.`;
    } catch (err) {
      status.textContent = err && err.message === 'cancelled' ? 'cancelled' : `could not save: ${err && err.message ? err.message : err}`;
    } finally { running = false; }
  });
  q('webm').addEventListener('click', async () => {
    if (running) return;
    running = true;
    try {
      const r = await recordWebm(viewer, { seconds: 10, fps: Number(q('fps').value), filename: 'spaceagora.webm', onProgress: (s, n) => { status.textContent = `recording ${s.toFixed(1)} of ${n} s`; } });
      status.textContent = `saved ${(r.bytes / 1e6).toFixed(1)} MB WebM`;
    } catch (err) {
      status.textContent = `could not record: ${err && err.message ? err.message : err}`;
    } finally { running = false; }
  });
  return {
    open() { box.hidden = false; refreshEstimate(); },
    close() { box.hidden = true; },
    get running() { return running; },
  };
}
