import re, sys, json, pathlib
REPO = pathlib.Path(__file__).resolve().parents[3]
src_html = pathlib.Path(sys.argv[1]).read_text()
out = pathlib.Path(sys.argv[2])
page_title = sys.argv[3] if len(sys.argv) > 3 else "AGORA Earth Orbit"
heading = sys.argv[4] if len(sys.argv) > 4 else "AGORA Earth · examples/AGORA_Earth.jl, 4-day propagation"
orbit_note = sys.argv[5] if len(sys.argv) > 5 else "200.6 km × 50,000 km, i 89.9°"
span_note = sys.argv[6] if len(sys.argv) > 6 else "4 days (≈6 orbits)"
foot = sys.argv[7] if len(sys.argv) > 7 else None
m = re.search(r"window\.SPACEAGORA_VIEWER = (.*?);\n</script>", src_html, re.S)
payload_json = m.group(1)
payload = json.loads(payload_json.replace("<\\/", "</"))
scene, frames = payload["scene"], payload["frames"]
planet = scene["planet"]
# Concatenate the viewer modules into one module scope: strip intra-viewer imports and export keywords.
order = ["data.js", "colormaps.js", "timeline.js", "globe.js", "atmosphere.js", "spacecraft.js", "lod.js", "ensemble.js", "paths.js", "references.js", "video.js", "ui.js", "main.js"]
parts = []
for name in order:
    s = (REPO / "viewer" / "src" / name).read_text()
    s = re.sub(r"^import [^\n]*from 'viewer/[^\n]*;\n", "", s, flags=re.M)
    s = re.sub(r"^import \* as THREE from 'three';\n", "", s, flags=re.M)
    s = re.sub(r"^import \{ OrbitControls \} from 'three/addons/controls/OrbitControls.js';\n", "", s, flags=re.M)
    s = re.sub(r"^import \{ STLLoader \} from 'three/addons/loaders/STLLoader.js';\n", "", s, flags=re.M)
    s = re.sub(r"^import \{ OBJLoader \} from 'three/addons/loaders/OBJLoader.js';\n", "", s, flags=re.M)
    s = re.sub(r"^import \{ GLTFLoader \} from 'three/addons/loaders/GLTFLoader.js';\n", "", s, flags=re.M)
    s = re.sub(r"^import \{ Muxer, ArrayBufferTarget \} from 'mp4-muxer';\n", "", s, flags=re.M)
    s = re.sub(r"^export (function|class|const|let)\b", r"\1", s, flags=re.M)
    if name == "ensemble.js":
        s = s.replace("function viridis(t) { const c = viridisRgb(t);", "function viridisColor(t) { const c = viridis(t);").replace("viridis(hi > lo ?", "viridisColor(hi > lo ?")
    parts.append(f"// ---- viewer/src/{name} ----\n{s}")
modules = "\n".join(parts)
assert "import " not in modules.replace("import.meta", ""), "unexpected import left in bundle"
mission_h = (frames_t := None) or None
epoch = scene["epoch"]["utc"]
count, sats, rows, stride = frames["count"], frames["sats"], frames["source_rows"], frames["stride_rows"]
sc = scene["spacecraft"][0]
importmap = json.dumps({"imports": {
    "three": "https://cdn.jsdelivr.net/npm/three@0.160.1/build/three.module.js",
    "three/addons/controls/OrbitControls.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/controls/OrbitControls.js",
    "three/addons/loaders/STLLoader.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/loaders/STLLoader.js",
    "three/addons/loaders/OBJLoader.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/loaders/OBJLoader.js",
    "three/addons/loaders/GLTFLoader.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/loaders/GLTFLoader.js",
    "three/addons/utils/BufferGeometryUtils.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/utils/BufferGeometryUtils.js",
    "mp4-muxer": "https://cdn.jsdelivr.net/npm/mp4-muxer@5.1.5/build/mp4-muxer.mjs",
}})
payload_js = payload_json  # already "</"-safe from the Julia bundler
html = f"""<meta charset="utf-8">
<title>{page_title}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
  :root {{
    --ground: #070b12; --panel: #10172380; --line: #24304380; --ink: #e2e8f0; --muted: #8aa0b8; --accent: #f2b950;
    --sans: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif; --mono: "IBM Plex Mono", "SFMono-Regular", Menlo, monospace;
  }}
  html, body {{ margin: 0; background: var(--ground); color: var(--ink); }}
  body {{ font-family: var(--sans); font-size: 14px; line-height: 1.45; padding-inline: 16px; padding-block: 12px 16px; overflow-x: hidden; }}
  .run {{ display: flex; flex-wrap: wrap; align-items: baseline; gap: 6px 22px; margin-bottom: 10px; }}
  .run h1 {{ font-family: var(--mono); font-weight: 500; font-size: 15px; letter-spacing: 0.04em; text-transform: uppercase; margin: 0; color: var(--accent); }}
  .run dl {{ display: flex; flex-wrap: wrap; gap: 4px 18px; margin: 0; }}
  .run dt {{ display: inline; color: var(--muted); font-size: 12px; letter-spacing: 0.06em; text-transform: uppercase; }}
  .run dd {{ display: inline; margin: 0 0 0 6px; font-family: var(--mono); font-variant-numeric: tabular-nums; }}
  .run .item {{ white-space: nowrap; }}
  #viewer {{ position: relative; width: 100%; height: calc(100vh - 92px); min-height: 460px; border: 1px solid var(--line); border-radius: 4px; overflow: hidden; background: var(--ground); }}
  #viewer canvas {{ display: block; }}
  #viewer .sa-ui, #viewer .sa-info {{ font-family: var(--mono); }}
  #viewer .sa-info h1 {{ font-family: var(--sans); }}
  #viewer .sa-ui button:focus-visible, #viewer .sa-ui select:focus-visible, #viewer .sa-ui input:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 1px; }}
  #sa-loading {{ position: absolute; inset: 0; display: flex; align-items: center; justify-content: center; color: var(--muted); font-family: var(--mono); pointer-events: none; }}
  .foot {{ margin-top: 10px; color: var(--muted); font-size: 12px; max-width: 70ch; }}
  @media (prefers-reduced-motion: reduce) {{ #viewer .sa-ui button {{ transition: none; }} }}
  @media (max-width: 480px) {{ #viewer {{ height: 70vh; }} }}
</style>
<script type="importmap">{importmap}</script>
<div class="run">
  <h1>{heading}</h1>
  <dl>
    <div class="item"><dt>Body</dt><dd>{planet["name"]}</dd></div>
    <div class="item"><dt>Epoch</dt><dd>{epoch}</dd></div>
    <div class="item"><dt>Span</dt><dd>{span_note}</dd></div>
    <div class="item"><dt>Orbit</dt><dd>{orbit_note}</dd></div>
    <div class="item"><dt>Rows</dt><dd>{rows} saved, {count} embedded (every {stride})</dd></div>
    <div class="item"><dt>Texture</dt><dd>{(payload["textures"].get(planet["texture"], {}) or {}).get("resolution", "none")}</dd></div>
  </dl>
</div>
<div id="viewer"><div id="sa-loading">Loading three.js and the scene…</div></div>
<p class="foot">{foot or "Drag to orbit, wheel to zoom, Space to pause. The trail selector sets how many orbits of history drag behind the spacecraft. Click the marker (or press F) to follow it and see the bus and panels up close; \"Planet-fixed\" holds the body still so the ground track drifts instead."}</p>
<script>
window.SPACEAGORA_VIEWER = {payload_js};
window.SPACEAGORA_VIEWER.manualStart = true;
</script>
<script type="module">
import * as THREE from 'three';
import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';
import {{ STLLoader }} from 'three/addons/loaders/STLLoader.js';
import {{ OBJLoader }} from 'three/addons/loaders/OBJLoader.js';
import {{ GLTFLoader }} from 'three/addons/loaders/GLTFLoader.js';
{modules}
const loading = document.getElementById('sa-loading');
try {{
  window.spaceagoraViewer = start(window.SPACEAGORA_VIEWER, document.getElementById('viewer'));
  loading.remove();
}} catch (err) {{
  loading.textContent = 'Viewer failed to start: ' + (err && err.message ? err.message : err);
  console.error(err);
}}
</script>
"""
out.write_text(html)
print("wrote", out, out.stat().st_size, "bytes; frames", count, "of", rows, "stride", stride, "sats", sats, "links", len(sc["links"]))
