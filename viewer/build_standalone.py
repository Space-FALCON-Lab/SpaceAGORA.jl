#!/usr/bin/env python3
"""Build the standalone viewer page: one HTML file that opens from disk, takes a
trajectory file through a form (or opens straight into embedded data) and
draws it with the SpaceAGORA viewer, with no Julia and no simulation.

  python3 viewer/build_standalone.py --out output/spaceagora_viewer.html
  python3 viewer/build_standalone.py --out page.html --data run.csv --planet venus --epoch 1993-05-26T00:00:07Z \
      --model data/models/magellan_nasa_3d_resources.glb --model-scale 1 --model-rotation 0,0,90 --title "Magellan"
  python3 viewer/build_standalone.py --out page.html --data a.csv b.csv c.csv --planet mars   # ensemble, a.csv nominal

Options: --textures {4k,8k,none} (default 4k, all five bodies embedded so any
body can be picked; --planet with --data embeds only that body), --units {m,km},
--dims "2,2,2", --mass 500, --reference ghost.csv --reference-name NAME,
--cdn (three.js from jsdelivr and the modules concatenated, the form the
claude.ai artifact host can show; the default page is fully offline).
"""
import argparse
import base64
import json
import pathlib
import re
import sys
import tomllib

VIEWER = pathlib.Path(__file__).resolve().parent
REPO = VIEWER.parent
TEXTURES = REPO / "data" / "textures"
MODULES = ["data.js", "colormaps.js", "timeline.js", "globe.js", "atmosphere.js", "spacecraft.js", "lod.js", "ensemble.js", "paths.js", "references.js", "video.js", "plots.js", "terrain.js", "ui.js", "main.js", "standalone.js"]
VENDOR = {
    "three": "three.module.js",
    "three/addons/controls/OrbitControls.js": "OrbitControls.js",
    "three/addons/loaders/STLLoader.js": "STLLoader.js",
    "three/addons/loaders/OBJLoader.js": "OBJLoader.js",
    "three/addons/loaders/GLTFLoader.js": "GLTFLoader.js",
    "three/addons/utils/BufferGeometryUtils.js": "BufferGeometryUtils.js",
    "mp4-muxer": "mp4-muxer.mjs",
}
CDN = {
    "three": "https://cdn.jsdelivr.net/npm/three@0.160.1/build/three.module.js",
    "three/addons/controls/OrbitControls.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/controls/OrbitControls.js",
    "three/addons/loaders/STLLoader.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/loaders/STLLoader.js",
    "three/addons/loaders/OBJLoader.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/loaders/OBJLoader.js",
    "three/addons/loaders/GLTFLoader.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/loaders/GLTFLoader.js",
    "three/addons/utils/BufferGeometryUtils.js": "https://cdn.jsdelivr.net/npm/three@0.160.1/examples/jsm/utils/BufferGeometryUtils.js",
    "mp4-muxer": "https://cdn.jsdelivr.net/npm/mp4-muxer@5.1.5/build/mp4-muxer.mjs",
}
MODEL_MIME = {".stl": "model/stl", ".obj": "model/obj", ".glb": "model/gltf-binary", ".gltf": "model/gltf+json"}


def data_url(path, mime):
    return f"data:{mime};base64," + base64.b64encode(pathlib.Path(path).read_bytes()).decode()


def js_data_url(path):
    return data_url(path, "text/javascript")


def script_safe(obj):
    return json.dumps(obj).replace("</", "<\\/")


def texture_table(tier, bodies=None):
    manifest = tomllib.loads((TEXTURES / "manifest.toml").read_text())
    out = {}
    for entry in manifest.get("texture", []):
        body = entry["body"].lower()
        if bodies and body not in bodies:
            continue
        res = entry.get("resolution", "4k").lower()
        path = TEXTURES / entry["file"]
        if not path.is_file():
            continue
        current = out.get(body)
        want = tier == res or (current is None and tier != res)
        if current is not None and current["resolution"] == tier:
            continue
        if want or current is None:
            mime = "image/jpeg" if path.suffix.lower() in (".jpg", ".jpeg") else "image/png"
            out[body] = {"url": data_url(path, mime), "lon_left_deg": float(entry.get("lon_left_deg", -180.0)), "resolution": res,
                         "width": int(entry.get("width", 0)), "height": int(entry.get("height", 0)), "source": entry.get("source", ""), "license": entry.get("license", "")}
    return out


def concatenated_modules():
    parts = []
    for name in MODULES:
        s = (VIEWER / "src" / name).read_text()
        s = re.sub(r"^import [^\n]*from 'viewer/[^\n]*;\n", "", s, flags=re.M)
        s = re.sub(r"^import \* as THREE from 'three';\n", "", s, flags=re.M)
        s = re.sub(r"^import \{ [A-Za-z]+ \} from 'three/addons/[^\n]*;\n", "", s, flags=re.M)
        s = re.sub(r"^import \{ Muxer, ArrayBufferTarget \} from 'mp4-muxer';\n", "", s, flags=re.M)
        s = re.sub(r"^export (function|class|const|let)\b", r"\1", s, flags=re.M)
        if name == "ensemble.js":
            s = s.replace("function viridis(t) { const c = viridisRgb(t);", "function viridisColor(t) { const c = viridis(t);").replace("viridis(hi > lo ?", "viridisColor(hi > lo ?")
        parts.append(f"// ---- viewer/src/{name} ----\n{s}")
    modules = "\n".join(parts)
    assert "import " not in modules.replace("import.meta", ""), "unexpected import left in the concatenated modules"
    return modules


def build_preset(args):
    if not args.data:
        return None
    preset = {"planet": args.planet, "epoch": args.epoch, "lengthUnit": args.units, "maxFrames": args.frames,
              "spacecraft": [{"dims_m": [float(x) for x in args.dims.split(",")], "mass_kg": args.mass}],
              "options": {"title": args.title or pathlib.Path(args.data[0]).stem, "trail_orbits": args.trail}, "models": [], "references": []}
    def table(path):
        text = pathlib.Path(path).read_text()
        return {"json": json.loads(text)} if path.lower().endswith(".json") else {"csv": text}
    if len(args.data) == 1:
        preset.update(table(args.data[0]))
    else:
        preset["samples"] = [dict(label=pathlib.Path(p).stem, nominal=(i == 0), **table(p)) for i, p in enumerate(args.data)]
    if args.model:
        ext = pathlib.Path(args.model).suffix.lower()
        if ext not in MODEL_MIME:
            sys.exit(f"unsupported model format {ext}")
        preset["models"].append({"id": 1, "filename": pathlib.Path(args.model).name, "dataUrl": data_url(args.model, MODEL_MIME[ext]),
                                 "scale": args.model_scale, "rotation_deg": [float(x) for x in args.model_rotation.split(",")],
                                 "articulations": json.loads(args.model_articulations) if args.model_articulations else []})
    if args.reference:
        preset["references"].append(dict(name=args.reference_name, lengthUnit=args.units, target=1, **table(args.reference)))
    return preset


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", nargs="*", help="trajectory file(s) to embed (CSV or JSON); several files form an ensemble, the first nominal")
    ap.add_argument("--planet", default="earth", choices=["earth", "mars", "venus", "moon", "titan"])
    ap.add_argument("--epoch", default="2000-01-01T12:00:00Z")
    ap.add_argument("--units", default="m", choices=["m", "km"])
    ap.add_argument("--dims", default="2,2,2")
    ap.add_argument("--mass", type=float, default=500.0)
    ap.add_argument("--frames", type=int, default=4000)
    ap.add_argument("--trail", type=float, default=1.0)
    ap.add_argument("--title", default="")
    ap.add_argument("--model")
    ap.add_argument("--model-scale", type=float, default=1.0)
    ap.add_argument("--model-rotation", default="0,0,0")
    ap.add_argument("--model-articulations", help="JSON list of {region:{min:[..],max:[..]}, axis:[..], angle_deg, pivot:[..]} in model units (as articulation_payload writes)")
    ap.add_argument("--reference")
    ap.add_argument("--reference-name", default="reference")
    ap.add_argument("--textures", default="4k", choices=["4k", "8k", "none"])
    ap.add_argument("--cdn", action="store_true")
    args = ap.parse_args()

    template = (VIEWER / "standalone.html").read_text()
    preset = build_preset(args)
    bodies = {args.planet} if preset else None
    textures = {} if args.textures == "none" else texture_table(args.textures, bodies)
    title = args.title or ("SpaceAGORA viewer" if not preset else f"SpaceAGORA viewer · {preset['options']['title']}")
    if args.cdn:
        importmap = {"imports": dict(CDN)}
        html = template.replace('<script type="importmap">__IMPORTMAP__</script>', '<script type="importmap">' + json.dumps(importmap) + "</script>")
        # the module script imports 'viewer/...' specifiers; replace it with the concatenated modules
        head_imports = ("import * as THREE from 'three';\nimport { OrbitControls } from 'three/addons/controls/OrbitControls.js';\n"
                        "import { STLLoader } from 'three/addons/loaders/STLLoader.js';\nimport { OBJLoader } from 'three/addons/loaders/OBJLoader.js';\n"
                        "import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';\nimport { Muxer, ArrayBufferTarget } from 'mp4-muxer';\n")
        html = html.replace("import { start } from 'viewer/main.js';\nimport { buildPayload, buildEnsemblePayload, parseCsv, fileToDataUrl, modelFormat } from 'viewer/standalone.js';\n",
                            head_imports + concatenated_modules() + "\n")
    else:
        imports = {spec: js_data_url(VIEWER / "vendor" / rel) for spec, rel in VENDOR.items()}
        for name in MODULES:
            imports["viewer/" + name] = js_data_url(VIEWER / "src" / name)
        html = template.replace('<script type="importmap">__IMPORTMAP__</script>', '<script type="importmap">' + json.dumps({"imports": imports}) + "</script>")
    html = html.replace("__TITLE__", title.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    html = html.replace("__TEXTURES__", script_safe(textures)).replace("__PRESET__", script_safe(preset))
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB; textures: {', '.join(sorted(textures)) or 'none'}; data: {'embedded' if preset else 'form'})")


if __name__ == "__main__":
    main()
