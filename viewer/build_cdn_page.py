#!/usr/bin/env python3
"""Turn an exported SpaceAGORA viewer page into its external-host variant.

  python3 viewer/build_cdn_page.py simulation_results_viewer.html artifact.html
  python3 viewer/build_cdn_page.py in.html out.html "Title" "Heading" "Orbit note" "Span note" "Footer text"

The default exported page is fully self-contained: three.js, the renderer
modules, the textures and the trajectory are all embedded as data URLs, and it
opens from disk. Some hosts (for example the claude.ai artifact host) refuse
`data:` script URLs. This tool rewrites the page so that the pinned three.js,
its loaders and mp4-muxer are fetched from cdn.jsdelivr.net, and the renderer
modules are concatenated into one inline module script. The trajectory,
textures and models stay embedded exactly as exported.

What the built page requests from the network, and nothing else:
  - https://cdn.jsdelivr.net/npm/...   three.js r160.1, its addons, mp4-muxer 5.1.5
                                        and, unless --no-pathtracer, the optional
                                        path-traced lighting packages
  - https://fonts.googleapis.com/...   only with --fonts (IBM Plex; the default
                                        page uses system fonts)
The page lists these hosts and the third-party notices in its footer. It does
not add a Content-Security-Policy and does not work around any host policy.

Metadata arguments are HTML-escaped. The renderer modules are taken from the
page itself (the same code the page was exported and checked with); pass
--modules repo to use the checked-out viewer/src instead.
"""
import argparse
import base64
import html
import importlib.util
import json
import pathlib
import re
import sys

VIEWER = pathlib.Path(__file__).resolve().parent

# Renderer modules in dependency order; standalone.js belongs to the form page only.
MODULE_ORDER = ["data.js", "colormaps.js", "timeline.js", "globe.js", "atmosphere.js", "spacecraft.js", "lod.js",
                "ensemble.js", "paths.js", "references.js", "groundtrack.js", "video.js", "plots.js", "terrain.js",
                "plumes.js", "dust.js", "lighting.js", "ui.js", "main.js"]
PATHTRACER_SPECIFIERS = ("three-gpu-pathtracer", "three-mesh-bvh", "three/examples/jsm/postprocessing/Pass.js")
PAYLOAD_MARKER = re.compile(r"\nwindow\.SPACEAGORA_VIEWER = (.*?);\n</script>", re.S)
IMPORTMAP_TAG = re.compile(r'<script type="importmap">(.*?)</script>', re.S)
DATA_JS_PREFIX = "data:text/javascript;base64,"
VIEWER_IMPORT = re.compile(r"^import [^\n]*from 'viewer/[^\n]*;\n", re.M)
VENDOR_IMPORT = re.compile(r"^import (\* as [A-Za-z_$][\w$]*|\{[^}]*\}) from '([^']+)';\n", re.M)
EXPORT_PREFIX = re.compile(r"^export (async function|function|class|const|let|var)\b", re.M)
TOPLEVEL_DECL = re.compile(r"^(?:async function|function|class|const|let|var)\s+([A-Za-z_$][\w$]*)", re.M)
STATIC_IMPORT = re.compile(r"^\s*import\s+(?!\(|\.meta)", re.M)


class BuildError(Exception):
    """A page that cannot be converted; the message says what to change."""


def _standalone_builder():
    """main's viewer/build_standalone.py, the single source of the pinned CDN map."""
    spec = importlib.util.spec_from_file_location("spaceagora_build_standalone", VIEWER / "build_standalone.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cdn_import_map(pathtracer=True):
    cdn = dict(_standalone_builder().CDN)
    if not pathtracer:
        for key in PATHTRACER_SPECIFIERS:
            cdn.pop(key, None)
    return cdn


def read_payload(page_html):
    match = PAYLOAD_MARKER.search(page_html)
    if match is None:
        raise BuildError("the input is not an exported SpaceAGORA viewer page: no 'window.SPACEAGORA_VIEWER = ...;' script found")
    text = match.group(1)
    try:
        payload = json.loads(text.replace("<\\/", "</"))
    except json.JSONDecodeError as error:
        raise BuildError(f"the embedded viewer payload is not valid JSON: {error}") from error
    if "</" in text:
        raise BuildError("the embedded payload contains a raw '</' sequence; the exporter must write it as '<\\/'")
    for key in ("scene", "frames", "textures"):
        if key not in payload:
            raise BuildError(f"the embedded payload has no '{key}' block")
    return payload, text


def page_modules(page_html):
    """{name: source} for the viewer/*.js modules embedded in the page's import map."""
    match = IMPORTMAP_TAG.search(page_html)
    if match is None:
        raise BuildError("the input page has no import map; use --modules repo to take the renderer from viewer/src")
    try:
        imports = json.loads(match.group(1))["imports"]
    except (json.JSONDecodeError, KeyError, TypeError) as error:
        raise BuildError(f"the page's import map is not readable: {error}") from error
    modules = {}
    for specifier, url in imports.items():
        if not specifier.startswith("viewer/"):
            continue
        if not url.startswith(DATA_JS_PREFIX):
            raise BuildError(f"module {specifier} is not embedded as a data URL; use --modules repo")
        modules[specifier[len("viewer/"):]] = base64.b64decode(url[len(DATA_JS_PREFIX):]).decode("utf-8")
    return modules


def repo_modules():
    return {name: (VIEWER / "src" / name).read_text() for name in MODULE_ORDER}


def assemble_modules(sources, cdn):
    """Concatenate the renderer modules into one module scope.

    Returns (vendor_imports, body). Intra-viewer imports and export keywords are
    removed; vendor imports are collected and re-emitted once per specifier so
    every name a module used stays bound (three, its addons and mp4-muxer).
    """
    missing = [name for name in MODULE_ORDER if name not in sources]
    if missing:
        raise BuildError("renderer modules missing from the page: " + ", ".join(missing))
    namespace = {}
    named = {}
    parts = []
    for name in MODULE_ORDER:
        text = sources[name]
        text = VIEWER_IMPORT.sub("", text)

        def collect(match):
            clause, specifier = match.group(1), match.group(2)
            if specifier not in cdn:
                raise BuildError(f"{name} imports '{specifier}', which has no pinned CDN location in viewer/build_standalone.py")
            if clause.startswith("* as "):
                alias = clause[len("* as "):]
                if namespace.setdefault(specifier, alias) != alias:
                    raise BuildError(f"'{specifier}' is imported under two namespace aliases")
            else:
                for item in clause.strip("{} \n").split(","):
                    item = item.strip()
                    if item:
                        named.setdefault(specifier, []).append(item) if item not in named.get(specifier, []) else None
            return ""

        text = VENDOR_IMPORT.sub(collect, text)
        text = EXPORT_PREFIX.sub(r"\1", text)
        if name == "ensemble.js":
            # ensemble.js imports colormaps' viridis under an alias and defines its own viridis(); in one scope they clash.
            text = text.replace("function viridis(t) { const c = viridisRgb(t);", "function viridisColor(t) { const c = viridis(t);")
            text = text.replace("viridis(hi > lo ?", "viridisColor(hi > lo ?")
        if STATIC_IMPORT.search(text):
            line = STATIC_IMPORT.search(text).group(0).strip()
            raise BuildError(f"{name} still contains a static import after rewriting ({line!r}); the tool does not know this form")
        parts.append(f"// ---- viewer/src/{name} ----\n{text}")
    declared = {}
    for part in parts:
        name = part.split("\n", 1)[0][len("// ---- viewer/src/"):-len(" ----")]
        for match in TOPLEVEL_DECL.finditer(part):
            ident = match.group(1)
            if ident in declared:
                raise BuildError(f"top-level name {ident!r} is declared in both {declared[ident]} and {name}; "
                                 "the CDN page concatenates the modules into one scope, so rename one of them")
            declared[ident] = name
    lines = []
    for specifier in sorted(set(namespace) | set(named), key=lambda s: (s != "three", s)):
        if specifier in namespace:
            lines.append(f"import * as {namespace[specifier]} from '{specifier}';")
        if specifier in named:
            lines.append("import { " + ", ".join(named[specifier]) + f" }} from '{specifier}';")
    return "\n".join(lines), "\n".join(parts)


def notices(payload, cdn, fonts):
    items = ["three.js r160 (MIT, three.js authors) and its loaders and controls, and mp4-muxer 5.1.5 (MIT, Vanilagy), "
             "loaded from cdn.jsdelivr.net at the pinned versions listed in the import map."]
    if any(key in cdn for key in PATHTRACER_SPECIFIERS):
        items.append("Optional path-traced lighting: three-gpu-pathtracer 0.0.23 and three-mesh-bvh 0.7.8 (MIT), fetched from "
                     "cdn.jsdelivr.net only when that lighting mode is selected.")
    for body, texture in sorted((payload.get("textures") or {}).items()):
        source = texture.get("source") or "unspecified source"
        licence = texture.get("license") or "license not recorded"
        items.append(f"{body.capitalize()} texture: {source} ({licence}); embedded in this page.")
    for model in (payload.get("models") or {}).values() if isinstance(payload.get("models"), dict) else (payload.get("models") or []):
        if isinstance(model, dict) and model.get("source"):
            items.append(f"Model {model['source']}: embedded display geometry; see the repository's data/models attribution.")
    hosts = ["cdn.jsdelivr.net"]
    if fonts:
        hosts += ["fonts.googleapis.com", "fonts.gstatic.com"]
    items.append("Network requests made by this page: " + ", ".join(hosts) + ". Everything else is embedded. "
                 "The default exported page makes no network requests at all.")
    return items


def build(page_html, title, heading, orbit_note, span_note, foot, modules="page", fonts=False, pathtracer=True):
    payload, payload_text = read_payload(page_html)
    scene, frames = payload["scene"], payload["frames"]
    planet = scene["planet"]
    cdn = cdn_import_map(pathtracer)
    sources = page_modules(page_html) if modules == "page" else repo_modules()
    vendor_imports, body = assemble_modules(sources, cdn)
    resolution = ((payload.get("textures") or {}).get(planet.get("texture"), {}) or {}).get("resolution", "none")
    esc = lambda value: html.escape(str(value), quote=True)
    default_foot = ("Drag to orbit, wheel to zoom, Space to pause. Click a spacecraft marker or visible model (or press F) "
                    "to follow it; Planet-fixed holds the body still so the ground track drifts instead.")
    notice_html = "".join(f"<li>{esc(item)}</li>" for item in notices(payload, cdn, fonts))
    font_links = ('<link rel="preconnect" href="https://fonts.googleapis.com">\n'
                  '<link rel="preconnect" href="https://fonts.gstatic.com/" crossorigin>\n'
                  '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap">\n') if fonts else ""
    sans = '"IBM Plex Sans", "Segoe UI", system-ui, sans-serif' if fonts else '"Segoe UI", system-ui, sans-serif'
    mono = '"IBM Plex Mono", "SFMono-Regular", Menlo, monospace' if fonts else '"SFMono-Regular", Menlo, Consolas, monospace'
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
{font_links}<style>
  :root {{
    --ground: #070b12; --panel: #10172380; --line: #24304380; --ink: #e2e8f0; --muted: #8aa0b8; --accent: #f2b950;
    --sans: {sans}; --mono: {mono};
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
  .foot, .notices {{ margin-top: 10px; color: var(--muted); font-size: 12px; max-width: 80ch; }}
  .notices ul {{ margin: 4px 0 0; padding-left: 18px; }}
  @media (prefers-reduced-motion: reduce) {{ #viewer .sa-ui button {{ transition: none; }} }}
  @media (max-width: 480px) {{ #viewer {{ height: 70vh; }} }}
</style>
<script type="importmap">{json.dumps({"imports": cdn})}</script>
</head>
<body>
<div class="run">
  <h1>{esc(heading)}</h1>
  <dl>
    <div class="item"><dt>Body</dt><dd>{esc(planet.get("name", "unknown"))}</dd></div>
    <div class="item"><dt>Epoch</dt><dd>{esc(scene.get("epoch", {}).get("utc", "unknown"))}</dd></div>
    <div class="item"><dt>Span</dt><dd>{esc(span_note)}</dd></div>
    <div class="item"><dt>Orbit</dt><dd>{esc(orbit_note)}</dd></div>
    <div class="item"><dt>Rows</dt><dd>{esc(frames.get("source_rows", "?"))} saved, {esc(frames.get("count", "?"))} embedded (every {esc(frames.get("stride_rows", "?"))})</dd></div>
    <div class="item"><dt>Texture</dt><dd>{esc(resolution)}</dd></div>
  </dl>
</div>
<div id="viewer"><div id="sa-loading">Loading three.js from cdn.jsdelivr.net and the scene…</div></div>
<script>
window.addEventListener('error', (e) => {{ const el = document.getElementById('sa-loading'); if (el) el.textContent = 'Viewer failed to start: ' + (e.message || e.error || 'script error'); }});
window.addEventListener('unhandledrejection', (e) => {{ const el = document.getElementById('sa-loading'); if (el) el.textContent = 'Viewer failed to start: ' + (e.reason && e.reason.message ? e.reason.message : e.reason); }});
</script>
<p class="foot">{esc(foot or default_foot)}</p>
<div class="notices"><strong>Third-party notices and network use</strong><ul>{notice_html}</ul></div>
<script>
window.SPACEAGORA_VIEWER = {payload_text};
window.SPACEAGORA_VIEWER.manualStart = true;
</script>
<script type="module">
{vendor_imports}
{body}
const loading = document.getElementById('sa-loading');
try {{
  window.spaceagoraViewer = start(window.SPACEAGORA_VIEWER, document.getElementById('viewer'));
  loading.remove();
}} catch (err) {{
  loading.textContent = 'Viewer failed to start: ' + (err && err.message ? err.message : err);
  console.error(err);
}}
</script>
</body>
</html>
"""


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", help="exported viewer page (simulation_results_viewer.html)")
    parser.add_argument("out", help="page to write")
    parser.add_argument("title", nargs="?", default="SpaceAGORA viewer")
    parser.add_argument("heading", nargs="?", default="SpaceAGORA run")
    parser.add_argument("orbit", nargs="?", default="see the run configuration")
    parser.add_argument("span", nargs="?", default="see the run configuration")
    parser.add_argument("foot", nargs="?", default=None)
    parser.add_argument("--modules", choices=("page", "repo"), default="page", help="renderer modules from the page (default) or from viewer/src")
    parser.add_argument("--fonts", action="store_true", help="load IBM Plex from Google Fonts (adds fonts.googleapis.com and fonts.gstatic.com requests)")
    parser.add_argument("--no-pathtracer", action="store_true", help="omit the optional path-traced lighting packages from the import map")
    args = parser.parse_args(argv)
    source = pathlib.Path(args.source)
    if not source.is_file():
        parser.exit(2, f"error: source page not found: {source}\n")
    try:
        page = build(source.read_text(), args.title, args.heading, args.orbit, args.span, args.foot,
                     modules=args.modules, fonts=args.fonts, pathtracer=not args.no_pathtracer)
    except (BuildError, UnicodeDecodeError, OSError) as error:
        parser.exit(2, f"error: {error}\n")
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB); renderer modules from {args.modules}; "
          f"network hosts: cdn.jsdelivr.net{', fonts.googleapis.com, fonts.gstatic.com' if args.fonts else ''}")


if __name__ == "__main__":
    main()
