"""Tests for viewer/build_cdn_page.py: page and module assembly, metadata escaping,
notices and the external hosts a built page may contact. No network access."""
import base64
import importlib.util
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = pathlib.Path(__file__).resolve().parents[2]
TOOL = REPO / "viewer" / "build_cdn_page.py"
spec = importlib.util.spec_from_file_location("build_cdn_page", TOOL)
cdn = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cdn)


def data_js(text):
    return "data:text/javascript;base64," + base64.b64encode(text.encode()).decode()


def fake_modules(marker="PAGE_MARKER_7c1d", clash=False):
    """Minimal renderer modules with the import forms the real ones use."""
    mods = {name: f"import * as THREE from 'three';\nimport {{ decodeFloat32 }} from 'viewer/data.js';\nexport function {name[:-3]}Fn() {{ return THREE && decodeFloat32; }}\n"
            for name in cdn.MODULE_ORDER}
    mods["data.js"] = f"export function decodeFloat32(b) {{ return b; }}\nexport const MARKER = '{marker}';\n"
    mods["video.js"] = "import { Muxer, ArrayBufferTarget } from 'mp4-muxer';\nimport { elapsedString } from 'viewer/timeline.js';\nexport function createVideoDialog() { return [Muxer, ArrayBufferTarget, elapsedString]; }\n"
    mods["lod.js"] = "import * as THREE from 'three';\nimport { STLLoader } from 'three/addons/loaders/STLLoader.js';\nimport { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';\nexport function createAssemblies() { return [STLLoader, GLTFLoader]; }\n"
    mods["main.js"] = "import * as THREE from 'three';\nimport { OrbitControls } from 'three/addons/controls/OrbitControls.js';\nimport { MARKER } from 'viewer/data.js';\nexport function start(payload, container) { container.textContent = MARKER; return { payload, OrbitControls }; }\n"
    mods["ensemble.js"] = "import { viridis as viridisRgb } from 'viewer/colormaps.js';\nexport function viridis(t) { const c = viridisRgb(t); return c; }\nexport function paint(hi, lo) { return viridis(hi > lo ? 1 : 0); }\n"
    mods["colormaps.js"] = "export function viridis(t) { return [t, t, t]; }\n"
    if clash:
        mods["plots.js"] = "export function createPlotPanel() {}\nexport function decodeFloat32(x) { return x; }\n"
    return mods


def fake_page(payload, modules, importmap_extra=None):
    imports = {"three": data_js("export const REVISION='160';"), "mp4-muxer": data_js("export class Muxer {}; export class ArrayBufferTarget {};")}
    imports.update({"viewer/" + name: data_js(text) for name, text in modules.items()})
    imports.update(importmap_extra or {})
    payload_text = json.dumps(payload).replace("</", "<\\/")
    return ('<!doctype html><html><head><meta charset="utf-8"><title>src</title>\n'
            f'<script type="importmap">{json.dumps({"imports": imports})}</script>\n</head><body>\n<div id="viewer"></div>\n'
            f'<script>\nwindow.SPACEAGORA_VIEWER = {payload_text};\n</script>\n'
            "<script type=\"module\">\nimport { start } from 'viewer/main.js';\nstart(window.SPACEAGORA_VIEWER, document.getElementById('viewer'));\n</script>\n</body></html>\n")


def payload_fixture():
    return {"scene": {"planet": {"name": "Earth", "texture": "earth"}, "epoch": {"utc": "2020-01-01T00:00:00.000Z"}, "spacecraft": [{"id": 1, "name": "craft </script> one"}]},
            "frames": {"count": 3, "sats": 1, "source_rows": 9, "stride_rows": 3},
            "textures": {"earth": {"resolution": "4k", "source": "NASA Blue Marble", "license": "public-domain (NASA)", "url": "data:image/jpeg;base64,AAAA"}},
            "models": {}, "references": [], "options": {"title": "fixture"}}


def module_script(page):
    match = re.search(r'<script type="module">\n(.*?)\n</script>', page, re.S)
    return match.group(1)


class Assembly(unittest.TestCase):
    def test_modules_come_from_the_page_and_vendor_imports_are_reemitted(self):
        page = cdn.build(fake_page(payload_fixture(), fake_modules()), "T", "H", "O", "S", None)
        script = module_script(page)
        self.assertIn("PAGE_MARKER_7c1d", script)
        head = script.split("// ---- viewer/src/data.js ----")[0]
        self.assertIn("import * as THREE from 'three';", head)
        self.assertIn("import { OrbitControls } from 'three/addons/controls/OrbitControls.js';", head)
        self.assertIn("import { STLLoader } from 'three/addons/loaders/STLLoader.js';", head)
        self.assertIn("import { Muxer, ArrayBufferTarget } from 'mp4-muxer';", head)
        body = script.split("// ---- viewer/src/data.js ----", 1)[1]
        self.assertIsNone(cdn.STATIC_IMPORT.search(body))
        self.assertNotIn("\nexport ", body)
        self.assertIn("function viridisColor(t) { const c = viridis(t);", body)
        self.assertIn("start(window.SPACEAGORA_VIEWER", script)
        importmap = json.loads(re.search(r'<script type="importmap">(.*?)</script>', page, re.S).group(1))["imports"]
        self.assertTrue(all(url.startswith("https://cdn.jsdelivr.net/npm/") for url in importmap.values()))
        self.assertNotIn("data:", json.dumps(importmap))
        self.assertIn("three-gpu-pathtracer", importmap)

    def test_no_pathtracer_and_repo_modules(self):
        page = cdn.build(fake_page(payload_fixture(), fake_modules()), "T", "H", "O", "S", None, pathtracer=False, modules="repo")
        importmap = json.loads(re.search(r'<script type="importmap">(.*?)</script>', page, re.S).group(1))["imports"]
        self.assertNotIn("three-gpu-pathtracer", importmap)
        self.assertNotIn("three-mesh-bvh", importmap)
        script = module_script(page)
        self.assertNotIn("PAGE_MARKER_7c1d", script)  # repo sources, not the page's
        for name in cdn.MODULE_ORDER:
            self.assertIn(f"// ---- viewer/src/{name} ----", script)
        self.assertNotIn("// ---- viewer/src/standalone.js ----", script)
        self.assertIn("function start(payload", script)
        self.assertIsNone(cdn.STATIC_IMPORT.search(script.split("// ---- viewer/src/data.js ----", 1)[1]))

    def test_duplicate_top_level_names_are_rejected(self):
        with self.assertRaises(cdn.BuildError) as raised:
            cdn.build(fake_page(payload_fixture(), fake_modules(clash=True)), "T", "H", "O", "S", None)
        self.assertIn("decodeFloat32", str(raised.exception))
        self.assertIn("data.js", str(raised.exception))
        self.assertIn("plots.js", str(raised.exception))

    def test_unknown_vendor_specifier_is_rejected(self):
        mods = fake_modules()
        mods["plots.js"] = "import { pad } from 'left-pad';\nexport function createPlotPanel() { return pad; }\n"
        with self.assertRaises(cdn.BuildError) as raised:
            cdn.build(fake_page(payload_fixture(), mods), "T", "H", "O", "S", None)
        self.assertIn("left-pad", str(raised.exception))

    def test_missing_module_and_missing_marker_are_rejected(self):
        mods = fake_modules()
        del mods["terrain.js"]
        with self.assertRaises(cdn.BuildError) as raised:
            cdn.build(fake_page(payload_fixture(), mods), "T", "H", "O", "S", None)
        self.assertIn("terrain.js", str(raised.exception))
        with self.assertRaises(cdn.BuildError):
            cdn.build("<html><body>no payload here</body></html>", "T", "H", "O", "S", None)

    def test_real_repo_modules_assemble_and_parse(self):
        page = cdn.build(fake_page(payload_fixture(), fake_modules()), "T", "H", "O", "S", None, modules="repo")
        script = module_script(page)
        node = shutil.which("node")
        if node is None:
            self.skipTest("node is not available for a syntax check")
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "assembled.mjs"
            path.write_text(script)
            result = subprocess.run([node, "--check", str(path)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


class Escaping(unittest.TestCase):
    def test_metadata_is_html_escaped(self):
        nasty = '<b onclick="x">&"quoted"</b>'
        page = cdn.build(fake_page(payload_fixture(), fake_modules()), nasty, nasty, nasty, nasty, "</script><script>alert(1)</script>")
        self.assertNotIn(nasty, page)
        self.assertIn("&lt;b onclick=&quot;x&quot;&gt;&amp;&quot;quoted&quot;&lt;/b&gt;", page)
        self.assertNotIn("</script><script>alert(1)</script>", page)
        self.assertIn("&lt;/script&gt;&lt;script&gt;alert(1)&lt;/script&gt;", page)

    def test_payload_text_is_kept_verbatim_and_script_safe(self):
        source = fake_page(payload_fixture(), fake_modules())
        page = cdn.build(source, "T", "H", "O", "S", None)
        source_payload = re.search(r"window\.SPACEAGORA_VIEWER = (.*?);\n</script>", source, re.S).group(1)
        built_payload = re.search(r"window\.SPACEAGORA_VIEWER = (.*?);\n", page, re.S).group(1)
        self.assertEqual(source_payload, built_payload)
        self.assertIn("craft <\\/script> one", built_payload)
        script_blocks = re.findall(r"<script[^>]*>(.*?)</script>", page, re.S)
        self.assertTrue(all("</script" not in block for block in script_blocks))

    def test_raw_close_tag_in_payload_is_rejected(self):
        source = fake_page(payload_fixture(), fake_modules()).replace("craft <\\/script> one", "craft </script> one")
        with self.assertRaises((cdn.BuildError, AttributeError)):
            cdn.build(source, "T", "H", "O", "S", None)


class NoticesAndHosts(unittest.TestCase):
    def external_hosts(self, page):
        return sorted(set(re.findall(r"https?://([A-Za-z0-9.-]+)[/\"]", page)))

    def test_default_page_lists_only_jsdelivr(self):
        page = cdn.build(fake_page(payload_fixture(), fake_modules()), "T", "H", "O", "S", None)
        self.assertEqual(self.external_hosts(page), ["cdn.jsdelivr.net"])
        self.assertIn("three.js r160 (MIT", page)
        self.assertIn("mp4-muxer 5.1.5 (MIT", page)
        self.assertIn("Earth texture: NASA Blue Marble (public-domain (NASA))", page)
        self.assertIn("Built-in library and font hosts: cdn.jsdelivr.net.", page)
        self.assertNotIn("fonts.googleapis.com", page)
        self.assertNotIn("Content-Security-Policy", page)

    def test_fonts_option_adds_the_font_hosts_and_says_so(self):
        page = cdn.build(fake_page(payload_fixture(), fake_modules()), "T", "H", "O", "S", None, fonts=True)
        self.assertEqual(self.external_hosts(page), ["cdn.jsdelivr.net", "fonts.googleapis.com", "fonts.gstatic.com"])
        self.assertIn("Built-in library and font hosts: cdn.jsdelivr.net, fonts.googleapis.com, fonts.gstatic.com.", page)

    def test_cdn_map_is_the_standalone_builders_pinned_map(self):
        standalone = cdn._standalone_builder()
        self.assertEqual(cdn.cdn_import_map(True), standalone.CDN)
        self.assertTrue(all("three@0.160.1" in url for key, url in standalone.CDN.items() if key.startswith("three") and "pathtracer" not in key and "bvh" not in key))


class CommandLine(unittest.TestCase):
    def test_positional_arguments_match_the_demo_drivers(self):
        with tempfile.TemporaryDirectory() as directory:
            src = pathlib.Path(directory) / "in.html"
            out = pathlib.Path(directory) / "sub" / "artifact.html"
            src.write_text(fake_page(payload_fixture(), fake_modules()))
            result = subprocess.run([sys.executable, str(TOOL), str(src), str(out), "Title <x>", "Heading", "Orbit", "Span", "Foot"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(out.is_file())
            self.assertIn("wrote", result.stdout)
            page = out.read_text()
            self.assertIn("<title>Title &lt;x&gt;</title>", page)
            self.assertIn("<h1>Heading</h1>", page)

    def test_bad_input_fails_with_a_message(self):
        with tempfile.TemporaryDirectory() as directory:
            bad = pathlib.Path(directory) / "bad.html"
            bad.write_text("<html><body>not an exported page</body></html>")
            result = subprocess.run([sys.executable, str(TOOL), str(bad), str(pathlib.Path(directory) / "out.html")], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not an exported SpaceAGORA viewer page", result.stderr)
            self.assertFalse((pathlib.Path(directory) / "out.html").exists())
            missing = subprocess.run([sys.executable, str(TOOL), str(pathlib.Path(directory) / "absent.html"), "x.html"], capture_output=True, text=True)
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("source page not found", missing.stderr)


if __name__ == "__main__":
    unittest.main()
