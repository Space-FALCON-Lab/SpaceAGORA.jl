"""inspect_surface.py <html> <out-dir>
Headless views of the lit ground: down-sun (opposition) and up-sun views of the
landing site, a wide view for terrain shadows, and the globe's terminator.
"""
import sys, json, asyncio, pathlib
from playwright.async_api import async_playwright

html, outdir = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2])
outdir.mkdir(parents=True, exist_ok=True)
SIZE = (1400, 900)

PLACE = """(([frac, azDeg, dkm, elevDeg, follow]) => {
  const v = window.spaceagoraViewer;
  v.timeline.playing = false;
  v.timeline.seekFraction(frac);
  v.state.select(0);
  v.state.setFollow(follow);
  v.renderAt(v.timeline.t);
  const s = v.lighting.sun;
  const d = s.position.clone().sub(s.target.position).normalize();
  const focus = v.lod.items[0].group.getWorldPosition(new v.camera.position.constructor());
  const up = focus.clone().sub(v.world.position).normalize();
  const sh = d.clone().addScaledVector(up, -d.dot(up)).normalize();     // the sun's azimuth, horizontal
  const across = up.clone().cross(sh).normalize();
  const el = elevDeg * Math.PI / 180, az = azDeg * Math.PI / 180;        // az = 0 puts the camera up-sun (phase near 0)
  const off = sh.clone().multiplyScalar(Math.cos(az) * Math.cos(el)).addScaledVector(across, Math.sin(az) * Math.cos(el))
    .addScaledVector(up, Math.sin(el)).multiplyScalar(dkm);
  v.camera.position.copy(focus).add(off);
  v.controls.target.copy(focus);
  v.camera.lookAt(focus);
  v.controls.update();
  v.renderAt(v.timeline.t);
  return {sunElevDeg: Math.asin(d.dot(up)) * 180 / Math.PI, phaseDeg: Math.acos(Math.max(-1, Math.min(1, d.dot(off.clone().normalize())))) * 180 / Math.PI,
          lighting: v.lighting.status, terrain: v.terrain.modelStatus, exposure: v.lighting.exposure};
})"""

GLOBE = """(([frac, dkm]) => {
  const v = window.spaceagoraViewer;
  v.timeline.playing = false;
  v.timeline.seekFraction(frac);
  v.state.setFollow(false);
  v.renderAt(v.timeline.t);
  const s = v.lighting.sun;
  const d = s.position.clone().sub(s.target.position).normalize();
  const c = v.world.position.clone();
  // 70 degrees off the sun, so the day side, the terminator and the night side are all in frame
  const any = Math.abs(d.z) < 0.9 ? new v.camera.position.constructor(0, 0, 1) : new v.camera.position.constructor(1, 0, 0);
  const perp = any.clone().addScaledVector(d, -any.dot(d)).normalize();
  const dir = d.clone().multiplyScalar(Math.cos(70 * Math.PI / 180)).addScaledVector(perp, Math.sin(70 * Math.PI / 180));
  v.camera.position.copy(c).addScaledVector(dir, dkm);
  v.controls.target.copy(c);
  v.camera.lookAt(c);
  v.controls.update();
  v.renderAt(v.timeline.t);
  return {lighting: v.lighting.status, exposure: v.lighting.exposure};
})"""

async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist", "--allow-file-access-from-files"])
        pg = await b.new_page(viewport={"width": SIZE[0], "height": SIZE[1]})
        logs = []
        pg.on("console", lambda m: logs.append(f"[{m.type}] {m.text}"))
        pg.on("pageerror", lambda e: logs.append(f"[pageerror] {e}"))
        await pg.goto(html.as_uri(), timeout=120000)
        await pg.wait_for_timeout(9000)
        for name, arg in (("downsun", [0.985, 0, 0.35, 22, True]),
                          ("upsun", [0.985, 180, 0.35, 22, True]),
                          ("crater_wide", [0.985, 0, 1.6, 12, True]),
                          ("crater_upsun_wide", [0.985, 180, 1.6, 12, True]),
                          ("lander_shadow", [0.995, 90, 0.05, 14, True]),
                          ("lander_shadow_wide", [0.995, 90, 0.25, 14, True]),
                          ("ground_level", [0.995, 180, 0.045, 2, True]),
                          ("ground_level_cross", [0.995, 90, 0.045, 2, True])):
            print(name, json.dumps(await pg.evaluate(PLACE, arg)))
            await pg.wait_for_timeout(1500)
            await pg.screenshot(path=str(outdir / f"{name}.png"))
        # the page's own follow camera, the view the demo opens on
        print("page_close", json.dumps(await pg.evaluate("""(() => {
          const v = window.spaceagoraViewer;
          v.timeline.playing = false; v.timeline.seekFraction(0.985); v.state.select(0); v.state.setFollow(true);
          return {terrain: v.terrain.modelStatus};
        })()""")))
        await pg.wait_for_timeout(2500)
        await pg.screenshot(path=str(outdir / "page_close.png"))
        print("globe", json.dumps(await pg.evaluate(GLOBE, [0.0, 9000])))
        await pg.wait_for_timeout(2000)
        await pg.screenshot(path=str(outdir / "globe.png"))
        await b.close()
        for l in logs[:25]: print("CONSOLE", l[:300])

asyncio.run(main())
