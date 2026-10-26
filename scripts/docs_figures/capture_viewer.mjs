// Captures screenshots of a saved SpaceAGORA viewer page in headless Chrome through the
// DevTools protocol. Called by viewer_screenshots.jl; needs Node.js 22 or newer (global
// WebSocket) and Google Chrome (set CHROME to its executable if it is not in the default place).
// Usage: node capture_viewer.mjs <steps.json>
// steps.json: {"url", "width", "height", "steps": [{"wait": ms} | {"eval": js} | {"click": [x, y]}
//              | {"wheel": [x, y, deltaY, count]} | {"drag": [x1, y1, x2, y2]} | {"shot": path}]}
import { spawn } from "node:child_process";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

const cfg = JSON.parse(readFileSync(process.argv[2], "utf8"));
const chromePath = process.env.CHROME ??
  ["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome", "/usr/bin/google-chrome", "/usr/bin/chromium"]
    .find((p) => existsSync(p));
if (!chromePath) throw new Error("Google Chrome not found; set CHROME to its executable.");
const profile = mkdtempSync(join(tmpdir(), "spaceagora-capture-"));
const chrome = spawn(chromePath, [
  "--headless=new", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--no-first-run",
  "--no-default-browser-check", "--hide-scrollbars", `--user-data-dir=${profile}`,
  "--remote-debugging-port=0", `--window-size=${cfg.width},${cfg.height}`, "about:blank",
], { stdio: "ignore" });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const exited = new Promise((r) => chrome.once("exit", r));
const shutdown = async () => {
  chrome.kill();
  await exited;
  rmSync(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
};

let target;
for (let i = 0; i < 100 && !target; i++) {
  await sleep(200);
  try {
    const port = readFileSync(join(profile, "DevToolsActivePort"), "utf8").split("\n")[0];
    const list = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json();
    target = list.find((t) => t.type === "page");
  } catch {}
}
if (!target) { await shutdown(); throw new Error("Chrome did not expose a page target."); }

const ws = new WebSocket(target.webSocketDebuggerUrl);
await new Promise((r) => ws.addEventListener("open", r, { once: true }));
let nextId = 1;
const pending = new Map();
const listeners = [];
ws.addEventListener("message", (ev) => {
  const msg = JSON.parse(ev.data);
  if (msg.id && pending.has(msg.id)) { pending.get(msg.id)(msg); pending.delete(msg.id); }
  else if (msg.method) listeners.forEach((f) => f(msg));
});
const send = (method, params = {}) => new Promise((resolve, reject) => {
  const id = nextId++;
  pending.set(id, (m) => (m.error ? reject(new Error(`${method}: ${m.error.message}`)) : resolve(m.result)));
  ws.send(JSON.stringify({ id, method, params }));
});
const exceptions = [];
listeners.push((m) => m.method === "Runtime.exceptionThrown" &&
  exceptions.push(m.params.exceptionDetails.exception?.description ?? m.params.exceptionDetails.text));

await send("Page.enable");
await send("Runtime.enable");
await send("Emulation.setDeviceMetricsOverride", { width: cfg.width, height: cfg.height, deviceScaleFactor: 1, mobile: false });
const loaded = new Promise((r) => listeners.push((m) => m.method === "Page.loadEventFired" && r()));
await send("Page.navigate", { url: cfg.url });
await loaded;

const mouse = (type, x, y, extra = {}) => send("Input.dispatchMouseEvent", { type, x, y, ...extra });
for (const step of cfg.steps) {
  if (step.wait) await sleep(step.wait);
  if (step.eval) {
    const r = await send("Runtime.evaluate", { expression: step.eval, awaitPromise: true, returnByValue: true });
    if (r.exceptionDetails) throw new Error(`step failed: ${r.exceptionDetails.text}`);
  }
  if (step.click) {
    const [x, y] = step.click;
    await mouse("mouseMoved", x, y);
    await mouse("mousePressed", x, y, { button: "left", buttons: 1, clickCount: 1 });
    await mouse("mouseReleased", x, y, { button: "left", buttons: 0, clickCount: 1 });
  }
  if (step.wheel) {
    const [x, y, dy, n] = step.wheel;
    for (let i = 0; i < (n ?? 1); i++) { await mouse("mouseWheel", x, y, { deltaX: 0, deltaY: dy }); await sleep(60); }
  }
  if (step.drag) {
    const [x1, y1, x2, y2] = step.drag;
    await mouse("mouseMoved", x1, y1);
    await mouse("mousePressed", x1, y1, { button: "left", buttons: 1, clickCount: 1 });
    for (let i = 1; i <= 20; i++) {
      await mouse("mouseMoved", x1 + (x2 - x1) * i / 20, y1 + (y2 - y1) * i / 20, { button: "left", buttons: 1 });
      await sleep(20);
    }
    await mouse("mouseReleased", x2, y2, { button: "left", buttons: 0, clickCount: 1 });
  }
  if (step.shot) {
    const jpeg = step.shot.endsWith(".jpg");
    const r = await send("Page.captureScreenshot", jpeg ? { format: "jpeg", quality: 88 } : { format: "png" });
    writeFileSync(step.shot, Buffer.from(r.data, "base64"));
    console.log("wrote", step.shot);
  }
}
ws.close();
await shutdown();
if (exceptions.length) console.log("page exceptions:\n  " + exceptions.join("\n  "));
