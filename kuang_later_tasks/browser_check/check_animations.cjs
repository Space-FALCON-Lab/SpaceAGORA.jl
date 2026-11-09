"use strict";

const assert = require("node:assert/strict");
const { createHash } = require("node:crypto");
const fs = require("node:fs");
const path = require("node:path");
const { fileURLToPath, pathToFileURL } = require("node:url");
process.env.PLAYWRIGHT_BROWSERS_PATH ??= path.join(__dirname, ".browsers");
const { chromium } = require("playwright");

const root = path.resolve(__dirname, "..");
const output = path.join(root, "results/backlog");
const gallery = path.join(root, "animations.html");
const validation = path.join(output, "animation_browser_validation.json");
const cases = [
  "producer_control", "production_competition", "busy_relay",
  "direct_route_control", "production_competition_refined", "busy_relay_refined"
];
const hash = file => createHash("sha256").update(fs.readFileSync(file)).digest("hex");
const writeEvidence = evidence => fs.writeFileSync(validation, JSON.stringify(evidence, null, 2) + "\n");

async function main() {
  writeEvidence({ status: "running", started_utc: new Date().toISOString() });
  const browser = await chromium.launch({ headless: true });
  try {
    const page = await browser.newPage({ viewport: { width: 1400, height: 1100 } });
    const errors = [];
    page.on("pageerror", error => errors.push(error.message));
    page.on("console", message => {
      if (message.type() === "error") errors.push(message.text());
    });
    await page.goto(pathToFileURL(gallery).href);
    assert.deepEqual(await page.locator("video").evaluateAll(videos =>
      videos.map(video => video.dataset.case)), cases);

    const links = await page.locator("a[href], source[src], video[poster]").evaluateAll(elements =>
      elements.map(element => new URL(
        element.getAttribute("href") ?? element.getAttribute("src") ?? element.getAttribute("poster"),
        document.baseURI).href));
    for (const link of links) assert.ok(fs.statSync(fileURLToPath(link)).isFile(), link);

    const videos = [];
    for (const name of cases) {
      const selector = `video[data-case="${name}"]`;
      const video = page.locator(selector);
      const directory = path.join(output, name, "animation");
      const frameMap = path.join(directory, "frames.csv");
      const rows = fs.readFileSync(frameMap, "utf8").trimEnd().split(/\r?\n/);
      assert.equal(rows[0].split(",")[0], "video_frame");
      const frames = rows.length - 1;
      assert.ok(rows.slice(1).every((row, i) => row.startsWith(`${i + 1},`)));
      await video.scrollIntoViewIfNeeded();
      await page.waitForFunction(selector => {
        const element = document.querySelector(selector);
        return element.error || element.readyState >= 1;
      }, selector);
      const metadata = await video.evaluate(element => ({
        duration_s: element.duration, width: element.videoWidth, height: element.videoHeight,
        source: element.currentSrc, error: element.error?.message ?? null
      }));
      assert.equal(metadata.error, null, `${name}: video metadata error`);
      assert.equal(metadata.source, pathToFileURL(path.join(directory, "animation.mp4")).href);
      assert.equal(metadata.width, 1600);
      assert.equal(metadata.height, 1080);
      assert.ok(Math.abs(metadata.duration_s - frames / 30) < 1 / 30, `${name}: duration mismatch`);

      await video.evaluate(async element => {
        element.muted = true;
        await element.play();
      });
      await page.waitForFunction(selector => {
        const element = document.querySelector(selector);
        return !element.paused && element.currentTime > 0.15 && element.readyState >= 2;
      }, selector);
      const playedSeconds = await video.evaluate(element => {
        element.pause();
        return element.currentTime;
      });
      const seeks = [metadata.duration_s / 2, metadata.duration_s - 1 / 60];
      for (const time of seeks) {
        await video.evaluate((element, time) => { element.currentTime = time; }, time);
        await page.waitForFunction(({ selector, time }) => {
          const element = document.querySelector(selector);
          return !element.seeking && element.readyState >= 2 &&
            Math.abs(element.currentTime - time) < 1 / 60;
        }, { selector, time });
      }
      const decodedFrames = await video.evaluate(element => element.getVideoPlaybackQuality().totalVideoFrames);
      assert.ok(decodedFrames > 0, `${name}: no decoded frames`);
      assert.equal(await video.evaluate(element => element.error), null);
      videos.push({
        case: name, frames, fps: 30, ...metadata, played_seconds: playedSeconds,
        seek_times_s: seeks, decoded_frames: decodedFrames,
        video_sha256: hash(path.join(directory, "animation.mp4")),
        frame_map_sha256: hash(frameMap)
      });
    }

    await page.screenshot({ path: path.join(output, "animation_gallery.png"), fullPage: true });
    await page.setViewportSize({ width: 480, height: 900 });
    assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth));
    assert.deepEqual(errors, []);
    const evidence = {
      status: "passed", checked_utc: new Date().toISOString(), browser: browser.version(),
      videos, local_links_checked: links.length, mobile_overflow_checked: true,
      page_errors: errors, gallery_sha256: hash(gallery), checker_sha256: hash(__filename)
    };
    writeEvidence(evidence);
    console.log(JSON.stringify(evidence, null, 2));
  } finally {
    await browser.close();
  }
}

main().catch(error => {
  writeEvidence({ status: "failed", checked_utc: new Date().toISOString(), error: error.stack });
  console.error(error);
  process.exitCode = 1;
});
