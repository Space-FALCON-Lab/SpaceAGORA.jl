"use strict";

const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { pathToFileURL } = require("node:url");
process.env.PLAYWRIGHT_BROWSERS_PATH ??= path.join(__dirname, ".browsers");
const { chromium } = require("playwright");

const output = path.resolve(__dirname, "../results/ground_stations");

async function main() {
  const browser = await chromium.launch({ headless: true });
  try {
    const page = await browser.newPage({ viewport: { width: 1200, height: 1100 } });
    const errors = [];
    page.on("pageerror", error => errors.push(error.message));
    await page.goto(pathToFileURL(path.join(output, "viewer.html")).href);
    assert.equal(await page.locator("#error").textContent(), "");
    assert.equal(await page.locator(".station").count(), 3);
    assert.equal(await page.locator(".satellite").count(), 40);
    const summary = await page.evaluate(() => {
      const data = JSON.parse(document.getElementById("ground-data").textContent);
      const slider = document.getElementById("time");
      const readKeys = selector => [...document.querySelectorAll(selector)].map(e => e.dataset.key).sort();
      const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
      let linksChecked = 0;
      for (let i = 0; i < data.frames.length; i++) {
        slider.value = String(i);
        slider.dispatchEvent(new Event("input", { bubbles: true }));
        const frame = data.frames[i];
        const active = [...frame.active].sort();
        const idle = frame.available.filter(key => !frame.active.includes(key)).sort();
        if (!same(readKeys(".active-link"), active) || !same(readKeys(".contact-link"), idle)) {
          throw new Error(`Rendered connections disagree with the recording at frame ${i}`);
        }
        for (const line of document.querySelectorAll(".active-link, .contact-link")) {
          const link = data.links.find(candidate => candidate.key === line.dataset.key);
          const satellite = document.querySelector(`[data-spacecraft="${link.spacecraft}"]`);
          const station = data.stations[link.station - 1].position_m;
          const expected = [Number(satellite.getAttribute("cx")), Number(satellite.getAttribute("cy")),
            500 + station[0] * 310 / data.radius_m, 400 - station[1] * 310 / data.radius_m];
          ["x1", "y1", "x2", "y2"].forEach((attribute, j) => {
            if (Math.abs(Number(line.getAttribute(attribute)) - expected[j]) > 1e-8) {
              throw new Error(`Wrong endpoint for ${line.dataset.key} at frame ${i}`);
            }
          });
          linksChecked++;
        }
        const cells = [...document.querySelectorAll("#delivery tr")];
        for (let s = 0; s < 3; s++) {
          if (cells[s].cells[2].textContent !== (frame.delivered_bits[s] / 1e6).toFixed(6)) {
            throw new Error(`Wrong delivery counter at frame ${i}, station ${s}`);
          }
        }
      }
      return {
        samples_checked: data.frames.length,
        rendered_links_checked: linksChecked,
        stations: data.stations.map(s => ({ name: s.name, longitude_deg: s.longitude_deg })),
        final_received_bits: data.frames.at(-1).delivered_bits
      };
    });
    await page.locator("#active").uncheck();
    assert.equal(await page.locator(".active-link").count(), 0);
    await page.locator("#contacts").uncheck();
    assert.equal(await page.locator(".contact-link").count(), 0);
    await page.locator("#active").check();
    await page.locator("#contacts").check();
    await page.locator("#time").fill("0");
    await page.locator("#time").dispatchEvent("input");
    await page.locator("#play").click();
    await page.waitForFunction(() => Number(document.getElementById("time").value) > 0);
    await page.locator("#play").click();
    assert.equal(await page.locator("#play").textContent(), "Play samples");
    assert.deepEqual(errors, []);
    await page.locator("#time").fill("0");
    await page.locator("#time").dispatchEvent("input");
    await page.screenshot({ path: path.join(output, "viewer.png"), fullPage: true });
    const evidence = {
      status: "passed", checked_utc: new Date().toISOString(),
      ...summary, toggles_checked: true, playback_checked: true,
      page_errors: errors,
      viewer_sha256: require("node:crypto").createHash("sha256")
        .update(fs.readFileSync(path.join(output, "viewer.html"))).digest("hex")
    };
    fs.writeFileSync(path.join(output, "viewer_validation.json"), JSON.stringify(evidence, null, 2) + "\n");
    fs.writeFileSync(path.resolve(output, "../../VIEWER_CHECK.md"), [
      "# Ground viewer browser validation", "",
      `Passed in local headless Chromium at ${evidence.checked_utc}.`,
      `Checked all **${summary.samples_checked} recorded samples** and **${summary.rendered_links_checked} rendered links**.`,
      "",
      "- Selected and available-idle links match the recorded sets at every sample.",
      "- Every rendered line ends at the correct spacecraft and ground-station coordinates.",
      "- All three displayed cumulative delivery counters match the recorded values.",
      "- Both visibility toggles and playback/pause work; there were no JavaScript page errors.",
      "- This validates the fixture's 2D diagnostic viewer, not the optional 3D renderer.",
      "",
      "[Open viewer](results/ground_stations/viewer.html) | [Screenshot](results/ground_stations/viewer.png) | " +
        "[Machine-readable evidence](results/ground_stations/viewer_validation.json)",
      "",
      `Validated viewer SHA-256: \`${evidence.viewer_sha256}\`. Rerun the browser check after regenerating the viewer.`,
      "", "See [README.md](README.md) for reproduction commands.", ""
    ].join("\n"));
    console.log(JSON.stringify(evidence, null, 2));
  } finally {
    await browser.close();
  }
}

main().catch(error => {
  fs.writeFileSync(path.join(output, "viewer_validation.json"), JSON.stringify({
    status: "failed", checked_utc: new Date().toISOString(), error: error.stack
  }, null, 2) + "\n");
  console.error(error);
  process.exitCode = 1;
});
