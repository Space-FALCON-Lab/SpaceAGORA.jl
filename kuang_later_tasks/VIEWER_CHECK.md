# Ground viewer browser validation

Passed in local headless Chromium at 2026-11-09T01:33:56.924Z.
Checked all **3924 recorded samples** and **53675 rendered links**.

- Selected and available-idle links match the recorded sets at every sample.
- Every rendered line ends at the correct spacecraft and ground-station coordinates.
- All three displayed cumulative delivery counters match the recorded values.
- Both visibility toggles and playback/pause work; there were no JavaScript page errors.
- This validates the fixture's 2D diagnostic viewer, not the optional 3D renderer.

[Open viewer](results/ground_stations/viewer.html) | [Screenshot](results/ground_stations/viewer.png) | [Machine-readable evidence](results/ground_stations/viewer_validation.json)

Validated viewer SHA-256: `fb4e5dfc11c4bafdd0cdc3985275fac6120fefce236b64207eb97041afc848ac`. Rerun the browser check after regenerating the viewer.

See [README.md](README.md) for reproduction commands.
