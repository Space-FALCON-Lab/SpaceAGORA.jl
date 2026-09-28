# CYGNSS reconstruction tools

These research tools reproduce recorded states and compare explicit torque hypotheses.
A successful run or a plausible animation does not establish flight accuracy. The
measured-wheel adequacy check currently fails with the retained input configuration;
keep that result visible. Scores and input-derived pages belong with the private
research evidence, not in this public guide.

## 1. Install and run the development tools

Start at the repository root after completing the normal package installation.
The core package keeps its own dependency environment. The two catalogue tools
(SatelliteToolboxSgp4 and SatelliteToolboxTle) have a separate, pinned environment:

```sh
julia --project=scripts/dev -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/dev/run.jl test_cygnss_environment.jl
python3 -m venv .venv-cygnss
.venv-cygnss/bin/python -m pip install -r scripts/dev/requirements-cygnss.txt
.venv-cygnss/bin/python scripts/dev/viewer_demos/test_cygnss_fetch.py
```

On Windows, use `.venv-cygnss\Scripts\python.exe` for the Python commands.
The Julia launcher keeps the root project first and adds the development project
behind it. It does not edit the root Project or Manifest, or bypass the stale-dependency
check. Use this launcher for the build and diagnostic scripts below.

Native GRAM is not required. Both examples use NRLMSISE-00 and space-weather
indices. SPICE supplies the clocks, Earth orientation and planetary positions.
For a fresh checkout, download these public kernels from the
[NASA NAIF archive](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/)
into a separate directory (about 125 MB):

```sh
export SPACEAGORA_SPICE_PATH=/absolute/path/to/cygnss-spice
mkdir -p "$SPACEAGORA_SPICE_PATH/lsk" "$SPACEAGORA_SPICE_PATH/pck" \
  "$SPACEAGORA_SPICE_PATH/spk/planets" "$SPACEAGORA_SPICE_PATH/tf"
curl --fail --location --output "$SPACEAGORA_SPICE_PATH/lsk/naif0012.tls" \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
curl --fail --location --output "$SPACEAGORA_SPICE_PATH/pck/pck00011.tpc" \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc
curl --fail --location --output "$SPACEAGORA_SPICE_PATH/pck/earth_latest_high_prec.bpc" \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
curl --fail --location --output "$SPACEAGORA_SPICE_PATH/spk/planets/de430.bsp" \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
curl --fail --location --output "$SPACEAGORA_SPICE_PATH/tf/earth_assoc_itrf93.tf" \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf
```

These are POSIX-shell commands. On Windows, create the same subdirectories and
save the linked files with those names, then set `$env:SPACEAGORA_SPICE_PATH`.
Keep the kernels with your reconstruction inputs: the Earth orientation file is
updated by NAIF, so record its hash when retaining a run. Missing kernels or
coverage outside their supported dates are errors, not zero rotations.

Set `SPACEAGORA_CYGNSS_DATA` to a data directory outside the public checkout and
`SPACEAGORA_VIEWER_DEMO_OUT` to a private output directory. Environment-variable
examples below use a POSIX shell; PowerShell uses `$env:NAME='value'`.

## 2. Inputs, clocks and frames

### 2.1 Seven-spacecraft navigation data

The public NASA source is [CYGNSS Level 1 version 3.2](https://podaac.jpl.nasa.gov/dataset/CYGNSS_L1_V3.2).
The reconstruction window is June 6 through June 9, 2025. Its roster is FM01,
FM02, FM03, FM04, FM05, FM07 and FM08. The downloader verifies all four days for
all seven spacecraft before writing a complete product. Catalogue elements alone
are not a substitute for missing navigation data.

The PVT positions and velocities are Earth-fixed (WGS84/ITRF93). The builder uses
the full SPICE state transform into J2000, including the derivative of the rotation
for velocity. SGP4 produces TEME; its output takes the separate TEME-to-J2000 path.
Positions use metres and velocities metres per second throughout the saved products.

NASA requires an Earthdata login for the protected arrays. Follow the
[Earthdata DAP access instructions](https://urs.earthdata.nasa.gov/documentation/for_users/data_access/dap_services)
and store a token in a local file. Set `EARTHDATA_TOKEN_FILE` to its path; otherwise
the downloader uses `~/.edl_token`. Never put a token in a command, repository,
issue, or shared log. Requests carrying it and their redirects are restricted to
NASA Earthdata HTTPS hosts.

If Python reports `CERTIFICATE_VERIFY_FAILED`, configure a trusted CA bundle for
the same virtual environment before retrying. On this macOS setup the following
resolved that error while keeping certificate and hostname verification enabled:

```sh
.venv-cygnss/bin/python -m pip install certifi
export SSL_CERT_FILE="$(.venv-cygnss/bin/python -m certifi)"
```

In PowerShell, use
`$env:SSL_CERT_FILE = & .venv-cygnss\Scripts\python.exe -m certifi` after installing
`certifi` with that interpreter. This is a certificate-store correction, not an
Earthdata credential change. Do not disable TLS verification.

```sh
export SPACEAGORA_CYGNSS_DATA=/absolute/path/to/cygnss-data
export EARTHDATA_TOKEN_FILE=/absolute/path/to/local-token-file
.venv-cygnss/bin/python scripts/dev/viewer_demos/fetch_cygnss_l1_states.py \
  --plan-only --out "$SPACEAGORA_CYGNSS_DATA/cygnss_l1_pvt_ecef_20250606_96hr.feather"
.venv-cygnss/bin/python scripts/dev/viewer_demos/fetch_cygnss_l1_states.py \
  --out "$SPACEAGORA_CYGNSS_DATA/cygnss_l1_pvt_ecef_20250606_96hr.feather" \
  --cache "$SPACEAGORA_CYGNSS_DATA/l1_opendap_cache"
julia --project=. scripts/dev/run.jl viewer_demos/build_cygnss_ics.jl
julia --project=. scripts/dev/run.jl viewer_demos/cygnss_constellation.jl
```

`--plan-only` requests public metadata, so it needs no token. The fetch requests
only the eight PVT/time/spacecraft variables, not the full science granules. Cache
entries are decoded before reuse and written atomically. A sidecar
`<output>.manifest.json` records the inventory, response hashes, row counts and
output hash. `data_complete=false` means the protected-data stage has not finished.
Keep that manifest beside the table. An incomplete roster is an error.

The builder writes `constellation_ics_20250606.json`, its catalogue counterpart,
`catalog_vs_flown_20250606.json`, and
`cygnss_constellation_tracks_20250606_96hr.feather` in the same data directory.
Inspect spacecraft coverage, frame-check results and fit residuals before accepting
the constellation page. The fitted initial states are scored on the same arc used
by the fit. This is an in-sample reconstruction, not an independent prediction test.

The default simulation spans all 96 hours for all seven spacecraft, after up to
four fitting runs of the same duration. It needs network access for the public
space-weather indices on first use. Under
`$SPACEAGORA_VIEWER_DEMO_OUT/cygnss_constellation/`, fitting products go in `fit/`,
the final results and scene go in `run/`, and the self-contained viewer is
`page/simulation_results_viewer.html`. Open that HTML file in a browser. The page
contains seven simulated spacecraft and seven translucent navigation references.
Both are sampled more sparsely for display and interpolated between samples.
Use the saved results and printed residuals for quantitative comparisons. A
reference shows "not covered" outside its available navigation interval, including
the final instant if the navigation record ends before the simulation.
The fit adjusts the initial velocity magnitude, so the fitted trajectories need
not start exactly on the original navigation velocities.

The output directory must be absent or empty. Choose a new
`SPACEAGORA_VIEWER_DEMO_OUT` for another run; the demo refuses to overwrite a
previous one. For a short installation check, set
`SPACEAGORA_DEMO_CYGNSS_HOURS=0.02` and `SPACEAGORA_DEMO_CYGNSS_FIT_SMA=0`, then
unset both before the full reconstruction. A short run does not validate the
96-hour example.


### 2.2 Private attitude and command data

The lab-held FM01 exports and mission mass properties have restricted distribution.
Ask the data owner for access through the existing private telemetry repository.
Do not copy them, screenshots of their plots, or generated pages into a public PR.
Required files in `SPACEAGORA_CYGNSS_DATA`:

- `cyg01_slew_adcs.feather`: original `t_rel`, absolute counter `t`, scalar-first
  `q_eci_0..3`, body-rate `w_eci_0..2`, wheel tachometer `Omega_rw_0..2` in rpm.
- `cyg01_slew_pv_eci.feather`: the matching packet grid and `r_eci_0..2`, `v_eci_0..2`.
- `cyg01_adcs_constants.toml`: body inertia and signed wheel axes/inertia,
  with their calibration and provenance retained privately.
- `cyg01_slew_commands.feather`: `t_rel`, `rwCmd_0..2`, `tqDmdCtrl_0..2`, extracted
  without resampling from `actCommands.rwCmd` and `acsLvlhPoint.tqDmdCtrl`.

The flight counter is read as SPICE ET (TDB seconds past J2000), not seconds added
to a UTC DateTime. Retain raw counter, GPS and vehicle-calendar columns to audit
this choice. UTC labels come through the leap-second kernel; the simulation epoch
comes directly from ET. Telemetry quaternions map to the engine's scalar-last
inertial-to-body convention as `(-q1,-q2,-q3,q0)`. Body rates use the verified
positive `w_eci` convention. Repeated navigation fixes are deduplicated and their
uniform cadence is inferred from position/speed; this does not remove uncertainty
in their absolute phase. Keep that distinction in along-track comparisons.

## 3. Two reconstruction paths and their limits

The measured-wheel diagnostic applies `-dH_w/dt - omega × H_w` from splined
tachometer speeds. It checks the drift of inertial `I*omega + H_w` against the
wheel momentum exchanged. Its fifth check remains a failing physical-adequacy
criterion when omitted external momentum is too large. Do not turn it into a skip
or adjust constants to make the check green.

The command reconstruction uses signed **per-wheel** commands, piecewise-linear
interpolation at the actual packet times, and the exact integral of that interpolation
to propagate wheel momentum from one measured initial value. It replaces the
measured-wheel effector and applies `-A*tau_w - omega × H_w` once. The engine still
supplies the body's `-omega × I*omega` term. Time outside the command record is
rejected. The public example assumes command gain one in N m; the private record
retains the unresolved gain/inertia and polarity qualifications.

`rwCmd` is the executed-wheel-command hypothesis to evaluate. `tqDmdCtrl` is a
control-demand comparison: later campaign analysis supersedes its early
interpretation as a body-frame torque, and demand is not executed torque during
the slew. Both runs omit the torque-rod model. Keep their failures and residuals;
command replay is not validation of a closed-loop controller or actuator model.
The older short-window replay score is not an acceptance target.

The existing CYGNSS-1 constants are preserved to reproduce the original diagnostic.
They include fitted values and disputed wheel-scale assumptions. In the command
replay, the initial wheel momentum uses the tachometer and the memo's wheel-inertia
scale, including its per-wheel multipliers; subsequent momentum changes use the
recorded commands at gain one. Those two ingredients use different momentum scales.
The measured-wheel and command-replay attitude scores therefore cannot be compared
as a matched test of which command channel represents the executed torque.

The demo records this limitation in `reconstruction_summary.toml`. No new fit,
wheel-polarity change or mass-properties replacement is performed by this tool.
A comparison using one consistent inertia and command scale must be labelled as a
separate scientific rebaseline, with its inputs resolved by the data owner first.

```sh
julia --project=. scripts/dev/run.jl viewer_demos/cygnss_slew_checks.jl
julia --project=. scripts/dev/run.jl viewer_demos/cygnss_slew.jl
```

The checks command returns nonzero when any convention or adequacy check fails.
The demo can still complete and save that failure in `reconstruction_summary.toml`.
Its `window` and `window_gravity_gradient` directories retain the measured-wheel
comparisons; `command_rwCmd` and `command_tqDmdCtrl` retain the command comparisons.
The private page is under `cygnss_slew/page/`. It explicitly labels the failed
wheels-only diagnostic. `SPACEAGORA_DEMO_CYGNSS_SLEW_FULL_HOUR=0` omits the long
extrapolation; `SPACEAGORA_DEMO_CYGNSS_COMMAND_REPLAY=0` runs the original diagnostic
without a command file. Neither setting changes the physical acceptance criterion.

## 4. Acceptance record

Public CI runs synthetic frame, loader, torque-accounting and download-contract
tests without credentials. Private flight runs must retain input hashes, exact
source revision, constants provenance, sampling choice, comparison window and
execution logs. A complete seven-spacecraft download and a successful page render
are separate checks. Current validation and scientific qualifications accompany the
private review handoff; absent data are recorded as unrun, never as passing.
