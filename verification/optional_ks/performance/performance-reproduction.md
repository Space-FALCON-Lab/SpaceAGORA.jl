# Reproduce the bounded performance experiment

The accepted measurement source is SpaceAGORA numerical source
c091191a1ae0dc602ae50025889d34d4b2800a96 and KS
05d969c24a6595c2752f7ff947af77c504053623, Julia1.12.7. The later KS public
planet docstring relocation keeps its own source identity; no measurements
are relabelled. Reproduce the numerical experiment on the recorded source,
not an unpinned branch tip. Repeated wall times need not reproduce exactly.

The companion archive contains this directory's original scripts, policy,
source-bound prerequisites, raw calibration/worker records, process monitors,
environments and reviews. Source Git bundles preserve review history separately.
For installation, use an ordinary Git/LFS SpaceAGORA checkout at the recorded
source and its verification/optional_ks/install.jl in a fresh environment.
The historical installer at c091 pins an earlier KS source. After that installer
finishes, explicitly select the measured05d source in the new benchmark-only
environment BEFORE preparing any gate or starting calibration:

```
JULIA_LOAD_PATH=@:@stdlib julia --startup-file=no --project=NEW_ENV -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/gfalcon2/KSPropagator.jl",rev="05d969c24a6595c2752f7ff947af77c504053623")); Pkg.instantiate(); using KSPropagator; println(pkgdir(KSPropagator))'
```

Use the printed directory as KS_INSTALLED_ROOT. This is a separate historical
benchmark environment; the final adapter installer may pin the later reviewed
documentation-only candidate. Normal authorized access to the private KS repository is required. The direct
SpaceAGORA Git-URL Pkg.add route failed; do not substitute it for the verified
checkout/installer route. Do not reuse historical absolute-path Manifest entries
on another machine: the installer resolves a fresh local environment and the
experiment freezes that actual Project/Manifest and resolved versions.

Work in a COPY of the extracted evidence directory, preserving the original raw
archive. Call it EVIDENCE_ROOT. Set JULIA_LOAD_PATH=@:@stdlib,
JULIA_NUM_THREADS=1 and OPENBLAS_NUM_THREADS=1. Supply the actual installed KS
root and SpaceAGORA checkout path to:

```
python3 prepare-performance-gate.py EVIDENCE_ROOT SPACEAGORA_ROOT KS_INSTALLED_ROOT
```

This helper verifies the qualified source hashes and prepares the gate BEFORE
execution. Its original absolute-path predecessor is retained as
prepare-performance-gate-original.py. The parameterized helper was added for
reproduction after the original experiment started; it was not retroactively
included in the measured experiment identity. Never run it or edit any bound
script, policy, input, source, Project/Manifest or prerequisite record while a
calibration or worker is in progress.

Run from EVIDENCE_ROOT using absolute paths for the environment and baseline
argument. Use new output names, never overwrite original attempts. The process
guardian is intentionally tied to the recorded machine's documented dormant
PID55147. On a new machine require no other Julia process; its absence is
accepted. Any unexpected Julia process stops/invalidates this experiment.
Do not kill other users' processes. Record a prospectively reviewed alternate
background policy if different conditions are unavoidable.

```
python3 run-performance-process.py NEW_CALIBRATION_MONITOR.json julia --startup-file=no --project=NEW_ENV performance-calibration.jl ABSOLUTE_EVIDENCE_ROOT/baseline-attempt1/results.toml NEW_CALIBRATION
mkdir NEW_TIMING
python3 run-performance-process.py NEW_TIMING/two_body-dp8-monitor.json julia --startup-file=no --project=NEW_ENV performance-worker.jl NEW_CALIBRATION/calibration.toml two_body dp8 NEW_TIMING/two_body-dp8
python3 run-performance-process.py NEW_TIMING/two_body-ks-monitor.json julia --startup-file=no --project=NEW_ENV performance-worker.jl NEW_CALIBRATION/calibration.toml two_body ks NEW_TIMING/two_body-ks
python3 run-performance-process.py NEW_TIMING/J2_plus_drag-ks-monitor.json julia --startup-file=no --project=NEW_ENV performance-worker.jl NEW_CALIBRATION/calibration.toml J2_plus_drag ks NEW_TIMING/J2_plus_drag-ks
python3 run-performance-process.py NEW_TIMING/J2_plus_drag-dp8-monitor.json julia --startup-file=no --project=NEW_ENV performance-worker.jl NEW_CALIBRATION/calibration.toml J2_plus_drag dp8 NEW_TIMING/J2_plus_drag-dp8
julia --startup-file=no performance-summary.jl NEW_TIMING NEW_SUMMARY.toml
```

Wait for terminal exit0 with a valid guardian after EACH command before the next.
The original package-load/first-use/warm call boundaries, accuracy calibration,
complete86400s endpoint, physical clock/mass/heat gates, all stable cross-pair
position/velocity matches and residual compilation gates remain mandatory.
Inspect every preserved candidate and failure. Three warm repetitions and this
sampled-error match support only this bounded observed-machine comparison.
Substantial recorded swap use, OS noise and the guardian's2s observation
interval preclude a general machine-independent speed claim.

## Preserved retries and timing-quality diagnosis

The actual experiment required a separate new J2KS output directory after its
first attempt was invalidated by a0.01s background CPU-counter advance. The
original hardcoded-layout summary remains unchanged. For the retained artifacts,
use Python3.11+ `performance-retry-summary.py performance-selection.json
NEW_SUMMARY.json`. The selection binds every exact timing/monitor artifact and
its original launch working directory; the prospective plan is separately hashed.
Do not rebuild it by choosing fast repetitions. `build-performance-selection.py`
is specific to the recorded paths and refuses overwrite; new experiments need
their own prospectively declared selection plan and complete raw files.

The J2KS retry passes physical/paired-error checks but has ten long monitor
gaps, including5034s. The additive analysis diagnoses those gaps and withholds
a J2speedratio. Stable-compilation flags do not establish continuous timing
coverage. No cause is assigned and no slow call is dropped. A fresh timing
study should prospectively require continuous observation and record any system
suspension/availability evidence, rather than reuse these wall times as speed
validation. All original measurements remain in the archive.
