# Verification data outside this repository

Two kinds of verification data are kept out of this public repository and
synced into gitignored paths under `data/telemetry/` on demand:

- **Simulator reference sets** (Basilisk and GMAT runs of the parity matrix,
  about half a gigabyte) and **public flight sets** (ESA OPS-SAT-1 re-entry
  telemetry) live in the lab-org repository
  `Space-FALCON-Lab/spaceagora-verification-data`. Lab members have access
  through the organisation. Its `manifest.toml` records provenance, the
  generator scripts, and per-file checksums; the fetch script verifies them.
- **Flight telemetry** provided under data-sharing agreements
  (`data/telemetry/CYGNSS/`, `data/telemetry/GRIFEX/`) lives in a separate
  access-restricted repository with per-person access granted by the PI. It
  must never be committed here. The CYGNSS preparation and comparison scripts
  live beside that data (`CYGNSS/tools/` there, `data/telemetry/CYGNSS/tools/`
  after a fetch) because their inputs and outputs are not public.

The public repository keeps only the small truth files that CI grades (Odyssey,
VEx, GMAT Earth aerobraking; under 3 MB) plus the verification harness itself.
`ci_tracked_feather_registry_gate` fails on any tracked feather that neither
`test/telemetry_benchmark_manifest.toml` nor `data/assets_manifest.toml` names.

## Getting the data

```bash
scripts/dev/fetch_private_telemetry.sh references   # simulator reference sets
scripts/dev/fetch_private_telemetry.sh telemetry    # flight telemetry (restricted)
scripts/dev/fetch_private_telemetry.sh              # both; a source you cannot reach is skipped
```

The references source can be narrowed to named datasets from the manifest, for
example `fetch_private_telemetry.sh references basilisk_full_1Ms`. Repository
locations and cache directories can be overridden with the
`SPACEAGORA_VERIFICATION_DATA_*` and `SPACEAGORA_PRIVATE_TELEMETRY_*`
environment variables documented at the top of the script.

## Behavior when data is absent

Code that consumes this data must degrade cleanly rather than fail:

- The Basilisk parity testsets and the CYGNSS testsets in
  `test/gmat_scenario_matrix.jl` skip (with a pointer to the fetch script)
  when their directories are missing.
- Study scripts under `benchmarks/studies/` that require it document the
  required files in their READMEs and error with a clear message.

Public CI never has this data; anything wired into CI must go through a
presence guard of this kind.

## Handling rules

- Never commit these files, attach them to issues/PRs, or copy them into other
  repositories. Flight telemetry additionally must not be re-shared: access is
  granted per person by the PI, and clones are not to be passed on.
- Keep quantitative results derived from flight telemetry within what has been
  cleared for publication.
