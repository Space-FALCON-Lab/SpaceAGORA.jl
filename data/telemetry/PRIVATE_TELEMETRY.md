# Verification data outside this repository

Two kinds of verification data are kept out of this public repository. When
you need them, a script copies them into Git-ignored folders under
`data/telemetry/`.

- **Simulator reference sets** (Basilisk and GMAT runs of the parity matrix,
  about half a gigabyte) and **public flight sets** (ESA OPS-SAT-1 re-entry
  telemetry) are in the lab repository
  `Space-FALCON-Lab/spaceagora-verification-data`, which lab members can read
  through the organisation. Its `manifest.toml` names each dataset, where it
  lives in that repository and where it goes here, and records its provenance.
  After copying a dataset, the fetch script checks the dataset's `SHA256SUMS`
  file when it has one.
- **Flight telemetry** received under data-sharing agreements
  (`data/telemetry/CYGNSS/`, `data/telemetry/GRIFEX/`) is in a separate
  restricted repository; the PI grants access person by person. It must never
  be committed here. The CYGNSS preparation and comparison scripts live beside
  that data (`CYGNSS/tools/` there, `data/telemetry/CYGNSS/tools/` after a
  fetch) because their inputs and outputs are not public.

The public repository keeps only the small truth files that CI grades
(Odyssey, VEx and the GMAT Earth aerobraking case, under 3 MB together) and the
verification code. `ci_tracked_feather_registry_gate` fails on any tracked
Feather file that is neither named in `test/telemetry_benchmark_manifest.toml`
nor inside a required in-repository asset of `data/assets_manifest.toml`. The
fetched lab folders do not count, so a file force-added under them still fails.

## Getting the data

```bash
scripts/dev/fetch_private_telemetry.sh references   # simulator reference sets
scripts/dev/fetch_private_telemetry.sh telemetry    # flight telemetry (restricted)
scripts/dev/fetch_private_telemetry.sh              # both; a source you cannot reach is skipped
```

You can limit the references to named datasets from the manifest, for example
`fetch_private_telemetry.sh references basilisk_full_1Ms`. The
`SPACEAGORA_VERIFICATION_DATA_*` and `SPACEAGORA_PRIVATE_TELEMETRY_*`
environment variables, listed at the top of the script, change the repositories
and the local cache folders.

## When the data is missing

Code that uses this data must skip cleanly rather than fail. The Basilisk
parity testsets and the CYGNSS testsets in `test/gmat_scenario_matrix.jl` are
skipped, with a pointer to the fetch script, when their files are missing; so
are the same testsets when `benchmarks/studies/scenario_matrix/run_cygnss_96hr.jl`
runs that file. Public CI never has this data, so anything wired into CI needs
a presence check like these.

## Handling rules

- Never commit these files, attach them to issues or pull requests, or copy
  them into other repositories. Do not pass flight telemetry on to anyone:
  the PI grants access to each person, and clones are not to be shared.
- Keep quantitative results derived from flight telemetry within what has been
  cleared for publication.
