# Private mission telemetry

Some verification datasets used by SpaceAGORA.jl are flight telemetry provided
under data-sharing agreements that do not permit public redistribution. This
public repository therefore follows a strict split:

- **Public (this repo):** the verification harness, scenario/study code,
  simulator-to-simulator reference data (GMAT/Basilisk examples), and public
  mission products.
- **Private (separate access-restricted repo):** raw and processed flight
  telemetry for `data/telemetry/CYGNSS/` and `data/telemetry/GRIFEX/`. Both
  paths are gitignored here (see `.gitignore`) and must never be committed.

## Getting the data

If you have been granted access to the private repo:

```bash
scripts/dev/fetch_private_telemetry.sh
```

This clones/updates a local cache and syncs mission directories into
`data/telemetry/`. The repo location can be overridden with
`SPACEAGORA_PRIVATE_TELEMETRY_REPO` (owner/name) and the cache location with
`SPACEAGORA_PRIVATE_TELEMETRY_CACHE`.

## Behavior when data is absent

Code that consumes private telemetry must degrade cleanly rather than fail:

- The CYGNSS testsets in `test/gmat_scenario_matrix.jl` skip (with a pointer to
  the fetch script) when the feather files are missing.
- Study scripts under `benchmarks/studies/` that require private data document
  the required files in their READMEs and error with a clear message.

Public CI never has this data; anything wired into CI must go through a
presence guard of this kind.

## Handling rules

- Never commit these files, attach them to issues/PRs, or copy them into other
  repositories - including private forks with broader access.
- Access is granted per person by the PI; do not re-share clones.
- Keep quantitative results derived from the data within what has been cleared
  for publication.
