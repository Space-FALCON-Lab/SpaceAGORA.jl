# Review record: documentation status and output cadence

## Scope

Correct the main documentation version label and the meaning of the retained
`num_steps_to_save` setting. Author verification, not independent review.

## Changed Files

- `docs/src/index.md`: identify the site as main and link release tags and migration guidance.
- `docs/src/user/upgrading.md`: identify 0.6.1 as unreleased.
- `docs/src/user/outputs.md`: describe output cadence and remove the incorrect buffering claim.
- `src/orchestration/simulation/config/run_settings.jl`: correct the field comment only.

## Findings

The setting is retained in configuration but does not control runtime buffering or
flushing. The active output cadence uses `data_rate`. No executable code changed.
GitHub's latest release was v0.6.0 when the correction was prepared.

## P1 Assessment

No P1 findings in this documentation correction.

## Tests Added/Updated

No behavior tests added for prose and a source comment. The documentation rendered
successfully with external-link checking disabled, and the changed output was
inspected. Strict builds, including a retry on 2026-10-23, failed with HTTP 429 at unchanged external links. The first attempt failed at
the viewer README, independently verified through GitHub. The retry also reported
the MuJoCo user guide and GRAMSuite grid-generation guide; these remain external
link-check failures rather than successful strict verification. The production build policy remains strict and unchanged.
`git diff --check` passes. Hosted documentation CI remains required before merge.

## Residual Risk

External-link rate limiting prevented a complete strict local build. This does
not claim release readiness or document deployment. The unreleased heading must
be updated when the release is actually tagged.

Unresolved P1: No
