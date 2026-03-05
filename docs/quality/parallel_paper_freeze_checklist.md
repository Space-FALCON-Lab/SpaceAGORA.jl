# Parallel Paper Freeze Checklist

This checklist is the repo-specific protocol to freeze an immutable paper checkpoint before changing accuracy-sensitive code.

## Naming Convention

- Checkpoint branch: `checkpoint/parallel-paper-YYYY-MM-DD`
- Tag: `parallel-paper-freeze-YYYY-MM-DD`
- Freeze artifact root: `output/freezes/parallel_paper_YYYY-MM-DD/`

## One-Command Freeze

From repo root:

```bash
./freeze_parallel_paper_checkpoint.sh \
  --freeze-date=2026-03-05 \
  --profile=full \
  --run-script=./run_parallel_paper_protocol.sh \
  --run-flags='passes_main=1; passes_sweeps=1; sweeps=r5_full_smart; cold+warm' \
  --run-outdir=output/performance/paper_protocol_20260304_224630
```

Notes:

- Use `--push-remote=0` to create local-only refs first.
- Use `--dry-run=1` to validate settings without creating archives.

## What Gets Captured

1. Immutable Git refs:
- checkpoint branch at exact commit
- annotated tag with purpose + run metadata

2. Reproducibility input archive:
- `Project.toml`, `Manifest.toml`
- `.AGORA/Project.toml`, `.AGORA/Manifest.toml`
- runtime scripts and performance harness files

3. Results archive:
- selected paper run outdir (`output/performance/paper_protocol_*`)

4. Metadata and figures archives:
- hardware/metadata CSV files found under selected run outdir
- figure files (`png`, `pdf`, `svg`) from `docs/` and run outdir

5. Integrity + inventory:
- `SHA256SUMS.txt`
- `FREEZE_MANIFEST.md` with commit/tag/submodules/machine/commands

## Post-Freeze Storage Policy

1. Git durability:
- push checkpoint branch and tag to `origin`

2. Artifact durability:
- copy `output/freezes/<freeze_id>/` to at least one external durable store (NAS/cloud/object store)

## Restore Verification (Required Once)

1. Fresh clone and checkout tag:

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git
cd SpaceAGORA.jl
git checkout parallel-paper-freeze-YYYY-MM-DD
```

2. Restore archives into clone and verify checksums:

```bash
cd output/freezes/parallel_paper_YYYY-MM-DD
sha256sum -c SHA256SUMS.txt
```

3. Confirm expected outputs:
- key ladder report markdown
- key summary CSV tables
- figures used in paper/docs

