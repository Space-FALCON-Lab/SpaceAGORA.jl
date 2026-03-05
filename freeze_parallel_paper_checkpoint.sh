#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

timestamp_utc() {
  date -u '+%Y-%m-%dT%H:%M:%SZ'
}

log() {
  printf '[freeze] %s\n' "$*"
}

die() {
  printf '[freeze][error] %s\n' "$*" >&2
  exit 1
}

usage() {
  cat <<'EOF'
Usage: ./freeze_parallel_paper_checkpoint.sh [options]

Options:
  --freeze-date=YYYY-MM-DD       Freeze date for naming (default: today).
  --freeze-id=ID                 Override freeze artifact folder name.
  --checkpoint-branch=NAME       Default: checkpoint/parallel-paper-<date>
  --tag=NAME                     Default: parallel-paper-freeze-<date>
  --purpose=TEXT                 Manifest/tag purpose text.
  --profile=NAME                 Runtime profile label (default: full).
  --run-script=PATH              Run script path to record in metadata.
  --run-flags=TEXT               Run flags text to record in metadata.
  --run-outdir=PATH              Run outputs root to archive (default: latest output/performance/paper_protocol_*)
  --run-log=PATH                 Run log path to archive (default: latest output/performance/overnight_*.log)
  --figure-dirs=CSV              Comma-separated dirs to scan for figures (default: docs,<run-outdir>)
  --archive-root=PATH            Artifact root folder (default: output/freezes)
  --allow-dirty=0|1              Allow dirty working tree (default: 1)
  --push-remote=0|1              Push checkpoint branch + tag (default: 1)
  --dry-run=0|1                  Show plan only (default: 0)
  --help                         Show this message.

Examples:
  ./freeze_parallel_paper_checkpoint.sh
  ./freeze_parallel_paper_checkpoint.sh --push-remote=0 --dry-run=1
  ./freeze_parallel_paper_checkpoint.sh --run-outdir=output/performance/paper_protocol_20260304_224630
EOF
}

need_cmd() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || die "Missing required command: $cmd"
}

first_existing_path() {
  local p
  for p in "$@"; do
    if [[ -n "$p" && -e "$p" ]]; then
      printf '%s\n' "$p"
      return 0
    fi
  done
  return 1
}

collect_existing_paths() {
  local out_file="$1"
  shift
  : > "$out_file"
  local p
  for p in "$@"; do
    if [[ -n "$p" && -e "$p" ]]; then
      printf '%s\n' "$p" >> "$out_file"
    fi
  done
}

collect_figure_paths() {
  local out_file="$1"
  shift
  : > "$out_file"
  local d
  for d in "$@"; do
    if [[ -d "$d" ]]; then
      find "$d" -type f \( -name '*.png' -o -name '*.pdf' -o -name '*.svg' \) >> "$out_file"
    fi
  done
  sort -u -o "$out_file" "$out_file"
}

collect_metadata_paths() {
  local out_file="$1"
  local run_outdir="$2"
  : > "$out_file"
  if [[ -d "$run_outdir" ]]; then
    find "$run_outdir" -type f | grep -E '(runtime_hardware_info|hardware_info|rung_metadata|rung_env|metadata)' > "$out_file" || true
  fi
  sort -u -o "$out_file" "$out_file"
}

append_archive_row() {
  local file="$1"
  local table_file="$2"
  local checksum
  local size
  checksum="$(sha256sum "$file" | awk '{print $1}')"
  size="$(stat -c '%s' "$file")"
  printf '| `%s` | %s | `%s` |\n' "$(basename "$file")" "$size" "$checksum" >> "$table_file"
}

need_cmd git
need_cmd tar
need_cmd sha256sum
need_cmd awk
need_cmd sed
need_cmd date
need_cmd find

[[ -f Project.toml ]] || die "Run from SpaceAGORA repo root (Project.toml missing)."
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || die "Not inside a git repo."

INVOCATION="$(printf '%q ' "$0" "$@")"

FREEZE_DATE="$(date +%F)"
FREEZE_ID=""
CHECKPOINT_BRANCH=""
TAG_NAME=""
PURPOSE="Parallel paper state freeze before accuracy-related code changes."
PROFILE="full"
RUN_SCRIPT="./run_parallel_paper_protocol.sh"
RUN_FLAGS="(defaults in tagged scripts)"
RUN_OUTDIR="$(ls -1dt output/performance/paper_protocol_* 2>/dev/null | head -n1 || true)"
RUN_LOG="$(ls -1t output/performance/overnight_*.log 2>/dev/null | head -n1 || true)"
FIGURE_DIRS_CSV=""
ARCHIVE_ROOT="output/freezes"
ALLOW_DIRTY="1"
PUSH_REMOTE="1"
DRY_RUN="0"

for arg in "$@"; do
  case "$arg" in
    --freeze-date=*)
      FREEZE_DATE="${arg#*=}"
      ;;
    --freeze-id=*)
      FREEZE_ID="${arg#*=}"
      ;;
    --checkpoint-branch=*)
      CHECKPOINT_BRANCH="${arg#*=}"
      ;;
    --tag=*)
      TAG_NAME="${arg#*=}"
      ;;
    --purpose=*)
      PURPOSE="${arg#*=}"
      ;;
    --profile=*)
      PROFILE="${arg#*=}"
      ;;
    --run-script=*)
      RUN_SCRIPT="${arg#*=}"
      ;;
    --run-flags=*)
      RUN_FLAGS="${arg#*=}"
      ;;
    --run-outdir=*)
      RUN_OUTDIR="${arg#*=}"
      ;;
    --run-log=*)
      RUN_LOG="${arg#*=}"
      ;;
    --figure-dirs=*)
      FIGURE_DIRS_CSV="${arg#*=}"
      ;;
    --archive-root=*)
      ARCHIVE_ROOT="${arg#*=}"
      ;;
    --allow-dirty=*)
      ALLOW_DIRTY="${arg#*=}"
      ;;
    --push-remote=*)
      PUSH_REMOTE="${arg#*=}"
      ;;
    --dry-run=*)
      DRY_RUN="${arg#*=}"
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $arg (see --help)"
      ;;
  esac
done

[[ -n "$FREEZE_ID" ]] || FREEZE_ID="parallel_paper_${FREEZE_DATE}"
[[ -n "$CHECKPOINT_BRANCH" ]] || CHECKPOINT_BRANCH="checkpoint/parallel-paper-${FREEZE_DATE}"
[[ -n "$TAG_NAME" ]] || TAG_NAME="parallel-paper-freeze-${FREEZE_DATE}"

if [[ -z "$RUN_OUTDIR" || ! -d "$RUN_OUTDIR" ]]; then
  die "Run outdir not found. Pass --run-outdir=<path>."
fi
if [[ -z "$RUN_LOG" || ! -f "$RUN_LOG" ]]; then
  log "Run log missing; continuing without run log archive."
  RUN_LOG=""
fi

if [[ "$ALLOW_DIRTY" != "1" && "$ALLOW_DIRTY" != "0" ]]; then
  die "--allow-dirty must be 0 or 1"
fi
if [[ "$PUSH_REMOTE" != "1" && "$PUSH_REMOTE" != "0" ]]; then
  die "--push-remote must be 0 or 1"
fi
if [[ "$DRY_RUN" != "1" && "$DRY_RUN" != "0" ]]; then
  die "--dry-run must be 0 or 1"
fi

WORKTREE_STATUS="$(git status --short)"
if [[ -n "$WORKTREE_STATUS" && "$ALLOW_DIRTY" == "0" ]]; then
  die "Working tree is dirty and --allow-dirty=0. Commit/stash first."
fi

HEAD_SHA="$(git rev-parse HEAD)"
CURRENT_BRANCH="$(git branch --show-current)"
REMOTE_URL="$(git remote get-url origin 2>/dev/null || echo "<no-origin>")"
JULIA_VERSION="$(julia --version 2>/dev/null | head -n1 || echo "julia_unavailable")"
SUBMODULE_STATUS="$(git submodule status --recursive || true)"
HOST_NAME="$(hostname)"
UNAME_LINE="$(uname -a)"
CPU_MODEL="$(lscpu 2>/dev/null | awk -F: '/Model name/{gsub(/^[ \t]+/, "", $2); print $2; exit}' || true)"
CPU_THREADS="$(nproc 2>/dev/null || echo "unknown")"
MEM_TOTAL="$(awk '/MemTotal/{printf "%.1f GiB", $2/1048576}' /proc/meminfo 2>/dev/null || echo "unknown")"

FREEZE_DIR="${ARCHIVE_ROOT}/${FREEZE_ID}"
INPUTS_ARCHIVE="${FREEZE_DIR}/${FREEZE_ID}_repro_inputs.tar.gz"
OUTPUTS_ARCHIVE="${FREEZE_DIR}/${FREEZE_ID}_outputs_$(basename "$RUN_OUTDIR").tar.gz"
METADATA_ARCHIVE="${FREEZE_DIR}/${FREEZE_ID}_hardware_metadata.tar.gz"
FIGURES_ARCHIVE="${FREEZE_DIR}/${FREEZE_ID}_figures.tar.gz"
LOGS_ARCHIVE="${FREEZE_DIR}/${FREEZE_ID}_logs.tar.gz"
SHA_FILE="${FREEZE_DIR}/SHA256SUMS.txt"
MANIFEST_FILE="${FREEZE_DIR}/FREEZE_MANIFEST.md"
ARCHIVE_TABLE="${FREEZE_DIR}/.archive_table.md"
TAG_MESSAGE_FILE="${FREEZE_DIR}/.tag_message.txt"
REPRO_INPUT_LIST="${FREEZE_DIR}/repro_inputs_filelist.txt"
METADATA_LIST="${FREEZE_DIR}/hardware_metadata_filelist.txt"
FIGURE_LIST="${FREEZE_DIR}/figure_filelist.txt"
LOG_LIST="${FREEZE_DIR}/log_filelist.txt"

if [[ "$DRY_RUN" == "1" ]]; then
  log "DRY RUN summary"
  log "HEAD: $HEAD_SHA (branch=$CURRENT_BRANCH)"
  log "Checkpoint branch: $CHECKPOINT_BRANCH"
  log "Tag: $TAG_NAME"
  log "Run outdir: $RUN_OUTDIR"
  log "Run log: ${RUN_LOG:-<none>}"
  log "Freeze dir: $FREEZE_DIR"
  log "Push remote: $PUSH_REMOTE"
  if [[ -n "$WORKTREE_STATUS" ]]; then
    log "Dirty worktree detected (allowed)."
  else
    log "Worktree clean."
  fi
  exit 0
fi

mkdir -p "$FREEZE_DIR"

log "Creating checkpoint branch ref: $CHECKPOINT_BRANCH -> $HEAD_SHA"
if git show-ref --verify --quiet "refs/heads/$CHECKPOINT_BRANCH"; then
  EXISTING_BRANCH_SHA="$(git rev-parse "$CHECKPOINT_BRANCH")"
  if [[ "$EXISTING_BRANCH_SHA" != "$HEAD_SHA" ]]; then
    die "Branch $CHECKPOINT_BRANCH already exists at $EXISTING_BRANCH_SHA (expected $HEAD_SHA)."
  fi
else
  git branch "$CHECKPOINT_BRANCH" "$HEAD_SHA"
fi

cat > "$TAG_MESSAGE_FILE" <<EOF
parallel paper freeze checkpoint

purpose: $PURPOSE
commit_sha: $HEAD_SHA
checkpoint_branch: $CHECKPOINT_BRANCH
freeze_date: $FREEZE_DATE
profile: $PROFILE
julia_version: $JULIA_VERSION
run_script: $RUN_SCRIPT
run_flags: $RUN_FLAGS
run_outdir: $RUN_OUTDIR
run_log: ${RUN_LOG:-<none>}
EOF

log "Creating/updating annotated tag metadata file"
if git rev-parse -q --verify "refs/tags/$TAG_NAME" >/dev/null 2>&1; then
  EXISTING_TAG_SHA="$(git rev-list -n1 "$TAG_NAME")"
  if [[ "$EXISTING_TAG_SHA" != "$HEAD_SHA" ]]; then
    die "Tag $TAG_NAME already exists at $EXISTING_TAG_SHA (expected $HEAD_SHA)."
  fi
  log "Tag $TAG_NAME already exists at desired commit; reusing."
else
  git tag -a "$TAG_NAME" "$HEAD_SHA" -F "$TAG_MESSAGE_FILE"
fi

log "Collecting reproducibility inputs"
collect_existing_paths "$REPRO_INPUT_LIST" \
  Project.toml \
  Manifest.toml \
  .AGORA/Project.toml \
  .AGORA/Manifest.toml \
  .gitmodules \
  run_parallel_paper_protocol.sh \
  test/performance_smart_parallel_ladder.jl \
  test/performance_runtime_analysis.jl \
  test/performance_paper_pipeline.jl \
  cutover_after_rung_completion.sh \
  progress_email_daemon.sh \
  configure_icloud_progress_email.sh

if [[ -s "$REPRO_INPUT_LIST" ]]; then
  tar -czf "$INPUTS_ARCHIVE" -T "$REPRO_INPUT_LIST"
else
  die "No reproducibility inputs found to archive."
fi

log "Archiving run outputs: $RUN_OUTDIR"
tar -czf "$OUTPUTS_ARCHIVE" -C "$(dirname "$RUN_OUTDIR")" "$(basename "$RUN_OUTDIR")"

log "Collecting hardware metadata files"
collect_metadata_paths "$METADATA_LIST" "$RUN_OUTDIR"
if [[ -s "$METADATA_LIST" ]]; then
  tar -czf "$METADATA_ARCHIVE" -T "$METADATA_LIST"
fi

if [[ -z "$FIGURE_DIRS_CSV" ]]; then
  FIGURE_DIRS_CSV="docs,$RUN_OUTDIR"
fi
IFS=',' read -r -a FIGURE_DIRS <<< "$FIGURE_DIRS_CSV"
log "Collecting figure files from: ${FIGURE_DIRS[*]}"
collect_figure_paths "$FIGURE_LIST" "${FIGURE_DIRS[@]}"
if [[ -s "$FIGURE_LIST" ]]; then
  tar -czf "$FIGURES_ARCHIVE" -T "$FIGURE_LIST"
fi

collect_existing_paths "$LOG_LIST" \
  "$RUN_LOG" \
  output/performance/progress_email_daemon.log \
  output/performance/progress_email_daemon.nohup.log \
  output/performance/hourly_monitor.log \
  output/performance/hourly_monitor_alerts.log
if [[ -s "$LOG_LIST" ]]; then
  tar -czf "$LOGS_ARCHIVE" -T "$LOG_LIST"
fi

: > "$SHA_FILE"
(
  cd "$FREEZE_DIR"
  for archive in "$INPUTS_ARCHIVE" "$OUTPUTS_ARCHIVE" "$METADATA_ARCHIVE" "$FIGURES_ARCHIVE" "$LOGS_ARCHIVE"; do
    if [[ -f "$archive" ]]; then
      sha256sum "$(basename "$archive")" >> "$(basename "$SHA_FILE")"
    fi
  done
)

: > "$ARCHIVE_TABLE"
printf '| Archive | Bytes | SHA256 |\n' >> "$ARCHIVE_TABLE"
printf '|---|---:|---|\n' >> "$ARCHIVE_TABLE"
for archive in "$INPUTS_ARCHIVE" "$OUTPUTS_ARCHIVE" "$METADATA_ARCHIVE" "$FIGURES_ARCHIVE" "$LOGS_ARCHIVE"; do
  if [[ -f "$archive" ]]; then
    append_archive_row "$archive" "$ARCHIVE_TABLE"
  fi
done

{
  echo "# FREEZE_MANIFEST"
  echo
  echo "Generated: $(timestamp_utc)"
  echo
  echo "## Freeze Identity"
  echo "- Freeze ID: \`$FREEZE_ID\`"
  echo "- Purpose: $PURPOSE"
  echo "- Checkpoint branch: \`$CHECKPOINT_BRANCH\`"
  echo "- Tag: \`$TAG_NAME\`"
  echo "- Commit: \`$HEAD_SHA\`"
  echo
  echo "## Repro Context"
  echo "- Profile: \`$PROFILE\`"
  echo "- Run script: \`$RUN_SCRIPT\`"
  echo "- Run flags: \`$RUN_FLAGS\`"
  echo "- Run outdir: \`$RUN_OUTDIR\`"
  echo "- Run log: \`${RUN_LOG:-<none>}\`"
  echo "- Freeze command: \`$INVOCATION\`"
  echo
  echo "## Repository State"
  echo "- Current branch at freeze time: \`$CURRENT_BRANCH\`"
  echo "- Remote origin: \`$REMOTE_URL\`"
  if [[ -n "$WORKTREE_STATUS" ]]; then
    echo "- Working tree: dirty (allowed by --allow-dirty=1)"
    echo
    echo '```text'
    echo "$WORKTREE_STATUS"
    echo '```'
  else
    echo "- Working tree: clean"
  fi
  echo
  echo "## Submodule SHAs"
  echo '```text'
  if [[ -n "$SUBMODULE_STATUS" ]]; then
    echo "$SUBMODULE_STATUS"
  else
    echo "<none>"
  fi
  echo '```'
  echo
  echo "## Machine Metadata"
  echo "- Hostname: \`$HOST_NAME\`"
  echo "- OS: \`$UNAME_LINE\`"
  echo "- CPU model: \`${CPU_MODEL:-unknown}\`"
  echo "- CPU threads: \`$CPU_THREADS\`"
  echo "- Memory: \`$MEM_TOTAL\`"
  echo "- Julia: \`$JULIA_VERSION\`"
  echo
  echo "## Archives"
  cat "$ARCHIVE_TABLE"
  echo
  echo "## Integrity"
  echo "- SHA256 file: \`$SHA_FILE\`"
} > "$MANIFEST_FILE"

if [[ "$PUSH_REMOTE" == "1" ]]; then
  log "Pushing checkpoint branch + tag to origin"
  git push origin "$CHECKPOINT_BRANCH"
  git push origin "$TAG_NAME"
else
  log "Skipping remote push (push-remote=0)"
fi

log "Freeze complete."
log "Manifest: $MANIFEST_FILE"
log "Checksums: $SHA_FILE"
