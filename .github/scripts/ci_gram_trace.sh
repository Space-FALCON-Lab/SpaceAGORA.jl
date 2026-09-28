#!/usr/bin/env bash
# TEMPORARY diagnostic: which files under data/GRAMSuite.jl a CI job touches.
# Prints paths only, never file contents.
#
#   start   install bpftrace and trace path syscalls system-wide
#   mark    after the checkout: set every tracked file's atime to the epoch, so
#           a later read (relatime) moves it forward
#   report  list the tracked files read (atime moved), and every path probed
#           by open/stat/access (bpftrace), marking the ones that do not exist
set -uo pipefail
cd "${GITHUB_WORKSPACE:-.}"
G=data/GRAMSuite.jl
T="${RUNNER_TEMP:-/tmp}/gram-trace"
mkdir -p "$T"

case "${1:-}" in
start)
  sudo apt-get install -y -qq bpftrace > "$T/apt.log" 2>&1 || { echo "gram_trace: bpftrace unavailable"; exit 0; }
  cat > "$T/probe.bt" <<'EOF'
tracepoint:syscalls:sys_enter_openat,
tracepoint:syscalls:sys_enter_newfstatat,
tracepoint:syscalls:sys_enter_statx,
tracepoint:syscalls:sys_enter_faccessat,
tracepoint:syscalls:sys_enter_faccessat2,
tracepoint:syscalls:sys_enter_readlinkat
{ printf("%s\t%s\n", comm, str(args.filename)); }
tracepoint:syscalls:sys_enter_access,
tracepoint:syscalls:sys_enter_newstat,
tracepoint:syscalls:sys_enter_newlstat,
tracepoint:syscalls:sys_enter_execve
{ printf("%s\t%s\n", comm, str(args.filename)); }
EOF
  sudo BPFTRACE_MAX_STRLEN=200 setsid bpftrace "$T/probe.bt" > "$T/bt.log" 2> "$T/bt.err" < /dev/null &
  echo $! > "$T/bt.pid"
  sleep 3
  echo "gram_trace: bpftrace started ($(wc -l < "$T/bt.log") lines so far); $(head -c 300 "$T/bt.err")"
  ;;
mark)
  git -C "$G" ls-files -z | (cd "$G" && xargs -0 touch -c -a -h -d @0) 2>/dev/null
  echo "gram_trace: atimes reset on $(git -C "$G" ls-files | wc -l) tracked files; mount: $(findmnt -no OPTIONS /)"
  ;;
report)
  sudo pkill -INT bpftrace 2>/dev/null; sleep 2
  job="${GRAM_TRACE_JOB:-${GITHUB_JOB:-job}}"
  if [ ! -d "$G/.git" ] && [ ! -f "$G/.git" ]; then echo "gram_trace[$job]: no GRAMSuite checkout"; exit 0; fi
  (cd "$G" && git ls-files -z | xargs -0 stat -c '%X %s %n' 2>/dev/null) | awk '$1 > 0' | cut -d' ' -f2- > "$T/read.txt"
  echo "gram_trace[$job]: tracked files read: $(wc -l < "$T/read.txt"), bytes $(awk '{s+=$1} END {print s+0}' "$T/read.txt")"
  echo "gram_trace[$job]: tracked files total: $(git -C "$G" ls-files | wc -l), bytes $(cd "$G" && git ls-files -z | xargs -0 stat -c %s 2>/dev/null | awk '{s+=$1} END {print s+0}')"
  echo "gram_trace[$job]: fetched git dir bytes $(du -sb .git/modules/GRAMSuite.jl 2>/dev/null | cut -f1)"
  sed "s/^/gram_read[$job] /" "$T/read.txt"
  if [ -s "$T/bt.log" ]; then
    root="$(pwd)/$G/"
    grep -F "GRAMSuite.jl" "$T/bt.log" | grep -vE '^(git|git-remote-http|git-lfs|tar|zstd|rm|find|xargs|touch|stat)\b' \
      | cut -f2 | sort -u > "$T/probed.txt"
    while IFS= read -r p; do
      rel="${p#"$root"}"
      if [ -e "$p" ]; then state=exists; else state=MISSING; fi
      printf 'gram_probe[%s] %s %s\n' "$job" "$state" "$rel"
    done < "$T/probed.txt"
    echo "gram_trace[$job]: distinct paths probed: $(wc -l < "$T/probed.txt"); bpftrace lines $(wc -l < "$T/bt.log")"
  else
    echo "gram_trace[$job]: no bpftrace log: $(head -c 500 "$T/bt.err" 2>/dev/null)"
  fi
  ;;
esac
