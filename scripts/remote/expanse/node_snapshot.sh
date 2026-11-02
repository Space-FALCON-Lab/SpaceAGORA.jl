#!/usr/bin/env bash
# Node state before one benchmark point, appended to the file given as $1. Read-only, no side effects on the child.
f="${1:?output file}"
{
  echo "== $(date -Iseconds) host=$(hostname)"
  awk '/^cpu MHz/ { s += $4; if (!n || $4 < lo) lo = $4; if ($4 > hi) hi = $4; n++ }
       END { if (n) printf "cpu MHz mean=%.1f min=%.1f max=%.1f over %d cpus\n", s / n, lo, hi, n }' /proc/cpuinfo
  echo "loadavg: $(cat /proc/loadavg)"
  echo "busiest (pcpu is lifetime average, not instantaneous):"
  ps -eo pid,psr,pcpu,args --sort=-pcpu | head -6 | cut -c1-160
} >> "$f" 2>&1 || true
