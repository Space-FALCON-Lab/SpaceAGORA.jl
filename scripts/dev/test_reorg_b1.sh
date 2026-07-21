#!/bin/bash
# Test-suite division, step B1 (PI-approved July 15): move gate scripts,
# runtime smokes, stress scripts, and coverage probes into subdirectories.
# The telemetry benchmark manifest does NOT move (it is a compiled-in default
# in src + the calibration subpackage + published example specs).
#
# Run from the repo root on a branch cut from main AFTER PRs #46 and #47 merge.
# All moves via git mv; all reference edits are exact-path and avoid the
# SpaceAGORACalibration.jl/test/ bare-name look-alikes.
set -euo pipefail

[ -f test/telemetry_benchmark_manifest.toml ] || { echo "run from repo root"; exit 1; }

# ---------- 1. moves ----------
mkdir -p test/gates test/probes

for f in test/ci_*_gate.jl; do git mv "$f" test/gates/; done

git mv test/ci_clean_depot_smoke.jl test/ci_threaded_smoke.jl \
       test/ci_no_gram_smoke.jl test/ci_examples_suite_smoke.jl \
       test/ci_examples_regression.jl test/smoke/

git mv test/ci_flake_guard.jl test/ci_nightly_montecarlo_stress.jl test/stress/

for f in test/coverage_debt_*_probes.jl; do
  git mv "$f" "test/probes/$(basename "$f" | sed 's/^coverage_debt_//')"
done
for f in test/coverage_*_probes.jl; do
  [ -e "$f" ] && git mv "$f" test/probes/
done

# ---------- 2. repo-root anchors inside moved files (depth 1 -> 2) ----------
for f in test/gates/*.jl test/smoke/ci_*.jl test/stress/ci_*.jl test/probes/*.jl; do
  # single pass: bump the repo-root anchor one level, whether bare or with a tail
  perl -pi -e 's/joinpath\(\@__DIR__, "\.\."([,)])/joinpath(\@__DIR__, "..", ".."$1/g' "$f"
done

# ---------- 3. runner shims that include the moved files ----------
# contracts runners: gates moved one level deeper under test/gates/
for r in test/contracts/pr_runtests.jl test/contracts/nightly_runtests.jl test/contracts/runtests.jl; do
  perl -pi -e 's/joinpath\(\@__DIR__, "\.\.", "(ci_[a-z0-9_]+_gate\.jl)"\)/joinpath(\@__DIR__, "..", "gates", "$1")/g' "$r"
done
# coverage runner: 3 gates + probes
perl -pi -e 's/joinpath\(\@__DIR__, "\.\.", "(ci_[a-z0-9_]+_gate\.jl)"\)/joinpath(\@__DIR__, "..", "gates", "$1")/g' test/coverage/runtests.jl
perl -pi -e 's/joinpath\(\@__DIR__, "\.\.", "(coverage_[a-z0-9_]+_probes\.jl)"\)/joinpath(\@__DIR__, "..", "probes", "$1")/g' test/coverage/runtests.jl
# smoke + stress shims: files moved INTO these dirs, ".." drops
perl -pi -e 's/joinpath\(\@__DIR__, "\.\.", "(ci_[a-z0-9_]+\.jl)"\)/joinpath(\@__DIR__, "$1")/g' test/smoke/runtests.jl test/stress/runtests.jl
# examples runner: regression smoke gains smoke/ segment
perl -pi -e 's/joinpath\(\@__DIR__, "\.\.", "\.\.", "ci_examples_regression\.jl"\)/joinpath(\@__DIR__, "..", "..", "smoke", "ci_examples_regression.jl")/' test/integration/examples/runtests.jl

# ---------- 4. suite drivers ----------
# 01: reads ci_clean_depot_smoke.jl as text by path
perl -pi -e 's{"ci_clean_depot_smoke\.jl"}{"smoke", "ci_clean_depot_smoke.jl"}g; s{test/ci_clean_depot_smoke\.jl}{test/smoke/ci_clean_depot_smoke.jl}g' test/suites/01_contract_and_api_tests.jl
# 02: spawns coverage_threaded_probes.jl by path
perl -pi -e 's{"coverage_threaded_probes\.jl"}{"probes", "coverage_threaded_probes.jl"}g; s{test/coverage_threaded_probes\.jl}{test/probes/coverage_threaded_probes.jl}g' test/suites/02_callbacks_parallel_and_smoke_tests.jl
# 05: includes 3 probes
perl -pi -e 's/joinpath\(\@__DIR__, "\.\.", "(coverage_[a-z0-9_]+_probes\.jl)"\)/joinpath(\@__DIR__, "..", "probes", "$1")/g' test/suites/05_thruster_control_and_quality_tests.jl
# 09: bare-filename list of the nine coverage_debt probes + subprocess path
perl -pi -e 's/"coverage_debt_([a-z0-9_]+_probes\.jl)"/"$1"/g' test/suites/09_coverage_debt_probe_drivers.jl
perl -pi -e 's{joinpath\((REPO_ROOT|[A-Za-z_]+), "test", ([A-Za-z_]+)\)}{joinpath($1, "test", "probes", $2)}g; s{joinpath\(\@__DIR__, "\.\.", ([a-z_]+)\)}{joinpath(\@__DIR__, "..", "probes", $1)}g' test/suites/09_coverage_debt_probe_drivers.jl

# ---------- 5. workflows (exact paths; " test/ci_" prefix cannot match the calibration copies) ----------
for w in .github/workflows/julia-ci.yml .github/workflows/nightly-stress.yml; do
  perl -pi -e 's{ test/(ci_[a-z0-9_]+_gate\.jl)}{ test/gates/$1}g' "$w"
  perl -pi -e 's{ test/(ci_(clean_depot_smoke|threaded_smoke|no_gram_smoke|examples_suite_smoke|examples_regression)\.jl)}{ test/smoke/$1}g' "$w"
done
# PR #46 paths filter: the two stress file entries are covered by test/stress/**
perl -0pi -e 's/      - "test\/ci_flake_guard\.jl"\n//; s/      - "test\/ci_nightly_montecarlo_stress\.jl"\n//' .github/workflows/nightly-stress.yml

# ---------- 6. three-way lockstep markers (gate <-> PR template <-> docs policy) ----------
perl -pi -e 's{test/ci_public_api_surface_gate\.jl}{test/gates/ci_public_api_surface_gate.jl}g; s{test/ci_hpc_extensibility_docs_gate\.jl}{test/gates/ci_hpc_extensibility_docs_gate.jl}g' \
  test/gates/ci_docs_contract_gate.jl .github/pull_request_template.md docs/src/documentation_policy.md

# ---------- 7. remaining docs + self-referential gate internals ----------
perl -pi -e 's{test/(ci_[a-z0-9_]+_gate\.jl)}{test/gates/$1}g; s{test/(ci_(clean_depot_smoke|threaded_smoke|examples_regression)\.jl)}{test/smoke/$1}g' \
  docs/quality/verification_contract.md docs/architecture/canonical_topology_contract.md
perl -pi -e 's{"test", "(ci_[a-z0-9_]+_gate\.jl)"}{"test", "gates", "$1"}g; s{test/ci_([a-z0-9_]+)_gate\.jl}{test/gates/ci_$1_gate.jl}g' \
  test/gates/ci_canonical_path_contract_gate.jl test/gates/ci_p1_findings_gate.jl

echo "B1 moves + rewrites applied. Now run the verification block."

# ---------- 8. late-found references (from dry-run leftover scan) ----------
# canonical-path gate's own exclusion rule for gate files follows them
perl -pi -e 's{startswith\(rel, joinpath\("test", "ci_"\)\)}{startswith(rel, joinpath("test", "gates", "ci_"))}' test/gates/ci_canonical_path_contract_gate.jl
# remaining docs pages naming gate paths
perl -pi -e 's{test/(ci_[a-z0-9_]+_gate\.jl)}{test/gates/$1}g' \
  docs/quality/api_naming_contract.md docs/architecture/src_completeness_contract.md docs/architecture/shim_window_manifest.md
# suite driver renamed truthfully (drives the probes; the debt is paid)
git mv test/suites/09_coverage_debt_probe_drivers.jl test/suites/09_probe_drivers.jl
perl -pi -e 's{"09_coverage_debt_probe_drivers\.jl"}{"09_probe_drivers.jl"}' test/integration/runtests.jl
perl -pi -e 's{Coverage Debt Probe Drivers}{Probe Drivers}; s{for the July 2026 coverage-debt register files;}{for the probe files under test/probes/ (the paid-off July 2026 coverage-debt register);}' test/suites/09_probe_drivers.jl
# self-exclusion idiom shared by several gates: they skip fellow gate files by
# relpath prefix; the prefix follows the gates into test/gates/
for f in $(grep -rln 'joinpath("test", "ci_")' test/ 2>/dev/null); do
  perl -pi -e 's{joinpath\("test", "ci_"\)}{joinpath("test", "gates", "ci_")}g' "$f"
done
# third repo-root anchor style: dirname chains gain one level
for f in test/gates/*.jl test/smoke/ci_*.jl test/stress/ci_*.jl test/probes/*.jl; do
  perl -pi -e 's{dirname\(dirname\(\@__FILE__\)\)}{dirname(dirname(dirname(\@__FILE__)))}g' "$f"
done
# suite 05 uses the REPO_ROOT form for its three probe includes
perl -pi -e 's{joinpath\(REPO_ROOT, "test", "(coverage_[a-z0-9_]+_probes\.jl)"\)}{joinpath(REPO_ROOT, "test", "probes", "$1")}g' test/suites/05_thruster_control_and_quality_tests.jl
