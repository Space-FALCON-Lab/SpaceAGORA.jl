#!/bin/sh
# Run the Basilisk driver against a MuJoCo-enabled Basilisk 2.12 build (see README.md).
#   BSK_DIST=<build>/dist3 BSK_PY=<python> sh run_bsk.sh <case.toml> <outdir> [duration_s]
exec env PYTHONPATH="${BSK_DIST:?set BSK_DIST}" "${BSK_PY:-python3}" "$(dirname "$0")/bsk_driver.py" "$@"
