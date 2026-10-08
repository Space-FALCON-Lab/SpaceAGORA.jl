"""Run the README's bounded CLI quickstart and check its fresh output artifacts."""

import argparse
import csv
import json
import math
import os
import signal
from pathlib import Path
import struct
import subprocess
import tempfile
import zlib


REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = "AGORA_Basic_Quickstart.jl"
SMOKE_SECONDS = 120.0
REQUIRED_COLUMNS = (
    "time", "sc1_altitude", "sc1_periapsis_altitude",
    "sc1_pos_1", "sc1_pos_2", "sc1_pos_3",
    "sc1_vel_1", "sc1_vel_2", "sc1_vel_3",
    "sc1_longitude_deg", "sc1_latitude_deg",
)
PLOTS = (
    "quickstart_altitude_speed.png",
    "quickstart_inertial_trajectory.png",
    "quickstart_3d_orbit.png",
    "quickstart_ground_track.png",
)


def png_dimensions(path):
    """Check PNG framing/CRCs and positive dimensions, not visual correctness."""
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError(f"Invalid PNG signature: {path.name}")
    offset = 8
    dimensions = None
    has_pixels = False
    while offset + 12 <= len(data):
        size = struct.unpack_from(">I", data, offset)[0]
        kind = data[offset + 4:offset + 8]
        end = offset + 12 + size
        if end > len(data):
            break
        payload = data[offset + 8:end - 4]
        crc = struct.unpack_from(">I", data, end - 4)[0]
        if zlib.crc32(kind + payload) != crc:
            raise ValueError(f"Invalid PNG checksum: {path.name}")
        if offset == 8:
            if kind != b"IHDR" or size != 13:
                raise ValueError(f"Missing PNG image header: {path.name}")
            dimensions = struct.unpack_from(">II", payload)
            if 0 in dimensions:
                raise ValueError(f"Empty PNG dimensions: {path.name}")
        if kind == b"IDAT" and size:
            has_pixels = True
        if kind == b"IEND":
            if size or end != len(data) or not has_pixels:
                raise ValueError(f"Incomplete PNG image: {path.name}")
            return dimensions
        offset = end
    raise ValueError(f"Truncated PNG image: {path.name}")


def verify_outputs(output_dir):
    """Validate chart inputs and all four promised files for this smoke run."""
    output_dir = Path(output_dir)
    csv_path = output_dir / "simulation_results.csv"
    times = []
    with csv_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        missing = sorted(set(REQUIRED_COLUMNS) - set(fields))
        if missing or len(fields) != len(set(fields)):
            raise ValueError(f"Invalid quickstart CSV columns; missing={missing}")
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Extra CSV fields on row {row_number}")
            for name in REQUIRED_COLUMNS:
                try:
                    value = float(row[name])
                except (TypeError, ValueError) as error:
                    raise ValueError(f"Invalid {name} on CSV row {row_number}") from error
                if not math.isfinite(value):
                    raise ValueError(f"Nonfinite {name} on CSV row {row_number}")
            times.append(float(row["time"]))
    if len(times) < 2:
        raise ValueError("Quickstart CSV must contain at least two samples")
    if not all(later > earlier for earlier, later in zip(times, times[1:])):
        raise ValueError("Quickstart CSV times must be strictly increasing")
    if not (math.isclose(times[0], 0.0, abs_tol=1e-6, rel_tol=0.0)
            and math.isclose(times[-1], SMOKE_SECONDS, abs_tol=1e-6, rel_tol=0.0)):
        raise ValueError("Quickstart CSV does not cover the 120-second smoke interval")
    plots = {name: png_dimensions(output_dir / "plots" / name) for name in PLOTS}
    return {"rows": len(times), "first_time": times[0], "last_time": times[-1], "plots": plots}


def run_smoke(julia="julia", *, repo_root=REPO_ROOT, timeout=1200):
    repo_root = Path(repo_root).resolve()
    # Absolute paths preserve the README command while isolating cwd/output.
    command = [julia, "--startup-file=no", "--compiled-modules=existing",
               f"--project={repo_root}", str(repo_root / "src/cli/main.jl"),
               "run", f"--example={EXAMPLE}", "--smoke"]
    env = os.environ.copy()
    # Let the documented --smoke switch enable smoke mode and CSV saving itself.
    for name in ("SPACEAGORA_CLI_OUTPUT_DIR", "SPACEAGORA_EXAMPLE_SMOKE",
                 "SPACEAGORA_EXAMPLE_SMOKE_RESULTS", "SPACEAGORA_VISUALIZATION"):
        env.pop(name, None)
    env.update(SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME=str(SMOKE_SECONDS),
               SPACEAGORA_EXAMPLE_PARALLEL="0", GKSwstype="100")
    # A unique, initially empty cwd prevents old repository outputs masking failure.
    with tempfile.TemporaryDirectory(prefix="spaceagora-readme-") as workdir:
        # The Julia CLI starts another Julia process. On timeout, stop the whole
        # session before removing its cwd. This runner targets POSIX CI hosts.
        with subprocess.Popen(command, cwd=workdir, env=env, start_new_session=True) as process:
            try:
                returncode = process.wait(timeout=timeout)
            except BaseException:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
                raise
            if returncode:
                raise subprocess.CalledProcessError(returncode, command)
        return verify_outputs(Path(workdir) / "output")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default="julia", help="Julia executable")
    args = parser.parse_args()
    print(json.dumps(run_smoke(args.julia), sort_keys=True))


if __name__ == "__main__":
    main()
