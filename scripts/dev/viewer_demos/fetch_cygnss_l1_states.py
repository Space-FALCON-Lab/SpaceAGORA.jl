#!/usr/bin/env python3
"""Fetch the flown CYGNSS spacecraft states for the 96-hour viewer window.

    python3 scripts/dev/viewer_demos/fetch_cygnss_l1_states.py

What it does, and why it does it this way:

  * The states live in the NASA CYGNSS Level 1 Science Data Record version 3.2
    (collection short name ``CYGNSS_L1_V3.2``, PO.DAAC/POCLOUD), one granule per
    spacecraft per day.  A granule is about 1.09 GB because it carries the full
    delay-Doppler map cube; the eight spacecraft-state variables in it are about
    0.5% of that.  Downloading whole granules for a four-day window and seven
    spacecraft would move roughly 30 GB to extract about 150 MB of useful bytes,
    so this script never downloads a granule.
  * It uses the OPeNDAP (Hyrax) service NASA publishes for the collection and a
    DAP4 constraint expression naming only the position, velocity and time
    variables.  That is variable subsetting at the server: the response carries
    only the requested arrays.
  * Every response is cached verbatim under ``l1_opendap_cache/`` in the output
    directory, so a second run moves no bytes at all.  The archive is a shared
    resource; one pass of 28 requests is the whole network cost.

Authentication is a NASA Earthdata Login bearer token read from ``~/.edl_token``
at the moment of use.  It is sent only to ``*.earthdata.nasa.gov`` and is never
written to a file, a log or this script's output.

The output is one Feather table of the Earth-fixed states, which
``build_cygnss_ics.jl`` reads and converts to J2000.  The frame conversion is
deliberately not done here: the repository already has a verified one in
``cygnss_ics.jl`` and a second implementation would be a second thing to be
wrong.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.feather as feather

REPO_ROOT = Path(__file__).resolve().parents[3]
TELEM_DIR = REPO_ROOT / "data" / "telemetry" / "CYGNSS"

CMR_GRANULES = "https://cmr.earthdata.nasa.gov/search/granules.json"
SHORT_NAME = "CYGNSS_L1_V3.2"
TOKEN_PATH = Path.home() / ".edl_token"

# The 96-hour viewer window: 2025-06-06T00:00:00Z to 2025-06-09T24:00:00Z.
WINDOW_START = "2025-06-06T00:00:00Z"
WINDOW_END = "2025-06-09T23:59:59Z"

# The eight spacecraft-state variables.  Position and velocity are the PVT
# (GPS navigation solution) pair, not the DDM-time interpolation of it, because
# the PVT pair is what the mirrored FM01 and FM04 tables use.
DAP4_VARIABLES = (
    "/pvt_timestamp_utc",
    "/sc_pos_x_pvt",
    "/sc_pos_y_pvt",
    "/sc_pos_z_pvt",
    "/sc_vel_x_pvt",
    "/sc_vel_y_pvt",
    "/sc_vel_z_pvt",
    "/spacecraft_num",
)

# Fill values declared by the granule's own metadata for these variables.
POS_FILL = -99999999
VEL_FILL = -9999

DAP4_ATOMIC = {
    "Int8": "i1", "UInt8": "u1", "Byte": "u1",
    "Int16": "i2", "UInt16": "u2",
    "Int32": "i4", "UInt32": "u4",
    "Int64": "i8", "UInt64": "u8",
    "Float32": "f4", "Float64": "f8",
}


# ---------------------------------------------------------------------------
# DAP4 response decoding
# ---------------------------------------------------------------------------

def dap4_dechunk(blob: bytes):
    """Split a DAP4 data response into its DMR text and its data bytes.

    The body is a sequence of chunks, each introduced by a four-byte header
    whose first byte is flags (0x01 last chunk, 0x02 error, 0x04 the payload is
    little-endian) and whose remaining three bytes are the chunk length, big
    endian.  The first chunk is the DMR; the rest are the serialized variables.
    """
    offset = 0
    dmr = None
    little = True
    parts = []
    while offset < len(blob):
        flags = blob[offset]
        length = int.from_bytes(blob[offset + 1:offset + 4], "big")
        offset += 4
        body = blob[offset:offset + length]
        offset += length
        if flags & 0x02:
            raise RuntimeError("DAP4 error chunk: " + body[:500].decode("latin-1"))
        if dmr is None:
            dmr = body.decode("latin-1")
            little = bool(flags & 0x04)
        else:
            parts.append(body)
        if flags & 0x01 and parts:
            break
    if dmr is None:
        raise RuntimeError("DAP4 response carried no DMR chunk")
    return dmr, b"".join(parts), little


def dap4_variable_layout(dmr: str):
    """Ordered (name, numpy code, shape) of every atomic variable in a DMR."""
    sizes = {m.group(1): int(m.group(2))
             for m in re.finditer(r'<Dimension name="([^"]+)" size="(\d+)"/>', dmr)}
    layout = []
    pattern = r'<(%s) name="([^"]+)">(.*?)</\1>' % "|".join(DAP4_ATOMIC)
    for m in re.finditer(pattern, dmr, re.S):
        kind, name, body = m.group(1), m.group(2), m.group(3)
        shape = []
        for d in re.finditer(r'<Dim (?:name="([^"]+)"|size="(\d+)")/>', body):
            shape.append(sizes[d.group(1).lstrip("/")] if d.group(1) else int(d.group(2)))
        layout.append((name, DAP4_ATOMIC[kind], tuple(shape)))
    return layout


def dap4_read(blob: bytes):
    """Decode a DAP4 data response into a name -> array mapping.

    Each variable is serialized in DMR order and followed by a four-byte CRC32
    checksum, which is skipped here: the transport is HTTPS over TCP and the
    arrays are validated downstream against their declared fill values and
    against the mirrored FM01 and FM04 tables.
    """
    dmr, data, little = dap4_dechunk(blob)
    order = "<" if little else ">"
    out = {}
    offset = 0
    for name, code, shape in dap4_variable_layout(dmr):
        count = int(np.prod(shape)) if shape else 1
        dtype = np.dtype(order + code)
        nbytes = count * dtype.itemsize
        values = np.frombuffer(data[offset:offset + nbytes], dtype=dtype)
        if values.size != count:
            raise RuntimeError(f"DAP4 response truncated in variable {name}")
        offset += nbytes + 4
        out[name] = values.reshape(shape).astype(dtype.newbyteorder("=")) if shape else values[0]
    return out, dmr


def dmr_attribute(dmr: str, name: str) -> str | None:
    m = re.search(r'<Attribute name="%s" type="String">\s*<Value>(.*?)</Value>' % re.escape(name),
                  dmr, re.S)
    return m.group(1).strip() if m else None


# ---------------------------------------------------------------------------
# Network
# ---------------------------------------------------------------------------

def _bearer_header() -> dict:
    """The Earthdata Login authorization header, built at the moment of use."""
    token = TOKEN_PATH.read_text().strip()
    if not token:
        raise SystemExit(f"{TOKEN_PATH} is empty; a NASA Earthdata Login token is required")
    return {"Authorization": "Bearer " + token}


def _get(url: str, headers: dict | None = None, timeout: int = 300) -> bytes:
    request = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def cmr_granules(start: str, end: str):
    """Every CYGNSS_L1_V3.2 granule overlapping the window, newest CMR page size."""
    query = urllib.parse.urlencode({
        "short_name": SHORT_NAME,
        "temporal": f"{start},{end}",
        "page_size": 200,
        "sort_key": "start_date",
    })
    feed = json.loads(_get(f"{CMR_GRANULES}?{query}"))["feed"]["entry"]
    granules = []
    for entry in feed:
        opendap = next((l["href"] for l in entry["links"]
                        if l.get("title", "").startswith("OPeNDAP")), None)
        if opendap is None:
            raise RuntimeError(f"granule {entry['title']} has no OPeNDAP link")
        granules.append({
            "title": entry["title"],
            "concept_id": entry["id"],
            "size_mb": float(entry["granule_size"]),
            "time_start": entry["time_start"],
            "opendap": opendap,
        })
    return granules


def fetch_granule_variables(granule: dict, cache_dir: Path, pause_s: float):
    """The eight state variables of one granule, from cache or from OPeNDAP."""
    cached = cache_dir / (granule["title"] + ".dap")
    if cached.exists():
        return cached.read_bytes(), 0
    url = granule["opendap"] + ".dap?" + urllib.parse.urlencode(
        {"dap4.ce": ";".join(DAP4_VARIABLES)})
    blob = _get(url, headers=_bearer_header())
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached.write_bytes(blob)
    time.sleep(pause_s)
    return blob, len(blob)


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def granule_rows(blob: bytes, title: str):
    """One granule's unique PVT epochs as a dict of columns.

    The sample dimension runs at the 2 Hz DDM rate while the navigation solution
    updates at 1 Hz, so consecutive samples repeat a PVT epoch; only the first
    of each repeat is kept.  ``pvt_timestamp_utc`` is seconds since the
    granule's ``time_coverage_start``, which the same DMR carries, so the
    absolute time is formed here and no assumption about the day boundary is
    needed.
    """
    fields, dmr = dap4_read(blob)
    coverage_start = dmr_attribute(dmr, "time_coverage_start")
    if coverage_start is None:
        raise RuntimeError(f"{title}: DMR has no time_coverage_start")
    epoch_unix = iso_to_unix(coverage_start)

    t_rel = fields["pvt_timestamp_utc"]
    keep = np.ones(t_rel.size, dtype=bool)
    keep[1:] = t_rel[1:] != t_rel[:-1]

    pos = np.stack([fields[f"sc_pos_{a}_pvt"] for a in "xyz"], axis=1)[keep]
    vel = np.stack([fields[f"sc_vel_{a}_pvt"] for a in "xyz"], axis=1)[keep]
    valid = (pos != POS_FILL).all(axis=1) & (vel != VEL_FILL).all(axis=1)

    return {
        "spacecraft_num": np.full(keep.sum(), int(fields["spacecraft_num"]), dtype=np.int8),
        "pvt_unix_seconds": epoch_unix + t_rel[keep],
        "sc_pos_x_pvt_m": pos[:, 0].astype(np.float64),
        "sc_pos_y_pvt_m": pos[:, 1].astype(np.float64),
        "sc_pos_z_pvt_m": pos[:, 2].astype(np.float64),
        "sc_vel_x_pvt_mps": vel[:, 0].astype(np.float64),
        "sc_vel_y_pvt_mps": vel[:, 1].astype(np.float64),
        "sc_vel_z_pvt_mps": vel[:, 2].astype(np.float64),
        "valid_fix": valid,
        "source_file": np.full(keep.sum(), title + ".nc", dtype=object),
        "time_coverage_start": np.full(keep.sum(), coverage_start, dtype=object),
    }


def iso_to_unix(stamp: str) -> float:
    """Seconds since 1970-01-01T00:00:00Z of an ISO-8601 UTC stamp, exactly.

    ``datetime`` truncates to microseconds and these stamps carry nanoseconds,
    so the sub-second part is parsed as a decimal string instead.  No leap
    second falls in 2025, so the UTC-to-Unix mapping over this window is a plain
    offset; the resulting number is a UTC-based count and is converted to
    ephemeris time by SPICE downstream.
    """
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})(\.\d+)?Z?$", stamp.strip())
    if m is None:
        raise ValueError(f"unrecognized UTC stamp: {stamp!r}")
    y, mo, d, hh, mm, ss = (int(m.group(i)) for i in range(1, 7))
    frac = float(m.group(7)) if m.group(7) else 0.0
    days = days_from_civil(y, mo, d)
    return days * 86400.0 + hh * 3600.0 + mm * 60.0 + ss + frac


def days_from_civil(y: int, m: int, d: int) -> int:
    """Days from 1970-01-01 to a proleptic-Gregorian date (Howard Hinnant's
    civil_from_days inverse, the standard branch-free form)."""
    y -= m <= 2
    era = (y if y >= 0 else y - 399) // 400
    yoe = y - era * 400
    doy = (153 * (m + (-3 if m > 2 else 9)) + 2) // 5 + d - 1
    doe = yoe * 365 + yoe // 4 - yoe // 100 + doy
    return era * 146097 + doe - 719468


# ---------------------------------------------------------------------------
# Self-test: the decoder, without the network
# ---------------------------------------------------------------------------

def _chunk(payload: bytes, flags: int) -> bytes:
    return bytes([flags]) + len(payload).to_bytes(3, "big") + payload


def self_test() -> int:
    """Check the DAP4 decoder without touching the archive.

    Two checks.  The first builds a chunked DAP4 response by hand from a known
    DMR and known arrays and reads it back, which pins the chunk framing, the
    variable ordering, the endianness flag and the per-variable checksum skip.
    The second, when the local files are there, decodes a cached response and
    compares it with the independently produced FM04 mirror table: that is the
    real check, because it is against data this pipeline had no part in making.
    """
    failures = 0

    dmr = ('<?xml version="1.0" encoding="ISO-8859-1"?>\n'
           '<Dataset xmlns="http://xml.opendap.org/ns/DAP/4.0#" name="t">\n'
           '    <Dimension name="sample" size="4"/>\n'
           '    <Float64 name="pvt_timestamp_utc"><Dim name="/sample"/>\n'
           '        <Attribute name="time_coverage_start" type="String">\n'
           '            <Value>2025-06-06T00:00:00.499261819Z</Value>\n'
           '        </Attribute>\n'
           '    </Float64>\n'
           '    <Int32 name="sc_pos_x_pvt"><Dim name="/sample"/></Int32>\n'
           '    <Int8 name="spacecraft_num"></Int8>\n'
           '</Dataset>\n\r\n')
    times = np.array([-0.5, 0.5, 1.5, 2.5], dtype="<f8")
    xs = np.array([1, 2, 3, POS_FILL], dtype="<i4")
    data = times.tobytes() + b"\0\0\0\0" + xs.tobytes() + b"\0\0\0\0" \
        + np.int8(4).tobytes() + b"\0\0\0\0"
    blob = _chunk(dmr.encode("latin-1"), 0x04) + _chunk(data[:9], 0x04) \
        + _chunk(data[9:], 0x05)
    got, back = dap4_read(blob)
    for name, want in (("pvt_timestamp_utc", times), ("sc_pos_x_pvt", xs)):
        if not np.array_equal(got[name], want):
            print(f"FAIL round trip of {name}: {got[name]} != {want}")
            failures += 1
    if int(got["spacecraft_num"]) != 4:
        print(f"FAIL round trip of spacecraft_num: {got['spacecraft_num']}")
        failures += 1
    if dmr_attribute(back, "time_coverage_start") != "2025-06-06T00:00:00.499261819Z":
        print("FAIL time_coverage_start not recovered from the DMR")
        failures += 1
    if iso_to_unix("2025-06-06T00:00:00.499261819Z") != 1749168000.0 + 0.499261819:
        print("FAIL iso_to_unix on the granule's own coverage stamp")
        failures += 1
    if iso_to_unix("1970-01-01T00:00:00Z") != 0.0:
        print("FAIL iso_to_unix at the Unix epoch")
        failures += 1

    cached = TELEM_DIR / "l1_opendap_cache" / \
        "cyg04.ddmi.s20250606-000000-e20250606-235959.l1.power-brcs.a32.d33.dap"
    mirror = TELEM_DIR / "cyg04_nasa_pvt_96hr.feather"
    if cached.exists() and mirror.exists():
        rows = granule_rows(cached.read_bytes(), cached.stem)
        table = feather.read_table(mirror)
        want = {c: np.asarray(table.column(c)) for c in
                ("pvt_unix_seconds", "sc_pos_x_pvt_m", "sc_vel_x_pvt_mps")}
        key = np.round(rows["pvt_unix_seconds"], 3)
        index = {k: i for i, k in enumerate(key)}
        checked = 0
        for j, k in enumerate(np.round(want["pvt_unix_seconds"], 3)):
            i = index.get(k)
            if i is None:
                continue
            checked += 1
            for column in ("sc_pos_x_pvt_m", "sc_vel_x_pvt_mps"):
                if rows[column][i] != want[column][j]:
                    print(f"FAIL {column} at {k}: {rows[column][i]} != {want[column][j]}")
                    failures += 1
                    break
        if checked < 80_000:
            print(f"FAIL only {checked} epochs of the mirror matched the decoded granule")
            failures += 1
        else:
            print(f"decoded granule agrees with the FM04 mirror on {checked} epochs")
    else:
        print("skipped the mirror comparison: cached response or FM04 mirror absent")

    print("self-test: " + ("ok" if failures == 0 else f"{failures} failures"))
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true",
                        help="check the DAP4 decoder and exit, without any network access")
    parser.add_argument("--out", default=str(TELEM_DIR / "cygnss_l1_pvt_ecef_20250606_96hr.feather"),
                        help="Feather table of the Earth-fixed states to write")
    parser.add_argument("--cache", default=str(TELEM_DIR / "l1_opendap_cache"),
                        help="directory holding the verbatim OPeNDAP responses")
    parser.add_argument("--pause-s", type=float, default=1.0,
                        help="delay after each request that actually hit the archive")
    args = parser.parse_args()
    if args.self_test:
        return self_test()

    print(f"CMR: {SHORT_NAME} granules over {WINDOW_START} .. {WINDOW_END}")
    granules = cmr_granules(WINDOW_START, WINDOW_END)
    total_mb = sum(g["size_mb"] for g in granules)
    spacecraft = sorted({g["title"].split(".")[0] for g in granules})
    print(f"  {len(granules)} granules, {len(spacecraft)} spacecraft ({', '.join(spacecraft)}), "
          f"{total_mb / 1024:.2f} GB if downloaded whole")

    cache_dir = Path(args.cache)
    columns: dict[str, list] = {}
    moved = 0
    reused = 0
    for granule in sorted(granules, key=lambda g: (g["title"].split(".")[0], g["time_start"])):
        blob, got = fetch_granule_variables(granule, cache_dir, args.pause_s)
        moved += got
        reused += 0 if got else len(blob)
        rows = granule_rows(blob, granule["title"])
        for key, value in rows.items():
            columns.setdefault(key, []).append(value)
        print(f"  {granule['title']:<62} {len(rows['pvt_unix_seconds']):>6} epochs  "
              f"{len(blob) / 1e6:6.2f} MB {'downloaded' if got else 'cached'}  "
              f"(granule is {granule['size_mb']:.0f} MB)")

    merged = {k: np.concatenate(v) for k, v in columns.items()}
    table = pa.table({k: pa.array(v.tolist() if v.dtype == object else v)
                      for k, v in merged.items()})
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    feather.write_feather(table, out, compression="lz4")

    print(f"\nbytes moved from the archive this run : {moved / 1e6:.2f} MB")
    print(f"bytes served from the local cache     : {reused / 1e6:.2f} MB")
    print(f"whole-granule equivalent              : {total_mb / 1024:.2f} GB "
          f"({100.0 * (moved + reused) / (total_mb * 1e6):.2f}% of it subset out)")
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB, {table.num_rows} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
