#!/usr/bin/env python3
"""Convert the raw CYGNSS FM4 48hr PVT CTAS export (xlsx) to a feather file.

Stage 1 of the CYGNSS 48hr data pipeline (see CYGNSS/README in the private
telemetry repo). Reads the SwRI CTAS query export
`query-45991-FM4_202506_08_PVT.xlsx` and writes an Arrow/feather file with the
ORIGINAL column names preserved plus nothing added or filtered:

  TIME OFFSET                              seconds from the first sample epoch
  ENG_PVT                                  UTC timestamp string of each sample
  OBS4.ENG_PVT.DDMI_PVT_SCPOS_{X,Y,Z} (m)  spacecraft position, WGS84 ECEF
  OBS4.ENG_PVT.DDMI_PVT_SCVEL_{X,Y,Z} (m/s) spacecraft velocity, WGS84 ECEF
  OBS4.ENG_PVT.DDMI_PVT_NUMSATS (numsats)  GPS satellites used in the fix
  OBS4.ENG_PVT.DDMI_PVT_GDOP (GDOP)        GDOP scaled x10 (per data provider)
  OBS4.ENG_PVT.DDMI_PVT_VALID (null)       PVT validity flag (2 == valid)

IMPORTANT: the GPS PVT coordinates are ECEF, NOT inertial. The verification
harness interprets telemetry columns as J2000 by default, so stage 2
(add_cygnss_48hr_eci_columns.jl) must be run afterwards to append
pos_ii_*/vel_ii_* J2000 columns via the repo's own SPICE ITRF93 convention.

Usage: python3 convert_cygnss_48hr_xlsx.py [input.xlsx] [output.feather]
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_IN = HERE / "CYGNSS" / "raw_email" / "query-45991-FM4_202506_08_PVT.xlsx"
DEFAULT_OUT = HERE / "CYGNSS" / "cygnss_48hr_raw_ecef.feather"


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_IN
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUT
    print(f"reading {src} ...")
    df = pd.read_excel(src, engine="openpyxl")
    df = df.sort_values("TIME OFFSET").reset_index(drop=True)

    t = df["TIME OFFSET"].to_numpy()
    dt = t[1:] - t[:-1]
    n_gaps = int((dt != 1).sum())
    print(f"rows={len(df)}  t=[{t[0]}..{t[-1]}] s  non-1s-steps={n_gaps}")
    print(f"epoch (first ENG_PVT stamp): {df['ENG_PVT'].iloc[0]}")
    valid = df["OBS4.ENG_PVT.DDMI_PVT_VALID (null)"]
    print(f"PVT_VALID==2: {(valid == 2).sum()}/{len(df)}")

    dst.parent.mkdir(parents=True, exist_ok=True)
    df.to_feather(dst)
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
