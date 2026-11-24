"""Join the two 100-orbit paper tables without rounding their values."""

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
KEY = ("target_altitude_km", "delta_i_deg", "helpers", "scheduler")
RULES = {"gve_sma": "a-max", "gve_inc": "i-max"}


def read_cases(filename):
    cases = {}
    with (HERE / filename).open(newline="") as file:
        for row in csv.DictReader(file):
            key = tuple(row[name] for name in KEY)
            if row["orbits"] != "100" or key in cases:
                raise ValueError(f"{filename}: unexpected window or duplicate case {key}")
            cases[key] = row
    return cases


def main():
    dv = read_cases("cumulative_dv_RT_N_branches_100orbits.csv")
    oe = read_cases("orbital_element_deltas_100orbits.csv")
    if dv.keys() != oe.keys() or len(dv) != 162:
        raise ValueError("Expected the same 162 unique cases in both paper tables")

    rows = []
    for key, v in dv.items():
        e = oe[key]
        rows.append((*key[:3], RULES[key[3]], v["dv_R_mps"], v["dv_T_mps"],
                     v["dv_N_positive_mps"], v["dv_N_negative_mps"], e["da_m"],
                     e["de"], e["di_deg"], e["draan_deg"] if e["raan_defined"] == "true" else ""))

    with (HERE / "combined_paper_data_100orbits.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("target altitude (km)", "inclination (deg)", "number of helpers", "rule",
                         "Delta v_R (m/s)", "Delta v_T (m/s)", "Delta v_N+ (m/s)",
                         "Delta v_N- (m/s)", "Delta a (m)", "Delta e", "Delta i (deg)",
                         "Delta Omega (deg)"))
        writer.writerows(rows)


if __name__ == "__main__":
    main()
