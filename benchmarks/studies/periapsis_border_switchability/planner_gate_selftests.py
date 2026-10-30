"""Python gate list of the periapsis-border study (Revision 22 item 6, Revision 23 item 7, Revision 24 item 6, Revision 25): runs every geometry/target selftest and the default-off bit-identity digest.
Usage: python3 -I planner_gate_selftests.py [--fast]   (--fast skips the digest and the long Revision 20 selftest). Exit code 0 = all pass. Unset SPACEAGORA_RA_TARGET_M for the digest."""
import os, pathlib, subprocess, sys

HERE = pathlib.Path(__file__).resolve().parent
GATES = [("planner_rev22_rendezvous.py", ["selftest"], False), ("planner_rev23_geometry.py", ["selftest"], False), ("planner_rev24_targets.py", ["selftest"], False), ("planner_select_geometry.py", ["selftest"], False),
         ("planner_rev20_ic.py", ["--selftest"], True), ("planner_p25_digest.py", ["--check", str(HERE / "planner_p25_digest.sha256")], True)]


def main(fast=False):
    bad = []
    for f, args, slow in GATES:
        if fast and slow: continue
        r = subprocess.run([sys.executable, "-I", str(HERE / f)] + args, capture_output=True, text=True, env={k: v for k, v in os.environ.items() if k != "SPACEAGORA_RA_TARGET_M"})
        print(("PASS " if r.returncode == 0 else "FAIL ") + f + " " + " ".join(args)); 
        if r.returncode != 0: bad.append(f); print(r.stdout[-400:], r.stderr[-600:])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main("--fast" in sys.argv))
