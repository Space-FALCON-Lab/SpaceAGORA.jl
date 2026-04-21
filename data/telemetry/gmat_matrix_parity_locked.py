import os
import subprocess
from pathlib import Path

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Keep this script and all outputs self-contained in the selected base directory.
BASE_DIR = Path(os.path.expanduser("~/Documents/Paper_Graphics/refactor_plotting/gmat_validation_plotting"))
OUTPUT_DIR = BASE_DIR / "gmat_results"
SCRIPT_DIR = BASE_DIR / "gmat_scripts"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SCRIPT_DIR.mkdir(parents=True, exist_ok=True)

GMAT_BIN = Path("../../../GMAT_2025/R2025a/bin/GMAT_Beta").expanduser().resolve()
if not GMAT_BIN.is_file():
    raise FileNotFoundError(f"GMAT executable not found: {GMAT_BIN}")

GMAT_ROOT = GMAT_BIN.parent.parent
GRAVITY_DIR = GMAT_ROOT / "data" / "gravity"

MAX_STEP_S = 10.0
ACCURACY = "1e-12"
MIN_STEP_S = 0.001
PROPAGATION_END_S = 1_000_000.0

# Explicit parity matrix: (degree, order)
# - J0: (0,0)
# - J2-only zonal: (2,0)
# - High-order full: (50,50)
HARMONIC_CASES = [(0, 0), (2, 0), (50, 50)]

THIRD_BODY_EFFECTS = {
    "Earth": "Sun, Luna",
    "Mars": "Sun",
    "Venus": "Sun",
    "Titan": "Saturn",
}

# Use absolute potential file paths to eliminate GMAT path resolution ambiguity.
POTENTIAL_FILENAMES = {
    "Earth": "JGM2.cof",
    "Mars": "Mars50c.cof",
    "Venus": "MGNP180U.cof",
    "Titan": None,
}

INITIAL_KEPLERIAN_STATES = {
    "Earth": {
        "epoch": "01 Jan 2026 12:00:00.000",
        "sma": 7000.0,
        "ecc": 0.00001,
        "inc": 28.5,
        "raan": 45.0,
        "aop": 0.0,
        "ta": 0.0,
    },
    "Mars": {
        "epoch": "01 Jan 2026 12:00:00.000",
        "sma": 4000.0,
        "ecc": 0.002,
        "inc": 45.0,
        "raan": 0.0,
        "aop": 0.0,
        "ta": 90.0,
    },
    "Venus": {
        "epoch": "01 Jan 2026 12:00:00.000",
        "sma": 7000.0,
        "ecc": 0.001,
        "inc": 30.0,
        "raan": 0.0,
        "aop": 0.0,
        "ta": 180.0,
    },
    "Titan": {
        "epoch": "01 Jan 2026 12:00:00.000",
        "sma": 6000.0,
        "ecc": 0.005,
        "inc": 20.0,
        "raan": 0.0,
        "aop": 0.0,
        "ta": 270.0,
    },
}

GMAT_TEMPLATE = """
%--------------------------------------------------------------------------
% Coordinate System: Planet-Centered Inertial
%--------------------------------------------------------------------------
Create CoordinateSystem PlanetInertial;
PlanetInertial.Origin = {body};
PlanetInertial.Axes   = BodyInertial;

%--------------------------------------------------------------------------
% Spacecraft Initial Conditions (Keplerian)
%--------------------------------------------------------------------------
Create Spacecraft Sat;
Sat.CoordinateSystem = PlanetInertial;
Sat.DateFormat = UTCGregorian;
Sat.Epoch = '{epoch}';
Sat.DisplayStateType = Keplerian;

Sat.SMA  = {sma};
Sat.ECC  = {ecc};
Sat.INC  = {inc};
Sat.RAAN = {raan};
Sat.AOP  = {aop};
Sat.TA   = {ta};

%--------------------------------------------------------------------------
% Force Model & Propagator
%--------------------------------------------------------------------------
Create ForceModel MultiFM;
MultiFM.CentralBody = {body};
MultiFM.PrimaryBodies = {{{body}}};
MultiFM.GravityField.{body}.Degree = {degree};
MultiFM.GravityField.{body}.Order = {order};
MultiFM.GravityField.{body}.StmLimit = 100;
{gravity_potential_line}
MultiFM.GravityField.{body}.TideModel = 'None';
MultiFM.AtmosphereModel.Drag = False;
MultiFM.SRP = False;
{point_masses_line}

Create Propagator DefaultProp;
DefaultProp.FM = MultiFM;
DefaultProp.Type = RungeKutta78;
DefaultProp.MaxStep = {max_step};
DefaultProp.Accuracy = {accuracy};
DefaultProp.MinStep = {min_step};

Create ReportFile DataLog;
DataLog.Filename = '{full_output_path}';
DataLog.Precision = 16;
DataLog.Add = {{
    Sat.ElapsedSecs,
    Sat.SMA, Sat.ECC, Sat.INC,
    Sat.PlanetInertial.X, Sat.PlanetInertial.Y, Sat.PlanetInertial.Z,
    Sat.PlanetInertial.VX, Sat.PlanetInertial.VY, Sat.PlanetInertial.VZ
}};
DataLog.WriteHeaders = On;
DataLog.ZeroFill = Off;
DataLog.FixedWidth = Off;
DataLog.Delimiter = ',';

BeginMissionSequence;
Propagate DefaultProp(Sat) {{Sat.ElapsedSecs = {end_s}}};
"""


def resolve_potential_path(body: str) -> Path | None:
    filename = POTENTIAL_FILENAMES.get(body)
    if filename is None:
        return None
    candidate = GRAVITY_DIR / filename
    if not candidate.is_file():
        raise FileNotFoundError(
            f"Potential file for {body} not found: {candidate}. "
            f"Update POTENTIAL_FILENAMES/GRAVITY_DIR to match your GMAT install."
        )
    return candidate


def build_force_model_lines(body: str, degree: int, order: int, third_bodies: str, potential_path: Path | None):
    if degree == 0 and order == 0:
        gravity_potential_line = "% No potential file needed for pure point-mass gravity"
    else:
        if potential_path is None:
            raise ValueError(
                f"{body} requested Degree/Order={degree}/{order} but no potential file is configured."
            )
        gravity_potential_line = f"MultiFM.GravityField.{body}.PotentialFile = '{potential_path.as_posix()}';"

    if third_bodies:
        point_masses_line = f"MultiFM.PointMasses = {{{third_bodies}}};"
    else:
        point_masses_line = "MultiFM.PointMasses = {};"

    return gravity_potential_line, point_masses_line


def run_case(body: str, state: dict, degree: int, order: int, third_body_on: bool):
    third_bodies = THIRD_BODY_EFFECTS[body] if third_body_on else ""
    potential_path = resolve_potential_path(body)

    # Titan has no potential file by default in this setup; skip harmonics safely.
    if potential_path is None and not (degree == 0 and order == 0):
        print(f"Skipping {body} J{degree} TB{third_body_on}: no potential file configured")
        return

    sim_id = f"Sim_{body}_J{degree}_TB{third_body_on}"
    output_path = (OUTPUT_DIR / f"{sim_id}.csv").resolve()
    script_path = (SCRIPT_DIR / f"{sim_id}.script").resolve()

    gravity_potential_line, point_masses_line = build_force_model_lines(
        body, degree, order, third_bodies, potential_path
    )

    script_content = GMAT_TEMPLATE.format(
        body=body,
        epoch=state["epoch"],
        sma=state["sma"],
        ecc=state["ecc"],
        inc=state["inc"],
        raan=state["raan"],
        aop=state["aop"],
        ta=state["ta"],
        degree=degree,
        order=order,
        gravity_potential_line=gravity_potential_line,
        point_masses_line=point_masses_line,
        max_step=MAX_STEP_S,
        accuracy=ACCURACY,
        min_step=MIN_STEP_S,
        end_s=PROPAGATION_END_S,
        full_output_path=output_path.as_posix(),
    )

    script_path.write_text(script_content, encoding="ascii")

    print(f"Running {sim_id} ...")
    subprocess.run(
        [GMAT_BIN.as_posix(), "--run", script_path.as_posix(), "--minimize", "--exit"],
        check=True,
    )


if __name__ == "__main__":
    print(f"GMAT executable : {GMAT_BIN}")
    print(f"GMAT root       : {GMAT_ROOT}")
    print(f"Gravity dir     : {GRAVITY_DIR}")
    print(f"Output dir      : {OUTPUT_DIR}")
    print(f"Script dir      : {SCRIPT_DIR}")

    for body, state in INITIAL_KEPLERIAN_STATES.items():
        for degree, order in HARMONIC_CASES:
            for third_body_on in (False, True):
                run_case(body, state, degree, order, third_body_on)

    print("Done.")
