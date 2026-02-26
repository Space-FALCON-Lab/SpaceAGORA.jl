from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from orchestrator.app.config import Settings
from orchestrator.app.errors import ExecutionError
from orchestrator.app.models import SimulationSpec


@dataclass
class ExecutionResult:
    success: bool
    run_dir: Path
    stdout_path: Path
    stderr_path: Path
    spec_snapshot_path: Path
    script_path: Path
    result_csv_path: Optional[Path]
    error_message: Optional[str] = None


class JuliaSimulationRunner:
    def __init__(self, settings: Settings) -> None:
        self._settings = settings

    def run(self, job_id: str, spec: SimulationSpec, wall_timeout_sec: float) -> ExecutionResult:
        run_dir = self._settings.worker_tmp_dir / job_id
        run_dir.mkdir(parents=True, exist_ok=True)

        spec_snapshot_path = run_dir / "simulation_spec.json"
        spec_snapshot_path.write_text(spec.model_dump_json(indent=2), encoding="utf-8")

        script_path = run_dir / "run_job.jl"
        script_path.write_text(self._render_julia_script(spec), encoding="utf-8")

        stdout_path = run_dir / "stdout.log"
        stderr_path = run_dir / "stderr.log"

        cmd = [
            "julia",
            f"--project={self._settings.repo_root / '.ABTS'}",
            str(script_path),
        ]

        with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            try:
                completed = subprocess.run(
                    cmd,
                    cwd=run_dir,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    timeout=wall_timeout_sec,
                    check=False,
                    text=True,
                )
            except subprocess.TimeoutExpired:
                message = f"simulation timed out after {wall_timeout_sec} seconds"
                return ExecutionResult(
                    success=False,
                    run_dir=run_dir,
                    stdout_path=stdout_path,
                    stderr_path=stderr_path,
                    spec_snapshot_path=spec_snapshot_path,
                    script_path=script_path,
                    result_csv_path=None,
                    error_message=message,
                )

        result_csv_path = run_dir / "simulation_results.csv"
        if completed.returncode != 0:
            message = f"julia exited with return code {completed.returncode}"
            return ExecutionResult(
                success=False,
                run_dir=run_dir,
                stdout_path=stdout_path,
                stderr_path=stderr_path,
                spec_snapshot_path=spec_snapshot_path,
                script_path=script_path,
                result_csv_path=result_csv_path if result_csv_path.exists() else None,
                error_message=message,
            )

        return ExecutionResult(
            success=True,
            run_dir=run_dir,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            spec_snapshot_path=spec_snapshot_path,
            script_path=script_path,
            result_csv_path=result_csv_path if result_csv_path.exists() else None,
        )

    def _render_julia_script(self, spec: SimulationSpec) -> str:
        repo = str(self._settings.repo_root).replace('"', '\\"')
        spice = spec.file_paths.spice.replace('"', '\\"')
        gram_data = spec.file_paths.gram_data.replace('"', '\\"')
        gram_py = spec.file_paths.gram_py.replace('"', '\\"')

        gravity_ctor = {
            "inverse_squared": "InverseSquaredGravityModel()",
            "inverse_squared_j2": "InverseSquaredJ2GravityModel()",
            "constant": "ConstantGravityModel()",
        }[spec.models.gravity_model]

        if spec.models.density_model == "none":
            density_ctor = "NoAtmosphereModel()"
        elif spec.models.density_model == "exponential":
            density_ctor = "ExponentialAtmosphereModel(planet.ρ_ref, planet.h_ref, planet.H)"
        elif spec.models.density_model == "gram":
            density_ctor = (
                f'GRAMAtmosphereModel(planet_name="{spec.planet.lower()}", '
                f'gram_root_directory="{gram_py}", gram_data_directory="{gram_data}", spice_directory="{spice}", initial_time=time)'
            )
        else:
            raise ExecutionError(f"Unsupported density model in script renderer: {spec.models.density_model}")

        if spec.models.aerodynamic_model == "mach_dependent":
            aero_ctor = "AerodynamicCoefficientfM()"
        else:
            aero_ctor = "AerodynamicCoefficientConstant()"

        include_aero = spec.models.density_model != "none"
        effectors_expr = (
            "(gravEffector, aeroEffector)" if include_aero else "(gravEffector,)"
        )

        srp_bool = "true" if spec.models.srp else "false"
        wind_bool = "true" if spec.models.wind else "false"
        keplerian_bool = "true" if spec.keplerian else "false"
        orientation_bool = "true" if spec.orientation_sim else "false"

        return f"""
include(\"{repo}/src/simulation_model/SimulationModel.jl\")
include(\"{repo}/src/simulation/Run.jl\")

using .SimulationModel
using StaticArrays

planet = {spec.planet}(\"\", \"{spice}\")
gravEffector = {gravity_ctor}
""" + (
            f"aeroEffector = {aero_ctor}\n"
            if include_aero
            else ""
        ) + f"""
dynamic_effectors = {effectors_expr}

main_bus = Link{{0}}(
    root=true,
    r=SVector{{3, Float64}}(0.0, 0.0, 0.0),
    q=SVector{{4, Float64}}(0.0, 0.0, 0.0, 1.0),
    ṙ=SVector{{3, Float64}}(0.0, 0.0, 0.0),
    dims=SVector{{3, Float64}}(2.2, 2.6, 1.7),
    ref_area=2.2 * 1.7,
    m=391.0
)

ic = InitialCondition(
    ra={spec.initial_conditions.ra_m},
    rp={spec.initial_conditions.rp_m},
    i={spec.initial_conditions.i_deg},
    ω={spec.initial_conditions.omega_deg},
    Ω={spec.initial_conditions.raan_deg},
    ν={spec.initial_conditions.true_anomaly_deg}
)

spacecraft = SpacecraftModel(root=main_bus, prop_mass=50.0, initial_condition=ic, id=1)
dynamics_model = DynamicsModel([spacecraft], dynamic_effectors)

time = InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0)
density_model = {density_ctor}
thermal_model = MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet)
environment_model = EnvironmentModel(planet=planet, EI={spec.entry_interface_km}, density_model=density_model, thermal_model=thermal_model, wind={wind_bool})

control_model = ControlModel(control_effectors=(), control_rates=Float64[])
guidance_model = GuidanceModel(guidance_effectors=(), guidance_rates=Float64[])
navigation_model = NavigationModel(navigation_effectors=(), navigation_rates=Float64[])

mission_configuration = MissionConfiguration(
    mission_type=\"Time\",
    keplerian={keplerian_bool},
    number_of_orbits=1,
    mission_time={spec.mission_time_sec},
    orientation_sim={orientation_bool},
    num_steps_to_save={spec.num_steps_to_save}
)

sim_settings = SimulationSettings(
    results={str(spec.outputs.results).lower()},
    verbose=false,
    results_directory=\"output\",
    generate_plots={str(spec.outputs.generate_plots).lower()},
    generate_filenames=false,
    normalize=true,
    save_csv={str(spec.outputs.save_csv).lower()}
)

args = SimulationConfiguration(
    dynamics_model=dynamics_model,
    initial_time=time,
    environment_model=environment_model,
    control_model=control_model,
    guidance_model=guidance_model,
    mission_configuration=mission_configuration,
    navigation_model=navigation_model,
    simulation_settings=sim_settings
)

run_analysis(args)
println(\"RUN_OK\")
"""
