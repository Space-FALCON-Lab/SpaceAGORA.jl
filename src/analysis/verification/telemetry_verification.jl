module TelemetryVerification

export VerificationRequest, VerificationResult, run_verification, run_verification_cli, run_study

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output")
const DEFAULT_MANIFEST_PATH = joinpath(REPO_ROOT, "test", "telemetry_benchmark_manifest.toml")
const STRICT_REL_ORBIT = 1e-7
const STRICT_ABS_ORBIT = 1e-9
const STRICT_REL_ATM = 1e-7
const STRICT_ABS_ATM = 1e-9
const STRICT_DT_ORBIT = 60.0
const STRICT_DT_ATM = 0.2
const TELEMETRY_SOLVER_MAXITERS_QUICK_DEFAULT = 5_000_000
const TELEMETRY_SOLVER_MAXITERS_FULL_DEFAULT = 20_000_000

using Statistics
using Dates
using Printf
using LinearAlgebra
using StaticArrays
using SPICE
using CSV
using DataFrames
using Arrow
using TOML

if isdefined(parentmodule(@__MODULE__), :SimulationEngine)
    using ..SimulationEngine
    const SimulationModel = SimulationEngine.SimulationModel
else
    error(
        "TelemetryVerification requires SimulationEngine in the parent module. " *
        "Load src/simulation/engine/simulation_engine.jl before including telemetry_verification.jl."
    )
end
using .SimulationModel
include(joinpath(REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")

include(joinpath(@__DIR__, "telemetry_verification", "types.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "manifest_parsing.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "example_support.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "scenario_builders.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "telemetry_loading.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "comparison_metrics.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "calibration.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "error_tables.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "reporting.jl"))
include(joinpath(@__DIR__, "telemetry_verification", "runner.jl"))

# Architecture contract: runner delegates simulation execution to SimulationEngine.run_simulation.
end # module TelemetryVerification
