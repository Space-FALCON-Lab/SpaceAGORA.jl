const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance")
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
const EARTH_GRAM_SURROGATE_FILE = joinpath(REPO_ROOT, "data", "GRAM_surrogate", "earth_surrogate.jls")
const PERF_BASELINE_SCENARIO = "single_j2"
const _PERF_POLICY_ENV_NAMES = (
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE",
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL",
    "SPACEAGORA_CONTROL_CALLBACK_PARALLEL",
    "SPACEAGORA_MULTIBODY_PARALLEL",
)
const _PERF_POLICY_ENV_BASELINE = Dict{String, Union{Nothing, String}}(
    name => (haskey(ENV, name) ? String(ENV[name]) : nothing)
    for name in _PERF_POLICY_ENV_NAMES
)
const _PERF_THREADS_BACKEND_WARNING_EMITTED = Ref(false)

using SpaceAGORA
using GRAMSuite

import SpaceAGORA.ParallelProfiles
using .ParallelProfiles

using CSV
using DataFrames
using Dates
using Distributed
using LinearAlgebra
using Random
using Sockets
using SPICE
using StaticArrays
using Statistics

if myid() == 1
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
    using Plots
end

import SpaceAGORA.SimulationModel
using .SimulationModel

import SpaceAGORA.SimulationEngine

# GRAMSuite exports InitialTime, GRAMAtmosphereModel, and GRAMAtmosphereModelSurrogate,
# which collide with the SpaceAGORA types brought in by `using .SimulationModel`.
# Explicit imports resolve the ambiguity in favour of the SpaceAGORA versions.
import SpaceAGORA.SimulationModel.SimConfig: InitialTime
import SpaceAGORA.SimulationModel.EnvironmentModels: GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end


include(joinpath(@__DIR__, "case_catalog.jl"))
include(joinpath(@__DIR__, "measurement.jl"))
include(joinpath(@__DIR__, "reporting.jl"))
include(joinpath(@__DIR__, "cli.jl"))
