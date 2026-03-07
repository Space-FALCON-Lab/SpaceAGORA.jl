const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

using SpaceAGORA

if !isdefined(@__MODULE__, :SimulationEngine)
    const SimulationEngine = SpaceAGORA.SimulationEngine
end
if !isdefined(@__MODULE__, :SimulationModel)
    const SimulationModel = SpaceAGORA.SimulationModel
end
using .SimulationModel

if !isdefined(@__MODULE__, :SM)
    const SM = SimulationModel
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SpaceAGORA.run_simulation
end
if !isdefined(@__MODULE__, :quat_mult)
    const quat_mult = SimulationModel.quat_mult
end
if !isdefined(@__MODULE__, :SPICE_PATH)
    const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
end

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft, run_and_report
