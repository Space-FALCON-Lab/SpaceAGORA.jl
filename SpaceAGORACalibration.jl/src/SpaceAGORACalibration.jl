module SpaceAGORACalibration

include("Spec.jl")
include("ParamSpace.jl")
include("Backend.jl")
include("Objective.jl")
include("GlobalSearch.jl")
include("LocalRefine.jl")
include("Robustness.jl")
include("Store.jl")
include("Report.jl")
include("Engine.jl")

using .Spec: CalibrationSpec, BudgetSpec, ParameterSpec, ParameterKind, load_spec, default_spec, validate_spec
using .ParamSpace: Candidate
using .Backend: AbstractBackend, MockBackend, CommandBackend, BackendEvaluation, evaluate_candidate
using .Engine: CalibrationResult, run_calibration

export CalibrationSpec, BudgetSpec, ParameterSpec, ParameterKind
export Candidate
export AbstractBackend, MockBackend, CommandBackend, BackendEvaluation, evaluate_candidate
export CalibrationResult, run_calibration
export load_spec, default_spec, validate_spec

end # module SpaceAGORACalibration
