module ComparisonLayout

export CODE_ROOT, LIBRARY_ROOT, DEFAULT_OUTPUT, VALIDATION_ROOT

const CODE_ROOT = @__DIR__
const SCITECH_ROOT = normpath(joinpath(CODE_ROOT, "..", ".."))
const LIBRARY_ROOT = joinpath(SCITECH_ROOT, "9_SpaceAGORA.jl-main_full integration")
const DEFAULT_OUTPUT = normpath(joinpath(CODE_ROOT, "..", "Output"))
const VALIDATION_ROOT = joinpath(SCITECH_ROOT, "Tests",
    "Output_3_rangeOnly_vs_positiveAlongTrack", "validation")

end
