const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Keep the reporting fixture isolated from other benchmark harness definitions
# in the shared contract runner. The child executes this same entry point.
const _COPY_OVERHEAD_INLINE = abspath(PROGRAM_FILE) == abspath(@__FILE__) ||
    get(ENV, "SPACEAGORA_COPY_OVERHEAD_GATE_INLINE", "") == "1"

if !_COPY_OVERHEAD_INLINE
    child = addenv(`$(Base.julia_cmd()) --startup-file=no --project=$(REPO_ROOT) $(@__FILE__)`,
        "SPACEAGORA_COPY_OVERHEAD_GATE_INLINE" => "1")
    success(child) || error("Runtime analysis copy-overhead gate failed in its isolated process")
else
    include(joinpath(@__DIR__, "runtime_analysis_copy_overhead_checks.jl"))
end
println("runtime_analysis_copy_overhead_gate_ok")
