const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

required = (
    joinpath("src", "io", "config", "io_config.jl"),
    joinpath("src", "io", "serialization", "io_serialization.jl"),
    joinpath("src", "io", "logging", "io_logging.jl"),
    joinpath("src", "io", "outputs", "io_outputs.jl"),
)

missing = String[]
for rel in required
    isfile(joinpath(REPO_ROOT, rel)) || push!(missing, rel)
end

if !isempty(missing)
    println("io_surface_missing_files:")
    for m in missing
        println("  - $m")
    end
    error("IO surface non-empty gate failed")
end

println("io_surface_nonempty_gate_ok")
