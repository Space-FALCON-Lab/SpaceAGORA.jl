const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const STRUCTURE_ROOT = joinpath(REPO_ROOT, "src", "vehicle", "structure")
const SPACECRAFT_ROOT = joinpath(REPO_ROOT, "src", "vehicle", "spacecraft")

required_structure_files = (
    joinpath("src", "vehicle", "structure", "structure_models.jl"),
    joinpath("src", "vehicle", "structure", "assembly_graph.jl"),
    joinpath("src", "vehicle", "structure", "mass_properties.jl"),
    joinpath("src", "vehicle", "structure", "geometry_properties.jl"),
)

for rel in required_structure_files
    isfile(joinpath(REPO_ROOT, rel)) ||
        error("Vehicle structure boundary violation: missing canonical structure owner file: $rel")
end

isfile(joinpath(REPO_ROOT, "src", "vehicle", "spacecraft", "spacecraft_analysis.jl")) &&
    error("Vehicle structure boundary violation: retired file still exists at src/vehicle/spacecraft/spacecraft_analysis.jl")

const STRUCTURE_OWNER_FUNCTIONS = (
    "traverse_bodies",
    "get_COM",
    "update_inertia_tensor!",
    "update_inertia_tensor",
    "get_inertia_tensor",
    "set_inertia_tensor!",
    "get_spacecraft_mass",
    "get_spacecraft_reference_area",
    "get_spacecraft_length",
    "get_SA_area",
    "get_SC_area",
    "get_normal_vector",
    "get_tangent_vector",
)

violations = String[]
for (root, _, files) in walkdir(SPACECRAFT_ROOT)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        for fn in STRUCTURE_OWNER_FUNCTIONS
            occursin(Regex("function\\s+" * fn * "\\b"), src) || continue
            push!(violations, "$rel: declares structure-owner function '$fn' in spacecraft integration folder")
        end
    end
end

if !isempty(violations)
    println("vehicle_structure_boundary_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Vehicle structure boundary gate failed")
end

println("vehicle_structure_boundary_gate_ok")
