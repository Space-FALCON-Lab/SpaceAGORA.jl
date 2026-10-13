using TOML

"""Select one existing scenario without changing any scientific setting."""
function select_telemetry_scenario(document::AbstractDict, name::AbstractString)
    scenarios = document["scenarios"]
    names = [entry["name"] for entry in scenarios]
    length(unique(names)) == length(names) || error("Duplicate telemetry scenario names")
    selected = filter(entry -> entry["name"] == name, scenarios)
    length(selected) == 1 || error("Unknown telemetry scenario: $name")
    return merge(document, Dict("scenarios" => selected))
end

function telemetry_scenario_main(args)
    length(args) == 3 || error("Expected input manifest, scenario name and output manifest")
    input, name, output = args
    abspath(input) != abspath(output) || error("Refuse to overwrite the source manifest")
    document = select_telemetry_scenario(TOML.parsefile(input), name)
    mkpath(dirname(output))
    open(io -> TOML.print(io, document), output, "w")
    TOML.parsefile(output) == document || error("Scenario manifest roundtrip changed values")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    telemetry_scenario_main(ARGS)
end
