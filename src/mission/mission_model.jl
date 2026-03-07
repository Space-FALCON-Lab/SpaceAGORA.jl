"""
    mission_def(mission)

Typed-only mission-model entrypoint. Legacy Dict-based compatibility parsing has been removed.
"""
function mission_def(mission)
    if mission isa AbstractDict
        throw(ArgumentError("Dict-based mission model definitions are not supported. Use typed mission/config objects."))
    end
    return mission
end
