@inline _mission_token(x) = lowercase(strip(replace(string(x), r"[_-]+" => " ")))

@inline function _mission_kind(x)::Symbol
    key = _mission_token(x)
    if key in ("drag passage",)
        return :drag_passage
    elseif key in ("entry",)
        return :entry
    elseif key in ("orbits", "orbit", "missionorbits")
        return :orbits
    elseif key in ("time", "missiontime")
        return :time
    elseif key in ("aerobraking campaign",)
        return :campaign
    end
    return :unknown
end

"""
    def_miss(args)

Typed-only mission-definition normalization entrypoint.
Dict-based mission configuration is no longer supported.
"""
function def_miss(args)
    if args isa AbstractDict
        throw(ArgumentError("Dict-based mission configuration is not supported. Use typed mission/config objects."))
    end
    return args
end
