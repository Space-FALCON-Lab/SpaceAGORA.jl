if !isdefined(@__MODULE__, :LEGACY_CONTROL_STATE_LOCK)
    const LEGACY_CONTROL_STATE_LOCK = ReentrantLock()
end

if !isdefined(@__MODULE__, :_legacy_get_cnf)
    @inline function _legacy_get_cnf(args=nothing; cnf=nothing)
        if cnf !== nothing
            return cnf
        end
        if args isa AbstractDict && haskey(args, :cnf)
            return args[:cnf]
        end
        if args isa NamedTuple && hasproperty(args, :cnf)
            return getproperty(args, :cnf)
        end
        if (@isdefined config) && isdefined(config, :cnf)
            return getproperty(config, :cnf)
        end
        throw(ArgumentError("Legacy control state `cnf` not found. Pass `cnf=` or args[:cnf]."))
    end
end

if !isdefined(@__MODULE__, :_legacy_get_solution)
    @inline function _legacy_get_solution(args=nothing; solution=nothing, cnf=nothing)
        if solution !== nothing
            return solution
        end
        if args isa AbstractDict && haskey(args, :solution)
            return args[:solution]
        end
        if args isa NamedTuple && hasproperty(args, :solution)
            return getproperty(args, :solution)
        end
        if cnf !== nothing && hasproperty(cnf, :solution)
            return getproperty(cnf, :solution)
        end
        if (@isdefined config) && isdefined(config, :solution)
            return getproperty(config, :solution)
        end
        throw(ArgumentError("Legacy solution state `solution` not found. Pass `solution=` or args[:solution]."))
    end
end
