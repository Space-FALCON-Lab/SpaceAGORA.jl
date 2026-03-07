if !isdefined(@__MODULE__, :CONTROL_BRIDGE_STATE_LOCK)
    const CONTROL_BRIDGE_STATE_LOCK = ReentrantLock()
end

if !isdefined(@__MODULE__, :_bridge_get_cnf)
    @inline function _bridge_get_cnf(args=nothing; cnf=nothing)
        if cnf !== nothing
            return cnf
        end
        if args isa NamedTuple && hasproperty(args, :cnf)
            return getproperty(args, :cnf)
        end
        if !(args isa NamedTuple) && hasproperty(args, :cnf)
            return getproperty(args, :cnf)
        end
        if (@isdefined config) && isdefined(config, :cnf)
            return getproperty(config, :cnf)
        end
        throw(ArgumentError("Control state `cnf` not found. Pass `cnf=` or typed args.cnf field."))
    end
end

if !isdefined(@__MODULE__, :_bridge_get_solution)
    @inline function _bridge_get_solution(args=nothing; solution=nothing, cnf=nothing)
        if solution !== nothing
            return solution
        end
        if args isa NamedTuple && hasproperty(args, :solution)
            return getproperty(args, :solution)
        end
        if !(args isa NamedTuple) && hasproperty(args, :solution)
            return getproperty(args, :solution)
        end
        if cnf !== nothing && hasproperty(cnf, :solution)
            return getproperty(cnf, :solution)
        end
        if (@isdefined config) && isdefined(config, :solution)
            return getproperty(config, :solution)
        end
        throw(ArgumentError("Solution state `solution` not found. Pass `solution=` or typed args.solution field."))
    end
end
