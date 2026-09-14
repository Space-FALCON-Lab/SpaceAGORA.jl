module IOLogging

@inline function io_log_enabled(args)::Bool
    if hasproperty(args, :simulation_settings) && hasproperty(args.simulation_settings, :verbose)
        return Bool(getproperty(args.simulation_settings, :verbose))
    end
    return false
end

@inline function io_log_info(args, msg::AbstractString)
    if io_log_enabled(args)
        @info msg
    end
    return nothing
end

@inline function io_log_warn(msg::AbstractString)
    @warn msg
    return nothing
end

export io_log_enabled
export io_log_info
export io_log_warn

end # module IOLogging
