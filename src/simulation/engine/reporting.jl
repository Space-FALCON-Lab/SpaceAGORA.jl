function _warn_compat_normalize_flag!(args)
    if !args.simulation_settings.normalize || !_typed_normalize_warning_enabled()
        return nothing
    end
    if _normalize_warning_emitted[]
        return nothing
    end
    _normalize_warning_emitted[] = true
    @warn "SimulationSettings.normalize=true is legacy-only in typed run_simulation; propagation is always SI-native (m, s, kg). Set normalize=false to silence this warning."
    return nothing
end

function _enforce_typed_normalize_policy!(args)
    if !args.simulation_settings.normalize
        return nothing
    end
    if _typed_allow_compat_normalize()
        _warn_compat_normalize_flag!(args)
        return nothing
    end
    throw(ArgumentError(
        "SimulationSettings.normalize=true is unsupported in typed run_simulation. " *
        "Typed propagation is SI-native (m, s, kg). Set normalize=false, or set " *
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE=1 only for legacy transition checks."
    ))
end

function _debug_print_nan_parameter_paths!(x, path::AbstractString="p")
    if x isa Number
        if isnan(x)
            println("NaN found in parameter: $path")
        end
        return nothing
    end

    if x isa Base.RefValue{<:Number}
        xv = x[]
        if isnan(xv)
            println("NaN found in parameter: $path[]")
        end
        return nothing
    end

    if x isa AbstractArray{<:Number}
        for (idx, xv) in pairs(x)
            if isnan(xv)
                println("NaN found in parameter: $path[$idx]")
            end
        end
        return nothing
    end

    # Skip generic arrays of non-numeric types to keep debug scans bounded.
    if x isa AbstractArray
        return nothing
    end

    T = typeof(x)
    if isstructtype(T)
        for field in fieldnames(T)
            val = getfield(x, field)
            _debug_print_nan_parameter_paths!(val, string(path, ".", field))
        end
    end
    return nothing
end

