
if !isdefined(@__MODULE__, :_legacy_get_cnf)
    @inline function _legacy_get_cnf(args=nothing; cnf=nothing)
        if cnf !== nothing
            return cnf
        end
        if args isa AbstractDict && haskey(args, :cnf)
            return args[:cnf]
        end
        if (@isdefined config) && isdefined(config, :cnf)
            return getproperty(config, :cnf)
        end
        throw(ArgumentError("Legacy control state `cnf` not found. Pass `cnf=` or args[:cnf]."))
    end
end

function no_maneuver(t0, thrust_mag, Δv, args, index_phase_aerobraking)
    thrust = 0
    return thrust
end

function abms(t0, thrust_mag, Δv, args, index_phase_aerobraking)
    if index_phase_aerobraking == 0
        thrust = thrust_mag
    else
        thrust = 0
    end

    return thrust
end

function deceleration_drag_passage(t0, thrust_mag, Δv, args, index_phase_aerobraking; cnf=nothing)
    cnf_state = _legacy_get_cnf(args; cnf=cnf)
    if cnf_state.drag_state == true
        thrust = thrust_mag
    else
        thrust = 0
    end

    return thrust
end
