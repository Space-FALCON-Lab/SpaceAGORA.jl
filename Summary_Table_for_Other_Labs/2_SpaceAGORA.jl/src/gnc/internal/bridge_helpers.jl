if !isdefined(@__MODULE__, :CONTROL_BRIDGE_STATE_LOCK)
    const CONTROL_BRIDGE_STATE_LOCK = ReentrantLock()
end

if !isdefined(@__MODULE__, :_bridge_verbose_enabled)
    @inline function _bridge_verbose_enabled(args=nothing)::Bool
        if get(ENV, "SPACEAGORA_DEBUG_LEGACY_CONTROL", "0") == "1"
            return true
        end
        if args !== nothing && hasproperty(args, :simulation_settings) && hasproperty(args.simulation_settings, :verbose)
            return Bool(getproperty(args.simulation_settings, :verbose))
        end
        if args !== nothing && hasproperty(args, :verbose)
            return Bool(getproperty(args, :verbose))
        end
        return false
    end
end

if !isdefined(@__MODULE__, :_bridge_required_field)
    @inline function _bridge_required_field(args, name::Symbol)
        if args !== nothing && hasproperty(args, name)
            return getproperty(args, name)
        end
        if args !== nothing && applicable(getindex, args, name)
            return getindex(args, name)
        end
        throw(ArgumentError("Required runtime field `$(name)` not found."))
    end
end

if !isdefined(@__MODULE__, :_bridge_optional_field)
    @inline function _bridge_optional_field(args, name::Symbol, default)
        if args !== nothing && hasproperty(args, name)
            return getproperty(args, name)
        end
        if args !== nothing && applicable(get, args, name, default)
            return get(args, name, default)
        end
        if args !== nothing && applicable(getindex, args, name)
            return getindex(args, name)
        end
        return default
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_topography_enabled)
    @inline function _bridge_aerobraking_topography_enabled(args)::Bool
        if args !== nothing && hasproperty(args, :environment_model) && hasproperty(args.environment_model, :topography)
            return Bool(getproperty(args.environment_model, :topography))
        end
        return _bridge_optional_field(args, :topography_model, "") == "Spherical Harmonics"
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_entry_interface_m)
    @inline function _bridge_aerobraking_entry_interface_m(args, mission=nothing)::Float64
        if args !== nothing && hasproperty(args, :environment_model) && hasproperty(args.environment_model, :EI)
            return Float64(getproperty(args.environment_model, :EI)) * 1e3
        end
        return Float64(_bridge_required_field(args, :EI)) * 1e3
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_exit_interface_m)
    @inline function _bridge_aerobraking_exit_interface_m(args, mission=nothing)::Float64
        default_m = _bridge_aerobraking_entry_interface_m(args, mission)
        return Float64(_bridge_optional_field(args, :AE, default_m / 1e3)) * 1e3
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_body_shape)
    @inline function _bridge_aerobraking_body_shape(args, mission=nothing)::String
        value = _bridge_optional_field(args, :body_shape, "Spacecraft")
        return String(value)
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_heat_load_solution)
    @inline _bridge_aerobraking_heat_load_solution(args)::Int = Int(_bridge_optional_field(args, :heat_load_sol, 0))
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_max_heat_rate)
    @inline function _bridge_aerobraking_max_heat_rate(args, mission=nothing)::Float64
        default_limit = if mission !== nothing && hasproperty(mission, :aerodynamics) && hasproperty(mission.aerodynamics, :heat_rate_limit)
            Float64(getproperty(mission.aerodynamics, :heat_rate_limit))
        else
            Inf
        end
        return Float64(_bridge_optional_field(args, :max_heat_rate, default_limit))
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_srp_enabled)
    @inline _bridge_aerobraking_srp_enabled(args)::Bool = Bool(_bridge_optional_field(args, :srp, false))
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_control_mode)
    @inline _bridge_aerobraking_control_mode(args)::Int = Int(_bridge_optional_field(args, :control_mode, 0))
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_struct_control_enabled)
    @inline _bridge_aerobraking_struct_control_enabled(args)::Bool = Bool(_bridge_optional_field(args, :struct_ctrl, false))
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_dry_mass)
    @inline function _bridge_aerobraking_dry_mass(args, mission=nothing)::Float64
        if args !== nothing && hasproperty(args, :dynamics_model) && hasproperty(args.dynamics_model, :spacecraft)
            spacecraft = getproperty(args.dynamics_model, :spacecraft)
            isempty(spacecraft) || return Float64(getproperty(first(spacecraft), :dry_mass))
        end
        if mission !== nothing && hasproperty(mission, :body) && hasproperty(mission.body, :dry_mass)
            return Float64(getproperty(mission.body, :dry_mass))
        end
        return Float64(_bridge_required_field(args, :dry_mass))
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_thrust_phi)
    @inline function _bridge_aerobraking_thrust_phi(args, mission=nothing)::Float64
        default_phi = if mission !== nothing && hasproperty(mission, :engines) && hasproperty(mission.engines, :ϕ)
            Float64(getproperty(mission.engines, :ϕ))
        else
            0.0
        end
        return Float64(_bridge_optional_field(args, :phi, default_phi))
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_control_in_loop)
    @inline _bridge_aerobraking_control_in_loop(args)::Bool = Bool(_bridge_optional_field(args, :control_in_loop, false))
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_integrator_name)
    @inline function _bridge_aerobraking_integrator_name(args)::String
        return String(_bridge_optional_field(args, :integrator, "Julia"))
    end
end

if !isdefined(@__MODULE__, :_bridge_aerobraking_drag_passage)
    @inline _bridge_aerobraking_drag_passage(args)::Bool = Bool(_bridge_optional_field(args, :drag_passage, false))
end

if !isdefined(@__MODULE__, :_make_aerobraking_runtime_settings)
    @inline function _make_aerobraking_runtime_settings(args, mission=nothing)
        return (
            topography_enabled=_bridge_aerobraking_topography_enabled(args),
            entry_interface_m=_bridge_aerobraking_entry_interface_m(args, mission),
            exit_interface_m=_bridge_aerobraking_exit_interface_m(args, mission),
            body_shape=_bridge_aerobraking_body_shape(args, mission),
            heat_load_solution=_bridge_aerobraking_heat_load_solution(args),
            max_heat_rate=_bridge_aerobraking_max_heat_rate(args, mission),
            srp_enabled=_bridge_aerobraking_srp_enabled(args),
            control_mode=_bridge_aerobraking_control_mode(args),
            struct_control_enabled=_bridge_aerobraking_struct_control_enabled(args),
            dry_mass=_bridge_aerobraking_dry_mass(args, mission),
            thrust_phi=_bridge_aerobraking_thrust_phi(args, mission),
            control_in_loop=_bridge_aerobraking_control_in_loop(args),
            integrator_name=_bridge_aerobraking_integrator_name(args),
            drag_passage=_bridge_aerobraking_drag_passage(args),
        )
    end
end

if !isdefined(@__MODULE__, :_make_aerobraking_runtime_context)
    @inline function _make_aerobraking_runtime_context(;
        mission,
        index_phase_aerobraking,
        ip,
        aerobraking_phase,
        t_prev,
        date_initial,
        time_0,
        args,
        initial_state,
        gram_atmosphere,
        gram,
        cnf=nothing,
        solution=nothing,
        )
        return (
            mission=mission,
            index_phase_aerobraking=index_phase_aerobraking,
            ip=ip,
            aerobraking_phase=aerobraking_phase,
            t_prev=t_prev,
            date_initial=date_initial,
            time_0=time_0,
            args=args,
            initial_state=initial_state,
            gram_atmosphere=gram_atmosphere,
            gram=gram,
            settings=_make_aerobraking_runtime_settings(args, mission),
            cnf=cnf,
            solution=solution,
        )
    end
end

if !isdefined(@__MODULE__, :_with_control_gain)
    @inline _with_control_gain(context::NamedTuple, control_gain) = (; context..., control_gain)
end

if !isdefined(@__MODULE__, :_with_time_switch)
    @inline _with_time_switch(context::NamedTuple, time_switch) = (; context..., time_switch)
end

if !isdefined(@__MODULE__, :_bridge_get_cnf)
    @inline function _bridge_get_cnf(args=nothing; cnf=nothing)
        if cnf !== nothing
            return cnf
        end
        if args !== nothing && hasproperty(args, :cnf)
            return getproperty(args, :cnf)
        end
        throw(ArgumentError("Control state `cnf` not found. Pass `cnf=` or typed args.cnf field."))
    end
end

if !isdefined(@__MODULE__, :_bridge_get_solution)
    @inline function _bridge_get_solution(args=nothing; solution=nothing, cnf=nothing)
        if solution !== nothing
            return solution
        end
        if args !== nothing && hasproperty(args, :solution)
            return getproperty(args, :solution)
        end
        if cnf !== nothing && hasproperty(cnf, :solution)
            return getproperty(cnf, :solution)
        end
        throw(ArgumentError("Solution state `solution` not found. Pass `solution=` or typed args.solution field."))
    end
end
