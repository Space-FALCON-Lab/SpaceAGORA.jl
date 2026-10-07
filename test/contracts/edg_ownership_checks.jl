module EDGOwnershipChecks
# These checks cover the typed route as well as the legacy boundary gate.
const OWNER = "src/gnc/guidance/aerobraking/typed_edg/"
const EXPECTED = Dict(
    "_edg_control_sat_state" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_control_pos_vel_mass" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_environment_state" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_in_drag_passage" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_ephemeris_time" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_planet_frame_lpi" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_targeting_prediction_environment" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_sample_prediction_atmosphere" => "src/gnc/guidance/aerobraking/typed_edg/services.jl",
    "_edg_heat_load_scale_height" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_total_ref_area" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_predict_mass" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_max_heat_load_for_links" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_weighted_aero_coefficients" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_heat_load_coefficients" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_eccentric_anomaly_from_true" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_mean_anomaly_from_true" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_drag_passage_duration" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_prediction_time_grid" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_closed_form_heat_load_trajectory" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_integrated_heat_load_trajectory" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_heat_load_lambdas" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_heat_load_alpha_profile" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_heat_load_track_env" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_constrained_heat_load_alpha_profile" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_profile_heat_rates" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_integrate_series" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_profile_heat_load" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_first_low_alpha_interval_indices" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_first_two_switch_alpha_profile" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_low_alpha_switch_window" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_balanced_tpbvp_heat_load_window" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_padded_heat_load_window" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_heat_load_profile_for_k" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_solve_heat_load_switches" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_edg_heat_load_low_drag_active" => "src/gnc/guidance/aerobraking/typed_edg/heat_load.jl",
    "_energy_depletion_heat_rate_calc" => "src/gnc/guidance/aerobraking/typed_edg/heat_rate.jl",
    "_energy_depletion_heatrate_root_alpha" => "src/gnc/guidance/aerobraking/typed_edg/heat_rate.jl",
    "_edg_maxwellian_heat_rate" => "src/gnc/guidance/aerobraking/typed_edg/heat_rate.jl",
    "_edg_heat_rate_alpha" => "src/gnc/guidance/aerobraking/typed_edg/heat_rate.jl",
    "_energy_depletion_struct_drag_area" => "src/gnc/guidance/aerobraking/typed_edg/structural_load.jl",
    "_energy_depletion_struct_load_root_alpha" => "src/gnc/guidance/aerobraking/typed_edg/structural_load.jl",
    "_edg_structural_alpha" => "src/gnc/guidance/aerobraking/typed_edg/structural_load.jl",
    "_edg_recompute_switches!" => "src/gnc/guidance/aerobraking/typed_edg/angle_decision.jl",
    "_edg_base_alpha" => "src/gnc/guidance/aerobraking/typed_edg/angle_decision.jl",
    "_edg_command_alpha!" => "src/gnc/guidance/aerobraking/typed_edg/angle_decision.jl",
    "_edg_orbit_metrics_from_rv" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_target_energy_from_apoapsis" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_targeting_constrained_alpha" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_targeting_prediction_time_grid" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_targeting_aero_acceleration" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_integrated_targeting_trajectory" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_integrated_max_energy_depletion_trajectory" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_predict_targeting_outcome" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_predict_max_energy_depletion_outcome" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_targeting_switch_outcomes" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_targeting_bracket_outcomes" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_targeting_outcome_with_heat_load" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_certify_targeting_candidates" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_disable_uncertified_targeting!" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_solve_targeting_switch" => "src/gnc/guidance/aerobraking/typed_edg/targeting.jl",
    "_edg_interpolate_bracket_value" => "src/gnc/guidance/aerobraking/typed_edg/guidance_decision.jl",
    "_edg_target_energy_from_reachable_bracket" => "src/gnc/guidance/aerobraking/typed_edg/guidance_decision.jl",
    "_edg_set_targeting_fallback!" => "src/gnc/guidance/aerobraking/typed_edg/guidance_decision.jl",
    "_edg_run_target_energy_bracketing!" => "src/gnc/guidance/aerobraking/typed_edg/guidance_decision.jl",
)

# Decision entry points have the same ownership rules as the retained kernels.
EXPECTED["control_decision!"] = OWNER * "angle_decision.jl"
EXPECTED["guidance_decision!"] = OWNER * "guidance_decision.jl"

const FORWARD_CALLS = Dict(
    "_edg_recompute_switches!" => :(EDGAlgorithms._edg_recompute_switches!(model.config, model.state, p, env, spacecraft, pos, vel, mass, heat_load_j_cm2, t, i)),
    "_edg_base_alpha" => :(EDGAlgorithms._edg_base_alpha(model.config, model.state, t, i)),
    "_edg_command_alpha!" => :(EDGAlgorithms._edg_command_alpha!(model.config, model.state, p, model.aoa_effector.controlled_panel_links, env, spacecraft, base_alpha, heat_load_j_cm2, heat_load_low_drag_active, i)),
    "_edg_heat_load_low_drag_active" => :(EDGAlgorithms._edg_heat_load_low_drag_active(model.config, model.state, t, i)),
    "_edg_interpolate_bracket_value" => :(_edg_algorithms()._edg_interpolate_bracket_value(exit_energy, energy_min, energy_max, value_at_min, value_at_max)),
    "_edg_target_energy_from_reachable_bracket" => :(_edg_algorithms()._edg_target_energy_from_reachable_bracket(planet, target_apoapsis_radius_m, energy_min, energy_max, periapsis_at_min, periapsis_at_max)),
    "_edg_set_targeting_fallback!" => :(_edg_algorithms()._edg_set_targeting_fallback!(config, state, i)),
    "_edg_run_target_energy_bracketing!" => :(_edg_algorithms()._edg_run_target_energy_bracketing!(model.config, model.state, u, p, t, i)),
)

const WRAPPER_FILES = Set((
    "src/gnc/control/targeting_control.jl",
    "src/gnc/control/heat_load_control.jl",
    "src/gnc/guidance/target_energy_bracketing.jl",
))

function source_map(root)
    result = Dict{String,String}()
    for (dir, _, files) in walkdir(joinpath(root, "src"))
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(dir, file)
            result[replace(relpath(path, root), '\\' => '/')] = read(path, String)
        end
    end
    return result
end


# Parse syntax without evaluating it. Parenthesized include arguments and
# indented, short-form, macro-wrapped or qualified methods remain visible.
function _walk(f, ex)
    ex isa Expr || return
    ex.head == :quote && return
    f(ex)
    foreach(arg -> _walk(f, arg), ex.args)
end

function _leaf_name(ex)
    ex isa Symbol && return String(ex)
    ex isa Expr && ex.head == :. && ex.args[end] isa QuoteNode &&
        return String(ex.args[end].value)
    return nothing
end

function _signature(ex)
    while ex isa Expr && ex.head in (:(::), :where)
        ex = ex.args[1]
    end
    return ex
end

function _definition(ex)
    ex.head in (:function, :(=)) || return nothing
    sig = _signature(ex.args[1])
    sig isa Expr && sig.head == :call || return nothing
    name = _leaf_name(sig.args[1])
    return name === nothing ? nothing : (name, sig.args[1], ex.args[2])
end

function _forwarder_only(body, name)
    haskey(FORWARD_CALLS, name) || return false
    statements = body isa Expr && body.head == :block ?
        filter(x -> !(x isa LineNumberNode), body.args) : [body]
    length(statements) in (1, 2) || return false
    call = first(statements)
    call isa Expr && call.head == :return && (call = only(call.args))
    call == FORWARD_CALLS[name] || return false
    return length(statements) == 1 ||
        (name == "_edg_recompute_switches!" && last(statements) == Expr(:return, :nothing))
end

function _include_calls(source)
    calls = Expr[]
    _walk(Meta.parseall(source)) do ex
        ex.head == :call && _leaf_name(ex.args[1]) == "include" && push!(calls, ex)
    end
    return calls
end

function _contains_control_literal(ex)
    ex isa String && return occursin("control", ex)
    ex isa Expr || return false
    return any(_contains_control_literal, ex.args)
end

has_control_include(source) =
    any(call -> any(_contains_control_literal, call.args[2:end]), _include_calls(source))

_syntax(ex) = ex isa Expr ?
    Expr(ex.head, (_syntax(arg) for arg in ex.args if !(arg isa LineNumberNode))...) : ex

function aggregator_violations(source)
    siblings = ("heat_rate.jl", "heat_load.jl", "structural_load.jl",
                "targeting.jl", "guidance_decision.jl", "angle_decision.jl")
    expected = [_syntax(Meta.parse("include(joinpath(@__DIR__, $(repr(file))))"))
                for file in siblings]
    actual = map(_syntax, _include_calls(source))
    return actual == expected ? String[] :
        ["$(OWNER)algorithms.jl: expected exactly the six ordered sibling includes"]
end

# Runtime checks complement the source inventory: aliases share the same
# function object, while a qualified extension carries its defining module.
function runtime_violations(algorithms, services, gnc_modules)
    errors = String[]
    for (name, path) in EXPECTED
        symbol = Symbol(name)
        owner = endswith(path, "/services.jl") ? services : algorithms
        isdefined(owner, symbol) || (push!(errors, "$owner: missing $name"); continue)
        owned = getfield(owner, symbol)
        all(method -> method.module === owner, methods(owned)) ||
            push!(errors, "$owner: foreign method on $name")
        for mod in gnc_modules
            isdefined(mod, symbol) || continue
            bound = getfield(mod, symbol)
            bound === owned && continue
            allowed = haskey(FORWARD_CALLS, name) &&
                nameof(mod) == (name in ("_edg_interpolate_bracket_value",
                    "_edg_target_energy_from_reachable_bracket", "_edg_set_targeting_fallback!",
                    "_edg_run_target_energy_bracketing!") ? :GuidanceHooks : :ControlHooks)
            allowed && all(method -> method.module === mod, methods(bound)) && continue
            push!(errors, "$mod: separate EDG binding $name")
        end
    end
    return sort!(errors)
end

function violations(sources)
    errors = String[]
    counts = Dict(name => 0 for name in keys(EXPECTED))
    for (path, source) in sources
        # Strip line comments for dependency checks; function bodies remain parsed
        # as complete top-level blocks so multiline signatures are covered.
        code = join((first(split(line, '#'; limit=2)) for line in split(source, '\n')), '\n')
        if startswith(path, OWNER)
            for forbidden in (r"\bControlHooks\b", r"\b_control_module\s*\(",
                              r"\b(?:rotate_link|_apply_solar_panel_aoa!)\s*\(",
                              r"\b(?:calcControlEffect!|calcGuidanceEffect!)\s*\(",
                              r"\.\s*(?:α|q|heat_loads)\s*(?:\.=|=(?!=))",
                              r"\b(?:AerobrakingEnergyDepletionControlModel|SolarPanelAngleOfAttackControlModel)\b")
                occursin(forbidden, code) && push!(errors, "$path: forbidden reverse dependency or mutation: $forbidden")
            end
        end
        if path == "src/gnc/guidance/target_energy_bracketing.jl"
            occursin(r"\b_control_module\s*\(|\bControlHooks\b", code) &&
                push!(errors, "$path: typed guidance depends on control")
        end
        _walk(Meta.parseall(source)) do ex
            definition = _definition(ex)
            definition === nothing && return
            name, callee, body = definition
            haskey(EXPECTED, name) || return
            if path == EXPECTED[name] && callee isa Symbol
                counts[name] += 1
            elseif path in WRAPPER_FILES && callee isa Symbol
                _forwarder_only(body, name) ||
                    push!(errors, "$path: $name must be an exact forwarding compatibility method")
            else
                push!(errors, "$path: duplicate EDG definition $name")
            end
        end
    end
    aggregator = OWNER * "algorithms.jl"
    haskey(sources, aggregator) ?
        append!(errors, aggregator_violations(sources[aggregator])) :
        push!(errors, "$aggregator: missing aggregator")
    for (name, count) in counts
        count == 1 || push!(errors, "$(EXPECTED[name]): expected one $name definition, found $count")
    end
    return sort!(errors)
end
end # module EDGOwnershipChecks
