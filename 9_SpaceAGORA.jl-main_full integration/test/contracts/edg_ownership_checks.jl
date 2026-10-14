module EDGOwnershipChecks
# These checks cover the typed route as well as the legacy boundary gate.
const OWNER = "src/models/gnc/guidance/aerobraking/typed_edg/"
const KERNELS = (
    "_edg_at_control_tick",
    "_edg_balanced_tpbvp_heat_load_window",
    "_edg_base_alpha",
    "_edg_candidate_constraint_sample",
    "_edg_candidate_energy_mapping",
    "_edg_certify_targeting_candidates",
    "_edg_chebyshev4",
    "_edg_chebyshev44",
    "_edg_chebyshev_t0_t4",
    "_edg_closed_form_correction_fit",
    "_edg_closed_form_density",
    "_edg_closed_form_density_fit",
    "_edg_closed_form_f1_f2",
    "_edg_closed_form_heat_load_trajectory",
    "_edg_command_alpha!",
    "_edg_constrained_heat_load_alpha_profile",
    "_edg_constraint_min_alpha",
    "_edg_constraint_status",
    "_edg_controlled_panel_area",
    "_edg_controlled_ref_area",
    "_edg_corrected_target_energy_from_apoapsis",
    "_edg_disable_uncertified_targeting!",
    "_edg_drag_passage_duration",
    "_edg_dynamic_spacecraft_aero_coefficients",
    "_edg_eccentric_anomaly_from_true",
    "_edg_exact_targeting_exit",
    "_edg_first_low_alpha_interval_indices",
    "_edg_first_two_switch_alpha_profile",
    "_edg_heat_budget_status",
    "_edg_heat_load_alpha_profile",
    "_edg_heat_load_coefficients",
    "_edg_heat_load_costates",
    "_edg_heat_load_lambdas",
    "_edg_heat_load_low_drag_active",
    "_edg_heat_load_profile_for_k",
    "_edg_heat_load_scale_height",
    "_edg_heat_load_security_required",
    "_edg_heat_load_shooting_profile",
    "_edg_heat_load_shooting_stage",
    "_edg_heat_load_track_env",
    "_edg_heat_rate_alpha",
    "_edg_integrate_heat_load_shooting",
    "_edg_integrate_series",
    "_edg_integrated_heat_load_trajectory",
    "_edg_integrated_max_energy_depletion_trajectory",
    "_edg_integrated_targeting_trajectory",
    "_edg_interpolate_bracket_value",
    "_edg_journal_plate_coefficients",
    "_edg_journal_spacecraft_aero_coefficients",
    "_edg_journal_targeting_plan",
    "_edg_low_alpha_switch_window",
    "_edg_max_heat_load_for_links",
    "_edg_maxwellian_heat_rate",
    "_edg_mean_anomaly_from_true",
    "_edg_orbit_metrics_from_rv",
    "_edg_ordered_feasible_prefix",
    "_edg_outbound_exit_index",
    "_edg_padded_heat_load_window",
    "_edg_pass_heat_load_for_links",
    "_edg_predict_mass",
    "_edg_predict_max_energy_depletion_outcome",
    "_edg_predict_targeting_outcome",
    "_edg_prediction_advance",
    "_edg_prediction_event",
    "_edg_prediction_gravity_acceleration",
    "_edg_prediction_step",
    "_edg_prediction_time_grid",
    "_edg_profile_heat_load",
    "_edg_profile_heat_rates",
    "_edg_recompute_second_heat_load_switch",
    "_edg_recompute_switches!",
    "_edg_remaining_closed_form_heat_load_trajectory",
    "_edg_remaining_closed_form_reference",
    "_edg_replan_due",
    "_edg_run_target_energy_bracketing!",
    "_edg_select_mapped_candidate",
    "_edg_set_targeting_fallback!",
    "_edg_solve_heat_load_switches",
    "_edg_solve_targeting_switch",
    "_edg_structural_alpha",
    "_edg_switch_due",
    "_edg_target_energy_from_apoapsis",
    "_edg_target_energy_from_reachable_bracket",
    "_edg_targeting_aero_acceleration",
    "_edg_targeting_bracket_outcomes",
    "_edg_targeting_candidate_certificate",
    "_edg_targeting_constrained_alpha",
    "_edg_targeting_control_period",
    "_edg_targeting_outcome_with_heat_load",
    "_edg_targeting_prediction_time_grid",
    "_edg_targeting_switch_outcomes",
    "_edg_targeting_switch_solution",
    "_edg_total_ref_area",
    "_edg_update_heat_load_security!",
    "_edg_vacuum_apoapsis_correction",
    "_edg_vacuum_drag_passage_exit",
    "_edg_weighted_aero_coefficients",
    "_energy_depletion_heat_rate_calc",
    "_energy_depletion_heatrate_root_alpha",
    "_energy_depletion_struct_drag_area",
    "_energy_depletion_struct_load_root_alpha",
    "control_decision!",
    "guidance_decision!",
)
const EXPECTED = Dict(name => OWNER * "algorithms.jl" for name in KERNELS)
EXPECTED["_edg_control_sat_state"] = OWNER * "services.jl"
EXPECTED["_edg_control_pos_vel_mass"] = OWNER * "services.jl"
EXPECTED["_edg_environment_state"] = OWNER * "services.jl"
EXPECTED["_edg_in_drag_passage"] = OWNER * "services.jl"
EXPECTED["_edg_ephemeris_time"] = OWNER * "services.jl"
EXPECTED["_edg_planet_frame_lpi"] = OWNER * "services.jl"
EXPECTED["_edg_targeting_prediction_environment"] = OWNER * "services.jl"
EXPECTED["_edg_sample_prediction_atmosphere"] = OWNER * "services.jl"
const RETIRED_KERNEL_FILES = ("heat_rate.jl", "heat_load.jl", "structural_load.jl",
    "targeting.jl", "targeting_plan.jl", "pass_heat.jl", "guidance_decision.jl", "angle_decision.jl")
const SERVICE_OVERLOADS = Set(("_edg_environment_state", "_edg_in_drag_passage", "_edg_ephemeris_time",
    "_edg_planet_frame_lpi", "_edg_targeting_prediction_environment", "_edg_sample_prediction_atmosphere"))
const FACADE_FILE = "src/models/gnc/guidance/adapters/edg/services.jl"
const FACADE_FORWARD = Dict(
    "_edg_environment_state" => :(_environment(u, p, t, i)),
    "_edg_in_drag_passage" => :(_in_passage(p, env)),
    "_edg_ephemeris_time" => :(_epoch(p, t)),
    "_edg_planet_frame_lpi" => :(_frame(p, t)),
    "_edg_targeting_prediction_environment" => :(_trajectory(p, r, v, t)),
    "_edg_sample_prediction_atmosphere" => :(_atmosphere(p, h, t)),
)

const FORWARD_CALLS = Dict(
    "_edg_recompute_switches!" => :(EDGAlgorithms._edg_recompute_switches!(model.config, model.state, p, env, spacecraft, pos, vel, mass, heat_load_j_cm2, t, i)),
    "_edg_base_alpha" => :(EDGAlgorithms._edg_base_alpha(model.config, model.state, t, i)),
    "_edg_command_alpha!" => :(EDGAlgorithms._edg_command_alpha!(model.config, model.state, p, model.aoa_effector.controlled_panel_links, env, spacecraft, base_alpha, heat_load_j_cm2, heat_load_low_drag_active, i; pass_heat_load_j_cm2=pass_heat_load_j_cm2)),
    "_edg_heat_load_low_drag_active" => :(EDGAlgorithms._edg_heat_load_low_drag_active(model.config, model.state, t, i)),
    "_edg_interpolate_bracket_value" => :(_edg_algorithms()._edg_interpolate_bracket_value(exit_energy, energy_min, energy_max, value_at_min, value_at_max)),
    "_edg_target_energy_from_reachable_bracket" => :(_edg_algorithms()._edg_target_energy_from_reachable_bracket(planet, target_apoapsis_radius_m, energy_min, energy_max, periapsis_at_min, periapsis_at_max)),
    "_edg_set_targeting_fallback!" => :(_edg_algorithms()._edg_set_targeting_fallback!(config, state, i)),
    "_edg_run_target_energy_bracketing!" => :(_edg_algorithms()._edg_run_target_energy_bracketing!(model.config, model.state, u, p, t, i)),
)

const WRAPPER_FILES = Set((
    "src/models/gnc/control/targeting_control.jl",
    "src/models/gnc/control/heat_load_control.jl",
    "src/models/gnc/guidance/aerobraking/typed_edg/guidance_models.jl",
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
    errors = String[]
    isempty(_include_calls(source)) || push!(errors, "$(OWNER)algorithms.jl: numerical source includes are forbidden")
    inventory = Any[]
    definitions = Dict{String,Int}()
    _walk(Meta.parseall(source)) do ex
        if ex.head == :const && ex.args[1] isa Expr && ex.args[1].head == :(=) &&
                ex.args[1].args[1] == :COMPATIBILITY_FUNCTIONS
            push!(inventory, _syntax(ex.args[1].args[2]))
        end
        definition = _definition(ex)
        definition === nothing && return
        name, _, _ = definition
        expected = Meta.parse("$name(args...; kwargs...) = _invoke_kernel(:$name, args...; kwargs...)")
        name in KERNELS && _syntax(ex) == _syntax(expected) ||
            push!(errors, "$(OWNER)algorithms.jl: $name must be an exact package forwarder")
        definitions[name] = get(definitions, name, 0) + 1
    end
    expected_inventory = Expr(:tuple, (QuoteNode(Symbol(name)) for name in KERNELS)...)
    inventory == [expected_inventory] || push!(errors, "$(OWNER)algorithms.jl: compatibility inventory must match all declared kernels")
    all(get(definitions, name, 0) == 1 for name in KERNELS) ||
        push!(errors, "$(OWNER)algorithms.jl: expected exactly one forwarder for every kernel")
    return errors
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
        # This relocated callback adapter was outside the numerical owner
        # before the move and retains its separate no-control rule below.
        if startswith(path, OWNER) && path != OWNER * "guidance_models.jl"
            for forbidden in (r"\bControlHooks\b", r"\b_control_module\s*\(",
                              r"\b(?:rotate_link|_apply_solar_panel_aoa!)\s*\(",
                              r"\b(?:calcControlEffect!|calcGuidanceEffect!)\s*\(",
                              r"\.\s*(?:α|q|heat_loads)\s*(?:\.=|=(?!=))",
                              r"\b(?:AerobrakingEnergyDepletionControlModel|SolarPanelAngleOfAttackControlModel)\b")
                occursin(forbidden, code) && push!(errors, "$path: forbidden reverse dependency or mutation: $forbidden")
            end
        end
        if path == "src/models/gnc/guidance/aerobraking/typed_edg/guidance_models.jl"
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
            elseif path == OWNER * "adapter_queries.jl" && name in SERVICE_OVERLOADS && callee isa Symbol
                counts[name] += 1
            elseif path == FACADE_FILE && haskey(FACADE_FORWARD, name) && callee isa Symbol
                _syntax(body isa Expr && body.head == :block ? only(filter(x -> !(x isa LineNumberNode), body.args)) : body) == FACADE_FORWARD[name] || push!(errors, "$path: $name must forward its capability query unchanged")
            elseif path in WRAPPER_FILES && callee isa Symbol
                _forwarder_only(body, name) ||
                    push!(errors, "$path: $name must be an exact forwarding compatibility method")
            else
                push!(errors, "$path: duplicate EDG definition $name")
            end
        end
    end
    services = OWNER * "services.jl"
    if haskey(sources, services)
        _syntax.(_include_calls(sources[services])) == [_syntax(:(include(joinpath(@__DIR__, "adapter_queries.jl"))))] ||
            push!(errors, "$services: only the reviewed adapter_queries.jl include is allowed")
    end
    aggregator = OWNER * "algorithms.jl"
    haskey(sources, aggregator) ?
        append!(errors, aggregator_violations(sources[aggregator])) :
        push!(errors, "$aggregator: missing aggregator")
    for (name, count) in counts
        expected_count = name in SERVICE_OVERLOADS ? 2 : 1
        count == expected_count || push!(errors, "$(EXPECTED[name]): expected $expected_count $name definitions, found $count")
    end
    for file in RETIRED_KERNEL_FILES
        haskey(sources, OWNER * file) && push!(errors, "$(OWNER)$file: numerical implementation belongs to SpaceAGORAEDG")
    end
    return sort!(errors)
end
end # module EDGOwnershipChecks
