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

const WRAPPER_FILES = Set((
    "src/gnc/control/targeting_control.jl",
    "src/gnc/control/heat_load_control.jl",
    "src/gnc/guidance/target_energy_bracketing.jl",
))

function source_map(root)
    result = Dict{String,String}()
    for owner in ("guidance", "control", "shared", "internal"), (dir, _, files) in walkdir(joinpath(root, "src", "gnc", owner))
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(dir, file)
            result[replace(relpath(path, root), '\\' => '/')] = read(path, String)
        end
    end
    return result
end


function _projection(ex)
    ex isa Symbol && return true
    ex isa Expr || return false
    return ex.head == :. && _projection(ex.args[1]) && ex.args[2] isa QuoteNode
end

function _forward_call(ex, name)
    ex isa Expr && ex.head == :return && (ex = only(ex.args))
    ex isa Expr && ex.head == :call || return false
    callee = ex.args[1]
    callee isa Expr && callee.head == :. || return false
    callee.args[2] == QuoteNode(Symbol(name)) || return false
    owner = callee.args[1]
    allowed_owner = owner == :EDGAlgorithms ||
        (owner isa Expr && owner.head == :call && owner.args == [:_edg_algorithms])
    return allowed_owner && all(_projection, ex.args[2:end])
end

function _forwarder_only(source, name)
    parsed = Meta.parse(source)
    parsed.head == :macrocall && (parsed = parsed.args[end])
    parsed.head == :function || return false
    body = filter(x -> !(x isa LineNumberNode), parsed.args[2].args)
    length(body) in (1,2) || return false
    _forward_call(first(body), name) || return false
    return length(body) == 1 || last(body) == Expr(:return, :nothing)
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
        for m in eachmatch(r"(?m)^(?:@inline )?function ([^\s(]+)[\s\S]*?^end\b", source)
            name = m.captures[1]
            haskey(EXPECTED, name) || continue
            if path == EXPECTED[name]
                counts[name] += 1
            elseif path in WRAPPER_FILES
                if !_forwarder_only(m.match, name)
                    push!(errors, "$path: $name must be a forwarding compatibility method")
                end
            else
                push!(errors, "$path: duplicate EDG definition $name")
            end
        end
    end
    for (name, count) in counts
        count == 1 || push!(errors, "$(EXPECTED[name]): expected one $name definition, found $count")
    end
    return sort!(errors)
end
end # module EDGOwnershipChecks
