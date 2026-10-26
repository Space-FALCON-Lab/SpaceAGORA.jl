# Accepted-step scheduling and history recording belong to simulation orchestration.
using ..InterLinkModels: InterLinkHistory, InterLinkKey, LinkContext, available_links, records_available
using ..SchedulingPolicies: schedule_interlinks!
using ..InterlinkTerminals: terminal_availability

export interlink_scheduler_callback

# Update the matching and its recorded history at initialization and accepted steps.
# The links score by their own `link_force`; `candidate_force` overrides it for
# callers that pass one (`candidate_force(target, partner, link)`).
function interlink_scheduler_callback(; candidate_force=nothing, initialization_active=nothing)
	function affect!(integrator; previous_active=nothing)
		args = integrator.p.args
		model = args.interlink_model
		planet = args.environment_model.planet
		# Terminal hardware the power budget sheds closes its links; the battery's
		# step-end update runs before this callback.
		selected = schedule_interlinks!(model, args.scheduling_policy_model,
			args.dynamics_model.spacecraft, integrator.u, planet.μ;
			candidate_force, previous_active, is_active=integrator.p.is_active,
			context=LinkContext(time_s=Float64(integrator.t), planet=planet),
			terminal_available=terminal_availability(model, args.dynamics_model.spacecraft))
		sample = InterLinkHistory(Float64(integrator.t), selected,
			records_available(model) ? available_links(model) : InterLinkKey[],
			_selected_scores(model, args.scheduling_policy_model, selected)...)
		if !isempty(model.history) && model.history[end].time == sample.time
			model.history[end] = sample
		else
			push!(model.history, sample)
		end
		DiffEqBase.derivative_discontinuity!(integrator, true)
		return nothing
	end
	initialize = (callback, state, time, integrator) -> affect!(integrator;
		previous_active=initialization_active === nothing ? nothing : initialization_active[]) # update link graph when the simulation is initialized
	return DiffEqBase.DiscreteCallback((state, time, integrator) -> true, affect!; # update link graph at each accepted simulation step
		initialize=initialize, save_positions=(false, false))
end

# The score that selected each link, read from the scheduler's workspace after the
# selection, and, when it keeps them, the weighted policy's term scores: one column
# per selected link.
function _selected_scores(model, scheduler, selected)
	work = scheduler.matching
	edges = [model.graph.edge_index_lookup[key] for key in selected]
	scores = Float64[work.score[edge] for edge in edges]
	terms = work.terms === nothing ? zeros(0, length(selected)) : work.terms[:, edges]
	return scores, terms
end

# The same scheduling, at initialization only: placed ahead of the run's other
# callbacks, the first selection is in place before they initialize, so the
# first saved row already shows what the links draw. The scheduler callback,
# after the others, repeats it and runs at every step. Both initialization
# passes score outflow from the matching carried into the segment, so the first
# pass cannot create new outflow credit for the second.
function interlink_scheduler_initialization(; candidate_force=nothing, initialization_active=nothing)
	scheduler = interlink_scheduler_callback(; candidate_force, initialization_active)
	function initialize(callback, state, time, integrator)
		# This pass precedes the normal planet-frame callback. A resumed run
		# starts at its checkpoint time, so stations must already be rotated to
		# that time before power, saving and other initializers read the selection.
		p = integrator.p
		initialization_active === nothing || (initialization_active[] = copy(p.args.interlink_model.graph.active))
		if has_ground_stations(p.args.interlink_model)
			p.args.environment_model.planet.L_PI .= _planet_lpi_at(p, Float64(time))
		end
		scheduler.initialize(callback, state, time, integrator)
	end
	return DiffEqBase.DiscreteCallback((state, time, integrator) -> false, integrator -> nothing;
		initialize=initialize, save_positions=(false, false))
end
