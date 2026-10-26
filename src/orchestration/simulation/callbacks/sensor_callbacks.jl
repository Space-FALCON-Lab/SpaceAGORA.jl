# Sensor sampling during a run (see Sensors, SensorMounting): each sensor takes
# a measurement from the accepted state at the run's start and every period
# after it, and keeps it for navigation effectors to read; the sensor columns.
# Sensors do not act on the dynamics.

# The field source of a magnetometer: the run's MagneticTorqueRodModel, or its
# default (the tilted dipole) when the run has none.
function _sensor_field_source(effectors)
    index = findfirst(e -> e isa MagneticTorqueRodModel, effectors)
    return index === nothing ? MagneticTorqueRodModel() : effectors[index]
end

# The direction to the Sun in the body frame and the lit fraction of its disk
# at inertial position `pos` and attitude `q`, as the power budget computes
# them; NaN when the run's ephemerides cannot place the Sun.
function _sensor_sun(p, pos::SVector{3, Float64}, q::SVector{4, Float64}, t::Float64)
    planet = p.args.environment_model.planet
    sun_dir = ephemerides_sun_ii(p.args.environment_model.ephemerides_model, planet, p.shared_buffers.et_start[] + t)
    sun_dir === nothing && return (SVector{3, Float64}(NaN, NaN, NaN), NaN)
    sun_ii = SVector{3, Float64}(sun_dir) * mean_solar_distance_m(planet)
    lit = Float64(eclipse_area_calc(pos, sun_ii, Float64(planet.Rp_e)))
    to_sun = sun_ii - pos
    return (rot(q) * (to_sun / norm(to_sun)), lit)
end

# The specific force on spacecraft `i` in the body frame: its acceleration at
# the accepted step (the solver's derivative, no new evaluation) less the
# acceleration of the run's gravitational effectors, each evaluated by the
# route the dynamics take.
function _sensor_specific_force(integrator, i::Int, q::SVector{4, Float64})
    p = integrator.p
    t = Float64(integrator.t)
    du = DiffEqBase.get_du(integrator)
    a = SVector{3, Float64}(du.sc[i].vel[1], du.sc[i].vel[2], du.sc[i].vel[3])
    sc_view = integrator.u.sc[i]
    spacecraft = p.args.dynamics_model.spacecraft[i]
    x = CallbackServices.build_state_sample(sc_view, spacecraft, hasproperty(sc_view, :q), i)
    gravity = SVector{3, Float64}(0.0, 0.0, 0.0)
    scopes = p.args.dynamics_model.effector_scopes
    if scopes === nothing
        for e in p.args.dynamics_model.dynamic_effectors
            is_gravitational(e) && (gravity += CallbackServices._effector_force_ii(e, sc_view, x, p, i, t))
        end
    else
        # Only the gravity applied to this spacecraft: a scoped effector outside
        # its scope applies none.
        for (j, e) in enumerate(p.args.dynamics_model.dynamic_effectors)
            is_gravitational(e) && _simulation_model_module.SpacecraftModels._effector_applies(scopes, j, i) &&
                (gravity += CallbackServices._effector_force_ii(e, sc_view, x, p, i, t))
        end
    end
    return rot(q) * (a - gravity / Float64(sc_view.mass))
end

# The body's angular acceleration at the accepted step (the solver's
# derivative); zero without attitude propagation.
function _sensor_angular_acceleration(integrator, i::Int)::SVector{3, Float64}
    du = DiffEqBase.get_du(integrator)
    hasproperty(du.sc[i], :ω) || return SVector{3, Float64}(0.0, 0.0, 0.0)
    return SVector{3, Float64}(du.sc[i].ω[1], du.sc[i].ω[2], du.sc[i].ω[3])
end

# The truth sensor `s` of spacecraft `i` measures at the accepted state.
function _sensor_truth(integrator, i::Int, s)
    req = sensor_requirements(s)
    wants_field = req.magnetic_field
    p = integrator.p
    sc_view = integrator.u.sc[i]
    pos = SVector{3, Float64}(sc_view.pos[1], sc_view.pos[2], sc_view.pos[3])
    vel = SVector{3, Float64}(sc_view.vel[1], sc_view.vel[2], sc_view.vel[3])
    q = hasproperty(sc_view, :q) ? SVector{4, Float64}(sc_view.q[1], sc_view.q[2], sc_view.q[3], sc_view.q[4]) :
        SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    ω = hasproperty(sc_view, :ω) ? SVector{3, Float64}(sc_view.ω[1], sc_view.ω[2], sc_view.ω[3]) : SVector{3, Float64}(0.0, 0.0, 0.0)
    nan3 = SVector{3, Float64}(NaN, NaN, NaN)
    field, sun, lit, force = nan3, nan3, NaN, nan3
    α, pos_pp, vel_pp, altitude, radius = nan3, nan3, nan3, NaN, NaN
    needs_force = sensor_needs_specific_force(s)
    if wants_field || needs_force
        # The field model and the environment sampling read the time from the
        # shared buffer the right-hand side keeps: set it to the sample's
        # time, then restore it.
        clock = p.shared_buffers.current_time
        previous = clock[]
        clock[] = Float64(integrator.t)
        wants_field && (field = rod_field_body(_sensor_field_source(p.args.dynamics_model.dynamic_effectors), pos, q, p))
        needs_force && (force = _sensor_specific_force(integrator, i, q))
        clock[] = previous
    end
    needs_force && (α = _sensor_angular_acceleration(integrator, i))
    if req.solar
        sun, lit = _sensor_sun(p, pos, q, Float64(integrator.t))
        radius = Float64(p.args.environment_model.planet.Rp_e)
    end
    if req.planet_frame
        frame = CallbackServices.sample_planet_frame(sc_view, p, i, Float64(integrator.t))
        pos_pp, vel_pp, altitude = SVector{3, Float64}(frame.pos_pp), SVector{3, Float64}(frame.vel_pp), Float64(frame.alt_m)
    end
    return SensorTruth(pos, vel, q, ω, Float64(sc_view.mass), Float64(integrator.t), field, sun, lit, force,
                       α, pos_pp, vel_pp, altitude, radius)
end

function _sample_sensors!(integrator, members::Vector{Tuple{Int, Int}}; first_sample::Bool=false)
    spacecraft = integrator.p.args.dynamics_model.spacecraft
    for (i, k) in members
        entry = mounted_sensors(spacecraft[i])[k]
        s = entry.sensor
        truth = _sensor_truth(integrator, i, s)
        ctx = actuator_context(entry.link)
        # The run's first sample: the sensor's own state resets first.
        first_sample && initialize_sensor!(s, truth, ctx, s.rng)
        s.latest = sensor_measure(s, truth, ctx, s.rng)
        s.latest_time_s = Float64(integrator.t)
    end
    return nothing
end

"""
    get_sensor_callbacks(num_sats, args) -> Tuple

For runs with sensors: one callback per distinct sampling period, which
samples its sensors at the run's start and then at every multiple of the
period after it, each a solver stop. The sampling times follow the run's
start, not the solve's, so a checkpointed run samples at the same times as a
run in one piece. A run resumed from a checkpoint file takes the sensors'
random streams, latest measurements, sampling times and gyro biases from the
file and keeps sampling on the same grid; a sensor kind's own state beyond
those fields (what `initialize_sensor!` resets) is in the file when the kind
implements `checkpoint_state`. Empty without sensors.
"""
function get_sensor_callbacks(num_sats::Int, args::SimulationConfiguration)
    spacecraft = args.dynamics_model.spacecraft
    groups = Dict{Float64, Vector{Tuple{Int, Int}}}()
    for i in eachindex(spacecraft), (k, entry) in enumerate(mounted_sensors(spacecraft[i]))
        push!(get!(groups, SensorMounting._check_period(entry.sensor), Tuple{Int, Int}[]), (i, k))
    end
    isempty(groups) && return ()
    sensor_sats = [i for i in eachindex(spacecraft) if !isempty(mounted_sensors(spacecraft[i]))]
    reset_done = Ref(false)
    callbacks = Any[]
    for period in sort!(collect(keys(groups)))
        members = groups[period]
        # Run state, kept across checkpoint segments: the run's start and the
        # number of the next sample, due at start + next * period.
        start = Ref(NaN)
        next = Ref(0)
        due() = start[] + next[] * period
        function sample_now!(integrator; first_sample::Bool=false)
            _sample_sensors!(integrator, members; first_sample=first_sample)
            next[] += 1
            while due() <= integrator.t
                next[] += 1
            end
            _maybe_add_control_tstop!(integrator, due())
            return nothing
        end
        function initialize!(cb, u, t, integrator)
            first_sample = isnan(start[])                # the run's first solve
            if first_sample && integrator.p.shared_buffers.callback_resume[]
                # A resumed run: the sensors' streams, measurements and times
                # come from the checkpoint; the samples keep the run's grid from
                # its start at 0, after the last one taken.
                reset_done[] = true
                start[] = 0.0
                taken = [integrator.p.args.dynamics_model.spacecraft[i] for (i, _) in members]
                times = [mounted_sensors(sc)[k].sensor.latest_time_s for (sc, (_, k)) in zip(taken, members)]
                finite = filter(isfinite, times)
                next[] = isempty(finite) ? 0 : round(Int, maximum(finite) / period) + 1
                first_sample = false
            elseif first_sample
                if !reset_done[]
                    for i in sensor_sats
                        reset_sensors!(integrator.p.args.dynamics_model.spacecraft[i])
                    end
                    reset_done[] = true
                end
                start[] = Float64(t)
                next[] = 0
            end
            # The run's first sample, or one due where the previous segment ended.
            Float64(t) >= due() ? sample_now!(integrator; first_sample=first_sample) : _maybe_add_control_tstop!(integrator, due())
            return nothing
        end
        affect! = integrator -> begin
            sample_now!(integrator)
            _mark_discontinuity!(integrator, false)        # sensors only observe
            return nothing
        end
        push!(callbacks, DiscreteCallback((u, t, integrator) -> t >= due(), affect!;
            initialize=initialize!, save_positions=(false, false)))
    end
    return Tuple(callbacks)
end

"""
    sensor_save_fields(args) -> Vector{SaveField}

The `sensor` field when a spacecraft carries sensors: per spacecraft every
sensor's latest measurement, sensors in mounting order and their components
in order (`sc{i}_sensor_{k}`), NaN before a sensor's first sample.
"""
function sensor_save_fields(args::SimulationConfiguration)
    any(sc -> !isempty(mounted_sensors(sc)), args.dynamics_model.spacecraft) || return SaveField[]
    num_sats = length(args.dynamics_model.spacecraft)
    getter = (u, t, integrator) -> [Float64[x for entry in mounted_sensors(integrator.p.args.dynamics_model.spacecraft[i])
                                             for x in entry.sensor.latest] for i in 1:num_sats]
    return SaveField[SaveField(:sensor, getter; per_satellite=true, column_prefix="sensor")]
end
