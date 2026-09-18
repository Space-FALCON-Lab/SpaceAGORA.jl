# Apollo-style powered descent guidance: the quadratic (E-guidance) law of the
# LM's P63 braking and P64 approach programs, a rate-of-descent P66 for the
# last meters, and the throttle logic of the descent engine. The guidance
# runs in a site-fixed Cartesian frame (x uprange along the approach, y
# crossrange, z up at the site) with velocities relative to the rotating
# body; the commanded thrust
# acceleration and the attitude that points the engine along it are handed
# to the descent control effector through the shared state.
#
# Klumpp, "Apollo Lunar Descent Guidance", Automatica 10 (1974): with
# time-to-go T and phase targets (r_T, v_T, a_T) the commanded acceleration is
#     a_c = a_T - 6 (v_T + v)/T + 12 (r_T - r)/T^2,
# and T is re-solved every cycle from the uprange cubic that adds a jerk
# target j_T: -j_T T^3/24 + a_T T^2/4 - (3 v_T + v) T/4 + (r_T - r) = 0.
# The braking targets sit beyond the high gate and the phase ends when the
# vehicle reaches the configured high-gate altitude.

"""
    DescentPhaseTargets(; name, r_T, v_T, a_T, jerk_uprange=0.0, t_go_initial_s, t_go_min_s, altitude_switch_m=-Inf)

Targets of one quadratic-guidance phase in the site frame (x uprange, y
crossrange, z up; meters, m/s, m/s²). The phase ends when the time-to-go
drops to `t_go_min_s` or the radial terrain clearance drops below `altitude_switch_m`.
"""
Base.@kwdef struct DescentPhaseTargets
    name::Symbol
    r_T::SVector{3, Float64}
    v_T::SVector{3, Float64}
    a_T::SVector{3, Float64}
    jerk_uprange::Float64 = 0.0
    t_go_initial_s::Float64
    t_go_min_s::Float64
    altitude_switch_m::Float64 = -Inf
end

"""
    ApolloDescentConfig(; reference_radius_m, site_lat_deg, site_lon_deg, approach_azimuth_deg, braking, approach, ...)

The explicit reference-sphere radius in metres and landing site (planetocentric degrees, east longitude; `site_height_m` above
the reference sphere, NaN reads it from the terrain model), the heading of
flight at the site (`approach_azimuth_deg`, clockwise from north; Apollo 11
flew west, 270), the braking (P63) and approach (P64) targets, the P66
rate-of-descent schedule, and the descent engine's throttle envelope:
`throttle_min`..`throttle_max` is the throttleable band and a demand above
it commands full thrust in this simplified engine model.
"""
Base.@kwdef struct ApolloDescentConfig
    reference_radius_m::Float64
    site_lat_deg::Float64
    site_lon_deg::Float64
    site_height_m::Float64 = NaN
    approach_azimuth_deg::Float64 = 270.0
    braking::DescentPhaseTargets
    approach::DescentPhaseTargets
    vertical_rate_mps::Float64 = -1.0
    vertical_rate_low_mps::Float64 = -0.5
    vertical_rate_low_altitude_m::Float64 = 10.0
    vertical_gain_v::Float64 = 0.6
    horizontal_gain_v::Float64 = 0.15
    horizontal_gain_r::Float64 = 0.01
    horizontal_accel_max_mps2::Float64 = 0.4
    ignition_trim_s::Float64 = 26.0
    trim_throttle::Float64 = 0.10
    max_thrust_n::Float64 = 45_040.0
    throttle_min::Float64 = 0.10
    throttle_max::Float64 = 0.60
    engine_cutoff_altitude_m::Float64 = 1.7
end

"""
    apollo11_descent_targets(; hg=(7900, 0, 2300), vg=(-152, 0, -46), a_T=(1.5, 0, 1.0), beyond_s=80, ...)

Braking and approach targets in the spirit of Apollo 11: the braking target
is the high gate (`hg`, `vg` meters and m/s in the site frame) extended
`beyond_s` seconds along the high-gate velocity, so the phase ends at the
high-gate altitude instead of at its target; the approach target is a
hover point `approach_altitude_m` above the site descending at 1 m/s.
"""
function apollo11_descent_targets(; hg=(7_900.0, 0.0, 2_300.0), vg=(-152.0, 0.0, -46.0), a_T=(1.5, 0.0, 1.0), beyond_s::Real=80.0,
                                  high_gate_altitude_m::Real=2_500.0, approach_altitude_m::Real=30.0, approach_t_go_s::Real=115.0)
    hgv = SVector{3, Float64}(hg...); vgv = SVector{3, Float64}(vg...)
    braking = DescentPhaseTargets(name=:braking, r_T=hgv + Float64(beyond_s) * vgv, v_T=vgv, a_T=SVector{3, Float64}(a_T...),
        t_go_initial_s=560.0, t_go_min_s=0.9 * Float64(beyond_s), altitude_switch_m=Float64(high_gate_altitude_m))
    approach = DescentPhaseTargets(name=:approach, r_T=SVector{3, Float64}(0.0, 0.0, Float64(approach_altitude_m)), v_T=SVector{3, Float64}(0.0, 0.0, -1.0),
        a_T=SVector{3, Float64}(0.0, 0.0, 0.0), t_go_initial_s=Float64(approach_t_go_s), t_go_min_s=4.0)
    return braking, approach
end

const DESCENT_PHASES = (:braking, :approach, :vertical, :landed)

"""
    ApolloDescentState(num_sats)

Per-spacecraft guidance state shared with the control effector: the phase,
time-to-go, the commanded thrust (N), its inertial direction, the commanded
attitude (scalar-last, body to inertial) and body rate, the radial terrain clearance
above the terrain, the site-frame position and velocity, the phase start
times and the touchdown record.
"""
mutable struct ApolloDescentState
    phase::Vector{Symbol}
    t_go_s::Vector{Float64}
    thrust_cmd_n::Vector{Float64}
    throttle::Vector{Float64}
    thrust_dir_i::Vector{SVector{3, Float64}}
    attitude_cmd::Vector{SVector{4, Float64}}
    rate_cmd::Vector{SVector{3, Float64}}
    radar_altitude_m::Vector{Float64}
    site_r_m::Vector{SVector{3, Float64}}
    site_v_mps::Vector{SVector{3, Float64}}
    demand_fraction::Vector{Float64}
    phase_start_s::Matrix{Float64}          # 4 x num_sats, NaN until entered
    last_update_s::Vector{Float64}
    site_frame::Vector{Union{Nothing, NamedTuple}}
    touchdown_s::Vector{Float64}
    touchdown_v_mps::Vector{SVector{3, Float64}}
    touchdown_miss_m::Vector{Float64}
end

function ApolloDescentState(num_sats::Integer)
    n = Int(num_sats)
    n >= 1 || throw(ArgumentError("ApolloDescentState needs at least one spacecraft"))
    z3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    return ApolloDescentState(fill(:braking, n), fill(NaN, n), zeros(n), zeros(n), fill(SVector{3, Float64}(0.0, 0.0, 1.0), n),
        fill(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), n), fill(z3, n), fill(NaN, n), fill(z3, n), fill(z3, n), zeros(n),
        fill(NaN, 4, n), fill(NaN, n), Vector{Union{Nothing, NamedTuple}}(nothing, n), fill(NaN, n), fill(z3, n), fill(NaN, n))
end

"""
    ApolloDescentGuidanceModel(config, state, terrain=NoTerrainModel())

Guidance effector: runs the quadratic guidance at the guidance rate, reads
the radial terrain clearance from the terrain model, and writes the thrust and
attitude commands into `state` for [`ApolloDescentControlModel`](@ref).
"""
# Both models use the same selected run indices and explicit terrain datum.
function _descent_indices(state::ApolloDescentState, indices)
    n = length(state.phase)
    n > 0 || throw(ArgumentError("descent state must contain spacecraft"))
    for name in fieldnames(ApolloDescentState)
        data = getfield(state, name)
        (name === :phase_start_s ? size(data) == (4, n) : length(data) == n) ||
            throw(ArgumentError("descent state field $name has the wrong size"))
    end
    requested = collect(indices)
    all(i -> i isa Integer && !(i isa Bool) && 1 <= i <= n, requested) ||
        throw(ArgumentError("spacecraft_indices must contain run indices within the descent state"))
    ids = Int.(requested)
    isempty(ids) && throw(ArgumentError("select at least one descent spacecraft"))
    length(unique(ids)) == length(ids) || throw(ArgumentError("duplicate descent spacecraft index"))
    return Tuple(ids)
end

function _validate_descent_config(c::ApolloDescentConfig, terrain::AbstractTerrainModel)
    isfinite(c.reference_radius_m) && c.reference_radius_m > 0 ||
        throw(ArgumentError("descent reference_radius_m must be finite and positive"))
    terrain_radius(terrain, c.site_lat_deg, c.site_lon_deg, c.reference_radius_m)
    (isnan(c.site_height_m) || isfinite(c.site_height_m)) || throw(ArgumentError("site_height_m must be finite or NaN"))
    isfinite(c.approach_azimuth_deg) || throw(ArgumentError("approach azimuth must be finite"))
    for name in (:vertical_rate_mps, :vertical_rate_low_mps)
        value = getfield(c, name)
        isfinite(value) && value <= 0 || throw(ArgumentError("$name must be finite and nonpositive"))
    end
    for name in (:vertical_rate_low_altitude_m, :vertical_gain_v, :horizontal_gain_v,
                 :horizontal_gain_r, :horizontal_accel_max_mps2, :ignition_trim_s,
                 :engine_cutoff_altitude_m)
        value = getfield(c, name)
        isfinite(value) && value >= 0 || throw(ArgumentError("$name must be finite and nonnegative"))
    end
    isfinite(c.max_thrust_n) && c.max_thrust_n > 0 || throw(ArgumentError("max_thrust_n must be finite and positive"))
    0 <= c.throttle_min <= c.throttle_max <= 1 || throw(ArgumentError("invalid descent throttle envelope"))
    isfinite(c.trim_throttle) && 0 <= c.trim_throttle <= 1 || throw(ArgumentError("trim_throttle must be within [0,1]"))
    for ph in (c.braking, c.approach)
        all(isfinite, ph.r_T) && all(isfinite, ph.v_T) && all(isfinite, ph.a_T) &&
            isfinite(ph.jerk_uprange) || throw(ArgumentError("phase targets must be finite"))
        isfinite(ph.t_go_initial_s) && isfinite(ph.t_go_min_s) &&
            1 <= ph.t_go_min_s < ph.t_go_initial_s || throw(ArgumentError("phase time bounds must be finite and ordered, with t_go_min_s at least one second"))
        !isnan(ph.altitude_switch_m) || throw(ArgumentError("phase altitude gate must not be NaN"))
    end
    return nothing
end

struct ApolloDescentGuidanceModel{T <: AbstractTerrainModel, N} <: AbstractGuidanceModel
    config::ApolloDescentConfig
    state::ApolloDescentState
    terrain::T
    spacecraft_indices::NTuple{N, Int}
end
function ApolloDescentGuidanceModel(config::ApolloDescentConfig, state::ApolloDescentState,
                                    terrain::AbstractTerrainModel=NoTerrainModel();
                                    spacecraft_indices=eachindex(state.phase))
    _validate_descent_config(config, terrain)
    return ApolloDescentGuidanceModel(config, state, terrain, _descent_indices(state, spacecraft_indices))
end

# ---- site frame ---------------------------------------------------------------

"Unit up, north and east vectors in planet-fixed axes at a latitude and east longitude (degrees)."
@inline function _site_axes(lat_deg::Float64, lon_deg::Float64)
    φ = deg2rad(lat_deg); λ = deg2rad(mod(lon_deg, 360.0))
    up = SVector{3, Float64}(cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ))
    east = SVector{3, Float64}(-sin(λ), cos(λ), 0.0)
    north = SVector{3, Float64}(-sin(φ) * cos(λ), -sin(φ) * sin(λ), cos(φ))
    return up, north, east
end

"""
    descent_site_frame(config, terrain, reference_radius_m) -> NamedTuple

Site position (planet-fixed, m) and the site-frame axes: `x` uprange (against
the approach heading), `y` crossrange (left of the approach), `z` up.
"""
function descent_site_frame(config::ApolloDescentConfig, terrain::AbstractTerrainModel, reference_radius_m::Real)
    _validate_descent_config(config, terrain)
    Float64(reference_radius_m) == config.reference_radius_m || throw(ArgumentError("site frame reference radius differs from descent configuration"))
    up, north, east = _site_axes(config.site_lat_deg, config.site_lon_deg)
    h = isnan(config.site_height_m) ? terrain_height(terrain, config.site_lat_deg, config.site_lon_deg) : config.site_height_m
    isfinite(config.reference_radius_m + h) && config.reference_radius_m + h > 0 ||
        throw(ArgumentError("site radius must be finite and positive"))
    az = deg2rad(mod(config.approach_azimuth_deg, 360.0))
    flight = normalize(cos(az) * north + sin(az) * east)
    x = -flight
    y = normalize(cross(up, x))
    x = normalize(cross(y, up))
    return (origin_p=(Float64(reference_radius_m) + h) * up, x=x, y=y, z=up, height_m=h, reference_radius_m=Float64(reference_radius_m))
end

@inline function _to_site(frame, v_p::SVector{3, Float64})::SVector{3, Float64}
    return SVector{3, Float64}(dot(v_p, frame.x), dot(v_p, frame.y), dot(v_p, frame.z))
end
@inline function _from_site(frame, v_s::SVector{3, Float64})::SVector{3, Float64}
    return v_s[1] * frame.x + v_s[2] * frame.y + v_s[3] * frame.z
end

# ---- the law ----------------------------------------------------------------------

"Klumpp's quadratic guidance: the commanded (frame-relative) acceleration at time-to-go `T`."
@inline function descent_accel_command(r::SVector{3, Float64}, v::SVector{3, Float64}, ph::DescentPhaseTargets, T::Float64)::SVector{3, Float64}
    return ph.a_T - 6.0 * (ph.v_T + v) / T + 12.0 * (ph.r_T - r) / T^2
end

"Time-to-go from the uprange cubic with the phase's jerk target, by Newton from the previous value."
function descent_time_to_go(r::SVector{3, Float64}, v::SVector{3, Float64}, ph::DescentPhaseTargets, T::Float64)::Float64
    j = ph.jerk_uprange; aT = ph.a_T[1]; vT = ph.v_T[1]; rT = ph.r_T[1]
    f(T) = -j * T^3 / 24 + aT * T^2 / 4 - (3 * vT + v[1]) * T / 4 + (rT - r[1])
    df(T) = -j * T^2 / 8 + aT * T / 2 - (3 * vT + v[1]) / 4
    T = max(T, 1.0)
    for _ in 1:12
        d = df(T)
        abs(d) < 1e-12 && break
        Tn = clamp(T - f(T) / d, 0.5 * T, 1.5 * T)
        if abs(Tn - T) < 1e-4
            T = Tn
            break
        end
        T = Tn
    end
    return max(T, 1.0)
end

"""
    descent_attitude_command(thrust_dir_i, up_i, forward_i=zero) -> SVector{4}

Attitude (scalar-last, body to inertial) that points the descent engine's
thrust axis, body -z, along `thrust_dir_i`, rolled so the body +x axis (the
windows) lies along the projection of `up_i + forward_i` on the plane
normal to the thrust: windows up while the engine fires retrograde during
braking, and facing the site along the direction of flight once the
vehicle pitches back to vertical. The sum keeps the roll reference well
conditioned through the whole descent; a zero `forward_i` uses up alone
and a fixed horizontal seed when the thrust is vertical.
"""
function descent_attitude_command(thrust_dir_i::SVector{3, Float64}, up_i::SVector{3, Float64}, forward_i::SVector{3, Float64}=SVector{3, Float64}(0.0, 0.0, 0.0))::SVector{4, Float64}
    all(isfinite, thrust_dir_i) && norm(thrust_dir_i) > eps(Float64) ||
        throw(ArgumentError("thrust direction must be finite and nonzero"))
    all(isfinite, up_i) && all(isfinite, forward_i) || throw(ArgumentError("attitude references must be finite"))
    z_b = -normalize(thrust_dir_i)
    ref = up_i + forward_i
    x_ref = ref - dot(ref, z_b) * z_b
    if norm(x_ref) < 1e-6
        # thrust along the reference: any perpendicular direction will do
        seed = abs(z_b[1]) < 0.9 ? SVector{3, Float64}(1.0, 0.0, 0.0) : SVector{3, Float64}(0.0, 1.0, 0.0)
        x_ref = seed - dot(seed, z_b) * z_b
    end
    x_b = normalize(x_ref)
    y_b = normalize(cross(z_b, x_b))
    # rows are the body axes in inertial coordinates: the passive inertial-to-body matrix rot(q)
    A = SMatrix{3, 3, Float64}(x_b[1], y_b[1], z_b[1], x_b[2], y_b[2], z_b[2], x_b[3], y_b[3], z_b[3])
    return quaternion_from_passive_dcm(A)
end

"Scalar-last quaternion q with rot(q) == A, for a passive (inertial-to-body) rotation matrix A (Shepperd's method)."
function quaternion_from_passive_dcm(A::SMatrix{3, 3, Float64})::SVector{4, Float64}
    return SVector{4, Float64}(dcm_to_quaternion(A))
end

"Throttle fraction of the descent engine for a thrust demand (N): the throttleable band, or full thrust above it."
@inline function descent_throttle(config::ApolloDescentConfig, demand_n::Float64)::Float64
    frac = demand_n / config.max_thrust_n
    frac > config.throttle_max && return 1.0
    return clamp(frac, config.throttle_min, config.throttle_max)
end

@inline _phase_index(phase::Symbol)::Int = phase === :braking ? 1 : phase === :approach ? 2 : phase === :vertical ? 3 : 4

@inline function _descent_ephemeris_time(p::ODEParams, t::Float64)::Float64
    if hasproperty(p, :shared_buffers) && hasproperty(p.shared_buffers, :et_start)
        return p.shared_buffers.et_start[] + t
    end
    return t
end

@inline function _descent_sat_state(u, i::Int)
    sc = hasproperty(u, :sc) ? u.sc[i] : u
    pos = SVector{3, Float64}(sc.pos)
    vel = SVector{3, Float64}(sc.vel)
    mass = hasproperty(sc, :mass) ? Float64(sc.mass) : NaN
    return pos, vel, mass
end

"""
    descent_environment(model, u, p, t, i) -> NamedTuple

Planet-fixed and site-frame state of spacecraft `i`, its radial terrain clearance
above the terrain, and the frame matrices the commands are expressed in.
"""
function descent_environment(model::ApolloDescentGuidanceModel, u, p::ODEParams, t::Float64, i::Int)
    pos, vel, mass = _descent_sat_state(u, i)
    planet = p.args.environment_model.planet
    ephem = p.args.environment_model.ephemerides_model
    et = _descent_ephemeris_time(p, t)
    r_p, v_p = r_intor_p!(pos, vel, planet, et, ephem)
    l_pi = planet_frame_lpi(planet, et, ephem)
    i in model.spacecraft_indices || throw(ArgumentError("spacecraft index $i is not selected for descent"))
    frame = model.state.site_frame[i]
    if frame === nothing
        frame = descent_site_frame(model.config, model.terrain, model.config.reference_radius_m)
        model.state.site_frame[i] = frame
    end
    lla = rtolatlongrad(r_p, planet)
    lat_deg = rad2deg(lla[2]); lon_deg = rad2deg(lla[3])
    ground_radius = terrain_radius(model.terrain, lat_deg, lon_deg, model.config.reference_radius_m)
    radar = norm(r_p) - ground_radius
    r_s = _to_site(frame, r_p - frame.origin_p)
    v_s = _to_site(frame, v_p)
    return (pos_i=pos, vel_i=vel, mass=mass, r_p=r_p, v_p=v_p, l_pi=l_pi, frame=frame, radar_m=radar, r_s=r_s, v_s=v_s,
        lat_deg=lat_deg, lon_deg=lon_deg, planet=planet)
end

function calcGuidanceEffect!(model::ApolloDescentGuidanceModel, u, p::ODEParams, t::Float64, i::Int64)
    state = model.state
    config = model.config
    (i in model.spacecraft_indices && p.is_active[i]) || return nothing
    env = descent_environment(model, u, p, t, i)
    r_s, v_s, m = env.r_s, env.v_s, env.mass
    isfinite(m) && m > 0 || throw(ArgumentError("descent guidance requires finite positive spacecraft mass"))
    if isnan(state.phase_start_s[1, i])
        state.phase_start_s[1, i] = t
        state.t_go_s[i] = config.braking.t_go_initial_s
        state.last_update_s[i] = t
    end
    dt = max(0.0, t - state.last_update_s[i])
    state.last_update_s[i] = t
    state.radar_altitude_m[i] = env.radar_m
    state.site_r_m[i] = r_s
    state.site_v_mps[i] = v_s
    phase = state.phase[i]

    # gravity and the rotating-frame terms in planet-fixed axes, expressed in the site frame
    planet = env.planet
    g_p = -planet.μ * env.r_p / norm(env.r_p)^3
    ω = SVector{3, Float64}(planet.ω)
    a_frame_p = -2.0 * cross(ω, env.v_p) - cross(ω, cross(ω, env.r_p))   # a_rel = a_thrust + g + a_frame
    g_s = _to_site(env.frame, g_p + a_frame_p)

    if phase === :landed
        state.thrust_cmd_n[i] = 0.0; state.throttle[i] = 0.0; state.demand_fraction[i] = 0.0
        return nothing
    end

    a_c = SVector{3, Float64}(0.0, 0.0, 0.0)
    # Quadratic phases: re-solve the time-to-go, hand over to the next phase when
    # it runs out or the altitude gate is passed, and evaluate the new phase at
    # once so a phase never commands from targets it has already passed.
    T_prev = state.t_go_s[i] - dt
    while phase === :braking || phase === :approach
        ph = phase === :braking ? config.braking : config.approach
        T = descent_time_to_go(r_s, v_s, ph, T_prev)
        if T <= ph.t_go_min_s || env.radar_m <= ph.altitude_switch_m
            phase = phase === :braking ? :approach : :vertical
            state.phase[i] = phase
            state.phase_start_s[_phase_index(phase), i] = t
            T_prev = phase === :approach ? config.approach.t_go_initial_s : NaN
            continue
        end
        state.t_go_s[i] = T
        a_c = descent_accel_command(r_s, v_s, ph, T)
        break
    end
    if phase === :vertical
        state.t_go_s[i] = NaN
        vz_cmd = env.radar_m > config.vertical_rate_low_altitude_m ? config.vertical_rate_mps : config.vertical_rate_low_mps
        ax = -config.horizontal_gain_v * v_s[1] - config.horizontal_gain_r * r_s[1]
        ay = -config.horizontal_gain_v * v_s[2] - config.horizontal_gain_r * r_s[2]
        ah = hypot(ax, ay)
        if ah > config.horizontal_accel_max_mps2
            ax *= config.horizontal_accel_max_mps2 / ah; ay *= config.horizontal_accel_max_mps2 / ah
        end
        a_c = SVector{3, Float64}(ax, ay, config.vertical_gain_v * (vz_cmd - v_s[3]))
        if env.radar_m <= config.engine_cutoff_altitude_m
            state.phase[i] = :landed
            state.phase_start_s[4, i] = t
            state.thrust_cmd_n[i] = 0.0; state.throttle[i] = 0.0; state.demand_fraction[i] = 0.0
            return nothing
        end
    end

    a_thrust_s = a_c - g_s
    demand = m * norm(a_thrust_s)
    state.demand_fraction[i] = demand / config.max_thrust_n
    throttle = t - state.phase_start_s[1, i] < config.ignition_trim_s ? config.trim_throttle : descent_throttle(config, demand)
    state.throttle[i] = throttle
    state.thrust_cmd_n[i] = throttle * config.max_thrust_n
    dir_s = norm(a_thrust_s) > 1e-9 ? a_thrust_s / norm(a_thrust_s) : SVector{3, Float64}(0.0, 0.0, 1.0)
    dir_p = _from_site(env.frame, dir_s)
    dir_i = normalize(SVector{3, Float64}(env.l_pi' * dir_p))
    state.thrust_dir_i[i] = dir_i
    up_i = normalize(env.pos_i)
    forward_i = normalize(SVector{3, Float64}(env.l_pi' * (-env.frame.x)))   # direction of flight at the site
    state.attitude_cmd[i] = descent_attitude_command(dir_i, up_i, forward_i)
    state.rate_cmd[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
    return nothing
end
