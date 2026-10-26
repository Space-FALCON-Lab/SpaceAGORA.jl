# One burn the scheduler refused to plan: when, why, and the command it was given.
const BurnPlanRejection = NamedTuple{(:time_s, :rejection, :delta_v_mps, :direction_rad, :source_orbit),
                                     Tuple{Float64, Symbol, Float64, Float64, Int64}}

"""
    BaseThrusterModel(; thrust, direction, Δv, start_burn_time, stop_burn_time, Isp,
                      actuation=:ideal, engine=1, pointing_frequency_rad_s=0.05, pointing_damping=0.7)

The orbit-maneuver burn scheduler: one entry per spacecraft in each vector.
A burn of `thrust` newtons and specific impulse `Isp` seconds runs between
`start_burn_time` and `stop_burn_time` (seconds into the run), prograde when
`direction` is 0 and retrograde when it is π. With a guidance maneuver
command, or a nonzero `Δv`, the scheduler plans the burn itself, symmetric
about the next apoapsis, after the spacecraft leaves the atmosphere. New requests
are processed on `control_rates`, after guidance due at the same time. A planned
burn keeps its window until completion, including across checkpoints; solver
initialization and atmosphere exits do not create or revise plans. Guidance
requests are identified by `(delta_v_mps, direction_rad, source_orbit)`, so an
unchanged request is accepted only once. A configured `Δv` is likewise used once
unless changed; a fresh run resets this record.

`actuation = :ideal` (the default) applies the burn as a force along the
velocity, whatever the attitude. In both modes the burn switches exactly at
the window's ends: the solver stops there and a callback switches the burn
between integration steps, so no step straddles an end. `actuation = :mounted` fires a thruster
mounted on the spacecraft instead: the `engine`-th mounted thruster (mounting
order) runs at the burn's thrust during the window, so the force follows the
vehicle's attitude and the propellant follows the engine's own `Isp`. From
the moment a burn is planned until it ends, the scheduler also points the
engine along the burn with the spacecraft's mounted reaction wheels, at each
control tick: a proportional-derivative law on the angle between the
engine's thrust direction and the burn direction, with natural frequency
`pointing_frequency_rad_s` and damping ratio `pointing_damping`, scaled by the
spacecraft's inertia and tracking the turn of the velocity along the orbit. A
spacecraft without mounted wheels is not pointed. Give the engine the
scheduler's thrust and specific impulse so the planned burn delivers its Δv.

A burn the scheduler cannot plan (it needs more propellant than the vehicle
carries, its command is not finite, the thrust or specific impulse is not
positive) is not scheduled, and the reason is kept: see
`burn_plan_rejections`. The scheduler also executes a plan that an
`AlgorithmGuidance` of the burn-plan family (`BurnPlanFamily`) installs for a
spacecraft: an explicit installation replaces the previous plan and takes its
window, direction, thrust and specific impulse. For a `BurnPlanFamily` binding,
configure exactly one `BaseThrusterModel`, set that spacecraft's configured `Δv`
to zero and leave its window unset, and do not combine it with the scheduler's
maneuver-guidance models. After setup, explicit installations may replace earlier
plans, including between control ticks. This is distinct from ordinary maneuver
requests, which wait for completion of the accepted burn.

Checkpoints written before control-clock request processing are refused for this
scheduler: their pending clock positions and consumed requests cannot be safely
inferred. Resume them with the version that wrote them, or start a new run.

During a run the scheduler plans in its `Δv`, `direction`, `start_burn_time` and
`stop_burn_time` vectors, which show the run's schedule afterwards. Each run
starts from the schedule configured before it: `initialize_control!` restores
every entry the scheduler wrote in the previous run, unless it was changed since,
and starts the burn state and the refusal records empty.
"""
struct BaseThrusterModel <: AbstractThrusterModel
    thrust::Vector{Float64}
    direction::Vector{Float64}
    Δv::Vector{Float64}
    start_burn_time::Vector{Float64}
    stop_burn_time::Vector{Float64}
    Isp::Vector{Float64}
    actuation::Symbol
    engine::Int
    pointing_frequency_rad_s::Float64
    pointing_damping::Float64
    # Per-spacecraft state: the window with its thrust (start, stop, thrust) and
    # whether the burn is on, flipped only at the window's ends; in `:mounted`
    # actuation also the levels of every mounted actuator (the wheels'
    # refreshed at each control tick, the engine's when the burn switches).
    mounted_levels::Vector{Vector{Float64}}
    burn_windows::Vector{NTuple{3, Float64}}
    burning::Vector{Bool}
    # The burns refused in the run, per spacecraft (`burn_plan_rejections`).
    plan_rejections::Vector{Vector{BurnPlanRejection}}
    # Per spacecraft (Δv, direction, start, stop): the schedule a run started
    # from (empty before the first run), and the scheduler's last write of it in
    # the run (`nothing` before any).
    configured_schedule::Vector{NTuple{4, Float64}}
    written_schedule::Vector{Union{Nothing, NTuple{4, Float64}}}
    # Accepted requests survive completion and resume. source_orbit separates
    # equal guidance requests on different orbits; the configured value prevents
    # an accepted guidance command falling back to the scheduler's own Δv write.
    consumed_command::Vector{Union{Nothing, Tuple{Float64, Float64, Int64}}}
    consumed_configured::Vector{Union{Nothing, NTuple{2, Float64}}}

    function BaseThrusterModel(
        thrust::Vector{Float64},
        direction::Vector{Float64},
        Δv::Vector{Float64},
        start_burn_time::Vector{Float64},
        stop_burn_time::Vector{Float64},
        Isp::Vector{Float64},
        actuation::Symbol,
        engine::Integer,
        pointing_frequency_rad_s::Real,
        pointing_damping::Real,
    )
        n = length(thrust)
        if length(direction) != n || length(Δv) != n || length(start_burn_time) != n || length(stop_burn_time) != n || length(Isp) != n
            throw(ArgumentError("BaseThrusterModel vector lengths must match; got (thrust=$(length(thrust)), direction=$(length(direction)), Δv=$(length(Δv)), start_burn_time=$(length(start_burn_time)), stop_burn_time=$(length(stop_burn_time)), Isp=$(length(Isp)))"))
        end
        (actuation === :ideal || actuation === :mounted) ||
            throw(ArgumentError("BaseThrusterModel actuation must be :ideal or :mounted, got $(repr(actuation))."))
        engine >= 1 || throw(ArgumentError("BaseThrusterModel engine must be a positive mounted-thruster index, got $(engine)."))
        (isfinite(pointing_frequency_rad_s) && pointing_frequency_rad_s > 0.0) ||
            throw(ArgumentError("BaseThrusterModel pointing_frequency_rad_s must be finite and positive, got $(pointing_frequency_rad_s)."))
        (isfinite(pointing_damping) && pointing_damping > 0.0) ||
            throw(ArgumentError("BaseThrusterModel pointing_damping must be finite and positive, got $(pointing_damping)."))

        return new(thrust, direction, Δv, start_burn_time, stop_burn_time, Isp, actuation, Int(engine),
                   Float64(pointing_frequency_rad_s), Float64(pointing_damping),
                   [Float64[] for _ in 1:n], fill((NaN, NaN, 0.0), n), fill(false, n),
                   [BurnPlanRejection[] for _ in 1:n], NTuple{4, Float64}[],
                   Union{Nothing, NTuple{4, Float64}}[nothing for _ in 1:n],
                   Union{Nothing, Tuple{Float64, Float64, Int64}}[nothing for _ in 1:n],
                   Union{Nothing, NTuple{2, Float64}}[nothing for _ in 1:n])
    end
end

BaseThrusterModel(thrust::Vector{Float64}, direction::Vector{Float64}, Δv::Vector{Float64},
                  start_burn_time::Vector{Float64}, stop_burn_time::Vector{Float64}, Isp::Vector{Float64}) =
    BaseThrusterModel(thrust, direction, Δv, start_burn_time, stop_burn_time, Isp, :ideal, 1, 0.05, 0.7)

BaseThrusterModel(; thrust, direction, Δv, start_burn_time, stop_burn_time, Isp,
                  actuation::Symbol=:ideal, engine::Integer=1, pointing_frequency_rad_s::Real=0.05, pointing_damping::Real=0.7) =
    BaseThrusterModel(thrust, direction, Δv, start_burn_time, stop_burn_time, Isp, actuation, engine,
                      pointing_frequency_rad_s, pointing_damping)

"""
    SixAxisThrusterModel

Fixed-direction six-thruster CubeSat actuator description. Directions and
locations are expressed in the spacecraft body frame. Control logic should use
this as hardware data and perform allocation in a control effector.
"""
struct SixAxisThrusterModel <: AbstractThrusterModel
    directions_body::SMatrix{3, 6, Float64}
    locations_body::SMatrix{3, 6, Float64}
    max_thrust_n::SVector{6, Float64}
    isp_s::SVector{6, Float64}
    min_firing_time_s::Float64
end

function SixAxisThrusterModel(;
    directions_body = SMatrix{3, 6, Float64}(
         1.0,  0.0,  0.0,
        -1.0,  0.0,  0.0,
         0.0,  1.0,  0.0,
         0.0, -1.0,  0.0,
         0.0,  0.0,  1.0,
         0.0,  0.0, -1.0,
    ),
    locations_body = SMatrix{3, 6, Float64}(zeros(3, 6)),
    max_thrust_n = SVector{6, Float64}(ones(6)),
    isp_s = SVector{6, Float64}(fill(60.0, 6)),
    min_firing_time_s::Float64 = 0.0,
)
    dirs = SMatrix{3, 6, Float64}(directions_body)
    locs = SMatrix{3, 6, Float64}(locations_body)
    max_thrust = SVector{6, Float64}(max_thrust_n)
    isp = SVector{6, Float64}(isp_s)

    cols = ntuple(j -> begin
        d = SVector{3, Float64}(dirs[:, j])
        d_norm = norm(d)
        d_norm <= eps(Float64) && throw(ArgumentError("SixAxisThrusterModel direction $j is zero."))
        d / d_norm
    end, 6)
    normalized_dirs = SMatrix{3, 6, Float64}(hcat(cols...))

    any(x -> !isfinite(x) || x < 0.0, max_thrust) &&
        throw(ArgumentError("SixAxisThrusterModel max_thrust_n values must be finite and nonnegative."))
    any(x -> !isfinite(x) || x <= 0.0, isp) &&
        throw(ArgumentError("SixAxisThrusterModel isp_s values must be finite and positive."))
    min_firing_time_s < 0.0 &&
        throw(ArgumentError("SixAxisThrusterModel min_firing_time_s must be nonnegative."))

    return SixAxisThrusterModel(normalized_dirs, locs, max_thrust, isp, min_firing_time_s)
end
