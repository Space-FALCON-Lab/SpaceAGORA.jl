using LinearAlgebra
using OrdinaryDiffEq
using Printf
using Statistics
using SpaceAGORA

const EARTH_KS_PARAMS = KSPropagationParams(
    Re=6_378_137.0,
    μ=3.986004418e14,
    J2=1.08262668e-3,
    Ω=7.292115e-5,
)

function perigee_state(a_m::Real, eccentricity::Real, inclination_deg::Real, params)
    radius = Float64(a_m) * (1.0 - Float64(eccentricity))
    speed = sqrt(params.μ * (2.0 / radius - 1.0 / Float64(a_m)))
    inclination = deg2rad(Float64(inclination_deg))
    return [radius, 0.0, 0.0, 0.0, speed * cos(inclination), speed * sin(inclination)]
end

function cartesian_rhs(state, params)
    position = @view state[1:3]
    velocity = @view state[4:6]
    radius = norm(position)
    acceleration = -(params.μ / radius^3) .* position .+
        ks_j2_acceleration_si(position, params)
    return vcat(velocity, acceleration)
end

function cartesian_rk4_step(state, params, step_s)
    k1 = cartesian_rhs(state, params)
    k2 = cartesian_rhs(state .+ 0.5 * step_s .* k1, params)
    k3 = cartesian_rhs(state .+ 0.5 * step_s .* k2, params)
    k4 = cartesian_rhs(state .+ step_s .* k3, params)
    return state .+ (step_s / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

function propagate_cartesian(initial_state, params, final_time_s, nominal_step_s)
    state = copy(initial_state)
    time_s = 0.0
    times = Float64[time_s]
    states = Vector{Vector{Float64}}([copy(state)])
    rhs_evaluations = 0
    while time_s < final_time_s
        step_s = min(nominal_step_s, final_time_s - time_s)
        state = cartesian_rk4_step(state, params, step_s)
        time_s += step_s
        rhs_evaluations += 4
        push!(times, time_s)
        push!(states, copy(state))
    end
    return (; times, states, rhs_evaluations)
end

function final_ks_step(state, params, nominal_delta_s, final_time_s)
    trial = ks_rk4_step(state, params, 0.0, nominal_delta_s)
    rhs_evaluations = 4
    trial[10] <= final_time_s && return trial, rhs_evaluations
    remaining_time_s = final_time_s - state[10]
    delta_s = nominal_delta_s * remaining_time_s / (trial[10] - state[10])
    for _ in 1:6
        trial = ks_rk4_step(state, params, 0.0, delta_s)
        rhs_evaluations += 4
        if abs(trial[10] - final_time_s) <= 1.0e-10
            trial[10] = final_time_s
            return trial, rhs_evaluations
        end
        delta_s *= remaining_time_s / (trial[10] - state[10])
    end
    trial[10] = final_time_s
    return trial, rhs_evaluations
end

function propagate_ks(initial_state, params, final_time_s, nominal_step_s, semimajor_axis_m)
    state = cartesian_to_ks_state(initial_state[1:3], initial_state[4:6], params)
    nominal_delta_s = nominal_step_s / semimajor_axis_m
    times = Float64[state[10]]
    states = Vector{Vector{Float64}}([copy(initial_state)])
    rhs_evaluations = 0
    while state[10] < final_time_s
        next_state, step_rhs_evaluations = final_ks_step(
            state, params, nominal_delta_s, final_time_s)
        rhs_evaluations += step_rhs_evaluations
        state = next_state
        cartesian = ks_state_to_cartesian(state)
        push!(times, state[10])
        push!(states, vcat(cartesian.position_ii_m, cartesian.velocity_ii_m))
    end
    return (; times, states, rhs_evaluations)
end

function specific_energy_j2(state, params)
    position = @view state[1:3]
    velocity = @view state[4:6]
    radius = norm(position)
    z_ratio = position[3] / radius
    potential = -params.μ / radius +
        params.μ * params.J2 * params.Re^2 / (2.0 * radius^3) * (3.0 * z_ratio^2 - 1.0)
    return 0.5 * dot(velocity, velocity) + potential
end

function accuracy_metrics(result, reference, params)
    position_errors = Float64[]
    velocity_errors = Float64[]
    energies = Float64[]
    angular_momentum_z = Float64[]
    for (time_s, state) in zip(result.times, result.states)
        truth = reference(time_s)
        push!(position_errors, norm(state[1:3] - truth[1:3]))
        push!(velocity_errors, norm(state[4:6] - truth[4:6]))
        push!(energies, specific_energy_j2(state, params))
        push!(angular_momentum_z, cross(state[1:3], state[4:6])[3])
    end
    energy_scale = max(abs(first(energies)), eps(Float64))
    hz_scale = max(abs(first(angular_momentum_z)), eps(Float64))
    return (
        endpoint_position_error_m=last(position_errors),
        endpoint_velocity_error_m_s=last(velocity_errors),
        maximum_position_error_m=maximum(position_errors),
        maximum_velocity_error_m_s=maximum(velocity_errors),
        maximum_relative_energy_drift=maximum(abs.(energies .- first(energies))) / energy_scale,
        maximum_relative_hz_drift=maximum(abs.(angular_momentum_z .- first(angular_momentum_z))) / hz_scale,
    )
end

function timed_median(run; repetitions=7, batch_size=50)
    run() # warm compilation and caches
    samples = [@timed begin
        result = run()
        for _ in 2:batch_size
            result = run()
        end
        result
    end for _ in 1:repetitions]
    middle = sortperm(getproperty.(samples, :time))[cld(repetitions, 2)]
    sample = samples[middle]
    return sample.value, sample.time / batch_size, sample.bytes / batch_size
end

function run_case(case, params, steps_s; repetitions=7, batch_size=50)
    initial_state = perigee_state(case.a_m, case.e, case.i_deg, params)
    period_s = 2.0 * pi * sqrt(case.a_m^3 / params.μ)
    final_time_s = case.orbits * period_s
    reference_problem = ODEProblem(
        (state, _, _) -> cartesian_rhs(state, params),
        initial_state,
        (0.0, final_time_s),
    )
    reference = solve(
        reference_problem,
        Vern9();
        abstol=1.0e-9,
        reltol=1.0e-13,
        dense=true,
        save_everystep=true,
    )
    tighter_reference = solve(
        reference_problem,
        Vern9();
        abstol=1.0e-10,
        reltol=3.0e-14,
        save_everystep=false,
    )
    reference_endpoint_position_difference_m = norm(
        reference(final_time_s)[1:3] - tighter_reference.u[end][1:3])
    reference_endpoint_velocity_difference_m_s = norm(
        reference(final_time_s)[4:6] - tighter_reference.u[end][4:6])

    rows = NamedTuple[]
    for step_s in steps_s
        cartesian, cartesian_time, cartesian_bytes = timed_median(
            () -> propagate_cartesian(initial_state, params, final_time_s, step_s);
            repetitions=repetitions, batch_size=batch_size,
        )
        ks, ks_time, ks_bytes = timed_median(
            () -> propagate_ks(initial_state, params, final_time_s, step_s, case.a_m);
            repetitions=repetitions, batch_size=batch_size,
        )
        for (method, result, runtime_s, bytes) in (
            (:Cartesian_RK4, cartesian, cartesian_time, cartesian_bytes),
            (:KS_RK4, ks, ks_time, ks_bytes),
        )
            metrics = accuracy_metrics(result, reference, params)
            push!(rows, merge(
                (; case=case.name, method, nominal_step_s=step_s,
                    physical_duration_s=final_time_s, accepted_steps=length(result.times) - 1,
                    rhs_evaluations=result.rhs_evaluations, runtime_s,
                    allocated_bytes=bytes),
                metrics,
            ))
        end
    end
    return (
        rows=rows,
        validation=(;
            case=case.name,
            reference_endpoint_position_difference_m,
            reference_endpoint_velocity_difference_m_s,
        ),
    )
end

function print_report(rows)
    println("| Case | Method | Nominal step (s) | Steps | RHS evals | Runtime (ms) | Max position error (m) | End position error (m) | Max velocity error (m/s) | Relative energy drift | Relative hₙ drift | Allocated (MiB) |")
    println("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows
        @printf("| %s | %s | %.1f | %d | %d | %.3f | %.6g | %.6g | %.6g | %.3e | %.3e | %.3f |\n",
            row.case, row.method, row.nominal_step_s, row.accepted_steps,
            row.rhs_evaluations, 1.0e3 * row.runtime_s,
            row.maximum_position_error_m, row.endpoint_position_error_m,
            row.maximum_velocity_error_m_s, row.maximum_relative_energy_drift,
            row.maximum_relative_hz_drift, row.allocated_bytes / 2.0^20)
    end
end

function print_reference_validation(validations)
    println("\nReference convergence check (baseline Vern9 versus tighter Vern9):")
    println("| Case | Endpoint position difference (m) | Endpoint velocity difference (m/s) |")
    println("|---|---:|---:|")
    for result in validations
        @printf("| %s | %.6g | %.6g |\n", result.case,
            result.reference_endpoint_position_difference_m,
            result.reference_endpoint_velocity_difference_m_s)
    end
end

quick = lowercase(get(ENV, "SPACEAGORA_PROPAGATOR_BENCHMARK_QUICK", "0")) in ("1", "true", "yes")
cases = (
    (name="near-circular LEO", a_m=7_000_000.0, e=0.01, i_deg=51.6, orbits=quick ? 1 : 10),
    (name="high-eccentricity Earth orbit", a_m=(6_678_137.0 + 42_164_000.0) / 2.0,
        e=(42_164_000.0 - 6_678_137.0) / (42_164_000.0 + 6_678_137.0),
        i_deg=63.4, orbits=quick ? 1 : 3),
)
steps_s = quick ? (120.0,) : (120.0, 60.0, 30.0)
repetitions = quick ? 1 : 7
batch_size = quick ? 1 : 50
results = collect(
    run_case(case, EARTH_KS_PARAMS, steps_s; repetitions, batch_size) for case in cases
)
rows = reduce(vcat, getproperty.(results, :rows))
print_report(rows)
print_reference_validation(getproperty.(results, :validation))
