using Test
using LinearAlgebra
using SpaceAGORA

@testset "reusable KS dynamics" begin
    params = KSPropagationParams(
        Re=6.378137e6,
        μ=3.986004418e14,
        J2=1.08262668e-3,
        Ω=7.292115e-5,
    )
    mpc_params = AerobrakingMPCParams(
        Re=params.Re,
        μ=params.μ,
        J2=params.J2,
        Ω=params.Ω,
    )
    @test ks_j2_acceleration_si([7.0e6, 0.0, 0.0], mpc_params) ≈
        ks_j2_acceleration_si([7.0e6, 0.0, 0.0], params)

    cases = (
        ([7.0e6, 0.0, 0.0], [0.0, 7.5e3, 250.0]),
        ([-7.0e6, 1.0e3, -2.0e3], [-10.0, -7.5e3, 50.0]),
        ([2.0e6, -6.8e6, 3.0e5], [7.1e3, 2.0e3, -400.0]),
    )
    for (position, velocity) in cases
        state = cartesian_to_ks_state(position, velocity, params; elapsed_time_s=42.0)
        cartesian = ks_state_to_cartesian(state)
        @test cartesian.position_ii_m ≈ position rtol=2.0e-14 atol=2.0e-7
        @test cartesian.velocity_ii_m ≈ velocity rtol=2.0e-14 atol=2.0e-9
        @test cartesian.elapsed_time_s == 42.0
        @test ks_position(state[1:4]) ≈ position rtol=2.0e-14 atol=2.0e-7
        @test ks_velocity(state[1:4], state[5:8]) ≈ velocity rtol=2.0e-14 atol=2.0e-9
    end

    circular_position = [7.0e6, 0.0, 0.0]
    circular_velocity = [0.0, sqrt(params.μ / circular_position[1]), 0.0]
    state = cartesian_to_ks_state(circular_position, circular_velocity, params)
    specific_energy = 0.5 * dot(circular_velocity, circular_velocity) -
        params.μ / norm(circular_position)
    @test state[9] ≈ -specific_energy rtol=3.0 * eps(Float64) atol=0.0
    @test state[9] ≈ ks_energy_parameter(specific_energy) rtol=3.0 * eps(Float64) atol=0.0
    @test specific_energy_from_ks(state[9]) ≈ specific_energy rtol=3.0 * eps(Float64) atol=0.0
    @test ks_state_to_cartesian(state).specific_energy_j_kg ≈
        specific_energy rtol=3.0 * eps(Float64) atol=0.0
    next_state = ks_rk4_step(state, params, 0.0, 1.0e-7)
    @test all(isfinite, next_state)
    @test next_state[10] > state[10]

    drag = ks_drag_acceleration_si(
        circular_position,
        circular_velocity,
        params,
        10.0;
        density_kg_m3=1.0e-8,
        drag_coefficient=2.2,
        mass_kg=500.0,
    )
    relative_velocity = circular_velocity - [0.0, params.Ω * circular_position[1], 0.0]
    @test dot(drag, relative_velocity) < 0.0
    @test all(isfinite, ks_j2_acceleration_si(circular_position, params))
    j2_jacobian = ks_j2_acceleration_jacobian_si(circular_position, params)
    @test size(j2_jacobian) == (3, 3)
    finite_difference_j2 = hcat((
        (ks_j2_acceleration_si(circular_position .+ 0.1 .* unit, params) -
         ks_j2_acceleration_si(circular_position .- 0.1 .* unit, params)) ./ 0.2
        for unit in eachcol(Matrix{Float64}(I, 3, 3))
    )...)
    @test j2_jacobian ≈ finite_difference_j2 rtol=2.0e-7 atol=1.0e-15
    @test size(ks_rhs_jacobian(state, params)) == (10, 10)
    @test size(ks_step_jacobian(state, params, 0.0, 1.0e-7)) == (10, 10)
    @test_throws ArgumentError cartesian_to_ks_state(zeros(3), circular_velocity, params)
    @test_throws ArgumentError ks_drag_acceleration_si(
        circular_position,
        circular_velocity,
        params,
        10.0;
        density_kg_m3=1.0e-8,
        drag_coefficient=2.2,
        mass_kg=0.0,
    )
end

@testset "h_KS=-ε convention preserves the former physical trajectory" begin
    params = KSPropagationParams(
        Re=6.378137e6,
        μ=3.986004418e14,
        J2=1.08262668e-3,
        Ω=7.292115e-5,
    )
    position = [7.0e6, 0.0, 2.0e5]
    velocity = [0.0, 7.4e3, 500.0]
    new_state = cartesian_to_ks_state(position, velocity, params)
    old_state = copy(new_state)
    old_state[9] = 2.0 * new_state[9]

    ks_dynamics = SpaceAGORA.SimulationModel.DynamicsKS
    new_rhs = ks_rhs(new_state, params)
    u = new_state[1:4]
    u_prime = new_state[5:8]
    radius = dot(u, u)
    acceleration = ks_j2_acceleration_si(ks_position(u), params)
    acceleration4 = [acceleration[1], acceleration[2], acceleration[3], 0.0]
    expected_u_prime_derivative = -0.5 .* new_state[9] .* u .+
        0.5 .* radius .* (transpose(ks_dynamics._ks_L(u)) * acceleration4)
    @test new_rhs[5:8] ≈ expected_u_prime_derivative rtol=2.0e-14 atol=2.0e-10
    @test new_rhs[9] ≈ -radius * dot(ks_velocity(u, u_prime), acceleration)
    @test sqrt(new_state[9] / 2.0) == sqrt(old_state[9] / 4.0)
    rhs_jacobian = ks_rhs_jacobian(
        new_state,
        params;
        relative_step=1.0e-4,
    )
    @test rhs_jacobian[1:4, 5:8] ≈ Matrix{Float64}(I, 4, 4) rtol=1.0e-9 atol=1.0e-12
    @test rhs_jacobian[5:8, 9] ≈ -0.5 .* u rtol=1.0e-8 atol=1.0e-8
    @test rhs_jacobian[10, 1:4] ≈ 2.0 .* u rtol=1.0e-6 atol=1.0e-6
    @test norm(rhs_jacobian[9, :]) > 0.0

    function old_convention_rhs(state)
        u = state[1:4]
        u_prime = state[5:8]
        h_old = state[9]
        radius = dot(u, u)
        rvec = ks_position(u)
        vvec = ks_velocity(u, u_prime)
        acceleration = ks_j2_acceleration_si(rvec, params)
        acceleration4 = [acceleration[1], acceleration[2], acceleration[3], 0.0]
        du_prime = -0.25 .* h_old .* u .+
            0.5 .* radius .* (transpose(ks_dynamics._ks_L(u)) * acceleration4)
        dh_old = -2.0 * radius * dot(vvec, acceleration)
        return vcat(u_prime, du_prime, dh_old, radius)
    end

    function old_convention_step(state, delta_s)
        k1 = old_convention_rhs(state)
        k2 = old_convention_rhs(state .+ 0.5 .* delta_s .* k1)
        k3 = old_convention_rhs(state .+ 0.5 .* delta_s .* k2)
        k4 = old_convention_rhs(state .+ delta_s .* k3)
        return state .+ (delta_s / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    end

    for _ in 1:100
        new_state = ks_rk4_step(new_state, params, 0.0, 1.0e-7)
        old_state = old_convention_step(old_state, 1.0e-7)
    end
    @test new_state[1:8] ≈ old_state[1:8] rtol=2.0e-14 atol=2.0e-12
    @test new_state[10] ≈ old_state[10] rtol=2.0e-14 atol=2.0e-12
    @test 2.0 * new_state[9] ≈ old_state[9] rtol=2.0e-14 atol=2.0e-8
    new_cartesian = ks_state_to_cartesian(new_state)
    old_cartesian = ks_state_to_cartesian(old_state)
    @test new_cartesian.position_ii_m ≈ old_cartesian.position_ii_m rtol=2.0e-14 atol=2.0e-8
    @test new_cartesian.velocity_ii_m ≈ old_cartesian.velocity_ii_m rtol=2.0e-14 atol=2.0e-10
end
