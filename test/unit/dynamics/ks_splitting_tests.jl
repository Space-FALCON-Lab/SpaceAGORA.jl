using Test
using LinearAlgebra
using SpaceAGORA
using SpaceAGORA.SimulationModel.DynamicsKS
using StaticArrays

# Earth, 300 x 3000 km altitude, i = 51.6 deg, starting at perigee.
function splitting_test_orbit(params)
    rp = params.Re + 300.0e3
    ra = params.Re + 3000.0e3
    a = 0.5 * (rp + ra)
    vp = sqrt(params.μ * (2.0 / rp - 1.0 / a))
    inclination = deg2rad(51.6)
    state = cartesian_to_ks_state([rp, 0.0, 0.0],
        [0.0, vp * cos(inclination), vp * sin(inclination)], params)
    return SVector{10,Float64}(state), 2pi * sqrt(a^3 / params.μ)
end

function splitting_steps(state, params, delta_s, count; order)
    kick_terms = ks_j2_kick_terms(state, params)
    for _ in 1:count
        state, kick_terms = ks_splitting_step(state, kick_terms, params, delta_s; order)
    end
    return state
end

@testset "KS splitting propagation" begin
    j2_params = KSPropagationParams(Re=6.378137e6, μ=3.986004418e14, J2=1.08262668e-3)
    kepler_params = KSPropagationParams(Re=j2_params.Re, μ=j2_params.μ)
    initial, period = splitting_test_orbit(j2_params)
    initial_position = ks_state_to_cartesian(initial).position_ii_m
    orbit_span = ks_splitting_step_size(initial, 1)

    @testset "a_p = 0 reproduces the Keplerian orbit" begin
        for order in (2, 4)
            final = splitting_steps(initial, kepler_params, orbit_span, 1; order)
            @test ks_state_to_cartesian(final).position_ii_m ≈ initial_position rtol=1.0e-13
            @test final[10] ≈ period rtol=1.0e-13
            @test final[9] == initial[9]
        end
    end

    # Both split flows preserve μ = h_KS |u|² + 2|u′|² exactly.
    @testset "KS energy identity under J2" begin
        states, _ = ks_splitting_propagate(initial, j2_params, 10period, 64; order=4)
        identity_error = maximum(states) do x
            u = SVector(x[1], x[2], x[3], x[4])
            u_prime = SVector(x[5], x[6], x[7], x[8])
            abs(x[9] * dot(u, u) + 2.0 * dot(u_prime, u_prime) - j2_params.μ) / j2_params.μ
        end
        @test identity_error < 1.0e-12
        @test states[end][9] != initial[9]
    end

    @testset "order $order: differences shrink by $(2^order) when N doubles" for order in (2, 4)
        positions = [ks_state_to_cartesian(splitting_steps(initial, j2_params,
            orbit_span / (64 * 2^k), 2 * 64 * 2^k; order)).position_ii_m for k in 0:2]
        ratio = norm(positions[1] - positions[2]) / norm(positions[2] - positions[3])
        @test ratio ≈ 2^order rtol=0.02
    end

    @testset "one or three perturbation evaluations per step" begin
        for (order, stages) in ((2, 1), (4, 3))
            states, work = ks_splitting_propagate(initial, j2_params, 2period, 16;
                order, save_history=false)
            @test length(states) == 1
            @test work.evaluations == 1 + stages * work.steps
            @test work.final_state[10] >= 2period
        end
    end

    @test_throws ArgumentError ks_splitting_step(initial,
        ks_j2_kick_terms(initial, j2_params), j2_params, orbit_span; order=3)
    @test_throws ArgumentError ks_oscillator_flow(initial[1:4], initial[5:8], -1.0, 1.0)
end
