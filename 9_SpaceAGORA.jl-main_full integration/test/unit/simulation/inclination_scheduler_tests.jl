module InclinationSchedulerTests
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel
using SpaceAGORA.SimulationModel.SchedulingPolicies: inclination_rate, matching_workspace!

const MU = 3.98600436233e14

inclination(h) = atan(hypot(h[1], h[2]), h[3])

@testset "Physical inclination scheduler" begin
    @testset "GVE equivalence and finite differences" begin
        for e in (0.0, 0.2), inc in (1e-8, deg2rad(1), 0.8, 1.8, pi - 1e-8),
            u in (0.0, 0.3, 2.1, 4.5)
            a = 7.4281366e6
            nu = u - 0.4
            p = a * (1 - e^2)
            radius = p / (1 + e * cos(nu))
            hnorm = sqrt(MU * p)
            normal = SVector(0.0, -sin(inc), cos(inc))
            radial = SVector(cos(u), sin(u) * cos(inc), sin(u) * sin(inc))
            transverse = cross(normal, radial)
            target = (pos=radius * radial,
                vel=(MU / hnorm) * (e * sin(nu) * radial + (p / radius) * transverse),
                mass=227.0)
            force = SVector(0.001, -0.002, 0.003)
            expected = radius * cos(u) / hnorm * dot(force / target.mass, normal)
            actual = @inferred inclination_rate(target, force, MU)
            @test actual ≈ expected rtol=1e-10 atol=1e-20
            @test inclination_rate(target, 2force, MU) ≈ 2actual
            @test inclination_rate(merge(target, (mass=454.0,)), force, MU) ≈ actual / 2
            @test inclination_rate(target, zero(force), MU) == 0.0
            if 1e-6 < inc < pi - 1e-6
                h = cross(target.pos, target.vel)
                dh = cross(target.pos, force / target.mass)
                dt = 0.1
                numerical = (inclination(h + dt * dh) - inclination(h - dt * dh)) / (2dt)
                @test actual ≈ numerical rtol=1e-4 atol=3e-15
            end
        end
    end

    @testset "Equatorial boundary and invalid states" begin
        for direction in (-1.0, 1.0), fz in (-1.0, 0.0, 1.0)
            target = (pos=SVector(7e6, 0.0, 0.0), vel=SVector(0.0, direction * 7500, 0.0), mass=100.0)
            force = SVector(0.0, 0.0, fz)
            rate = inclination_rate(target, force, MU)
            @test rate ≈ direction * abs(fz) / (100 * 7500)
            h = cross(target.pos, target.vel)
            dt = 0.01
            @test rate ≈ (inclination(h + dt * cross(target.pos, force / target.mass)) -
                inclination(h)) / dt atol=5e-14
        end
        target = (pos=SVector(7e6, 0.0, 0.0), vel=SVector(0.0, 7500.0, 0.0), mass=100.0)
        for mass in (0.0, -1.0, Inf, NaN)
            @test_throws ArgumentError inclination_rate(merge(target, (; mass)), SVector(0.0, 0.0, 1.0), MU)
        end
        @test_throws ArgumentError inclination_rate(merge(target, (vel=zero(target.vel),)), zero(target.pos), MU)
        @test_throws ArgumentError inclination_rate(target, SVector(NaN, 0.0, 0.0), MU)
    end

    @testset "Target selection, actual forces, matching, and SMA regression" begin
        sc = [SpacecraftModel(n_terminal=1) for _ in 1:4]
        inc = deg2rad(1.0)
        velocity = SVector(0.0, 7500cos(inc), 7500sin(inc))
        pos = SVector(7e6, 0.0, 0.0)
        state = (sc=[
            (pos=pos, vel=velocity, mass=100.0),
            (pos=pos - SVector(0.0, 0.0, 1000.0), vel=velocity, mass=100.0),
            (pos=pos - SVector(0.0, 0.0, 2000.0), vel=velocity, mass=100.0),
            (pos=pos + SVector(0.0, 0.0, 1000.0), vel=velocity, mass=100.0)],)
        model = ConstellationInterLinkModel()
        for (j, power) in ((2, 10_000.0), (3, 20_000.0), (4, 30_000.0))
            register_allowed_link!(model, sc, (1, 1), (j, 1);
                parameters=InterLinkParameters(P=power, B=100.0, range=200e3))
        end
        update_availability!(model, sc, state)
        policy = SchedulingPolicyModel(:gve_inc; target_idx=1)
        score_candidates!(model, policy, state, MU; candidate_force=force_on_endpoint)
        scores = copy(matching_workspace!(model, policy).score)
        @test scores[1] > 0
        @test scores[2] ≈ 2scores[1]
        @test scores[3] ≈ -3scores[1]
        select_interlinks!(model, policy)
        @test model.graph.active == [false, true, false]
        for (kwargs, targets) in (((target_idx=2,), [2]),
            ((target_idx=nothing,), [1, 2, 3, 4]),
            ((target_idx=99, target_idxs=[1, 2]), [1, 2]))
            selected = SchedulingPolicyModel(:gve_inc; kwargs...)
            score_candidates!(model, selected, state, MU; candidate_force=force_on_endpoint)
            for j in 2:4
                expected = 0.0
                for (t, partner) in ((1, j), (j, 1))
                    t in targets || continue
                    force = force_on_endpoint(state.sc[t], state.sc[partner], model.graph.edge_properties[j - 1])
                    expected += inclination_rate(state.sc[t], force, MU)
                end
                @test matching_workspace!(model, selected).score[j - 1] ≈ expected atol=1e-20
            end
        end
        update_availability!(model, sc, state; is_active=[true, true, false, true])
        score_candidates!(model, policy, state, MU; candidate_force=force_on_endpoint)
        @test matching_workspace!(model, policy).score[2] == 0
        select_interlinks!(model, policy)
        @test model.graph.active == [true, false, false]
        negative = SchedulingPolicyModel(:gve_inc; target_idx=2)
        score_candidates!(model, negative, state, MU; candidate_force=force_on_endpoint)
        select_interlinks!(model, negative)
        @test !any(model.graph.active)
        for invalid in (SchedulingPolicyModel(:gve_inc; target_idx=5),
            SchedulingPolicyModel(:gve_inc; target_idxs=[1, 5]))
            @test_throws ArgumentError score_candidates!(model, invalid, state, MU; candidate_force=force_on_endpoint)
        end
        @test_throws ArgumentError score_candidates!(model, policy, state, NaN; candidate_force=force_on_endpoint)
        @test_throws ArgumentError SchedulingPolicyModel(:gve_inc; target_idxs=Int[])
        @test_throws ArgumentError SchedulingPolicyModel(:gve_inc; target_idxs=[1, 1])
        for kwargs in ((target_idx=1,), (target_idx=nothing,), (target_idxs=[1, 2],))
            sma = SchedulingPolicyModel(:gve_sma; kwargs...)
            score_candidates!(model, sma, state, MU; candidate_force=force_on_endpoint)
            targets = haskey(kwargs, :target_idxs) ? kwargs.target_idxs :
                kwargs.target_idx === nothing ? collect(1:4) : [kwargs.target_idx]
            for j in 2:4
                expected = 0.0
                if model.graph.available[j - 1]
                    for (t, partner) in ((1, j), (j, 1))
                        t in targets || continue
                        target = state.sc[t]
                        force = force_on_endpoint(target, state.sc[partner], model.graph.edge_properties[j - 1])
                        a = inv(2 / norm(target.pos) - dot(target.vel, target.vel) / MU)
                        expected += 2a^2 / MU * dot(target.vel, force) / target.mass
                    end
                end
                @test matching_workspace!(model, sma).score[j - 1] ≈ expected
            end
        end
    end
end
end
