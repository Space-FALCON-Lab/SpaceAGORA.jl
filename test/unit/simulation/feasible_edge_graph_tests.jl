module FeasibleEdgeGraphTests

# The feasible-edge graph of a constellation (`feasible_edge_graph`): each
# constraint of an edge (line of sight, range against power, field of regard at
# both ends, terminal count, duty cycle), at and across its boundary.

using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3
const SPEED = sqrt(PLANET.μ / RADIUS)
const OMNI = OpticalTerminal(boresight=[1.0, 0.0, 0.0], field_of_regard_rad=π, frame=:inertial)
const CONTEXT = LinkContext(time_s=0.0, planet=PLANET)

sample(pos; vel=SVector(0.0, SPEED, 0.0), q=nothing, index=0) =
    StateSample(SVector{3, Float64}(pos), SVector{3, Float64}(vel), 100.0; q_ib=q, spacecraft_index=index)

# Spacecraft on the y axis of a circular orbit's tangent, `offsets` metres from
# (RADIUS, 0, 0), each with an omnidirectional terminal unless given terminals.
along_track(offsets) = [sample((RADIUS, y, 0.0); index=i) for (i, y) in enumerate(offsets)]
omni(count) = [[OMNI] for _ in 1:count]

beam_link(; range_m=1e6, occultation=false, beam=nothing, power_w=1000.0) =
    LaserBeamLink(; power_w, magnification=100.0, range_m, occultation, beam)

@testset "Feasible-edge graph" begin
    @testset "Line of sight past the planet" begin
        # Opposite sides of the planet: the segment passes through it.
        samples = [sample((RADIUS, 0.0, 0.0)), sample((-RADIUS, 0.0, 0.0))]
        blocked = feasible_edge_graph(EdgeFeasibilityModel(beam_link(range_m=3e7, occultation=true), omni(2)), samples;
            context=CONTEXT)
        @test blocked.pairs == [(1, 2)]
        @test !blocked.line_of_sight[1] && blocked.in_range[1] && !blocked.feasible[1]
        open_sky = feasible_edge_graph(EdgeFeasibilityModel(beam_link(range_m=3e7), omni(2)), samples; context=CONTEXT)
        @test open_sky.line_of_sight[1] && open_sky.feasible[1]
        # Without a planet nothing blocks.
        @test feasible_edge_graph(EdgeFeasibilityModel(beam_link(range_m=3e7, occultation=true), omni(2)), samples).feasible[1]

        # A chord whose closest point is 50 km above the surface: clear below a
        # 50 km ray altitude, blocked at it (touching counts as blocked) and above.
        height = PLANET.Rp_e + 50e3
        chord = [sample((height, -1000e3, 0.0)), sample((height, 1000e3, 0.0))]
        clear_at(margin) = feasible_edge_graph(EdgeFeasibilityModel(
            beam_link(range_m=3e6, occultation=(; min_ray_altitude_m=margin)), omni(2)), chord; context=CONTEXT).line_of_sight[1]
        @test clear_at(0.0)
        @test clear_at(50e3 - 1e-3)
        @test !clear_at(50e3)
        @test !clear_at(50e3 + 1e-3)

        # Line of sight and range agree with the link's own availability in a run.
        link = beam_link(range_m=2e6, occultation=(; min_ray_altitude_m=80e3))
        model = EdgeFeasibilityModel(link, omni(2))
        for angle in range(0.0, 2π; length=13), separation in (1e5, 1.9e6, 2.1e6, 9e6)
            a = sample((RADIUS, 0.0, 0.0))
            b = sample(RADIUS .* (cos(angle), sin(angle), 0.0) .+ (0.0, 0.0, separation / 4))
            graph = feasible_edge_graph(model, [a, b]; context=CONTEXT)
            @test (graph.line_of_sight[1] && graph.in_range[1]) == link_available(link, nothing, nothing, a, b, CONTEXT)
        end
    end

    @testset "Range against power" begin
        # The link's own range, inclusive, with the separation exact.
        model = EdgeFeasibilityModel(beam_link(range_m=100e3), omni(2))
        in_range(offset) = feasible_edge_graph(model, along_track((0.0, offset))).in_range[1]
        @test in_range(100e3)
        @test !in_range(nextfloat(100e3))
        @test !in_range(0.0)
        graph = feasible_edge_graph(model, along_track((0.0, 60e3)))
        @test graph.range_m[1] == 60e3 && graph.max_range_m[1] == 100e3
        @test graph.direction[1] == SVector(0.0, 1.0, 0.0)

        # A function of power, capped by the link's range; each pair's power.
        rule(power) = 50e3 * sqrt(power / 1000.0)
        powered = EdgeFeasibilityModel(beam_link(range_m=1e6), omni(3); power_w=[4000.0, 1000.0, 16000.0], max_range=rule)
        graph = feasible_edge_graph(powered, along_track((0.0, 90e3, 180e3)))
        @test graph.pairs == [(1, 2), (1, 3), (2, 3)]
        @test graph.max_range_m ≈ [100e3, 50e3, 200e3]
        @test graph.in_range == [true, false, true]
        by_pair = EdgeFeasibilityModel(beam_link(range_m=150e3), omni(3); power_w=(i, j) -> 1000.0 * (i + j)^2,
            max_range=rule)
        @test feasible_edge_graph(by_pair, along_track((0.0, 90e3, 180e3))).max_range_m ≈ [150e3, 150e3, 150e3]
        @test power_limited_range_m(nothing, beam_link(), 10.0) == Inf

        # The least force at range: the link's own force law at ρ_max(P).
        beam = GaussianBeam(wavelength_m=1.064e-6, waist_radius_m=0.02, aperture_radius_m=0.05)
        link = beam_link(range_m=1e7, beam=beam, power_w=1000.0)
        limit = MinimumForceRange(min_force_n=1e-4)
        reach = power_limited_range_m(limit, link, 1000.0)
        @test 0 < reach < Inf
        force_at(distance) = norm(first(link_force(link, sample((0.0, 0.0, 0.0)), sample((distance, 0.0, 0.0)), CONTEXT)))
        @test force_at(reach * (1 - 1e-9)) >= 1e-4
        @test force_at(reach * (1 + 1e-9)) < 1e-4
        @test power_limited_range_m(limit, link, 4000.0) > reach
        @test power_limited_range_m(limit, link, 0.0) == 0.0
        @test power_limited_range_m(MinimumForceRange(min_force_n=0.0), link, 0.0) == Inf
        @test power_limited_range_m(MinimumForceRange(min_force_n=0.0), link, 10.0) == Inf
        # More force than the near field gives, or than the waist captures: no range.
        near_field = 100.0 * 1000.0 / 299_792_458.0
        @test power_limited_range_m(MinimumForceRange(min_force_n=near_field), link, 1000.0) == 0.0
        small_aperture = beam_link(range_m=1e7, beam=GaussianBeam(wavelength_m=1.064e-6, waist_radius_m=0.05,
            aperture_radius_m=0.01))
        waist_share = 1 - exp(-2 * 0.01^2 / 0.05^2)
        @test power_limited_range_m(MinimumForceRange(min_force_n=1.01 * waist_share * near_field), small_aperture,
            1000.0) == 0.0
        @test power_limited_range_m(MinimumForceRange(min_force_n=0.99 * waist_share * near_field), small_aperture,
            1000.0) > 0.0
        # Without a beam the force does not fall with range.
        @test power_limited_range_m(MinimumForceRange(min_force_n=near_field), beam_link(), 1000.0) == Inf
        @test power_limited_range_m(MinimumForceRange(min_force_n=nextfloat(near_field)), beam_link(), 1000.0) == 0.0
        @test_throws ArgumentError power_limited_range_m(limit, link, -1.0)
        @test_throws ArgumentError MinimumForceRange(min_force_n=-1.0)

        # In the graph, ρ_max(P) decides across its boundary.
        limited = EdgeFeasibilityModel(link, omni(2); max_range=limit)
        @test feasible_edge_graph(limited, along_track((0.0, reach * (1 - 1e-6)))).feasible[1]
        @test !feasible_edge_graph(limited, along_track((0.0, reach * (1 + 1e-6)))).feasible[1]
    end

    @testset "Field of regard at both ends" begin
        ahead = OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=π / 4, frame=:inertial)
        behind = OpticalTerminal(boresight=[0.0, -1.0, 0.0], field_of_regard_rad=π / 4, frame=:inertial)
        # Spacecraft 2 is ahead (+y) of spacecraft 1: 1 must look ahead, 2 behind.
        graph(terminals) = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), terminals), along_track((0.0, 50e3)))
        @test graph([[ahead], [behind]]).field_of_regard[1] == (true, true)
        @test graph([[ahead], [behind]]).feasible[1]
        @test graph([[ahead], [ahead]]).field_of_regard[1] == (true, false)
        @test !graph([[ahead], [ahead]]).feasible[1]
        @test graph([[behind], [behind]]).field_of_regard[1] == (false, true)
        @test graph([OpticalTerminal[], [behind]]).field_of_regard[1] == (false, true)
        @test graph([[behind, ahead], [behind]]).terminal_options[1] == ([2], [1])

        # The cone's edge: a sight just inside the half angle is seen, just outside is not.
        cone = π / 6
        seen(angle) = feasible_edge_graph(EdgeFeasibilityModel(beam_link(),
            [[OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=cone, frame=:inertial)], [OMNI]]),
            [sample((RADIUS, 0.0, 0.0)), sample((RADIUS + 50e3 * sin(angle), 50e3 * cos(angle), 0.0))]).field_of_regard[1][1]
        @test seen(cone - 1e-9)
        @test !seen(cone + 1e-9)
        @test seen(-(cone - 1e-9))
        # A field of regard of π sees every direction, the anti-boresight included.
        anti = OpticalTerminal(boresight=[0.0, -1.0, 0.0], field_of_regard_rad=π, frame=:inertial)
        @test graph([[anti], [anti]]).field_of_regard[1] == (true, true)

        # The local orbital frame: at (0, R, 0) moving along −x, along-track is −x.
        along = OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=0.1)
        @test along.frame === :rtn
        rtn_samples = [sample((0.0, RADIUS, 0.0); vel=(-SPEED, 0.0, 0.0)), sample((-50e3, RADIUS, 0.0); vel=(-SPEED, 0.0, 0.0))]
        rtn_graph = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[along], [OMNI]]), rtn_samples)
        @test rtn_graph.field_of_regard[1] == (true, true)
        radial = OpticalTerminal(boresight=[1.0, 0.0, 0.0], field_of_regard_rad=0.1)
        @test !feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[radial], [OMNI]]), rtn_samples).field_of_regard[1][1]
        @test_throws ArgumentError feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[along], [OMNI]]),
            [sample((0.0, RADIUS, 0.0); vel=(0.0, 1.0, 0.0)), sample((-50e3, RADIUS, 0.0))])

        # The body frame, through the attitude: identity as inertial; a half turn
        # about z points a +y boresight along −y.
        body_ahead = OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=π / 4, frame=:body)
        body_graph(q) = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[body_ahead], [OMNI]]),
            [sample((RADIUS, 0.0, 0.0); q=q), sample((RADIUS, 50e3, 0.0))])
        @test body_graph(SVector(0.0, 0.0, 0.0, 1.0)).field_of_regard[1][1]
        @test !body_graph(SVector(0.0, 0.0, 1.0, 0.0)).field_of_regard[1][1]
        @test_throws ArgumentError body_graph(nothing)

        # Under a general attitude, a :body terminal sees what the link's own
        # terminal pointing sees for mounted hardware with that boresight. A half
        # turn cannot tell a rotation from its transpose; this attitude can.
        attitude = normalize(SVector(0.3, -0.5, 0.2, 0.78))
        cone_half_angle = 0.6
        pointed = LaserBeamLink(power_w=1.0, magnification=1.0, range_m=1e6,
            pointing=TerminalPointing(field_of_regard_rad=cone_half_angle))
        function mounted(boresight, index)
            bus = Link(root=true, m=100.0)
            craft = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
                initial_condition=CartesianInitialCondition(SVector(RADIUS, 0.0, 0.0), SVector(0.0, SPEED, 0.0)),
                id=index, n_terminal=0)
            add_actuator!(craft, bus, LaserTerminal(boresight=boresight, max_power_w=1.0))
            return craft, only(mounted_link_terminals(craft)).terminal
        end
        own, own_terminal = mounted([0.2, 0.9, -0.3], 1)
        mine = OpticalTerminal(boresight=first(terminal_pose(own, own_terminal)), field_of_regard_rad=cone_half_angle,
            frame=:body)
        agreed = Bool[]
        for azimuth in range(0.0, 2π; length=9)[1:end-1], polar in range(0.1, π - 0.1; length=9)
            d = SVector(sin(polar) * cos(azimuth), sin(polar) * sin(azimuth), cos(polar))
            partner, partner_terminal = mounted(-d, 2)       # looks straight back, at the identity attitude
            a = StateSample(SVector(RADIUS, 0.0, 0.0), SVector(0.0, SPEED, 0.0), 100.0, attitude, SVector(0.0, 0.0, 0.0),
                own, 1)
            b = StateSample(SVector(RADIUS, 0.0, 0.0) + 50e3 * d, SVector(0.0, SPEED, 0.0), 100.0,
                SVector(0.0, 0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0), partner, 2)
            expected = link_available(pointed, own_terminal, partner_terminal, a, b, LinkContext())
            seen_here = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[mine], [OMNI]]), [a, b]).field_of_regard[1][1]
            @test seen_here == expected
            push!(agreed, expected)
        end
        @test any(agreed) && !all(agreed)

        # The local orbital frame is the basis of the RPO frame utilities.
        r, v = SVector(0.6, 0.8, 0.1) .* RADIUS, SVector(-4.0e3, 3.5e3, 4.5e3)
        basis = SpaceAGORA.SimulationModel.FrameTransforms.rtn_dcm_from_inertial(r, v)
        for sight in (SVector(1.0, 0.0, 0.0), normalize(SVector(0.3, -0.4, 0.8)))
            @test SpaceAGORA.SimulationModel.InterLinkModels._rtn_sight(sample(r; vel=v), sight) ≈ basis' * sight
        end

        @test_throws ArgumentError OpticalTerminal(boresight=[0.0, 0.0, 0.0], field_of_regard_rad=0.1)
        @test_throws ArgumentError OpticalTerminal(boresight=[1.0, 0.0, 0.0], field_of_regard_rad=0.0)
        @test_throws ArgumentError OpticalTerminal(boresight=[1.0, 0.0, 0.0], field_of_regard_rad=3.2)
        @test_throws ArgumentError OpticalTerminal(boresight=[1.0, 0.0, 0.0], field_of_regard_rad=0.1, frame=:lvlh)
        @test OpticalTerminal(boresight=[0.0, 2.0, 0.0], field_of_regard_rad=0.1).boresight == SVector(0.0, 1.0, 0.0)
    end

    @testset "Terminal count" begin
        # Spacecraft 1 between 2 (ahead) and 3 (behind), each pair in range.
        samples = along_track((0.0, 50e3, -50e3))
        one_terminal = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), omni(3); pairs=[(1, 2), (3, 1)]), samples)
        @test one_terminal.pairs == [(1, 2), (1, 3)]
        @test feasible_pairs(one_terminal) == [(1, 2), (1, 3)]
        @test one_terminal.terminal_count == [1, 1, 1]
        @test admissible_edge_set(one_terminal, [(1, 2)])
        @test admissible_edge_set(one_terminal, [(3, 1)])
        @test !admissible_edge_set(one_terminal, [(1, 2), (1, 3)])
        @test terminal_assignment(one_terminal, [(1, 2), (1, 3)]) === nothing
        two_terminals = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[OMNI, OMNI], [OMNI], [OMNI]];
            pairs=[(1, 2), (1, 3)]), samples)
        @test admissible_edge_set(two_terminals, [(1, 2), (1, 3)])
        @test !admissible_edge_set(two_terminals, [(1, 2), (1, 3)]; max_active=1)
        @test admissible_edge_set(two_terminals, Tuple{Int, Int}[]; max_active=0)
        assignment = terminal_assignment(two_terminals, [(1, 2), (1, 3)])
        @test sort([assignment[(1, 2)][1], assignment[(1, 3)][1]]) == [1, 2]
        @test assignment[(1, 2)][2] == 1 && assignment[(1, 3)][2] == 1

        # The first terminal sees both partners, the second only spacecraft 2:
        # the assignment moves (1, 2) to the second terminal to free the first.
        ahead = OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=π / 4, frame=:inertial)
        shared = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[OMNI, ahead], [OMNI], [OMNI]];
            pairs=[(1, 2), (1, 3)]), samples)
        @test shared.terminal_options[1] == ([1, 2], [1])
        @test shared.terminal_options[2] == ([1], [1])
        assignment = terminal_assignment(shared, [(1, 2), (1, 3)])
        @test assignment[(1, 2)] == (2, 1) && assignment[(1, 3)] == (1, 1)
        both_ahead = feasible_edge_graph(EdgeFeasibilityModel(beam_link(), [[ahead, ahead], [OMNI], [OMNI]];
            pairs=[(1, 2), (1, 3)]), samples)
        @test !both_ahead.feasible[2]
        @test !admissible_edge_set(both_ahead, [(1, 3)])
        @test_throws ArgumentError admissible_edge_set(one_terminal, [(2, 3)])
        @test_throws ArgumentError terminal_assignment(one_terminal, [(2, 3)])
    end

    @testset "Duty cycle" begin
        samples = along_track((0.0, 50e3, -50e3))
        terminals = [[OMNI, OMNI], [OMNI, OMNI], [OMNI, OMNI]]
        limited(limit) = EdgeFeasibilityModel(beam_link(), terminals; duty_cycle=limit)
        # Two link-steps per spacecraft in any three steps.
        spacecraft_limit = limited(DutyCycleLimit(on_steps=2, window_steps=3))
        graph = feasible_edge_graph(spacecraft_limit, samples; previous=[[(1, 2)], [(1, 3)]])
        @test graph.pairs == [(1, 2), (1, 3), (2, 3)]
        @test graph.duty_cycle == [false, false, true]
        @test graph.feasible == [false, false, true]
        @test graph.line_of_sight == [true, true, true]
        # The window's edge: a link-step three steps back has left a window of three.
        one_per_window = limited(DutyCycleLimit(on_steps=1, window_steps=3))
        @test feasible_edge_graph(one_per_window, samples; previous=[[(1, 2)], [], []]).duty_cycle[1]
        @test !feasible_edge_graph(one_per_window, samples; previous=[[(1, 2)], []]).duty_cycle[1]
        @test !feasible_edge_graph(one_per_window, samples; previous=[[(2, 1)]]).duty_cycle[1]
        @test all(feasible_edge_graph(one_per_window, samples).duty_cycle)
        # Per link: another link of the same spacecraft does not count.
        edge_limit = limited(DutyCycleLimit(on_steps=1, window_steps=2, scope=:edge))
        @test feasible_edge_graph(edge_limit, samples; previous=[[(1, 2)]]).duty_cycle == [false, true, true]
        # Links active together share a spacecraft's budget.
        tight = feasible_edge_graph(limited(DutyCycleLimit(on_steps=1, window_steps=1)), samples)
        @test all(tight.feasible)
        @test admissible_edge_set(tight, [(1, 2)])
        @test !admissible_edge_set(tight, [(1, 2), (1, 3)])
        @test admissible_edge_set(feasible_edge_graph(limited(DutyCycleLimit(on_steps=2, window_steps=1)), samples),
            [(1, 2), (1, 3)])
        @test !admissible_edge_set(feasible_edge_graph(spacecraft_limit, samples), [(1, 2)]; previous=[[(1, 2)], [(1, 3)]])
        @test admissible_edge_set(feasible_edge_graph(spacecraft_limit, samples), [(2, 3)]; previous=[[(1, 2)], [(1, 3)]])
        @test tight.duty_limit == DutyCycleLimit(on_steps=1, window_steps=1)
        @test DutyCycleLimit(on_steps=0, window_steps=1).on_steps == 0
        @test !feasible_edge_graph(limited(DutyCycleLimit(on_steps=0, window_steps=4)), samples).feasible[1]
        @test_throws ArgumentError DutyCycleLimit(on_steps=-1, window_steps=2)
        @test_throws ArgumentError DutyCycleLimit(on_steps=1, window_steps=0)
        @test_throws ArgumentError DutyCycleLimit(on_steps=1, window_steps=2, scope=:terminal)
    end

    @testset "Model, states and sample times" begin
        @test EdgeFeasibilityModel(beam_link(), omni(3)).pairs == [(1, 2), (1, 3), (2, 3)]
        @test EdgeFeasibilityModel(beam_link(power_w=7.0), omni(2)).power_w == [7.0]
        @test EdgeFeasibilityModel(beam_link(), omni(2); power_w=3.0).power_w == [3.0]
        momentum = EdgeFeasibilityModel(LaserMomentumLink(power_w=5.0, magnification=2.0, range_m=10e3), omni(2))
        @test momentum.link isa LaserBeamLink && link_range_m(momentum.link) == 10e3 && !momentum.link.occultation
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(2); pairs=[(1, 3)])
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(2); pairs=[(2, 2)])
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(2); pairs=[(1, 2), (2, 1)])
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(3); power_w=[1.0])
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(2); power_w=-1.0)
        # A number of metres is a cap that does not depend on the power; other values are refused.
        capped = EdgeFeasibilityModel(beam_link(range_m=1e6), omni(2); max_range=70e3)
        @test feasible_edge_graph(capped, along_track((0.0, 60e3))).max_range_m == [70e3]
        @test !feasible_edge_graph(capped, along_track((0.0, 80e3))).feasible[1]
        @test power_limited_range_m(70e3, beam_link(), 1.0) == 70e3
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(2); max_range=-1.0)
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(2); max_range=NaN)
        @test_throws ArgumentError EdgeFeasibilityModel(beam_link(), omni(2); max_range="far")
        @test_throws ArgumentError EdgeFeasibilityModel(LaserBeamLink(power_w=1.0, magnification=1.0, range_m=1.0,
            pointing=TerminalPointing(field_of_regard_rad=0.5)), omni(2))
        @test_throws ArgumentError EdgeFeasibilityModel(PowerBeamLink(power_w=1.0, range_m=1.0,
            wavelength_m=1e-6, waist_radius_m=0.1), omni(2))
        @test_throws ArgumentError feasible_edge_graph(EdgeFeasibilityModel(beam_link(), omni(3)), along_track((0.0, 1e3)))

        # A run state (spacecraft `i` at `state.sc[i]`), and a propagated
        # constellation sampled at several times.
        model = EdgeFeasibilityModel(beam_link(range_m=100e3), omni(2))
        state_at(t) = (sc=[(pos=SVector(RADIUS, 0.0, 0.0), vel=SVector(0.0, SPEED, 0.0), mass=100.0),
            (pos=SVector(RADIUS, 1e3 * t, 0.0), vel=SVector(0.0, SPEED, 0.0), mass=100.0)],)
        @test feasible_edge_graph(model, state_at(50.0); context=LinkContext(time_s=50.0)).feasible == [true]
        graphs = feasible_edge_graphs(model, state_at, [0.0, 50.0, 150.0]; planet=PLANET)
        @test [g.time_s for g in graphs] == [0.0, 50.0, 150.0]
        @test [only(g.feasible) for g in graphs] == [false, true, false]
        @test [only(g.range_m) for g in graphs] == [0.0, 50e3, 150e3]
    end
end

end
