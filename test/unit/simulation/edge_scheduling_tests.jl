module EdgeSchedulingTests

# Finite-horizon controllability with link-force inputs and the edge schedules
# that raise it, on small cases with known answers: the linear-model
# interface, the barycentre's uncontrollability under internal forces, the
# double-integrator Gramian in closed form, the metrics, the greedy and
# constant schedules, and minimum input energy.

using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const OMNI = OpticalTerminal(boresight=[1.0, 0.0, 0.0], field_of_regard_rad=π, frame=:inertial)
const ORIGIN = 7.0e6

sample(pos) = StateSample(SVector{3, Float64}(pos), SVector(0.0, 7.5e3, 0.0), 100.0)

# The feasible-edge graph of spacecraft at `positions`, each with `terminals`
# omnidirectional terminals, every pair in range and in sight.
function geometry(positions; terminals=1, duty_cycle=nothing)
    link = LaserBeamLink(power_w=1.0, magnification=1.0, range_m=1e7)
    model = EdgeFeasibilityModel(link, [fill(OMNI, terminals) for _ in positions]; duty_cycle)
    return feasible_edge_graph(model, [sample(p) for p in positions])
end

# Free-space double integrators, `[r_i; v_i]` per spacecraft.
function free_space(count)
    Ac = zeros(6 * count, 6 * count)
    for i in 1:count
        Ac[6 * (i - 1) .+ (1:3), 6 * (i - 1) .+ (4:6)] .= Matrix(1.0I, 3, 3)
    end
    return Ac
end

# A static two-coordinate window: A = I, and one column per edge, the same at every step.
const TOY_EDGES = [(1, 2), (1, 3), (2, 3)]
toy(steps; columns=([3.0, 0.0], [0.0, 2.0], [2.5, 0.0])) =
    EdgeLinearModel([EdgeLinearStep(Matrix(1.0I, 2, 2), collect(columns), TOY_EDGES) for _ in 1:steps])

@testset "Edge scheduling for controllability" begin
    @testset "Linear-model interface" begin
        step = EdgeLinearStep(Matrix(1.0I, 2, 2), [1.0 2.0; 3.0 4.0], [(2, 1), (1, 3)])
        @test step.edges == [(1, 2), (1, 3)]
        @test step.columns == [1:1, 2:2]
        blocks = EdgeLinearStep(Matrix(1.0I, 2, 2), [[1.0, 3.0], [2.0 5.0; 4.0 6.0]], [(1, 2), (2, 3)])
        @test blocks.B == [1.0 2.0 5.0; 3.0 4.0 6.0]
        @test blocks.columns == [1:1, 2:3]
        wide = EdgeLinearStep(Matrix(1.0I, 2, 2), [1.0 2.0 5.0; 3.0 4.0 6.0], [(1, 2), (2, 3)]; columns=[2:3, 1:1])
        @test wide.columns == [2:3, 1:1]
        @test_throws DimensionMismatch EdgeLinearStep(zeros(2, 3), zeros(2, 1), [(1, 2)])
        @test_throws DimensionMismatch EdgeLinearStep(zeros(2, 2), zeros(3, 1), [(1, 2)])
        @test_throws DimensionMismatch EdgeLinearStep(zeros(2, 2), zeros(2, 2), [(1, 2)])
        @test_throws ArgumentError EdgeLinearStep(zeros(2, 2), zeros(2, 2), [(1, 2), (2, 1)])
        @test_throws ArgumentError EdgeLinearStep(zeros(2, 2), zeros(2, 1), [(1, 1)])
        @test_throws ArgumentError EdgeLinearStep(zeros(2, 2), zeros(2, 2), [(1, 2), (2, 3)]; columns=[1:2, 2:2])
        @test_throws ArgumentError EdgeLinearStep(zeros(2, 2), zeros(2, 2), [(1, 2)]; columns=[3:3])
        @test_throws ArgumentError EdgeLinearStep(zeros(2, 2), [NaN; 0.0;;], [(1, 2)])
        @test_throws ArgumentError EdgeLinearStep(zeros(2, 2), [[1.0, 0.0]], [(1, 2), (2, 3)])

        @test toy(3).projection == Matrix(1.0I, 2, 2)
        @test toy(3).edges == TOY_EDGES
        @test EdgeLinearModel([zeros(2, 2), zeros(2, 2)], [ones(2, 1), ones(2, 1)], [(1, 2)], [1.0 0.0]).projection == [1.0 0.0]
        @test_throws ArgumentError EdgeLinearModel(EdgeLinearStep[], zeros(1, 2))
        @test_throws ArgumentError EdgeLinearModel(EdgeLinearStep[])
        @test_throws ArgumentError EdgeLinearModel([EdgeLinearStep(zeros(2, 2), ones(2, 1), [(1, 2)]),
            EdgeLinearStep(zeros(2, 2), ones(2, 1), [(1, 3)])])
        @test_throws DimensionMismatch EdgeLinearModel([EdgeLinearStep(zeros(2, 2), ones(2, 1), [(1, 2)]),
            EdgeLinearStep(zeros(3, 3), ones(3, 1), [(1, 2)])])
        @test_throws DimensionMismatch EdgeLinearModel([EdgeLinearStep(zeros(2, 2), ones(2, 1), [(1, 2)])], zeros(1, 3))
        @test_throws ArgumentError EdgeLinearModel([zeros(2, 2)], [ones(2, 1), ones(2, 1)], [(1, 2)])
    end

    @testset "Zero-order hold and barycentric coordinates" begin
        A, B = zoh_discretize([0.0 1.0; 0.0 0.0], [0.0, 1.0], 2.0)
        @test A ≈ [1.0 2.0; 0.0 1.0]
        @test B ≈ reshape([2.0, 2.0], 2, 1)
        A, B = zoh_discretize([-1.0 0.0; 0.0 0.0], [1.0 0.0; 0.0 1.0], 0.5)
        @test A ≈ [exp(-0.5) 0.0; 0.0 1.0]
        @test B ≈ [1 - exp(-0.5) 0.0; 0.0 0.5]
        @test_throws ArgumentError zoh_discretize(zeros(2, 2), zeros(2), 0.0)
        @test_throws DimensionMismatch zoh_discretize(zeros(2, 3), zeros(2), 1.0)
        @test_throws DimensionMismatch zoh_discretize(zeros(2, 2), zeros(3), 1.0)

        masses = [2.0, 3.0, 5.0]
        x = vcat([fill(Float64(i), 6) .* (1:6) for i in 1:3]...)
        chief = barycentric_transform(masses)
        @test chief.barycenter_rows == 1:6 && chief.relative_rows == 7:18
        @test chief.T * x ≈ vcat((2.0 * 1 + 3.0 * 2 + 5.0 * 3) / 10 .* (1:6), 1.0 .* (1:6), 2.0 .* (1:6))
        centred = barycentric_transform(masses; relative=:barycenter)
        cm = (2.0 * 1 + 3.0 * 2 + 5.0 * 3) / 10
        @test centred.T * x ≈ vcat(cm .* (1:6), (2 - cm) .* (1:6), (3 - cm) .* (1:6))
        @test rank(chief.T) == 18 && rank(centred.T) == 18
        @test barycentric_transform([1.0, 1.0]; block_size=1).T == [0.5 0.5; -1.0 1.0]
        @test_throws ArgumentError barycentric_transform([1.0, 0.0])
        @test_throws ArgumentError barycentric_transform(Float64[])
        @test_throws ArgumentError barycentric_transform([1.0]; relative=:jacobi)
        @test_throws ArgumentError barycentric_transform([1.0]; block_size=0)
    end

    @testset "Two-spacecraft chain: the barycentre is uncontrollable" begin
        m1, m2, h, N = 3.0, 5.0, 10.0, 6
        graph = geometry([(ORIGIN, 0.0, 0.0), (ORIGIN, 1e3, 0.0)])
        @test graph.direction[1] == SVector(0.0, 1.0, 0.0)
        Bc = edge_input_matrix(graph, [m1, m2])
        @test size(Bc) == (12, 1)
        @test Bc[4:6, 1] ≈ [0.0, -1 / m1, 0.0] && Bc[10:12, 1] ≈ [0.0, 1 / m2, 0.0]
        @test iszero(Bc[[1:3; 7:9], 1])
        # The column pushes the ends apart, as the laser link's own force does.
        a, b = sample((ORIGIN, 0.0, 0.0)), sample((ORIGIN, 1e3, 0.0))
        force_a, force_b = link_force(LaserMomentumLink(), a, b, LinkContext())
        @test normalize(Bc[4:6, 1]) ≈ normalize(force_a) && normalize(Bc[10:12, 1]) ≈ normalize(force_b)
        @test_throws ArgumentError edge_input_matrix(graph, [m1, m2]; edges=[(1, 3)])
        @test_throws ArgumentError edge_input_matrix(graph, [m1])
        @test_throws ArgumentError edge_input_matrix(graph, [m1, -m2])

        A, B = zoh_discretize(free_space(2), Bc, h)
        transform = barycentric_transform([m1, m2])
        barycenter = transform.T[transform.barycenter_rows, :]
        relative = transform.T[transform.relative_rows, :]
        # Proposition 1: internal forces leave the barycentre's rows of B empty.
        @test norm(barycenter * B) <= 1e-15 * norm(B)

        full = EdgeLinearModel(fill(EdgeLinearStep(A, B, [(1, 2)]), N))
        W = controllability_gramian(full)
        @test norm(barycenter * W * barycenter') <= 1e-14 * opnorm(W)
        @test abs(controllability_metric(SmallestEigenvalueMetric(), W)) <= 1e-12 * opnorm(W)
        # A barycentre target is out of reach.
        target = barycenter' * [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        @test !minimum_input_energy(full, nothing, target).reachable

        # In relative coordinates, only motion along the line of centres responds:
        # a double integrator of gain g = 1/m1 + 1/m2, in closed form.
        model = EdgeLinearModel(fill(EdgeLinearStep(A, B, [(1, 2)]), N), relative)
        W = controllability_gramian(model)
        g = 1 / m1 + 1 / m2
        half_steps = sum((i + 0.5)^2 for i in 0:N-1)
        closed = g^2 .* [h^4 * half_steps h^3 * N^2 / 2; h^3 * N^2 / 2 h^2 * N]
        @test W[[2, 5], [2, 5]] ≈ closed rtol = 1e-12
        @test norm(W[[1, 3, 4, 6], :]) <= 1e-12 * opnorm(W)
        @test abs(controllability_metric(SmallestEigenvalueMetric(), W)) <= 1e-12 * opnorm(W)

        along_line = EdgeLinearModel(fill(EdgeLinearStep(A, B, [(1, 2)]), N), relative[[2, 5], :])
        W2 = controllability_gramian(along_line)
        @test W2 ≈ closed rtol = 1e-12
        @test controllability_metric(SmallestEigenvalueMetric(), W2) ≈ eigmin(Symmetric(closed)) rtol = 1e-10

        # Minimum energy to open the separation by 2 m and stop.
        shift = [2.0, 0.0]
        energy = minimum_input_energy(along_line, nothing, shift)
        @test energy.energy ≈ dot(shift, closed \ shift) rtol = 1e-10
        @test energy.reachable && energy.residual <= 1e-12
        @test sort(collect(keys(energy.inputs[1]))) == [(1, 2)]
        reach = reachability_matrix(along_line).matrix
        @test reach * [only(energy.inputs[k][(1, 2)]) for k in 1:N] ≈ shift
        @test energy.energy ≈ sum(only(energy.inputs[k][(1, 2)])^2 for k in 1:N)
        # Across the line of centres is out of reach.
        sideways = minimum_input_energy(model, nothing, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        @test !sideways.reachable && sideways.residual ≈ 1.0
        # From a drifting start, the inputs correct the free response; none needed to keep it.
        drift = zeros(12)
        drift[11] = 0.01
        free = relative[[2, 5], :] * A^N * drift
        @test minimum_input_energy(along_line, nothing, free; initial_state=drift).energy <= 1e-20
        @test minimum_input_energy(along_line, nothing, shift + free; initial_state=drift).energy ≈ energy.energy rtol = 1e-10
        @test minimum_input_energy(along_line, nothing, [0.0, 0.0]).energy == 0.0
        @test_throws DimensionMismatch minimum_input_energy(along_line, nothing, [1.0])
        @test_throws DimensionMismatch minimum_input_energy(along_line, nothing, shift; initial_state=zeros(3))
        # No active link: nothing is reachable but the origin.
        silent = minimum_input_energy(along_line, [Tuple{Int, Int}[] for _ in 1:N], shift)
        @test !silent.reachable && silent.energy == 0.0
        # No projected coordinate: nothing to reach, and zero inputs.
        empty_output = EdgeLinearModel(fill(EdgeLinearStep(A, B, [(1, 2)]), N), zeros(0, 12))
        nothing_to_do = minimum_input_energy(empty_output, nothing, Float64[])
        @test nothing_to_do.reachable && nothing_to_do.energy == 0.0
        @test all(nothing_to_do.inputs[k][(1, 2)] == [0.0] for k in 1:N)
    end

    @testset "Three spacecraft in a plane: shape, not rotation" begin
        # Central link forces conserve the barycentre's motion and the angular
        # momentum about it. In free space, of the eight in-plane relative
        # coordinates, the six of the triangle's shape respond; the rigid
        # rotation (an angle and its rate) does not.
        positions = [(ORIGIN, 0.0, 0.0), (ORIGIN, 1e3, 0.0), (ORIGIN + 800.0, 500.0, 0.0)]
        graph = geometry(positions)
        masses = [3.0, 4.0, 5.0]
        A, B = zoh_discretize(free_space(3), edge_input_matrix(graph, masses), 5.0)
        transform = barycentric_transform(masses)
        relative = transform.T[transform.relative_rows, :]
        W = controllability_gramian(EdgeLinearModel(fill(EdgeLinearStep(A, B, graph.pairs), 4), relative))
        in_plane = [1, 2, 4, 5, 7, 8, 10, 11]
        @test rank(W[in_plane, in_plane]; rtol=1e-9) == 6
        @test norm(W[[3, 6, 9, 12], :]) <= 1e-12 * opnorm(W)
        full = controllability_gramian(EdgeLinearModel(fill(EdgeLinearStep(A, B, graph.pairs), 4)))
        barycenter = transform.T[transform.barycenter_rows, :]
        @test norm(barycenter * full * barycenter') <= 1e-14 * opnorm(full)
        centre = sum(masses .* [collect(p) for p in positions]) / sum(masses)
        spin_rate, spin_angle = zeros(18), zeros(18)
        for (i, p) in enumerate(positions)
            ρ = collect(p) - centre
            spin_rate[6 * (i - 1) .+ (4:5)] .= masses[i] .* [-ρ[2], ρ[1]]
            spin_angle[6 * (i - 1) .+ (1:2)] .= masses[i] .* [-ρ[2], ρ[1]]
        end
        @test dot(spin_rate, full, spin_rate) <= 1e-14 * opnorm(full) * dot(spin_rate, spin_rate)
        @test dot(spin_angle, full, spin_angle) <= 1e-14 * opnorm(full) * dot(spin_angle, spin_angle)
    end

    @testset "Reachability columns are the terminal gains" begin
        n, N = 4, 3
        As = [Matrix(1.0I, n, n) .+ 0.1 .* [sin(k + 2i + 3j) for i in 1:n, j in 1:n] for k in 1:N]
        Bs = [[cos(k * i + j) for i in 1:n, j in 1:2] for k in 1:N]
        P = [1.0 0.0 2.0 0.0; 0.0 1.0 0.0 -1.0]
        model = EdgeLinearModel(As, Bs, [(1, 2), (2, 3)], P)
        Φ(k) = reduce(*, reverse(As[k+1:N]); init=Matrix(1.0I, n, n))
        expected = reduce(hcat, [P * Φ(k) * Bs[k] for k in 1:N])
        reach = reachability_matrix(model)
        @test reach.matrix ≈ expected rtol = 1e-12
        @test reach.columns == [(1, (1, 2), 1), (1, (2, 3), 1), (2, (1, 2), 1), (2, (2, 3), 1), (3, (1, 2), 1), (3, (2, 3), 1)]
        chosen = reachability_matrix(model; active=[[(2, 3)], Tuple{Int, Int}[], [(1, 2), (3, 2)]])
        @test chosen.matrix ≈ expected[:, [2, 5, 6]] rtol = 1e-12
        @test controllability_gramian(model) ≈ expected * expected' rtol = 1e-12
        @test_throws ArgumentError reachability_matrix(model; active=[[(1, 2)]])
        @test_throws ArgumentError reachability_matrix(model; active=[[(1, 4)], [], []])

        # An edge with a vector input (three force components): three columns.
        vector_step = EdgeLinearStep(Matrix(1.0I, 3, 3), [Matrix(1.0I, 3, 3)], [(1, 2)])
        vector_model = EdgeLinearModel(fill(vector_step, 2))
        @test reachability_matrix(vector_model).columns == [(1, (1, 2), 1), (1, (1, 2), 2), (1, (1, 2), 3),
            (2, (1, 2), 1), (2, (1, 2), 2), (2, (1, 2), 3)]
        energy = minimum_input_energy(vector_model, nothing, [2.0, 0.0, 0.0])
        @test energy.inputs[1][(1, 2)] ≈ [1.0, 0.0, 0.0] && energy.inputs[2][(1, 2)] ≈ [1.0, 0.0, 0.0]
        @test energy.energy ≈ 2.0
    end

    @testset "Metrics" begin
        W = Diagonal([1.0, 2.0, 3.0])
        @test controllability_metric(LargestEigenvalueMetric(), W) ≈ 3.0
        @test controllability_metric(SmallestEigenvalueMetric(), W) ≈ 1.0
        @test controllability_metric(LogDetMetric(epsilon=0.5), W) ≈ log(1.5) + log(2.5) + log(3.5)
        @test controllability_metric(LogDetMetric(0.5), zeros(3, 3)) ≈ 3 * log(0.5)
        @test controllability_metric(DirectionEnergyMetric([0.0, 1.0, 0.0]), W) ≈ 2.0
        @test controllability_metric(DirectionEnergyMetric([1.0, 1.0, 0.0]), W) ≈ 3.0
        all_metrics = gramian_metrics(W; epsilon=0.5, direction=[0.0, 0.0, 2.0])
        @test all_metrics.largest_eigenvalue ≈ 3.0 && all_metrics.smallest_eigenvalue ≈ 1.0
        @test all_metrics.log_det ≈ log(1.5) + log(2.5) + log(3.5)
        @test all_metrics.direction_energy ≈ 12.0
        @test all_metrics.eigenvalues ≈ [1.0, 2.0, 3.0]
        @test abs.(all_metrics.dominant_mode) ≈ [0.0, 0.0, 1.0]
        @test isnan(gramian_metrics(W; epsilon=0.5).direction_energy)
        @test_throws ArgumentError LogDetMetric(epsilon=0.0)
        @test_throws ArgumentError DirectionEnergyMetric([0.0, 0.0])
        @test_throws DimensionMismatch controllability_metric(DirectionEnergyMetric([1.0, 0.0]), W)
    end

    @testset "Greedy and constant schedules" begin
        epsilon = 1e-3
        model = toy(2)
        # Greedy λ₁ piles both steps onto the strongest direction.
        largest = schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric()); max_active=1)
        @test largest.active == [[(1, 2)], [(1, 2)]]
        @test largest.gramian ≈ [18.0 0.0; 0.0 0.0]
        @test largest.terminals == [Dict{Tuple{Int, Int}, NTuple{2, Int}}(), Dict{Tuple{Int, Int}, NTuple{2, Int}}()]
        # Greedy log det and greedy λ_min balance the two directions.
        volume = schedule_edges(model, GreedyEdgeSchedule(LogDetMetric(epsilon=epsilon)); max_active=1)
        @test volume.active == [[(1, 2)], [(1, 3)]]
        @test volume.gramian ≈ [9.0 0.0; 0.0 4.0]
        hardest = GreedyEdgeSchedule(SmallestEigenvalueMetric(); tie_break=LogDetMetric(epsilon=epsilon))
        @test schedule_edges(model, hardest; max_active=1).active == [[(1, 2)], [(1, 3)]]
        @test schedule_edges(model, GreedyEdgeSchedule(SmallestEigenvalueMetric()); max_active=1).active ==
            [[(1, 2)], [(1, 3)]]
        @test schedule_edges(model, GreedyEdgeSchedule(DirectionEnergyMetric([0.0, 1.0])); max_active=1).active ==
            [[(1, 3)], [(1, 3)]]
        # On this small window the greedy volume schedule is the best of all schedules.
        logdet_of(set) = controllability_metric(LogDetMetric(epsilon), controllability_gramian(model; active=set))
        best = maximum(logdet_of([[a], [b]]) for a in TOY_EDGES, b in TOY_EDGES)
        @test logdet_of(volume.active) ≈ best
        @test controllability_metric(SmallestEigenvalueMetric(), volume.gramian) >
            controllability_metric(SmallestEigenvalueMetric(), largest.gramian)

        # Every admissible slot is filled.
        filled = schedule_edges(model, GreedyEdgeSchedule(LogDetMetric(epsilon=epsilon)); max_active=3)
        @test filled.active == [TOY_EDGES, TOY_EDGES]
        @test schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric()); max_active=0).active == [[], []]
        @test_throws ArgumentError schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric()); max_active=-1)

        # Feasibility by step.
        feasible = trues(3, 2)
        feasible[1, 1] = false
        masked = schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric()); max_active=1, feasible)
        @test masked.active == [[(2, 3)], [(1, 2)]]
        @test_throws DimensionMismatch schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric());
            max_active=1, feasible=trues(2, 2))

        # A duty cycle of one step in two per link, and the steps before the window.
        duty = DutyCycleLimit(on_steps=1, window_steps=2, scope=:edge)
        cycled = schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric()); max_active=1, duty_cycle=duty)
        @test cycled.active == [[(1, 2)], [(2, 3)]]
        after = schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric()); max_active=1, duty_cycle=duty,
            previous=[[(1, 2)]])
        @test after.active == [[(2, 3)], [(1, 2)]]

        # A constant set: given, or chosen by a metric.
        given = schedule_edges(model, ConstantEdgeSet(edges=[(3, 1)]); max_active=1)
        @test given.active == [[(1, 3)], [(1, 3)]]
        @test given.gramian ≈ [0.0 0.0; 0.0 8.0]
        @test schedule_edges(model, ConstantEdgeSet(edges=[(1, 3)]); max_active=1, feasible=[true true; false true; true true]).active ==
            [[], [(1, 3)]]
        chosen = schedule_edges(model, ConstantEdgeSet(metric=LogDetMetric(epsilon=epsilon)); max_active=1)
        @test chosen.active == [[(1, 2)], [(1, 2)]]
        pair = schedule_edges(model, ConstantEdgeSet(metric=LogDetMetric(epsilon=epsilon)); max_active=2)
        @test pair.active == [[(1, 2), (1, 3)], [(1, 2), (1, 3)]]
        @test schedule_edges(model, ConstantEdgeSet(metric=SmallestEigenvalueMetric(), tie_break=LogDetMetric(epsilon=epsilon));
            max_active=2).active == [[(1, 2), (1, 3)], [(1, 2), (1, 3)]]
        @test_throws ArgumentError ConstantEdgeSet()
        @test_throws ArgumentError ConstantEdgeSet(edges=[(1, 2)], metric=LargestEigenvalueMetric())
        @test_throws ArgumentError ConstantEdgeSet(edges=[(1, 2), (2, 1)])
        @test_throws ArgumentError ConstantEdgeSet(metric=1.0)
        @test_throws ArgumentError GreedyEdgeSchedule(LargestEigenvalueMetric(); tie_break=1.0)
        @test_throws ArgumentError schedule_edges(model, ConstantEdgeSet(edges=[(1, 2), (1, 3)]); max_active=1)
        @test_throws ArgumentError schedule_edges(model, ConstantEdgeSet(edges=[(1, 4)]); max_active=1)
    end

    @testset "Schedules on feasible-edge graphs: terminals and duty cycle" begin
        positions = [(ORIGIN, 0.0, 0.0), (ORIGIN, 1e3, 0.0), (ORIGIN + 800.0, 500.0, 0.0)]
        model = toy(2)
        # One terminal each: in a triangle, one link per step whatever max_active allows.
        single = [geometry(positions) for _ in 1:2]
        volume = GreedyEdgeSchedule(LogDetMetric(epsilon=1e-3))
        schedule = schedule_edges(model, volume; max_active=3, graphs=single)
        @test schedule.active == [[(1, 2)], [(1, 3)]]
        @test all(admissible_edge_set(single[k], schedule.active[k]; max_active=3) for k in 1:2)
        @test schedule.terminals == [Dict((1, 2) => (1, 1)), Dict((1, 3) => (1, 1))]
        # Two terminals each: every link at every step.
        double = [geometry(positions; terminals=2) for _ in 1:2]
        @test schedule_edges(model, volume; max_active=3, graphs=double).active == [TOY_EDGES, TOY_EDGES]
        # The graphs' duty cycle applies unless another is given.
        limited = [geometry(positions; terminals=2, duty_cycle=DutyCycleLimit(on_steps=1, window_steps=2)) for _ in 1:2]
        duty_schedule = schedule_edges(model, volume; max_active=3, graphs=limited)
        @test sum(length, duty_schedule.active) == 1
        @test sum(length, schedule_edges(model, volume; max_active=3, graphs=limited,
            duty_cycle=DutyCycleLimit(on_steps=4, window_steps=2)).active) == 6
        # A step where a link is out of range.
        far = [geometry(positions), geometry([(ORIGIN, 0.0, 0.0), (ORIGIN, 2e7, 0.0), (ORIGIN + 800.0, 500.0, 0.0)])]
        far_schedule = schedule_edges(model, GreedyEdgeSchedule(LargestEigenvalueMetric()); max_active=1, graphs=far)
        @test far_schedule.active == [[(1, 2)], [(1, 3)]]
        @test_throws ArgumentError schedule_edges(model, volume; max_active=1, graphs=single[1:1])
        @test_throws ArgumentError schedule_edges(model, volume; max_active=1, graphs=single, feasible=trues(3, 2))
        two_craft = [geometry(positions[1:2]) for _ in 1:2]
        @test_throws ArgumentError schedule_edges(model, volume; max_active=1, graphs=two_craft)
    end
end

end
