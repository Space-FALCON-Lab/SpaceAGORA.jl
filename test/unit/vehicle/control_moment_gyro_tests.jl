module ControlMomentGyroTests
using Test, StaticArrays, LinearAlgebra, SpaceAGORA
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel

mutable struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

function far_spacecraft(; ω0=SVector(0.0, 0.0, 0.0), pos=SVector(4.0e9, 0.0, 0.0))
    root = SM.Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(pos, SVector(0.0, 0.0, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=ω0)
    return SM.SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia, initial_condition=ic, id=1), root
end

function cmg_run(sc, levels; mission_time=20.0, solver=nothing, orientation_sim=true, EI_km=1.0)
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=orientation_sim,
        keplerian=false, EI_km=EI_km, verbose=false, results=false,
        solver_config=solver === nothing ? nothing : SM.SolverConfig(solver_mode=solver))
    effectors = isempty(levels) ? () : (HoldLevels(collect(Float64, levels)),)
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=effectors, control_rates=fill(1.0, length(effectors))))
    return run_simulation(args; isolate_state=false, return_results=true, return_solution=true)
end

# The vehicle's angular momentum, body plus rotors (CMGs and wheels), in the inertial frame.
function inertial_momentum(sc, x)
    q = SVector{4}(x.q)
    assembly = sc.root.rw_assembly
    n = length(assembly.cmgs)
    rotor(j, c) = assembly.cmg_inner_of[j] == 0 ? SM.cmg_momentum(c, x.cmg_angles[j]) :
        SM.cmg_momentum(assembly.cmg_inner[assembly.cmg_inner_of[j]], x.cmg_angles[j], x.cmg_angles[n + assembly.cmg_inner_of[j]])
    h = sum(rotor(j, c) for (j, c) in enumerate(assembly.cmgs))
    assembly.n_wheels > 0 && (h += SVector{3}(assembly.J_rw * x.h_wheels))
    return SM.rot(q)' * (sc.inertia_tensor * SVector{3}(x.ω) + h)
end

gyro(; kw...) = ControlMomentGyro(; rotor_momentum_nms=2.0, gimbal_axis=SVector(0.0, 0.0, 1.0), spin_axis=SVector(1.0, 0.0, 0.0), kw...)

@testset "Control moment gyros" begin
    @testset "the record and mounting" begin
        @test_throws ArgumentError gyro(rotor_momentum_nms=0.0)
        @test_throws ArgumentError ControlMomentGyro(rotor_momentum_nms=1.0, gimbal_axis=SVector(0.0, 0.0, 1.0), spin_axis=SVector(1.0, 0.0, 1.0))
        @test_throws ArgumentError gyro(max_gimbal_rate_rad_s=0.0)
        @test_throws ArgumentError gyro(min_gimbal_angle_rad=0.5, max_gimbal_angle_rad=0.4)
        @test_throws ArgumentError gyro(max_gimbal_angle_rad=NaN)
        @test_throws ArgumentError gyro(gimbal_angle_rad=0.2, max_gimbal_angle_rad=0.1)   # starts beyond its stop
        positional = ControlMomentGyro(2.0, SVector(0.0, 0.0, 1.0), SVector(1.0, 0.0, 0.0), Inf, 0.0)
        @test positional.min_gimbal_angle_rad == -Inf && positional.max_gimbal_angle_rad == Inf
        sc, root = far_spacecraft()
        panel = SM.Link(m=1.0, r=MVector(0.0, 1.0, 0.0))
        push!(sc.links, panel)
        @test_throws ArgumentError add_actuator!(sc, panel, gyro())                      # the root link only
        c = add_actuator!(sc, root, gyro(max_gimbal_rate_rad_s=0.5))
        @test_throws ArgumentError add_actuator!(sc, root, c)                            # mounted twice
        @test SM.mounted_cmg_count(sc) == 1 && [e.actuator for e in mounted_actuators(sc)] == [c]
        @test actuator_limits(c) == (-0.5, 0.5) && actuator_saturate(c, 2.0) == 0.5
        @test all(iszero, allocation_matrix(sc))                                         # its torque depends on its state
    end

    @testset "the rotor's momentum and torque axis" begin
        c = gyro()
        @test cmg_momentum(c, 0.0) ≈ SVector(2.0, 0.0, 0.0)
        @test cmg_momentum(c, π / 2) ≈ SVector(0.0, 2.0, 0.0) atol=1e-15
        @test cmg_torque_axis(c, 0.0) ≈ SVector(0.0, 2.0, 0.0)                           # g × h
        stop = gyro(min_gimbal_angle_rad=-0.2, max_gimbal_angle_rad=0.3)
        @test cmg_angle_saturate(stop, 0.1, 0.3) == 0.0 && cmg_angle_saturate(stop, -0.1, 0.3) == -0.1
        @test cmg_angle_saturate(stop, -0.1, -0.25) == 0.0 && cmg_angle_saturate(stop, 0.1, 0.0) == 0.1
        @test cmg_angle_saturate(c, 5.0, 1.0e6) == 5.0                                    # no stops by default
    end

    @testset "steering: exact away from singularities, bounded at one, within rate limits" begin
        sc, root = far_spacecraft()
        for (g, s) in ((SVector(0.0, 0.0, 1.0), SVector(1.0, 0.0, 0.0)), (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)),
                       (SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)))
            add_actuator!(sc, root, ControlMomentGyro(rotor_momentum_nms=2.0, gimbal_axis=g, spin_axis=s))
        end
        τ = SVector(0.01, -0.02, 0.03)
        rates = cmg_steering(sc, zeros(3), τ)
        made = -sum(rates[j] * cmg_torque_axis(c, 0.0) for (j, c) in enumerate(sc.root.rw_assembly.cmgs))
        @test made ≈ τ rtol=1e-5
        pair, proot = far_spacecraft()
        for _ in 1:2
            add_actuator!(pair, proot, gyro())                                           # parallel: torque only along y
        end
        rates = cmg_steering(pair, zeros(2), SVector(1.0, 0.0, 0.0))
        @test all(isfinite, rates) && maximum(abs, rates) < 1e-12                        # an unreachable request gives no motion
        rates = cmg_steering(pair, zeros(2), SVector(0.0, -0.4, 0.0))
        @test all(isfinite, rates) && rates[1] ≈ rates[2] && rates[1] > 0.0
        limited, lroot = far_spacecraft()
        add_actuator!(limited, lroot, gyro(max_gimbal_rate_rad_s=0.1))
        add_actuator!(limited, lroot, ControlMomentGyro(rotor_momentum_nms=2.0, gimbal_axis=SVector(1.0, 0.0, 0.0), spin_axis=SVector(0.0, 1.0, 0.0), max_gimbal_rate_rad_s=0.1))
        free = cmg_steering(limited, zeros(2), SVector(0.0, -1.0, -0.5); λ0=0.0)
        @test maximum(abs, free) ≈ 0.1 && free[1] / free[2] ≈ 2.0                        # scaled together
        # A CMG at its stop is left out when the torque would drive it further.
        stops, sroot = far_spacecraft()
        add_actuator!(stops, sroot, gyro(gimbal_angle_rad=0.3, max_gimbal_angle_rad=0.3))
        add_actuator!(stops, sroot, ControlMomentGyro(rotor_momentum_nms=2.0, gimbal_axis=SVector(1.0, 0.0, 0.0), spin_axis=SVector(0.0, 1.0, 0.0)))
        a1 = cmg_torque_axis(stops.root.rw_assembly.cmgs[1], 0.3)
        outward = cmg_steering(stops, [0.3, 0.0], -0.1 * a1 / norm(a1) + SVector(0.0, 0.0, 0.05); λ0=0.0)
        @test outward[1] == 0.0 && outward[2] ≈ -0.025                                   # the second carries what it can
        inward = cmg_steering(stops, [0.3, 0.0], 0.1 * a1 / norm(a1); λ0=0.0)
        @test inward[1] < 0.0                                                            # away from the stop
    end

    @testset "a run integrates the gimbal angle and conserves the angular momentum" begin
        sc, root = far_spacecraft(ω0=SVector(0.001, 0.0, 0.0))
        c = add_actuator!(sc, root, gyro(gimbal_angle_rad=0.1))
        res = cmg_run(sc, [0.02]; mission_time=20.0)
        x0, x = res.solution.u[1].sc[1], res.solution.u[end].sc[1]
        @test x.cmg_angles[1] ≈ 0.1 + 0.02 * 20.0 rtol=1e-9
        @test res.table.sc1_cmg_angle_1[end] ≈ x.cmg_angles[1] && res.table.sc1_actuator_level_1[end] == 0.02
        @test inertial_momentum(sc, x) ≈ inertial_momentum(sc, x0) rtol=1e-8
        # The body took up the rotor's momentum change: ω changed by about -(Δh)/I.
        Δh = cmg_momentum(c, x.cmg_angles[1]) - cmg_momentum(c, 0.1)
        @test norm(SVector{3}(x.ω) - SVector(0.001, 0.0, 0.0)) > 0.5 * norm(Δh) / sc.inertia_tensor[1, 1]
    end

    @testset "a gimbal stop holds the angle, with the torque consistent" begin
        sc, root = far_spacecraft(ω0=SVector(0.001, 0.0, 0.0))
        add_actuator!(sc, root, gyro(gimbal_angle_rad=0.1, max_gimbal_angle_rad=0.3))
        res = cmg_run(sc, [0.02]; mission_time=20.0)                                     # reaches the stop at 10 s
        x0, x = res.solution.u[1].sc[1], res.solution.u[end].sc[1]
        @test x.cmg_angles[1] == 0.3                                                    # on the stop, not past it
        @test inertial_momentum(sc, x) ≈ inertial_momentum(sc, x0) rtol=1e-8
        @test res.table.sc1_actuator_level_1[1] == 0.02 && res.table.sc1_actuator_level_1[end] == 0.0
    end

    @testset "the lower stop holds, and a command back into the range leaves a stop" begin
        sc, root = far_spacecraft(ω0=SVector(0.001, 0.0, 0.0))
        add_actuator!(sc, root, gyro(gimbal_angle_rad=0.1, min_gimbal_angle_rad=-0.05))
        res = cmg_run(sc, [-0.02]; mission_time=20.0)                                    # reaches the stop at 7.5 s
        x0, x = res.solution.u[1].sc[1], res.solution.u[end].sc[1]
        @test x.cmg_angles[1] == -0.05 && all(a -> a >= -0.05, (u.sc[1].cmg_angles[1] for u in res.solution.u))
        @test inertial_momentum(sc, x) ≈ inertial_momentum(sc, x0) rtol=1e-8
        back, broot = far_spacecraft()
        add_actuator!(back, broot, gyro(gimbal_angle_rad=0.3, max_gimbal_angle_rad=0.3))   # starts on its stop
        @test cmg_run(back, [-0.02]; mission_time=5.0).solution.u[end].sc[1].cmg_angles[1] ≈ 0.3 - 0.02 * 5.0 rtol=1e-9
    end

    @testset "the split solver modes stop the gimbal exactly once" begin
        # The stop event sees the full state in split_imex; in multirate the
        # gimbal moves only in the fast subsolves, each starting from the
        # clamped state once the stop is reached.
        for solver in (:split_imex, :multirate)
            sc, root = far_spacecraft(ω0=SVector(0.001, 0.0, 0.0), pos=SVector(7.0e6, 0.0, 0.0))
            add_actuator!(sc, root, gyro(gimbal_angle_rad=0.1, max_gimbal_angle_rad=0.31))
            res = cmg_run(sc, [0.02]; mission_time=20.0, solver=solver, EI_km=2000.0)   # the stop at 10.5 s, between ticks
            x0, x = res.solution.u[1].sc[1], res.solution.u[end].sc[1]
            @testset let solver = solver
                @test x.cmg_angles[1] == 0.31
                @test all(a -> a <= 0.31, (u.sc[1].cmg_angles[1] for u in res.solution.u))
                @test inertial_momentum(sc, x) ≈ inertial_momentum(sc, x0) rtol=1e-6
            end
        end
    end

    @testset "a held gimbal still couples through the gyroscopic term" begin
        sc, root = far_spacecraft(ω0=SVector(0.0, 0.0, 0.01))
        add_actuator!(sc, root, gyro())
        res = cmg_run(sc, Float64[]; mission_time=50.0)
        x0, x = res.solution.u[1].sc[1], res.solution.u[end].sc[1]
        @test x.cmg_angles[1] == 0.0
        @test inertial_momentum(sc, x) ≈ inertial_momentum(sc, x0) rtol=1e-8
        @test norm(x.ω[1:2]) > 1e-4                                                     # the rotor's momentum tips the spin axis
    end

    @testset "a spacecraft reused without controllers keeps its CMG and wheel still" begin
        # The first run commands the CMG and a mounted wheel; the second, on the
        # same records with no control effector, must not reuse those rates.
        sc, root = far_spacecraft()
        add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0)))
        add_actuator!(sc, root, gyro(gimbal_angle_rad=0.1))
        cmg_run(sc, [0.01, 0.02]; mission_time=5.0)
        x = cmg_run(sc, Float64[]; mission_time=5.0).solution.u[end].sc[1]
        @test x.cmg_angles[1] == 0.1 && x.h_wheels[1] == 0.0 && norm(SVector{3}(x.ω)) == 0.0
    end

    @testset "a failed CMG holds its gimbal" begin
        sc, root = far_spacecraft()
        c = add_actuator!(sc, root, gyro(gimbal_angle_rad=0.3))
        add_fault!(sc, ActuatorFault(c; kind=:failed, start_s=5.0))
        res = cmg_run(sc, [0.05]; mission_time=10.0)
        @test res.solution.u[end].sc[1].cmg_angles[1] ≈ 0.3 + 0.05 * 5.0 rtol=1e-9
    end

    @testset "the split solver modes integrate each CMG and mounted-wheel rate once" begin
        # Inside the entry interface, so the implicit component of split_imex
        # runs too; multirate's slow component runs at every step.
        for solver in (:split_imex, :multirate)
            sc, root = far_spacecraft(ω0=SVector(0.001, 0.0, 0.0), pos=SVector(7.0e6, 0.0, 0.0))
            add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0)))
            add_actuator!(sc, root, gyro(gimbal_angle_rad=0.1))
            res = cmg_run(sc, [0.01, 0.02]; mission_time=20.0, solver=solver, EI_km=2000.0)
            x0, x = res.solution.u[1].sc[1], res.solution.u[end].sc[1]
            @test x.cmg_angles[1] ≈ 0.1 + 0.02 * 20.0 rtol=1e-6
            @test x.h_wheels[1] ≈ -0.01 * 20.0 rtol=1e-6                                # the wheel takes up minus the body torque
            @test inertial_momentum(sc, x) ≈ inertial_momentum(sc, x0) rtol=1e-6
        end
    end

    @testset "runs without attitude refuse CMGs" begin
        sc, root = far_spacecraft()
        add_actuator!(sc, root, gyro())
        @test_throws ArgumentError cmg_run(sc, Float64[]; mission_time=1.0, orientation_sim=false)
    end

    @testset "double-gimbal CMGs: the inner gimbal record and its mounting" begin
        sc, root = far_spacecraft()
        c = add_actuator!(sc, root, gyro())                                       # g_o = z, s = x
        @test_throws ArgumentError CMGInnerGimbal(c; axis=SVector(1.0, 0.0, 0.0))     # along the spin axis
        @test_throws ArgumentError CMGInnerGimbal(c; axis=SVector(0.0, 0.0, 1.0))     # along the outer gimbal axis
        @test_throws ArgumentError CMGInnerGimbal(c; axis=SVector(0.0, 1.0, 0.0), gimbal_angle_rad=0.3, max_gimbal_angle_rad=0.2)
        loose = CMGInnerGimbal(gyro(); axis=SVector(0.0, 1.0, 0.0))
        @test_throws ArgumentError add_actuator!(sc, root, loose)                 # its CMG is not mounted
        panel = SM.Link(m=1.0, r=MVector(0.0, 1.0, 0.0)); push!(sc.links, panel)
        inner = CMGInnerGimbal(c; axis=SVector(0.0, 1.0, 0.0), max_gimbal_rate_rad_s=0.4)
        @test_throws ArgumentError add_actuator!(sc, panel, inner)                # the root link only
        @test add_actuator!(sc, root, inner) === inner
        @test_throws ArgumentError add_actuator!(sc, root, CMGInnerGimbal(c; axis=SVector(0.0, 1.0, 0.0)))   # one per CMG
        other = add_actuator!(sc, root, gyro())                                   # mounted after the inner gimbal
        @test [e.actuator for e in mounted_actuators(sc)] == [c, other, inner]   # inner gimbals after every CMG
        @test sc.root.rw_assembly.cmg_inner_of == [1, 0] && length(sc.root.rw_assembly.cmg_rates) == 3
        @test actuator_limits(inner) == (-0.4, 0.4)
    end

    @testset "double-gimbal CMGs: the rotor's momentum and its two torque axes" begin
        c = gyro()                                                                # h0 = 2, g_o = z, s = x
        inner = CMGInnerGimbal(c; axis=SVector(0.0, 1.0, 0.0))
        @test cmg_momentum(inner, 0.3, 0.0) ≈ cmg_momentum(c, 0.3)                # no inner angle: the single gimbal
        @test cmg_momentum(inner, 0.0, π / 2) ≈ SVector(0.0, 0.0, -2.0) atol=1e-15 # x turned about y to -z
        α, β, ε = 0.4, -0.3, 1e-7
        outer_axis, inner_axis = cmg_torque_axes(inner, α, β)
        @test outer_axis ≈ (cmg_momentum(inner, α + ε, β) - cmg_momentum(inner, α - ε, β)) / (2ε) rtol=1e-6
        @test inner_axis ≈ (cmg_momentum(inner, α, β + ε) - cmg_momentum(inner, α, β - ε)) / (2ε) rtol=1e-6
    end

    @testset "double-gimbal CMGs: a run integrates both angles and conserves the angular momentum" begin
        sc, root = far_spacecraft(ω0=SVector(0.001, 0.0, 0.0))
        c = add_actuator!(sc, root, gyro(gimbal_angle_rad=0.1))
        add_actuator!(sc, root, CMGInnerGimbal(c; axis=SVector(0.0, 1.0, 0.0), gimbal_angle_rad=-0.2))
        res = cmg_run(sc, [0.02, 0.01]; mission_time=20.0)
        x0, x = res.solution.u[1].sc[1], res.solution.u[end].sc[1]
        @test x.cmg_angles[1] ≈ 0.1 + 0.02 * 20.0 rtol=1e-9
        @test x.cmg_angles[2] ≈ -0.2 + 0.01 * 20.0 atol=1e-12
        @test res.table.sc1_cmg_angle_2[end] ≈ x.cmg_angles[2] && res.table.sc1_actuator_level_2[end] == 0.01
        @test inertial_momentum(sc, x) ≈ inertial_momentum(sc, x0) rtol=1e-8
        @test norm(SVector{3}(x.ω) - SVector(0.001, 0.0, 0.0)) > 1e-5             # the body took up the rotor's change
    end

    @testset "double-gimbal CMGs: steering is exact with two of them, and an inner stop holds" begin
        sc, root = far_spacecraft()
        for (g, s, i) in ((SVector(0.0, 0.0, 1.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)),
                          (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)))
            c = add_actuator!(sc, root, ControlMomentGyro(rotor_momentum_nms=2.0, gimbal_axis=g, spin_axis=s, gimbal_angle_rad=0.2))
            add_actuator!(sc, root, CMGInnerGimbal(c; axis=i, gimbal_angle_rad=0.1))
        end
        a = sc.root.rw_assembly
        angles = [0.2, 0.2, 0.1, 0.1]                                             # outer, outer, inner, inner
        τ = SVector(0.01, -0.02, 0.03)
        rates = cmg_steering(sc, angles, τ; λ0=0.0)
        columns = [cmg_torque_axes(a.cmg_inner[1], 0.2, 0.1)..., cmg_torque_axes(a.cmg_inner[2], 0.2, 0.1)...]
        made = -(rates[1] * columns[1] + rates[3] * columns[2] + rates[2] * columns[3] + rates[4] * columns[4])
        @test made ≈ τ rtol=1e-9
        held, hroot = far_spacecraft()
        c = add_actuator!(held, hroot, gyro())
        add_actuator!(held, hroot, CMGInnerGimbal(c; axis=SVector(0.0, 1.0, 0.0), max_gimbal_angle_rad=0.05))
        x = cmg_run(held, [0.0, 0.01]; mission_time=10.0).solution.u[end].sc[1]   # the inner stop at 5 s
        @test x.cmg_angles[2] == 0.05
    end
end

end # module ControlMomentGyroTests
