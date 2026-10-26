module ActuatorFaultTests
using Test, StaticArrays, LinearAlgebra, SpaceAGORA, TOML
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const G0 = 9.80665

mutable struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

# A controller whose command grows with the time of its tick.
mutable struct RampLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
    slope::Float64
end
SpaceAGORA.actuation_mode(::RampLevels) = :mounted
SpaceAGORA.actuator_levels(c::RampLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(c::RampLevels, u, p, t::Float64, i::Int64) = (c.levels[1] = c.slope * t; nothing)

function far_spacecraft(; mass=100.0, prop_mass=10.0)
    root = SM.Link(root=true, m=mass, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    return SM.SpacecraftModel(links=[root], root=root, prop_mass=prop_mass, inertia_tensor=root.inertia, initial_condition=ic, id=1), root
end

function held_run(sc, levels; mission_time=5.0, rate=1.0, effector=HoldLevels(collect(Float64, levels)), settings=(;), kw...)
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
    extra = (;)
    if !isempty(settings)
        st = base.simulation_settings
        fields = NamedTuple{fieldnames(typeof(st))}(map(n -> getfield(st, n), fieldnames(typeof(st))))
        extra = (; simulation_settings=SM.SimulationSettings(; merge(fields, settings)...))
    end
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(effector,), control_rates=[rate]), extra...)
    return run_simulation(args; isolate_state=false, return_results=true, return_solution=true, kw...)
end

@testset "Actuator faults" begin
    @testset "the record and mounting" begin
        jet = Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0))
        @test_throws ArgumentError ActuatorFault(jet; kind=:melted)
        @test_throws ArgumentError ActuatorFault(jet; kind=:failed, start_s=5.0, stop_s=5.0)
        @test_throws ArgumentError ActuatorFault(jet; kind=:degraded, value=1.5)
        @test_throws ArgumentError ActuatorFault(jet; kind=:biased)                 # needs a value
        sc, root = far_spacecraft()
        @test_throws ArgumentError add_fault!(sc, ActuatorFault(jet; kind=:failed))  # not mounted
        add_actuator!(sc, root, jet)
        f = add_fault!(sc, ActuatorFault(jet; kind=:failed, start_s=1.0))
        @test SM.has_faults(sc) && sc.root.faults.faults == [f]
    end

    @testset "each kind changes the level as documented" begin
        sc, root = far_spacecraft()
        panel = SM.Link(m=1.0, r=MVector(0.0, 1.0, 0.0))
        push!(sc.links, panel)
        jet = add_actuator!(sc, root, Thruster(max_thrust=2.0, direction=MVector(1.0, 0.0, 0.0)))
        wheel = add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0), max_torque_nm=0.1))
        rod = add_actuator!(sc, root, Magnetorquer(axis=SVector(1.0, 0.0, 0.0), max_dipole_am2=1.0))
        hinge = add_actuator!(sc, panel, PanelHinge(axis=SVector(0.0, 1.0, 0.0), angle_rad=0.3))
        faults = [add_fault!(sc, ActuatorFault(jet; kind=:degraded, value=0.5)),
                  add_fault!(sc, ActuatorFault(wheel; kind=:failed)),
                  add_fault!(sc, ActuatorFault(rod; kind=:biased, value=0.25)),
                  add_fault!(sc, ActuatorFault(hinge; kind=:failed))]
        SM.resolve_faults!(sc)
        @test [f.index for f in faults] == [1, 2, 3, 4]
        levels = [3.0, 0.05, 0.5, 1.0]
        @test SM.faulted_levels(sc, levels) === levels                               # none active yet
        @test SM.refresh_faults!(sc, levels, 0.0)
        out = SM.faulted_levels(sc, levels)
        @test out[1] == 1.0                    # 3 N on a 2 N jet at half: within the halved limit of 1 N
        @test out[2] == 0.0                    # failed wheel
        @test out[3] == 0.75                   # biased rod
        @test out[4] == 0.3                    # a failed hinge holds its angle
        @test levels == [3.0, 0.05, 0.5, 1.0]  # the command is untouched
        stuck_sc, stuck_root = far_spacecraft()
        w = add_actuator!(stuck_sc, stuck_root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)))
        add_fault!(stuck_sc, ActuatorFault(w; kind=:stuck, start_s=2.0, stop_s=4.0))
        SM.resolve_faults!(stuck_sc)
        @test !SM.refresh_faults!(stuck_sc, [0.01], 1.0)
        @test SM.refresh_faults!(stuck_sc, [0.02], 2.0) && SM.faulted_levels(stuck_sc, [0.07]) == [0.02]   # held at onset
        @test SM.refresh_faults!(stuck_sc, [0.07], 4.0) && SM.faulted_levels(stuck_sc, [0.07]) == [0.07]   # over at stop
        @test SM.fault_edge(stuck_sc, 2.0) && SM.fault_edge(stuck_sc, 4.0) && !SM.fault_edge(stuck_sc, 3.0)
    end

    @testset "a wheel that fails between ticks stops at the fault time exactly" begin
        function wheel_run(fail)
            sc, root = far_spacecraft()
            wheel = add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)))
            fail && add_fault!(sc, ActuatorFault(wheel; kind=:failed, start_s=2.5))
            return held_run(sc, [0.01]), sc.inertia_tensor
        end
        (failed, inertia), (healthy, _) = wheel_run(true), wheel_run(false)
        x = failed.solution.u[end].sc[1]
        @test x.ω[1] ≈ 0.01 * 2.5 / inertia[1, 1] rtol=1e-9                           # torque for 2.5 s only
        @test x.h_wheels[1] ≈ -0.01 * 2.5 rtol=1e-9
        @test healthy.solution.u[end].sc[1].h_wheels[1] ≈ -0.01 * 5.0 rtol=1e-9
        col = failed.table.sc1_fault_active_1
        @test all(i -> col[i] == (failed.table.time[i] >= 2.5 ? 1.0 : 0.0), eachindex(col))
        @test all(i -> failed.table.sc1_actuator_level_1[i] == (failed.table.time[i] >= 2.5 ? 0.0 : 0.01), eachindex(col))
    end

    @testset "a degraded thruster burns half its propellant in the window; a stuck hinge stays" begin
        sc, root = far_spacecraft(mass=500.0)
        jet = add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0), Isp=200.0))
        add_fault!(sc, ActuatorFault(jet; kind=:degraded, value=0.5, start_s=1.5, stop_s=3.5))
        res = held_run(sc, [1.0])
        used = 510.0 - res.solution.u[end].sc[1].mass
        @test used ≈ (3.0 * 1.0 + 2.0 * 0.5) / (200.0 * G0) rtol=1e-9                  # 3 s full, 2 s at half
        hsc, hroot = far_spacecraft()
        panel = SM.Link(m=2.0, r=MVector(0.0, 1.5, 0.0), dims=MVector(0.02, 2.0, 1.0))
        push!(hsc.links, panel)
        hinge = add_actuator!(hsc, panel, PanelHinge(axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.1))
        add_fault!(hsc, ActuatorFault(hinge; kind=:stuck, start_s=2.5))
        held_run(hsc, [1.0])
        @test 0.0 < hinge.angle_rad < 0.25 + 1e-12                                     # moved until the fault, then held
        stuck_at = hinge.angle_rad
        @test panel.q ≈ [0.0, sin(stuck_at / 2), 0.0, cos(stuck_at / 2)]
    end

    @testset "a stuck actuator keeps its onset level across checkpoint segments" begin
        # The command ramps by 0.001 N m per second at 10 Hz ticks; the wheel sticks
        # at 0.55 s with the command of the 0.5 s tick. The 1 s segments must not
        # take a later command.
        function ramp_run(; settings=(;), mission_time=5.0)
            sc, root = far_spacecraft()
            wheel = add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)))
            add_fault!(sc, ActuatorFault(wheel; kind=:stuck, start_s=0.55))
            return held_run(sc, [0.0]; rate=0.1, effector=RampLevels([0.0], 0.001), settings=settings,
                mission_time=mission_time)   # a tick every 0.1 s
        end
        after(r) = [r.table.sc1_actuator_level_1[i] for i in eachindex(r.table.time) if r.table.time[i] >= 0.55]
        dir = mktempdir()
        straight = ramp_run()
        split = ramp_run(settings=(checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=dir))
        @test !isempty(after(straight)) && all(l -> isapprox(l, 0.0005; rtol=1e-9), after(straight))
        @test !isempty(after(split)) && all(l -> isapprox(l, 0.0005; rtol=1e-9), after(split))
        @test split.solution.u[end].sc[1].h_wheels[1] ≈ straight.solution.u[end].sc[1].h_wheels[1] rtol=1e-9
        # A run resumed from the file keeps the onset level, which the file holds.
        resumed_dir = mktempdir()
        first_leg = ramp_run(settings=(checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=resumed_dir),
            mission_time=2.0)
        resumed = ramp_run(settings=(checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=resumed_dir,
            resume_from_checkpoint=true, results=true, results_directory=resumed_dir))
        # It resumed at 2 s rather than starting again, and its table holds the whole
        # run: the first leg's rows before 2 s (the row at 2 s is saved again after
        # the ticks due then), and its own after them.
        @test TOML.parsefile(joinpath(resumed_dir, "simulation_results.manifest.toml"))["resumed_from_s"] == 2.0
        @test resumed.table.time[1] == 0.0
        @test isequal(resumed.table[resumed.table.time .< 2.0, :], first_leg.table[first_leg.table.time .< 2.0, :])
        resumed_rows = [l for (t, l) in zip(resumed.table.time, resumed.table.sc1_actuator_level_1) if t >= 2.0 - 1e-9]
        @test !isempty(resumed_rows) && all(l -> isapprox(l, 0.0005; rtol=1e-9), resumed_rows)
        @test resumed.solution.u[end].sc[1].h_wheels[1] ≈ straight.solution.u[end].sc[1].h_wheels[1] rtol=1e-9
        # The fault moved to the other wheel of the same kind refuses the file.
        function two_wheel_run(faulted; settings, mission_time=5.0)
            sc, root = far_spacecraft()
            wheels = [add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0))),
                      add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(0.0, 1.0, 0.0)))]
            add_fault!(sc, ActuatorFault(wheels[faulted]; kind=:stuck, start_s=0.55))
            return held_run(sc, [0.001, 0.001]; rate=0.1, settings=settings, mission_time=mission_time)
        end
        two_dir = mktempdir()
        two_wheel_run(1; settings=(checkpoint_enabled=true, checkpoint_interval_s=1.0, checkpoint_directory=two_dir), mission_time=2.0)
        @test_throws ArgumentError two_wheel_run(2; settings=(resume_from_checkpoint=true, checkpoint_directory=two_dir))
        same_wheel = two_wheel_run(1; settings=(resume_from_checkpoint=true, checkpoint_directory=two_dir,
            results=true, results_directory=two_dir)).table
        @test TOML.parsefile(joinpath(two_dir, "simulation_results.manifest.toml"))["resumed_from_s"] == 2.0
        @test same_wheel.time[1] == 0.0 && same_wheel.time[end] ≈ 5.0              # the table holds the whole run
    end

    @testset "a thruster with firing dynamics stops at its fault time, not at the next tick" begin
        # Ticks every second from 1 s schedule the jet fully on (about 1 N, a fast
        # ramp); it fails at 2.5 s, between the ticks at 2 s and 3 s.
        function firing_run(fail)
            sc, root = far_spacecraft(mass=500.0)
            jet = add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0), Isp=200.0,
                min_firing_time=0.05, cutoff_frequency=1.0e6, firing_dynamics=true))
            fail && add_fault!(sc, ActuatorFault(jet; kind=:failed, start_s=2.5))
            return held_run(sc, [1.0])
        end
        failed, healthy = firing_run(true), firing_run(false)
        used(r) = 510.0 - r.solution.u[end].sc[1].mass
        @test used(failed) ≈ 1.5 / (200.0 * G0) rtol=1e-4                             # fired from 1 s to 2.5 s
        @test used(healthy) ≈ 4.0 / (200.0 * G0) rtol=1e-4
        t, col = failed.table.time, failed.table.sc1_actuator_level_1
        @test all(i -> t[i] < 2.5 || col[i] == 0.0, eachindex(t))
    end

    @testset "the thruster-level columns see the faults" begin
        sc, root = far_spacecraft()
        jet = add_actuator!(sc, root, Thruster(max_thrust=2.0, direction=MVector(1.0, 0.0, 0.0)))
        add_fault!(sc, ActuatorFault(jet; kind=:degraded, value=0.5, start_s=1.0))
        SM.ActuatorFaults.resolve_faults!(sc)
        c = HoldLevels([2.0])
        @test SM.SimulationCallbacks._effector_thruster_levels(c, 1, sc) == [1.0]
        SM.ActuatorFaults.refresh_faults!(sc, c.levels, 1.0)
        @test SM.SimulationCallbacks._effector_thruster_levels(c, 1, sc) == [0.5]
    end

    @testset "the power budget sheds by the faulted load from the fault's start" begin
        # About 12 W of arrays against a 5 W load and a laser drawing 10 W: the empty
        # battery sheds the laser until a fault halves its level at 2 s; the load
        # (10 W) then fits the arrays, and the laser is restored at once.
        sc, root = far_spacecraft()
        ephem = SM.SimpleEphemeridesModel()
        epoch = SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
        sun = SM.ephemerides_sun_ii(ephem, SM.Earth(), SM.ephemerides_epoch_seconds(ephem, epoch))
        add_power!(sc, root, SolarArray(normal=sun, area_m2=0.03, efficiency=0.3))
        add_power!(sc, root, Battery(capacity_wh=5.0, initial_soc=0.0, recovery_soc=0.5))
        add_power!(sc, root, PowerLoad(power_w=5.0))
        laser = add_actuator!(sc, root, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=2.0, efficiency=0.2))
        add_fault!(sc, ActuatorFault(laser; kind=:degraded, value=0.5, start_s=2.0))
        res = held_run(sc, [2.0]; mission_time=20.0)
        gen = res.table.sc1_power_generated_w[1]
        @test 10.0 < gen < 15.0
        @test res.solution.u[end].sc[1].battery_wh ≈ ((gen - 5.0) * 2.0 + (gen - 10.0) * 18.0) / 3600 rtol=1e-6
    end

    @testset "faults are opt-in" begin
        sc, root = far_spacecraft()
        add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)))
        @test !any(n -> occursin("fault", n), names(held_run(sc, [0.01]).table))
    end
end

end # module ActuatorFaultTests
