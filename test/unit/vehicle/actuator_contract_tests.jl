module ActuatorContractTests
using Test, StaticArrays, LinearAlgebra, Random, ComponentArrays, SpaceAGORA
using DiffEqCallbacks: PresetTimeCallback
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const CH = SM.ControlHooks
const G0 = 9.80665
const ZERO3 = SVector(0.0, 0.0, 0.0)

# A user-defined actuator kind, written the way a user would: a cold-gas puck
# that pushes along a fixed axis of its link.
struct ColdGasPuck <: AbstractActuator
    axis::SVector{3, Float64}
    location::SVector{3, Float64}
    max_force_n::Float64
    isp_s::Float64
end
SpaceAGORA.actuator_limits(a::ColdGasPuck) = (0.0, a.max_force_n)
function SpaceAGORA.actuator_wrench(a::ColdGasPuck, level, ctx::ActuatorContext)
    force = ctx.link_to_body * (level * a.axis)
    arm = ctx.link_to_body * a.location + ctx.link_r_body
    return force, cross(arm, force)
end
SpaceAGORA.actuator_consumption(a::ColdGasPuck, level, ctx) = (mass_rate_kg_s = level / (a.isp_s * G0), power_w = 0.0)

# A control effector that commands mounted actuators and writes no force,
# mass-flow or firing-level hook of its own.
mutable struct HoldLevels <: AbstractControlEffectorModel
    sat_idx::Int
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == c.sat_idx ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

function readout_params(spacecraft::Vector{SM.SpacecraftModel}, effectors::Tuple)
    planet = SM.Earth()
    args = SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=1.0,
            orientation_sim=true, num_steps_to_save=1),
        environment_model=SM.EnvironmentModel(planet=planet, EI=120.0,
            density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=SM.DynamicsModel(spacecraft, (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=effectors, control_rates=fill(1.0, length(effectors))),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
    return SM.ODEParams(n_sats=length(spacecraft), args=args)
end

function readout(p, i, state)
    f, tq, rw = MVector{3, Float64}(ZERO3), MVector{3, Float64}(ZERO3), MVector{3, Float64}(ZERO3)
    rate = SE._accumulate_control_effectors!(f, tq, rw, state, p, i, 2.0, false)
    return SVector{3, Float64}(f), SVector{3, Float64}(tq), SVector{3, Float64}(rw), rate
end

body_state(q; mass=500.0, ω=ZERO3) = ComponentVector(mass=mass, q=collect(q), ω=collect(ω))
spacecraft_with(root; id=1, links=[root]) = SM.SpacecraftModel(links=links, root=root, prop_mass=10.0, inertia_tensor=root.inertia, id=id)

include(joinpath(@__DIR__, "..", "..", "..", "templates", "actuator_template.jl"))

function wheel_spin_up(max_momentum)
    inertia = SMatrix{3, 3, Float64}(2.0, 0, 0, 0, 3.0, 0, 0, 0, 4.0)
    root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), ZERO3; q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1.0, inertia_tensor=inertia, initial_condition=ic, id=4)
    for axis in (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        add_actuator!(sc, root, ReactionWheel(spin_axis=axis, max_torque_nm=0.02, max_momentum_nms=max_momentum))
    end
    control = HoldLevels(1, [0.01, 0.0, 0.0])                # 0.01 N m on the body about +x
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=5.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(control,), control_rates=[1.0]))
    return run_simulation(args; return_solution=true).u[end].sc[1], inertia
end


function pyramid_run(levels)
    s2 = 1 / sqrt(2)
    axes = (SVector(s2, 0.0, s2), SVector(-s2, 0.0, s2), SVector(0.0, s2, s2), SVector(0.0, -s2, s2))
    inertia = SMatrix{3, 3, Float64}(2.0, 0, 0, 0, 3.0, 0, 0, 0, 4.0)
    root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), ZERO3; q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1.0, inertia_tensor=inertia, initial_condition=ic, id=7)
    foreach(a -> add_actuator!(sc, root, ReactionWheel(spin_axis=a)), axes)
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=5.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(HoldLevels(1, collect(levels)),), control_rates=[1.0]))
    return run_simulation(args; return_solution=true).u[end].sc[1], inertia, hcat(collect.(axes)...)
end

# A user-defined kind that asks for sampled environment and records what its
# context carried at its last evaluation. It puts no wrench on the vehicle.
struct EnvironmentRecorder <: AbstractActuator
    seen::Base.RefValue{Any}
    requirements::ActuatorRequirements
end
SpaceAGORA.actuator_limits(::EnvironmentRecorder) = (-Inf, Inf)
SpaceAGORA.actuator_requirements(a::EnvironmentRecorder) = a.requirements
function SpaceAGORA.actuator_wrench(a::EnvironmentRecorder, level, ctx::ActuatorContext)
    a.seen[] = ctx.environment
    return ZERO3, ZERO3
end

# A spacecraft far from Earth, attitude at the identity, for mounted runs.
function far_spacecraft(; id, mass=100.0, prop_mass=1.0)
    root = SM.Link(root=true, m=mass, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), ZERO3; q=SVector(0.0, 0.0, 0.0, 1.0))
    return SM.SpacecraftModel(links=[root], root=root, prop_mass=prop_mass, inertia_tensor=root.inertia, initial_condition=ic, id=id), root
end

# A run whose only controller holds `levels` for the mounted actuators of `sc`.
function held_run(sc, levels; mission_time=2.0, rate=1.0, kw...)
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(HoldLevels(1, collect(Float64, levels)),), control_rates=[rate]))
    return run_simulation(args; isolate_state=false, kw...)
end

# A 110 kg spacecraft in a 7000 km circular orbit, attitude at the identity
# (body x radial, so an engine along body x points 90 degrees off prograde),
# with a 1 N, 300 s engine along body x and, optionally, three wheels.
function orbit_spacecraft(; id, wheels=true)
    root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
    speed = sqrt(SM.Earth().μ / 7.0e6)
    ic = SM.CartesianInitialCondition(SVector(7.0e6, 0.0, 0.0), SVector(0.0, speed, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=10.0, inertia_tensor=root.inertia, initial_condition=ic, id=id)
    add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0), Isp=300.0))
    if wheels
        for axis in (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
            add_actuator!(sc, root, ReactionWheel(spin_axis=axis))
        end
    end
    return sc
end

# A prograde 1 N burn from 300.5 s to 400.5 s, applied ideally or through the
# mounted engine; the window's ends fall between the 1 s control ticks.
function burn_run(sc, actuation; mission_time=450.0, scene=false, sample_period=nothing)
    model = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[300.5],
        stop_burn_time=[400.5], Isp=[300.0], actuation=actuation)
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredGravityModel(),), ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(model,), control_rates=[1.0]))
    if sample_period !== nothing
        mission = args.mission_configuration
        fields = NamedTuple{fieldnames(typeof(mission))}(Tuple(getfield(mission, f) for f in fieldnames(typeof(mission))))
        args = SM.SimConfig._with_configuration(args;
            mission_configuration=SM.MissionConfiguration(; merge(fields, (data_rate=sample_period,))...))
    end
    # The scene adds the `thruster_level` field to the returned table.
    args = scene ? SM.SceneVisualization.with_visualization_scene(args, true) : args
    return run_simulation(args; return_results=true, return_solution=true)
end

@testset "Actuator contract" begin
    @testset "a thruster's wrench, limits and propellant follow from its record" begin
        root = SM.Link(root=true, m=10.0)
        jet = Thruster(max_thrust=2.0, location=MVector(1.0, 0.0, 0.0), direction=MVector(0.0, 0.0, 2.0), Isp=200.0)
        ctx = actuator_context(root)
        force, torque = actuator_wrench(jet, 1.5, ctx)
        @test force == SVector(0.0, 0.0, 1.5)                 # along the (normalized) thrust direction
        @test torque == SVector(0.0, -1.5, 0.0)               # (1, 0, 0) x (0, 0, 1.5)
        @test actuator_consumption(jet, 1.5, ctx).mass_rate_kg_s == 1.5 / (200.0 * G0)
        @test actuator_consumption(jet, 1.5, ctx).power_w == 0.0
        @test actuator_limits(jet) == (0.0, 2.0)
        @test actuator_saturate(jet, 5.0) == 2.0
        @test actuator_saturate(jet, -1.0) == 0.0
        @test actuator_saturate(jet, NaN) == 0.0
        @test actuator_fraction(jet, 1.0) == 0.5
        @test actuator_unit_wrench(jet, ctx) == actuator_wrench(jet, 1.0, ctx)
        @test actuator_consumption(Thruster(Isp=0.0), 1.0, ctx).mass_rate_kg_s == 0.0
        @test actuator_wrench(Thruster(direction=MVector(0.0, 0.0, 0.0)), 1.0, ctx) == (ZERO3, ZERO3)

        # On a rotated, offset link the torque is the column of the existing
        # thruster torque Jacobian (ThrusterHooks), so both use one link-to-body map.
        panel = SM.Link(r=MVector(0.0, 1.0, 0.5), q=MVector(0.0, 0.0, sin(pi / 8), cos(pi / 8)))
        side = Thruster(max_thrust=1.0, location=MVector(0.3, -0.2, 0.1), direction=MVector(0.0, 1.0, 1.0))
        R = Matrix(SM.rotate_to_body(panel))
        d = normalize(Vector(side.direction))
        _, panel_torque = actuator_wrench(side, 1.0, actuator_context(panel))
        @test panel_torque ≈ cross(R * Vector(side.location) + Vector(panel.r), R * d) atol=1e-14
    end

    @testset "add_actuator! mounts thrusters and user kinds in mounting order" begin
        root = SM.Link(root=true, m=10.0)
        pod = SM.Link(m=1.0)
        sc = spacecraft_with(root; links=[root, pod])
        puck = ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0)
        @test add_actuator!(sc, pod, puck) === puck
        a = add_actuator!(sc, root, Thruster(max_thrust=1.0, location=MVector(0.0, 0.0, 1.0), direction=MVector(3.0, 0.0, 0.0)))
        b = add_actuator!(sc, pod, Thruster(max_thrust=2.0, direction=MVector(0.0, 1.0, 0.0)))
        c = add_actuator!(sc, root, Thruster(max_thrust=3.0, location=MVector(1.0, 0.0, 0.0), direction=MVector(0.0, 0.0, 1.0)))
        @test a.direction == [1.0, 0.0, 0.0]                  # normalized on mounting
        @test !hasfield(SM.Link, :J_thruster) && !hasfield(SM.SpacecraftModel, :n_thrusters)   # the mounted hardware is the record
        @test SM.mounted_thruster_count(sc) == 3 && SM.mounted_actuator_count(sc) == 4
        order = [entry.actuator for entry in mounted_actuators(sc)]
        @test order[1] === a && order[2] === c && order[3] === b && order[4] === puck   # links in order, thrusters first
        @test pod.custom_actuators == [puck]
        B = allocation_matrix(sc)
        @test size(B) == (6, 4)
        @test B[:, 1] == vcat(SVector(1.0, 0.0, 0.0), cross(SVector(0.0, 0.0, 1.0), SVector(1.0, 0.0, 0.0)))
        @test B[1:3, 4] == SVector(1.0, 0.0, 0.0)
        # Links carry an empty list of user kinds by default.
        @test isempty(SM.Link().custom_actuators)
    end

    @testset "six_axis_thrusters reproduce the six-axis model" begin
        six = SM.SixAxisThrusterModel(locations_body=SMatrix{3, 6, Float64}(0.1 .* randn(MersenneTwister(1), 3, 6)),
            max_thrust_n=SVector{6}(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), isp_s=SVector{6}(fill(70.0, 6)), min_firing_time_s=0.02)
        jets = six_axis_thrusters(six)
        @test length(jets) == 6
        @test all(j -> Vector(jets[j].direction) == six.directions_body[:, j] && Vector(jets[j].location) == six.locations_body[:, j] &&
            jets[j].max_thrust == six.max_thrust_n[j] && jets[j].Isp == 70.0 && jets[j].min_firing_time == 0.02, 1:6)
        sc = spacecraft_with(SM.Link(root=true))
        foreach(j -> add_actuator!(sc, sc.root, j), jets)
        forces = SVector{6}(0.3, 0.0, 0.2, 0.0, 0.0, 0.7)
        f_rpo, tq_rpo = SM.rpo_thruster_wrench_body(forces, six)
        f, tq, mdot = mounted_wrench(sc, forces)
        @test f ≈ f_rpo atol=1e-15
        @test tq ≈ tq_rpo atol=1e-15
        @test mdot ≈ sum(forces) / (70.0 * G0)
    end

    @testset "allocation reproduces the existing allocators exactly" begin
        rng = MersenneTwister(20261015)
        six = SM.SixAxisThrusterModel(max_thrust_n=SVector{6}(0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
        sc = spacecraft_with(SM.Link(root=true))
        foreach(j -> add_actuator!(sc, sc.root, j), six_axis_thrusters(six))
        B = allocation_matrix(sc)
        mismatches = 0
        for _ in 1:2000
            command = SVector{3}(randn(rng, 3))
            m = 0.5 * rand(rng)
            a = allocate_actuators(B, command, 0.0, collect(six.max_thrust_n); method=:axis_clamp, minimum_fraction=m)
            mismatches += a == collect(SM.rpo_allocate_six_axis_thrusters(command, six; minimum_throttle_fraction=m)) ? 0 : 1
        end
        @test mismatches == 0

        # The descent controller's jet levels: the least-norm allocation of a
        # torque over the jets' full-thrust torque arms, clipped to 0..1.
        jet(r, thrust_dir, rating) = Thruster(location=MVector{3, Float64}(r), direction=MVector{3, Float64}(thrust_dir), max_thrust=rating)
        root = SM.Link(root=true, m=500.0, thrusters=[jet((0.0, 1.0, 0.0), (0.0, 0.0, 2.0), 10.0), jet((0.0, 0.0, 0.0), (0.0, 0.0, -1.0), 100.0)])
        pod = SM.Link(thrusters=[jet((0.0, 0.0, 1.0), (3.0, 0.0, 0.0), 20.0), jet((1.0, 0.0, 0.0), (0.0, 4.0, 0.0), 30.0),
                                 jet((0.5, 0.5, 0.0), (1.0, -1.0, 1.0), 15.0)])
        craft = SM.SpacecraftModel(links=[root, pod], root=root, id=97)
        layout = CH.descent_thruster_layout(craft)
        levels = zeros(5)
        mismatches = 0
        for _ in 1:2000
            torque = SVector{3}(30.0 .* randn(rng, 3))
            CH.descent_thruster_levels!(levels, layout, 0.0, torque)
            mismatches += allocate_actuators(layout.torque_arms_nm, torque, 0.0, 1.0) == levels[layout.jets] ? 0 : 1
        end
        @test mismatches == 0

        # Jets on a rotated, offset link: the descent layout's torque arms are
        # the mounted jets' torques at full thrust, so its levels act as allocated.
        tilted = SM.Link(r=MVector(0.4, -0.3, 1.2), q=MVector(normalize([0.1, 0.3, -0.2, 0.9])...),
            thrusters=[jet((0.5, 0.0, 0.0), (0.0, 1.0, 0.0), 12.0), jet((0.0, 0.5, 0.0), (0.0, 0.0, 1.0), 14.0),
                       jet((0.0, 0.0, 0.5), (1.0, 0.0, 0.0), 16.0)])
        engine_link = SM.Link(root=true, thrusters=[jet((0.0, 0.0, 0.0), (0.0, 0.0, -1.0), 100.0)])
        posed = SM.SpacecraftModel(links=[engine_link, tilted], root=engine_link, id=98)
        posed_layout = CH.descent_thruster_layout(posed)
        B_posed = allocation_matrix(posed)
        ratings = [12.0, 14.0, 16.0]
        @test posed_layout.jets == [2, 3, 4]
        @test posed_layout.torque_arms_nm ≈ B_posed[4:6, 2:4] .* ratings' atol=1e-12

        # A feasible wrench is reproduced by the least-norm levels.
        full = allocation_matrix(spacecraft_with(SM.Link(root=true; thrusters=six_axis_thrusters(SM.SixAxisThrusterModel(max_thrust_n=SVector{6}(fill(5.0, 6)))))))
        target = vcat(SVector(0.2, -0.1, 0.3), ZERO3)
        @test full * allocate_actuators(full, target, -Inf, Inf) ≈ target
        floor_levels = allocate_actuators(full, target, 0.05, 5.0; method=:axis_clamp)
        @test all(>=(0.05), floor_levels) && all(<=(5.0), floor_levels)       # both bounds hold, idle jets included
        @test_throws ArgumentError allocate_actuators(full, target, 0.0, 1.0; method=:nearest)
        @test_throws ArgumentError allocate_actuators(full[:, [2, 1, 3, 4, 5, 6]], target, 0.0, 1.0; method=:axis_clamp)
        @test_throws DimensionMismatch allocate_actuators(full, target, zeros(3), 1.0)
    end

    @testset "a controller in mounted mode gets its wrench from the mounted actuators" begin
        root = SM.Link(root=true, m=500.0)
        sc = spacecraft_with(root)
        add_actuator!(sc, root, Thruster(max_thrust=2.0, location=MVector(0.0, 1.0, 0.0), direction=MVector(0.0, 0.0, 1.0), Isp=100.0))
        add_actuator!(sc, root, ColdGasPuck(SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.5), 0.5, 50.0))
        control = HoldLevels(1, [1.5, 9.0])                    # the puck's 9 N saturates at 0.5 N
        p = readout_params(SM.SpacecraftModel[sc], (control,))
        q = normalize(SVector(0.1, -0.2, 0.3, 0.9))
        force, torque, wheel, rate = readout(p, 1, body_state(q))
        f_body, tq_body, mdot = mounted_wrench(sc, [1.5, 0.5])
        @test f_body == SVector(0.5, 0.0, 1.5)
        @test tq_body == cross(SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.5)) + cross(SVector(0.0, 0.0, 0.5), SVector(0.5, 0.0, 0.0))
        @test force ≈ SM.rot(q)' * f_body                      # inertial force from the current attitude
        @test torque == tq_body
        @test wheel == ZERO3
        @test rate ≈ -(1.5 / (100.0 * G0) + 0.5 / (50.0 * G0))
        @test SpaceAGORA.actuation_mode(control) === :mounted
        @test SpaceAGORA.actuation_mode(SM.MagneticMomentumManagerModel()) === :ideal

        control.sat_idx = 2                                    # drives no spacecraft of this run
        @test readout(p, 1, body_state(q)) == (ZERO3, ZERO3, ZERO3, 0.0)
        control.sat_idx = 1
        control.levels = [1.0]
        @test_throws ArgumentError readout(p, 1, body_state(q))   # one level for two actuators
        control.levels = [1.0, 0.0]
        @test_throws ArgumentError readout(p, 1, ComponentVector(mass=500.0))   # no attitude state
    end

    @testset "mounted RPO actuation drives the chaser's mounted jets" begin
        six = SM.SixAxisThrusterModel(max_thrust_n=SVector{6}(fill(1.0, 6)), isp_s=SVector{6}(fill(65.0, 6)))
        root = SM.Link(root=true, m=12.0)
        chaser = spacecraft_with(root)
        foreach(j -> add_actuator!(chaser, root, j), six_axis_thrusters(six))
        forces = SVector{6}(0.4, 0.0, 0.0, 0.25, 0.1, 0.0)
        rw = SVector(1e-4, -2e-4, 3e-4)
        function held_model(actuation)
            model = SM.RPOMPCControlModel(chaser_index=1, target_index=2, thrusters=six, actuation=actuation)
            f_body, tq_body = SM.rpo_thruster_wrench_body(forces, six)
            model.held = SM.ControlHooks.RPOHeldActuation(force_ii=f_body, torque_body=tq_body + rw,
                thruster_forces_n=forces, rw_torque_body=rw)          # held at an identity attitude
            actuation === :mounted && CH._rpo_fill_mounted_levels!(model, chaser, forces)   # what a control tick does
            return model
        end
        ideal, mounted = held_model(:ideal), held_model(:mounted)
        p_ideal = readout_params(SM.SpacecraftModel[chaser], (ideal,))
        p_mounted = readout_params(SM.SpacecraftModel[chaser], (mounted,))
        identity = SVector(0.0, 0.0, 0.0, 1.0)
        fi, ti, wi, ri = readout(p_ideal, 1, body_state(identity; mass=12.0))
        fm, tm, wm, rm = readout(p_mounted, 1, body_state(identity; mass=12.0))
        @test fm ≈ fi atol=1e-15
        @test tm ≈ ti atol=1e-15                               # jets' torque plus the wheels'
        @test wm == wi == rw
        @test rm ≈ ri
        # Between ticks the vehicle turns: the ideal force stays where the tick
        # put it, the mounted jets turn with the body.
        q = normalize(SVector(0.0, 0.0, sin(0.2), cos(0.2)))
        @test readout(p_ideal, 1, body_state(q; mass=12.0))[1] == fi
        @test readout(p_mounted, 1, body_state(q; mass=12.0))[1] ≈ SM.rot(q)' * fm atol=1e-15
        @test SpaceAGORA.actuator_levels(mounted, 2) === nothing
        @test SpaceAGORA.actuator_levels(SM.RPOMPCControlModel(actuation=:mounted), 1) === nothing   # before the first tick
        @test_throws ArgumentError SpaceAGORA.actuation_mode(SM.RPOMPCControlModel(actuation=:bogus))
        # Other hardware on the chaser stays idle rather than failing the evaluation.
        add_actuator!(chaser, root, ColdGasPuck(SVector(0.0, 1.0, 0.0), ZERO3, 0.5, 50.0))
        CH._rpo_fill_mounted_levels!(mounted, chaser, forces)
        @test SpaceAGORA.actuator_levels(mounted, 1) == vcat(collect(forces), 0.0)
        @test readout(p_mounted, 1, body_state(identity; mass=12.0))[1] ≈ fm atol=1e-15
        @test_throws ArgumentError CH._rpo_fill_mounted_levels!(mounted, spacecraft_with(SM.Link(root=true)), forces)
    end

    @testset "mounted controllers fill the thruster_level columns" begin
        root = SM.Link(root=true, m=50.0)
        sc = spacecraft_with(root)
        add_actuator!(sc, root, Thruster(max_thrust=4.0, direction=MVector(1.0, 0.0, 0.0)))
        add_actuator!(sc, root, Thruster(max_thrust=2.0, direction=MVector(-1.0, 0.0, 0.0)))
        control = HoldLevels(1, [1.0, 3.0])
        p = readout_params(SM.SpacecraftModel[sc], (control,))
        @test SM.SimulationCallbacks.thruster_level_counts(p.args) == [2]
        levels = SM.SimulationCallbacks._save_thruster_levels(1, [2], nothing, 0.0, (p=p,))
        @test levels == [[0.25, 1.0]]                          # 1 of 4 N, and 3 N saturated at 2 N
        @test SM.mounted_thruster_fractions(sc, [1.0]) == [0.25, 0.0]
    end

    @testset "the actuator template mounts and drives through the engine" begin
        T = ExampleActuatorTemplate
        root = SM.Link(root=true, m=200.0)
        sc = spacecraft_with(root)
        nozzle = add_actuator!(sc, root, T.ExampleColdGasNozzle(SVector(1.0, 0.0, 0.0), SVector(0.0, 0.2, 0.0), 0.5, 60.0))
        @test root.custom_actuators == [nozzle]
        control = T.ExampleNozzleController(1, [0.25])
        p = readout_params(SM.SpacecraftModel[sc], (control,))
        force, torque, _, rate = readout(p, 1, body_state(SVector(0.0, 0.0, 0.0, 1.0); mass=210.0))
        @test force == SVector(0.25, 0.0, 0.0)
        @test torque == cross(SVector(0.0, 0.2, 0.0), SVector(0.25, 0.0, 0.0))
        @test rate == -0.25 / (60.0 * 9.80665)
    end

    @testset "a reaction wheel mounts on the root link and its torque and momentum are limited" begin
        root = SM.Link(root=true, m=50.0)
        sc = spacecraft_with(root)
        wheel = ReactionWheel(spin_axis=SVector(0.0, 0.0, 2.0), max_torque_nm=0.01, max_momentum_nms=0.05)
        @test wheel.spin_axis == SVector(0.0, 0.0, 1.0)
        @test add_actuator!(sc, root, wheel) === wheel
        @test root.rw_assembly.n_wheels == 1 && root.rw_assembly.J_rw == reshape([0.0, 0.0, 1.0], 3, 1)
        @test root.rw_assembly.h_wheels == [0.0] && SM.mounted_wheel_count(sc) == 1 && !hasfield(SM.SpacecraftModel, :n_reaction_wheels)
        @test_throws ArgumentError add_actuator!(sc, SM.Link(m=1.0), ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)))
        @test actuator_limits(wheel) == (-0.01, 0.01)
        @test actuator_saturate(wheel, 1.0) == 0.01 && actuator_saturate(wheel, -1.0) == -0.01
        @test actuator_wrench(wheel, 0.005, actuator_context(root)) == (ZERO3, SVector(0.0, 0.0, 0.005))
        @test SM.wheel_momentum_saturate(wheel, 0.005, -0.05) == 0.0       # at the limit, storing more (dh/dt = -level)
        @test SM.wheel_momentum_saturate(wheel, -0.005, -0.05) == -0.005   # unloading is allowed
        @test SM.wheel_momentum_saturate(wheel, 0.005, -0.04) == 0.005
        @test_throws ArgumentError ReactionWheel(spin_axis=ZERO3)
        @test_throws ArgumentError ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0), max_torque_nm=0.0)
        @test isinf(ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)).max_momentum_nms)   # unlimited by default
        # Mounting order: thrusters, then wheels, then user-defined kinds.
        puck = add_actuator!(sc, root, ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0))
        jet = add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(0.0, 1.0, 0.0)))
        order = [entry.actuator for entry in mounted_actuators(sc)]
        @test order[1] === jet && order[2] === wheel && order[3] === puck
        force, torque, _, wheel_torque = mounted_wrench(sc, [0.0, 0.004, 0.0]; wheel_momentum=[0.0])
        @test force == ZERO3 && wheel_torque == SVector(0.0, 0.0, 0.004) && torque == wheel_torque
        @test mounted_wrench(sc, [0.0, 0.004, 0.0]; wheel_momentum=[-0.05])[4] == ZERO3      # saturated
        @test allocation_matrix(sc)[:, 2] == [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    end

    @testset "mounted wheels spin the body up and store the opposite momentum" begin
        final, inertia = wheel_spin_up(Inf)
        @test final.ω[1] ≈ 0.01 * 5.0 / 2.0 rtol=1e-8              # τ T / I_xx
        @test final.h_wheels ≈ [-0.05, 0.0, 0.0] atol=1e-12
        @test inertia * collect(final.ω) + collect(final.h_wheels) ≈ zeros(3) atol=1e-12   # angular momentum conserved
        limited, _ = wheel_spin_up(0.03)
        @test limited.h_wheels[1] ≈ -0.03 atol=2e-4                 # stops storing at the limit
        @test limited.ω[1] ≈ 0.03 / 2.0 rtol=1e-2
    end

    @testset "the momentum manager reads the integrated wheel momentum when the vehicle carries wheels" begin
        root = SM.Link(root=true, m=50.0)
        sc = spacecraft_with(root)
        for axis in (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
            add_actuator!(sc, root, ReactionWheel(spin_axis=axis))
        end
        manager = SM.MagneticMomentumManagerModel(spacecraft_index=1, mu_gain=0.01, h_wheels_0=SVector(9.0, 9.0, 9.0),
            commanded_torque=(t, r, v, q, w) -> SVector(0.0, 0.0, 0.0), b_field_ii=(t, r) -> SVector(0.0, 0.0, 3.0e-5))
        p = readout_params(SM.SpacecraftModel[sc], (manager,))
        u = (sc=[(pos=[7.0e6, 0.0, 0.0], vel=[0.0, 7.5e3, 0.0], q=[0.0, 0.0, 0.0, 1.0], ω=zeros(3), h_wheels=[0.02, -0.01, 0.0])],)
        CH.update_control!(manager, u, p, 0.0, 1)
        @test manager.h_wheels == SVector(0.02, -0.01, 0.0)       # not the configured h_wheels_0
        bare = SM.MagneticMomentumManagerModel(spacecraft_index=1, h_wheels_0=SVector(9.0, 9.0, 9.0),
            commanded_torque=(t, r, v, q, w) -> SVector(0.0, 0.0, 0.0), b_field_ii=(t, r) -> SVector(0.0, 0.0, 3.0e-5))
        CH.update_control!(bare, (sc=[(pos=[7.0e6, 0.0, 0.0], vel=[0.0, 7.5e3, 0.0], q=[0.0, 0.0, 0.0, 1.0], ω=zeros(3))],),
            readout_params(SM.SpacecraftModel[spacecraft_with(SM.Link(root=true, m=50.0))], (bare,)), 0.0, 1)
        @test bare.h_wheels == SVector(9.0, 9.0, 9.0)               # no wheels: its own accumulator
    end

    @testset "a torque rod's torque is its dipole crossed with the field" begin
        root = SM.Link(root=true, m=20.0)
        rod = Magnetorquer(axis=SVector(0.0, 0.0, 2.0), max_dipole_am2=1.5)
        @test rod.axis == SVector(0.0, 0.0, 1.0)
        @test actuator_limits(rod) == (-1.5, 1.5) && actuator_saturate(rod, -4.0) == -1.5
        field = SVector(2.0e-5, -1.0e-5, 3.0e-5)
        @test actuator_wrench(rod, 1.0, actuator_context(root; magnetic_field_body=field)) == (ZERO3, cross(SVector(0.0, 0.0, 1.0), field))
        @test actuator_wrench(rod, 1.0, actuator_context(root)) == (ZERO3, ZERO3)      # no field, no torque
        @test_throws ArgumentError Magnetorquer(axis=ZERO3)
        @test_throws ArgumentError Magnetorquer(axis=SVector(1.0, 0.0, 0.0), max_dipole_am2=-1.0)
        sc = spacecraft_with(root)
        @test add_actuator!(sc, root, rod) === rod && root.magnetorquers == [rod]
        wheel = add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)))
        puck = add_actuator!(sc, root, ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0))
        jet = add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(0.0, 1.0, 0.0)))
        @test [e.actuator for e in mounted_actuators(sc)] == AbstractActuator[jet, wheel, rod, puck]   # thrusters, wheels, rods, user kinds
        @test SM.mounted_magnetorquer_count(sc) == 1 && SM.mounted_actuator_count(sc) == 4
        @test mounted_wrench(sc, [0.0, 0.0, 1.0, 0.0]; magnetic_field_body=field)[2] == cross(SVector(0.0, 0.0, 1.0), field)
        @test isempty(SM.Link().magnetorquers)
    end

    @testset "a mounted rod acts in the run's field, as an equal permanent dipole does" begin
        q = normalize(SVector(0.1, 0.2, 0.3, 0.9))
        x = ComponentVector(pos=[7.0e6, 1.0e6, 2.0e5], vel=[0.0, 7.5e3, 0.0], mass=50.0, q=collect(q), ω=zeros(3))
        passive_root = SM.Link(root=true, m=50.0)
        push!(passive_root.magnets, SM.Magnet(m=MVector(0.3, -0.2, 0.8)))
        passive = spacecraft_with(passive_root)
        _, torque_passive = SE._evaluate_dynamic_effector(SM.MagneticTorqueRodModel(), x,
            SE.build_state_sample(x, passive, true, 1), readout_params(SM.SpacecraftModel[passive], ()), 1, 0.0)
        rod_root = SM.Link(root=true, m=50.0)
        rods = spacecraft_with(rod_root)
        for axis in (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
            add_actuator!(rods, rod_root, Magnetorquer(axis=axis))
        end
        _, torque_rods, _, _ = readout(readout_params(SM.SpacecraftModel[rods], (HoldLevels(1, [0.3, -0.2, 0.8]),)), 1, x)
        @test torque_rods ≈ torque_passive rtol=1e-12
        @test norm(torque_rods) > 0.0
    end

    @testset "the momentum manager drives mounted torque rods with its dipole command" begin
        root = SM.Link(root=true, m=50.0)
        sc = spacecraft_with(root)
        for axis in (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
            add_actuator!(sc, root, Magnetorquer(axis=axis))
        end
        estimate(t, r) = SVector(1.0e-5, -2.0e-5, 3.0e-5)
        manager = SM.MagneticMomentumManagerModel(spacecraft_index=1, mu_gain=0.01, h_wheels_0=SVector(0.02, 0.0, 0.01),
            commanded_torque=(t, r, v, q, w) -> ZERO3, b_field_ii=estimate, actuation=:mounted)
        p = readout_params(SM.SpacecraftModel[sc], (manager,))
        u = (sc=[(pos=[7.0e6, 0.0, 0.0], vel=[0.0, 7.5e3, 0.0], q=[0.0, 0.0, 0.0, 1.0], ω=zeros(3))],)
        @test SpaceAGORA.actuator_levels(manager, 1) === nothing                  # before the first tick
        CH.update_control!(manager, u, p, 0.0, 1)
        @test SpaceAGORA.actuation_mode(manager) === :mounted
        @test SpaceAGORA.actuator_levels(manager, 1) ≈ collect(manager.held_dipole_am2) atol=1e-15   # orthogonal rods carry its components
        x = ComponentVector(pos=[7.0e6, 0.0, 0.0], vel=[0.0, 7.5e3, 0.0], mass=50.0, q=[0.0, 0.0, 0.0, 1.0], ω=zeros(3))
        _, torque, _, _ = readout(p, 1, x)
        truth = SM.rod_field_body(SM.MagneticTorqueRodModel(), SVector(7.0e6, 0.0, 0.0), SVector(0.0, 0.0, 0.0, 1.0), p)
        @test torque ≈ cross(manager.held_dipole_am2, truth) rtol=1e-12          # the run's field, not the estimate
        bound = SpaceAGORA.bind_spacecraft(manager, 2)
        @test bound.mounted_levels !== manager.mounted_levels && isempty(bound.mounted_levels)
        rodless = SM.MagneticMomentumManagerModel(spacecraft_index=1, h_wheels_0=SVector(0.02, 0.0, 0.0),
            commanded_torque=(t, r, v, q, w) -> ZERO3, b_field_ii=estimate, actuation=:mounted)
        @test_throws ArgumentError CH.update_control!(rodless, u,
            readout_params(SM.SpacecraftModel[spacecraft_with(SM.Link(root=true, m=50.0))], (rodless,)), 0.0, 1)
        @test_throws ArgumentError SpaceAGORA.actuation_mode(SM.MagneticMomentumManagerModel(actuation=:bogus))
    end

    @testset "a laser terminal recoils with its transmitted power over c" begin
        c = 299_792_458.0
        root = SM.Link(root=true, m=30.0)
        laser = LaserTerminal(boresight=SVector(2.0, 0.0, 0.0), location=SVector(0.0, 0.5, 0.0), max_power_w=1.0e4,
            magnification=50.0, range_m=1.0e5, efficiency=0.4)
        @test laser.boresight == SVector(1.0, 0.0, 0.0)
        force, torque = actuator_wrench(laser, 3.0e3, actuator_context(root))
        @test force == SVector(-3.0e3 / c, 0.0, 0.0)                         # opposite the beam
        @test torque == cross(SVector(0.0, 0.5, 0.0), force)
        @test actuator_consumption(laser, 3.0e3, actuator_context(root)) == (mass_rate_kg_s = 0.0, power_w = 3.0e3 / 0.4)
        @test actuator_limits(laser) == (0.0, 1.0e4) && actuator_fraction(laser, 2.5e3) == 0.25
        @test_throws ArgumentError LaserTerminal(boresight=ZERO3, max_power_w=1.0)
        @test_throws ArgumentError LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=1.0, efficiency=1.5)
        link_parameters = SM.InterLinkParameters(laser)
        @test (link_parameters.P, link_parameters.B, link_parameters.range) == (1.0e4, 50.0, 1.0e5)
        mixed = spacecraft_with(SM.Link(root=true, m=30.0))
        rod = add_actuator!(mixed, mixed.root, Magnetorquer(axis=SVector(0.0, 0.0, 1.0)))
        puck = add_actuator!(mixed, mixed.root, ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0))
        beam = add_actuator!(mixed, mixed.root, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=10.0))
        @test mixed.root.laser_terminals == [beam] && SM.mounted_laser_count(mixed) == 1
        @test [e.actuator for e in mounted_actuators(mixed)] == AbstractActuator[rod, beam, puck]   # rods, lasers, user kinds
        sc = spacecraft_with(root)
        add_actuator!(sc, root, laser)
        p = readout_params(SM.SpacecraftModel[sc], (HoldLevels(1, [2.0e4]),))   # saturates at 1e4 W
        q = normalize(SVector(0.0, 0.0, sin(0.3), cos(0.3)))
        f_ii, tq_body, _, rate = readout(p, 1, body_state(q))
        @test f_ii ≈ SM.rot(q)' * SVector(-1.0e4 / c, 0.0, 0.0) rtol=1e-14
        @test tq_body == cross(SVector(0.0, 0.5, 0.0), SVector(-1.0e4 / c, 0.0, 0.0))
        @test rate == 0.0                                                  # no propellant
    end

    @testset "a hinge turns its link within its angle and rate limits at control ticks" begin
        root = SM.Link(root=true, m=50.0)
        panel = SM.Link(m=5.0, r=MVector(0.0, 1.5, 0.0))
        sc = SM.SpacecraftModel(links=[root, panel], root=root, prop_mass=1.0, inertia_tensor=root.inertia, id=5)
        hinge = PanelHinge(axis=SVector(0.0, 2.0, 0.0), min_angle_rad=-1.0, max_angle_rad=1.0, max_rate_rad_s=0.1, angle_rad=0.2)
        @test hinge.axis == SVector(0.0, 1.0, 0.0)
        @test add_actuator!(sc, panel, hinge) === hinge
        @test panel.q ≈ [0.0, sin(0.1), 0.0, cos(0.1)]                       # the link takes the hinge's angle
        @test_throws ArgumentError add_actuator!(sc, root, PanelHinge(axis=SVector(1.0, 0.0, 0.0)))   # not the root
        @test_throws ArgumentError add_actuator!(sc, panel, PanelHinge(axis=SVector(1.0, 0.0, 0.0)))  # one per link
        @test_throws ArgumentError PanelHinge(axis=SVector(1.0, 0.0, 0.0), min_angle_rad=1.0, max_angle_rad=0.0)
        @test actuator_wrench(hinge, 0.5, actuator_context(panel)) == (ZERO3, ZERO3)
        @test actuator_limits(hinge) == (-1.0, 1.0)
        @test SM.hinge_step(hinge, 5.0, 0.0) == 0.2                         # first tick: holds, rate-limited
        SM.apply_hinge_levels!(sc, [5.0], 0.0)
        @test hinge.angle_rad == 0.2 && hinge.last_update_s == 0.0
        SM.apply_hinge_levels!(sc, [5.0], 2.0)
        @test hinge.angle_rad ≈ 0.4                                          # 0.1 rad/s for 2 s toward the limit
        SM.apply_hinge_levels!(sc, [5.0], 100.0)
        @test hinge.angle_rad == 1.0                                          # the angle limit
        @test panel.q ≈ [0.0, sin(0.5), 0.0, cos(0.5)]
        @test SM.hinge_step(PanelHinge(axis=SVector(1.0, 0.0, 0.0)), 0.7, 0.0) == 0.7   # unlimited: at once
        @test [e.actuator for e in mounted_actuators(sc)] == AbstractActuator[hinge]
        # After the commanding effector's tick the engine applies its hinge levels.
        control = HoldLevels(1, [0.3])
        p = readout_params(SM.SpacecraftModel[sc], (control,))
        SM.SimulationCallbacks._apply_mounted_hinges!(control, p, 1, 110.0)
        @test hinge.angle_rad ≈ 0.3
        @test SM.SimulationCallbacks._apply_mounted_hinges!(SM.MagneticMomentumManagerModel(), p, 1, 120.0) === nothing
        @test hinge.angle_rad ≈ 0.3                                          # ideal effectors leave hinges alone
    end

    @testset "a mounted run moves the hinge at the controller's ticks" begin
        root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
        panel = SM.Link(m=5.0, r=MVector(0.0, 1.5, 0.0), dims=MVector(0.02, 2.0, 1.0))
        ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), ZERO3; q=SVector(0.0, 0.0, 0.0, 1.0))
        sc = SM.SpacecraftModel(links=[root, panel], root=root, prop_mass=1.0, inertia_tensor=root.inertia, initial_condition=ic, id=6)
        hinge = add_actuator!(sc, panel, PanelHinge(axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.1))
        base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=3.0,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
            keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM.SimConfig._with_configuration(base;
            control_model=SM.ControlModel(control_effectors=(HoldLevels(1, [0.5]),), control_rates=[1.0]))
        run_simulation(args; isolate_state=false)
        @test 0.1 - 1e-12 <= hinge.angle_rad <= 0.3 + 1e-12                  # 0.1 rad/s from the first tick
        @test hinge.angle_rad < 0.5                                          # still on its way
        @test panel.q ≈ [0.0, sin(hinge.angle_rad / 2), 0.0, cos(hinge.angle_rad / 2)]
    end

    @testset "redundant mounted wheels each integrate their own command" begin
        # A pure null-space command: the four torques cancel on the body, yet
        # every wheel spins up at minus its own level.
        null_cmd = [0.01, 0.01, -0.01, -0.01]
        final, inertia, J = pyramid_run(null_cmd)
        @test collect(final.h_wheels) ≈ -5.0 .* null_cmd atol=1e-12
        @test collect(final.ω) ≈ zeros(3) atol=1e-14
        # A mixed command: a body torque about x plus a null-space part.
        mixed = [0.02, 0.0, -0.01, -0.01]
        final, inertia, J = pyramid_run(mixed)
        @test collect(final.h_wheels) ≈ -5.0 .* mixed atol=1e-12
        @test collect(final.ω) ≈ inertia \ (J * mixed .* 5.0) rtol=1e-8
        @test inertia * collect(final.ω) + J * collect(final.h_wheels) ≈ zeros(3) atol=1e-12   # angular momentum conserved
    end

    @testset "actuators mount only on the spacecraft's own links" begin
        root = SM.Link(root=true, m=10.0)
        sc = spacecraft_with(root)
        stray = SM.Link(m=1.0)                                       # not in sc.links
        for actuator in (Thruster(direction=MVector(1.0, 0.0, 0.0)), Magnetorquer(axis=SVector(1.0, 0.0, 0.0)),
                         LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=1.0), PanelHinge(axis=SVector(1.0, 0.0, 0.0)),
                         ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0))
            @test_throws ArgumentError add_actuator!(sc, stray, actuator)
        end
        @test isempty(stray.thrusters) && isempty(stray.magnetorquers) && isempty(stray.laser_terminals) &&
              isempty(stray.hinges) && isempty(stray.custom_actuators)
        @test stray.q == [0.0, 0.0, 0.0, 1.0]                           # a refused hinge leaves the link alone
        @test SM.mounted_actuator_count(sc) == 0
    end

    @testset "a mounted run moves the vehicle and uses propellant as its jets say" begin
        root = SM.Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
        ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
        sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=10.0, inertia_tensor=root.inertia, initial_condition=ic, id=3)
        add_actuator!(sc, root, ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0))   # through the centre: no torque
        thrust_n, duration_s = 0.1, 100.0
        control = HoldLevels(1, [thrust_n])
        base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=duration_s,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
            keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM.SimConfig._with_configuration(base;
            control_model=SM.ControlModel(control_effectors=(control,), control_rates=[10.0]))
        result = run_simulation(args; return_solution=true)
        m0 = sc.dry_mass + sc.prop_mass
        mdot = thrust_n / (50.0 * G0)
        final = result.u[end].sc[1]
        @test result.t[end] ≈ duration_s
        @test final.mass ≈ m0 - mdot * duration_s rtol=1e-9
        @test final.vel[1] ≈ (thrust_n / mdot) * log(m0 / (m0 - mdot * duration_s)) rtol=1e-6   # rocket equation
        @test abs(final.vel[2]) < 1e-12 && abs(final.vel[3]) < 1e-12
    end

    @testset "actuator requirements combine and reach the actuators' context" begin
        @test actuator_requirements(ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0)) == ActuatorRequirements()
        @test actuator_requirements(Magnetorquer(axis=SVector(1.0, 0.0, 0.0))) == ActuatorRequirements(magnetic_field=true)
        @test (ActuatorRequirements(planet_frame=true) | ActuatorRequirements(magnetic_field=true)) ==
              ActuatorRequirements(planet_frame=true, magnetic_field=true)
        sc, root = far_spacecraft(id=8)
        seen = Ref{Any}(:never)
        add_actuator!(sc, root, EnvironmentRecorder(seen, ActuatorRequirements(planet_frame=true, atmosphere=true)))
        add_actuator!(sc, root, Magnetorquer(axis=SVector(1.0, 0.0, 0.0)))
        @test SM.mounted_requirements(sc) == ActuatorRequirements(planet_frame=true, atmosphere=true, magnetic_field=true)
        held_run(sc, [0.0, 0.0])
        env = seen[]
        @test env isa EnvironmentSample
        @test env.planet_frame isa PlanetFrameSample
        @test env.planet_frame.alt_m ≈ 4.0e9 - SM.Earth().Rp_e rtol=1e-3
        @test env.atmosphere isa AtmosphereSample && env.atmosphere.rho_kg_m3 == 0.0      # the run has no atmosphere
        @test env.solar === nothing                                                      # not asked for
        # A kind that asks for nothing gets no environment.
        sc2, root2 = far_spacecraft(id=9)
        plain = Ref{Any}(:never)
        add_actuator!(sc2, root2, EnvironmentRecorder(plain, ActuatorRequirements()))
        held_run(sc2, [0.0])
        @test plain[] === nothing
    end

    @testset "a thruster with firing dynamics delivers what its control tick scheduled" begin
        root = SM.Link(root=true, m=50.0)
        sc = spacecraft_with(root)
        jet = add_actuator!(sc, root, Thruster(max_thrust=2.0, direction=MVector(1.0, 0.0, 0.0),
            min_firing_time=0.2, cutoff_frequency=1.0e6, firing_dynamics=true))
        plain = add_actuator!(sc, root, Thruster(max_thrust=2.0, direction=MVector(0.0, 1.0, 0.0)))
        control = HoldLevels(1, [1.0, 1.0])
        p = readout_params(SM.SpacecraftModel[sc], (control,))
        state = body_state(SVector(0.0, 0.0, 0.0, 1.0); mass=60.0)
        @test readout(p, 1, state)[1] == SVector(0.0, 1.0, 0.0)          # the jet is idle until its first tick
        SM.SimulationCallbacks._apply_mounted_thruster_firing!(control, p, 1, 0.0, 1.0)
        @test jet.thrust ≈ 1.0 rtol=1e-5                                  # on at 2 N for half the period, fast ramp
        @test plain.thrust == 0.0                                         # a thruster without firing dynamics is untouched
        @test readout(p, 1, state)[1] == SVector(jet.thrust, 1.0, 0.0)
        @test SM.mounted_thruster_fractions(sc, [1.0, 1.0]) == [jet.thrust / 2.0, 0.5]
        @test applied_actuator_levels(sc, [1.0, 1.0]) == [jet.thrust, 1.0]
        control.levels .= [0.2, 1.0]                                      # 0.1 s on: under the 0.2 s minimum
        SM.SimulationCallbacks._apply_mounted_thruster_firing!(control, p, 1, 1.0, 1.0)
        @test jet.thrust < 1e-6
        held = jet.thrust
        @test SM.SimulationCallbacks._apply_mounted_thruster_firing!(SM.MagneticMomentumManagerModel(), p, 1, 2.0, 1.0) === nothing
        @test jet.thrust == held                                          # ideal effectors leave the jets alone
        # In a run, a command under the minimum on-time never fires; without firing dynamics it acts at once.
        function firing_run(dynamics, level)
            sc, root = far_spacecraft(id=10, mass=500.0, prop_mass=10.0)
            add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0), Isp=200.0,
                min_firing_time=0.5, cutoff_frequency=1.0e6, firing_dynamics=dynamics))
            return held_run(sc, [level]; mission_time=10.0, return_solution=true).u[end].sc[1], sc.dry_mass + sc.prop_mass
        end
        idle, m0 = firing_run(true, 0.1)
        direct, _ = firing_run(false, 0.1)
        @test collect(idle.vel) == zeros(3) && idle.mass == m0
        @test direct.vel[1] > 0.0 && direct.mass < m0
        # At full command the ticks fire the jet: at most the first period is lost.
        fired, _ = firing_run(true, 1.0)
        full, _ = firing_run(false, 1.0)
        @test 0.85 * full.vel[1] <= fired.vel[1] <= full.vel[1] * (1 + 1e-9)
        @test (m0 - fired.mass) <= (m0 - full.mass) * (1 + 1e-9)
    end

    @testset "mounted runs write every actuator's applied level and the actuators' power" begin
        root = SM.Link(root=true, m=50.0)
        sc = spacecraft_with(root)
        add_actuator!(sc, root, Thruster(max_thrust=2.0, direction=MVector(1.0, 0.0, 0.0)))
        add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0), max_torque_nm=0.01, max_momentum_nms=0.1))
        add_actuator!(sc, root, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=10.0, efficiency=0.25))
        add_actuator!(sc, root, ColdGasPuck(SVector(0.0, 1.0, 0.0), ZERO3, 0.5, 50.0))
        commands = [3.0, 0.05, 4.0, 1.0]
        @test applied_actuator_levels(sc, commands) == [2.0, 0.01, 4.0, 0.5]                 # each within its limits
        @test applied_actuator_levels(sc, commands; wheel_momentum=[-0.1]) == [2.0, 0.0, 4.0, 0.5]   # a full wheel stops
        @test mounted_power(sc, commands) == 16.0                                           # 4 W optical at 25 %
        p = readout_params(SM.SpacecraftModel[sc], (HoldLevels(1, commands),))
        counts = SM.SimulationCallbacks.actuator_level_counts(p.args)
        @test counts == [4]
        @test SM.SimulationCallbacks._save_actuator_levels(1, counts, nothing, 0.0, (p=p,)) == [[2.0, 0.01, 4.0, 0.5]]
        @test SM.SimulationCallbacks._save_actuator_power(1, counts, nothing, 0.0, (p=p,)) == [16.0]
        names(args) = [field.name for field in SM.SimulationCallbacks.default_save_fields(args)]
        @test :actuator_level in names(p.args) && :actuator_power in names(p.args)
        ideal = readout_params(SM.SpacecraftModel[sc], (SM.MagneticMomentumManagerModel(),))
        @test SM.SimulationCallbacks.actuator_level_counts(ideal.args) == [0]
        @test !(:actuator_level in names(ideal.args)) && !(:actuator_power in names(ideal.args))
        # In a run, the columns carry the held levels.
        sc2, root2 = far_spacecraft(id=11)
        add_actuator!(sc2, root2, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=10.0, efficiency=0.5))
        add_actuator!(sc2, root2, ColdGasPuck(SVector(1.0, 0.0, 0.0), ZERO3, 0.5, 50.0))
        table = held_run(sc2, [4.0, 1.0]; mission_time=3.0, return_results=true).table
        @test size(table, 1) > 1
        @test all(==(4.0), table.sc1_actuator_level_1)
        @test all(==(0.5), table.sc1_actuator_level_2)
        @test all(==(8.0), table.sc1_actuator_power)
    end

    @testset "wheels and torque rods draw power" begin
        wheel = ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0), spin_inertia_kgm2=0.02, base_power_w=1.5,
            motor_efficiency=0.8, regenerative_efficiency=0.25)
        @test SM.wheel_power(wheel, -0.01, 2.0) ≈ 1.5 + 0.01 * 100.0 / 0.8      # speeding up: 1 W into the wheel at 100 rad/s
        @test SM.wheel_power(wheel, 0.01, 2.0) ≈ 1.5 - 0.25 * 1.0               # slowing down: a quarter of 1 W returns
        @test SM.wheel_power(wheel, 0.01, NaN) == 1.5                            # stored momentum unknown: base power
        @test SM.wheel_power(ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)), 0.01, 2.0) == 0.0   # no power model by default
        ctx = actuator_context(SM.Link(root=true))
        @test actuator_consumption(wheel, 0.01, ctx) == (mass_rate_kg_s = 0.0, power_w = 1.5)
        rod = Magnetorquer(axis=SVector(1.0, 0.0, 0.0), max_dipole_am2=2.0, power_at_max_dipole_w=0.8)
        @test actuator_consumption(rod, 1.0, ctx).power_w ≈ 0.2 && actuator_consumption(rod, -2.0, ctx).power_w ≈ 0.8
        @test actuator_consumption(Magnetorquer(axis=SVector(1.0, 0.0, 0.0)), 5.0, ctx).power_w == 0.0
        @test_throws ArgumentError Magnetorquer(axis=SVector(1.0, 0.0, 0.0), power_at_max_dipole_w=1.0)   # needs a finite dipole limit
        @test_throws ArgumentError ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0), motor_efficiency=0.0)
        @test_throws ArgumentError ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0), spin_inertia_kgm2=-1.0)
        @test ReactionWheel(SVector(1.0, 0.0, 0.0), 0.1, 0.2).max_momentum_nms == 0.2 && Magnetorquer(SVector(1.0, 0.0, 0.0), 3.0).max_dipole_am2 == 3.0
        root = SM.Link(root=true, m=50.0)
        sc = spacecraft_with(root)
        add_actuator!(sc, root, wheel)
        add_actuator!(sc, root, rod)
        @test mounted_power(sc, [-0.01, 1.0]; wheel_momentum=[2.0]) ≈ 1.5 + 0.01 * 100.0 / 0.8 + 0.2
        @test mounted_power(sc, [-0.01, 1.0]) ≈ 1.5 + 0.2
        # In a run the column follows the spin-up: h = -0.002 t, so 1 W plus 0.0004 t.
        sc2, root2 = far_spacecraft(id=12)
        add_actuator!(sc2, root2, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0), spin_inertia_kgm2=0.01, base_power_w=1.0))
        table = held_run(sc2, [0.002]; mission_time=5.0, return_results=true).table
        @test table.sc1_actuator_power ≈ 1.0 .+ 0.0004 .* table.time rtol=1e-9
    end

    @testset "energy-depletion panels move through a mounted hinge" begin
        function panel_craft(hinge)
            root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
            panel = SM.Link(m=5.0, r=MVector(0.0, 1.5, 0.0), dims=MVector(0.02, 2.0, 1.0))
            sc = SM.SpacecraftModel(links=[root, panel], root=root, prop_mass=1.0, inertia_tensor=root.inertia, id=1)
            hinge === nothing || add_actuator!(sc, panel, hinge)
            return sc, panel
        end
        aoa = SM.SolarPanelAngleOfAttackControlModel(controlled_panel_links=(2,))
        apply!(sc, alpha, t) = SM.ControlHooks._apply_solar_panel_aoa!(aoa, sc, alpha, t)
        plain, plain_panel = panel_craft(nothing)
        hinged, hinged_panel = panel_craft(PanelHinge(axis=SVector(0.0, 1.5, 0.0)))   # the axis the panels turn about, no limits
        for (t, alpha) in ((0.0, 0.3), (1.0, 1.1), (2.0, 0.05))
            apply!(plain, alpha, t)
            apply!(hinged, alpha, t)
            @test hinged_panel.q == plain_panel.q && hinged_panel.α === plain_panel.α === alpha   # bit for bit
        end
        limited, limited_panel = panel_craft(PanelHinge(axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.1,
            min_angle_rad=0.0, max_angle_rad=π / 2, angle_rad=π / 2 - 0.3))
        apply!(limited, 0.3, 10.0)                       # the first tick holds the angle and records the time
        apply!(limited, 1.2, 12.0)                       # commanded to π/2 - 1.2; 0.1 rad/s for 2 s
        @test limited_panel.hinges[1].angle_rad ≈ π / 2 - 0.3 - 0.2
        @test limited_panel.α ≈ 0.5                      # the angle of attack it reached, not the command
        @test limited_panel.q ≈ [0.0, sin((π / 2 - 0.5) / 2), 0.0, cos((π / 2 - 0.5) / 2)]
        SM.ControlHooks._apply_solar_panel_aoa!(aoa, plain, 0.7)          # without a hinge the tick time is not needed
        @test plain_panel.α == 0.7
        @test_throws ArgumentError SM.ControlHooks._apply_solar_panel_aoa!(aoa, limited, 0.7)
    end

    @testset "a mounted burn points its engine with the wheels and fires it in the window" begin
        @test_throws ArgumentError SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[0.0],
            stop_burn_time=[1.0], Isp=[300.0], actuation=:bogus)
        legacy = SM.BaseThrusterModel([1.0], [0.0], [0.0], [0.0], [1.0], [300.0])
        @test SpaceAGORA.actuation_mode(legacy) === :ideal && SpaceAGORA.actuator_levels(legacy, 1) === nothing
        mounted = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[0.0],
            stop_burn_time=[1.0], Isp=[300.0], actuation=:mounted)
        @test SM.ControlHooks.control_thruster_levels(mounted, 1) === nothing      # the mounted engine reports its own level
        @test SM.ControlHooks._thruster_levels_at(mounted, 1, 0.5) === nothing
        # The engine switches only at the window's ends, through the flag the edge callback flips.
        flag_model = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[300.0],
            stop_burn_time=[400.0], Isp=[300.0], actuation=:mounted)
        p_flag = readout_params(SM.SpacecraftModel[orbit_spacecraft(id=25)], (flag_model,))
        for (t, on) in ((299.999, false), (300.0, true), (350.0, true), (400.0, false))
            SM.ControlHooks.refresh_burn!(flag_model, p_flag, 1, t)
            @test flag_model.burning[1] == on && SpaceAGORA.actuator_levels(flag_model, 1)[1] == (on ? 1.0 : 0.0)
        end
        @test SM.ControlHooks.burn_edge(flag_model, 300.0) && SM.ControlHooks.burn_edge(flag_model, 400.0)
        @test !SM.ControlHooks.burn_edge(flag_model, 350.0) && !SM.ControlHooks.burn_edge(legacy, 0.0)

        ideal = burn_run(orbit_spacecraft(id=21), :ideal)
        pointed = burn_run(orbit_spacecraft(id=22), :mounted)
        unpointed = burn_run(orbit_spacecraft(id=23; wheels=false), :mounted)
        speed(result) = norm(result.solution.u[end].sc[1].vel)
        # The engine runs exactly the window: 1 N for 100 s at 300 s of specific impulse; so does the ideal burn.
        @test 110.0 - pointed.solution.u[end].sc[1].mass ≈ 100.0 / (300.0 * G0) rtol=1e-6
        @test 110.0 - ideal.solution.u[end].sc[1].mass ≈ 100.0 / (300.0 * G0) rtol=1e-6
        # A burn already on when the run starts: the first step thrusts too.
        at_start = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[0.0],
            stop_burn_time=[100.0], Isp=[300.0])
        start_base = make_example_config(planet=SM.Earth(), spacecraft=orbit_spacecraft(id=26), mission_time=150.0,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(SM.InverseSquaredGravityModel(),), ephemerides_model=SM.SimpleEphemeridesModel(),
            orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
        start_args = SM.SimConfig._with_configuration(start_base;
            control_model=SM.ControlModel(control_effectors=(at_start,), control_rates=[1.0]))
        @test 110.0 - run_simulation(start_args; return_solution=true).u[end].sc[1].mass ≈ 100.0 / (300.0 * G0) rtol=1e-6
        ideal_mass_at(t) = ideal.solution(t).sc[1].mass
        @test ideal_mass_at(300.5) == ideal_mass_at(300.0) && ideal_mass_at(300.9) < ideal_mass_at(300.5) - 0.3 / (300.0 * G0)
        @test isapprox(ideal_mass_at(449.0), ideal_mass_at(400.5); atol=1e-12)
        # ... and lights and stops at the window's ends, between the control ticks.
        mass_at(t) = pointed.solution(t).sc[1].mass
        mdot = 1.0 / (300.0 * G0)
        @test mass_at(300.5) == mass_at(300.0) && mass_at(300.9) < mass_at(300.5) - 0.3 * mdot
        @test isapprox(mass_at(400.9), mass_at(400.5); atol=1e-12) && mass_at(400.5) < mass_at(400.1) - 0.3 * mdot
        # Pointed: the speed gained matches the ideal burn's to a small fraction of the burn's Δv.
        dv_burn = 1.0 / 110.0 * 100.0
        @test abs(speed(pointed) - speed(ideal)) < 0.01 * dv_burn
        # Not pointed: the engine fires 90 degrees off, and the orbit does not speed up the same way.
        @test abs(speed(unpointed) - speed(ideal)) > 0.2 * dv_burn
        # At ignition the engine (body +x) is on the velocity.
        k = findfirst(t -> t >= 300.5, pointed.solution.t)
        state = pointed.solution.u[k].sc[1]
        engine_ii = SM.rot(SVector{4}(state.q))' * SVector(1.0, 0.0, 0.0)
        @test acos(clamp(dot(engine_ii, normalize(SVector{3}(state.vel))), -1.0, 1.0)) < deg2rad(0.5)
        # The engine's level is the burn's thrust in the window and zero outside it.
        table = pointed.table
        inside = (table.time .> 300.5) .& (table.time .< 400.5)
        @test any(inside) && all(==(1.0), table.sc1_actuator_level_1[inside])
        @test all(==(0.0), table.sc1_actuator_level_1[.!inside .& ((table.time .< 300.5) .| (table.time .> 400.5))])
        # The ideal burn's saved firing level is the burn's state at each sample's
        # own time: 1 exactly in [start, stop), though the scheduler itself ticks
        # on the 1 s control grid. Samples are saved after the callbacks of
        # the step they fall in; here the steps that end at the window's ends each
        # hold a sample, which must still read the level at its own time.
        scene_run = burn_run(orbit_spacecraft(id=27), :ideal; scene=true, sample_period=0.25)
        sample_times, step_ends = scene_run.table.time, scene_run.solution.t
        step_start(edge) = maximum(filter(<(edge), step_ends))
        @test 300.5 in step_ends && 400.5 in step_ends
        @test any(step_start(300.5) .< sample_times .< 300.5) && any(step_start(400.5) .< sample_times .< 400.5)
        @test scene_run.table.sc1_thruster_level_1 == Float64.(300.5 .<= sample_times .< 400.5)
    end

    @testset "the saved burn level and actuator faults" begin
        # A mounted burn with a failed engine reports what the fault leaves (the
        # mounted fractions); an ideal burn does not fly the link thrusters, so a
        # fault there leaves its level, read at the sample's time, unchanged.
        burn_model(actuation) = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[300.0],
            stop_burn_time=[400.0], Isp=[300.0], actuation=actuation)
        faulted(id) = (sc = orbit_spacecraft(id=id);
            add_fault!(sc, ActuatorFault(sc.root.thrusters[1]; kind=:failed, start_s=0.0));
            SM.ActuatorFaults.resolve_faults!(sc); sc)
        mounted_sc, mounted = faulted(28), burn_model(:mounted)
        SM.ControlHooks.refresh_burn!(mounted, readout_params(SM.SpacecraftModel[mounted_sc], (mounted,)), 1, 350.0)
        @test mounted.burning[1] && SpaceAGORA.actuator_levels(mounted, 1)[1] == 1.0
        SM.ActuatorFaults.refresh_faults!(mounted_sc, mounted.mounted_levels[1], 350.0)
        @test SM.SimulationCallbacks._effector_thruster_levels(mounted, 1, mounted_sc, 350.0) == [0.0]
        ideal_sc, ideal = faulted(29), burn_model(:ideal)
        SM.ControlHooks.refresh_burn!(ideal, readout_params(SM.SpacecraftModel[ideal_sc], (ideal,)), 1, 0.0)
        @test SM.SimulationCallbacks._effector_thruster_levels(ideal, 1, ideal_sc, 350.0) == [1.0]
        @test SM.SimulationCallbacks._effector_thruster_levels(ideal, 1, ideal_sc, 400.0) == [0.0]
    end

    @testset "a mounted burn releases the engine and wheels of an inactive spacecraft" begin
        # Two spacecraft burn from 300.5 s to 400.5 s; the second becomes inactive at
        # 350 s, mid-burn, and the run ends at 380 s, with both windows open. From then
        # on its mounted tick is released (`skip_inactive`): engine and wheels off.
        model = SM.BaseThrusterModel(thrust=[1.0, 1.0], direction=[0.0, 0.0], Δv=[0.0, 0.0],
            start_burn_time=[300.5, 300.5], stop_burn_time=[400.5, 400.5], Isp=[300.0, 300.0], actuation=:mounted)
        sc1, sc2 = orbit_spacecraft(id=41), orbit_spacecraft(id=42)
        base = make_example_config(planet=SM.Earth(), spacecraft=sc1, mission_time=380.0,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(SM.InverseSquaredGravityModel(),), ephemerides_model=SM.SimpleEphemeridesModel(),
            orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM.SimConfig._with_configuration(base;
            dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sc1, sc2], (SM.InverseSquaredGravityModel(),)),
            control_model=SM.ControlModel(control_effectors=(model,), control_rates=[1.0]))
        deactivate = PresetTimeCallback([350.0], integrator -> (integrator.p.is_active[2] = false);
            save_positions=(false, false))
        sol = run_simulation(args; isolate_state=false, return_solution=true, extra_callbacks=(deactivate,))
        @test sol.prob.p.is_active == [true, false] && sol.t[end] == 380.0
        @test model.burning == [true, true]                    # the release leaves the schedule
        # The active spacecraft's tick goes on: engine at the burn's thrust, wheels pointing it.
        on = model.mounted_levels[1]
        @test length(on) == 4 && on[1] == 1.0 && any(!iszero, on[2:4])
        @test 110.0 - sol.u[end].sc[1].mass ≈ 79.5 / (300.0 * G0) rtol=1e-6
        # The inactive one's tick is released: engine and wheels off. Its state stays
        # where the deactivation left it, after 49.5 s of the burn.
        off = model.mounted_levels[2]
        @test length(off) == 4 && all(iszero, off)
        @test 110.0 - sol.u[end].sc[2].mass ≈ 49.5 / (300.0 * G0) rtol=1e-6
        # The release on its own: the mounted engine and wheels off, the schedule kept.
        @test release_command!(model, 380.0, 1) === nothing && all(iszero, model.mounted_levels[1]) && model.burning[1]
        @test release_command!(model, 380.0, 3) === nothing
        ideal = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[300.5],
            stop_burn_time=[400.5], Isp=[300.0])
        @test release_command!(ideal, 0.0, 1) === nothing && all(isempty, ideal.mounted_levels)
    end

    @testset "no burn edge or scheduling pass undoes a release" begin
        # The spacecraft becomes inactive at 0.5 s, before its window opens at 1.5 s,
        # and active again at 2.5 s; the scheduler ticks every second. The ignition
        # edge at 1.5 s records the burn but leaves the released engine off, as does
        # the tick at 2 s. After the reactivation the next tick, at 3 s, fires it.
        model = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[1.5],
            stop_burn_time=[4.5], Isp=[300.0], actuation=:mounted)
        base = make_example_config(planet=SM.Earth(), spacecraft=orbit_spacecraft(id=43), mission_time=4.0,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(SM.InverseSquaredGravityModel(),), ephemerides_model=SM.SimpleEphemeridesModel(),
            orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM.SimConfig._with_configuration(base;
            control_model=SM.ControlModel(control_effectors=(model,), control_rates=[1.0]))
        switch = PresetTimeCallback([0.5, 2.5], integrator -> (integrator.p.is_active[1] = integrator.t >= 2.5);
            save_positions=(false, false))
        # After the engine's own callbacks at these times: the commanded levels.
        levels = Dict{Float64, Vector{Float64}}()
        record = PresetTimeCallback([1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
            integrator -> (levels[integrator.t] = copy(model.mounted_levels[1])); save_positions=(false, false))
        sol = run_simulation(args; isolate_state=false, return_solution=true, extra_callbacks=(switch, record))
        @test sol.prob.p.is_active == [true] && sol.t[end] == 4.0
        @test all(t -> all(iszero, levels[t]), (1.0, 1.5, 2.0, 2.5))     # released across the edge and the tick
        @test levels[3.0][1] == 1.0 && levels[3.5][1] == 1.0                # fires once active again
        # The engine ran from 3 s to the end of the run, one second of the window.
        @test 110.0 - sol.u[end].sc[1].mass ≈ 1.0 / (300.0 * G0) rtol=1e-6

        # Each path on its own: scheduling (at the run's start, or when the spacecraft
        # leaves the atmosphere) and a burn edge keep the schedule, and leave an inactive
        # spacecraft's engine off and fire an active one's.
        direct = SM.BaseThrusterModel(thrust=[1.0], direction=[0.0], Δv=[0.0], start_burn_time=[1.5],
            stop_burn_time=[4.5], Isp=[300.0], actuation=:mounted)
        p = readout_params(SM.SpacecraftModel[orbit_spacecraft(id=44)], (direct,))
        u = SE.build_initial_conditions(p.args)
        edge! = SM.SimulationCallbacks._burn_edge_callback(direct, 1).affect!
        for active in (false, true)
            p.is_active[1] = active
            SM.SimulationCallbacks._schedule_thruster_control!((u=u, p=p, t=1.6), direct, 1)
            @test direct.burning[1] && direct.mounted_levels[1][1] == (active ? 1.0 : 0.0)
            fill!(direct.mounted_levels[1], 0.0)
            edge!((u=u, p=p, t=1.5))
            @test direct.burning[1] && direct.mounted_levels[1][1] == (active ? 1.0 : 0.0)
        end
    end

    @testset "two burn runs at once match two in sequence" begin
        # The scheduler keeps its run state in the model, so concurrent runs of
        # one configuration trace and fly exactly as sequential ones do.
        mktempdir() do dir
            traced(runs, path) = withenv(runs, "SPACEAGORA_MANEUVER_TRACE_CSV" => path)
            sequential_csv, concurrent_csv = joinpath(dir, "sequential.csv"), joinpath(dir, "concurrent.csv")
            sequential = traced(sequential_csv) do
                [burn_run(orbit_spacecraft(id=31), :ideal; scene=true) for _ in 1:2]
            end
            concurrent = traced(concurrent_csv) do
                fetch.([Threads.@spawn(burn_run(orbit_spacecraft(id=31), :ideal; scene=true)) for _ in 1:2])
            end
            # Concurrent first writes may each add the header line; the rows match.
            rows(path) = sort(filter(!startswith("event,"), readlines(path)))
            @test length(rows(sequential_csv)) == 6      # each run's burn_start, burn_end and next-tick cleanup
            @test rows(concurrent_csv) == rows(sequential_csv)
            reference = sequential[1]
            @test all(isequal(r.solution.u[end], reference.solution.u[end]) && isequal(r.table, reference.table)
                for r in (sequential[2], concurrent...))
        end
    end
end

end # module ActuatorContractTests
