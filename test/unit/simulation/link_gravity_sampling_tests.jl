module LinkGravitySamplingTests

using Test
using SpaceAGORA
using StaticArrays
using ComponentArrays
using LinearAlgebra

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const CRA = SM.ClothRobotArmDynamics

# A position-only test effector that makes a wrong spacecraft handle observable.
struct SpacecraftGravityProbe end
SM.gravity_backbone_structure(::SpacecraftGravityProbe) = :position_only_static_gravity
SM.gravity_backbone_acceleration_ii(::SpacecraftGravityProbe, x, env, t::Float64) =
    SVector(x.spacecraft.gravity_tag, 1.0e-6 * x.pos_ii[2], 0.0)

struct NongravityProbe end
SM.gravity_backbone_acceleration_ii(::NongravityProbe, x, env, t::Float64) =
    error("A non-gravity effector must not be evaluated as link gravity")

function one_link_plan()
    arm = SM.default_cloth_arm_model(
        link_lengths_m=(0.2,), link_radii_m=(0.01,), link_masses_kg=(1.0,),
        joint_axes=((0.0, 0.0, 1.0),),
    )
    base = SM.ClothArmBasePose([0.0, 0.0, 0.0])
    tip = SVector(0.2, 0.0, 0.0)
    return SM.RobotArmPlan(
        arm, base, [0.0, 1.0], zeros(1, 2), zeros(1, 2), zeros(1, 2),
        repeat(reshape(collect(tip), 3, 1), 1, 2), [0.0], [0.0], tip, 0.0,
        :cloth_quintic,
    )
end

function sampling_params(effectors, plan)
    controls = ntuple(2) do i
        SM.RobotArmControlEffector(
            spacecraft_index=i, plan=plan,
            k_translation_n_m=0.0, c_translation_n_s_m=0.0,
            k_rotation_n_m_rad=0.0, c_rotation_n_m_s_rad=0.0,
        )
    end
    return (
        args=(
            environment_model=(planet=make_no_gram_planet(:earth), ephemerides_model=SM.SimpleEphemeridesModel()),
            dynamics_model=(spacecraft=[(gravity_tag=1.0,), (gravity_tag=3.0,)], dynamic_effectors=effectors),
            control_model=(control_effectors=controls,),
        ),
        shared_buffers=(
            et_start=Ref(0.0),
            spice_runtime_counters=(planet_pxform_runtime_calls=Threads.Atomic{Int64}(0),),
            planet_frame_ephemeris_cache=Ref{Any}(nothing),
            robot_arm_present=Ref{Union{Nothing, Bool}}(nothing),
        ),
    )
end

@testset "Arm and articulated gravity preserve spacecraft context" begin
    plan = one_link_plan()
    p = sampling_params((SM.InverseSquaredGravityModel(), SpacecraftGravityProbe(), NongravityProbe()), plan)
    sc = ComponentVector(merge(
        (pos=[6.778e6, 1.0e5, 0.0], vel=[0.0, 7670.0, 0.0],
         q=[0.0, 0.0, 0.0, 1.0], ω=zeros(3)),
        CRA.coupled_cloth_robot_arm_state_shape(plan),
    ))
    CRA.initialize_coupled_cloth_robot_arm_state!(sc, plan)
    before = copy(sc)
    r = SVector{3, Float64}(sc.arm_r[:, 1])
    @test r != SVector{3, Float64}(sc.pos)
    planet = p.args.environment_model.planet
    point_mass = -planet.μ * r / norm(r)^3

    for sat_idx in 1:2
        expected = point_mass + SVector(p.args.dynamics_model.spacecraft[sat_idx].gravity_tag, 1.0e-6 * r[2], 0.0)
        @test SE._gravity_only_acceleration_ii(p, r, 0.0; sat_idx=sat_idx) ≈ expected
        @test SE._ArmLinkGravity(p, 0.0, sat_idx)(r) ≈ expected

        # Exercise the engine call site as well as the callable: forgetting to pass
        # sat_idx when constructing the callable must fail on spacecraft 2.
        du = zero(sc)
        forces, torques = zeros(3), zeros(3)
        SE._apply_coupled_robot_arm_rhs!(du, sc, p, sat_idx, 0.0, forces, torques)
        @test du.arm_v[:, 1] ≈ expected
        @test forces == zeros(3)
        @test torques == zeros(3)
    end
    @test SE._gravity_only_acceleration_ii(p, r, 0.0) ==
        SE._gravity_only_acceleration_ii(p, r, 0.0; sat_idx=1)
    @test sc == before
    @test p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls[] == 0

    empty_p = sampling_params((), plan)
    @test SE._ArmLinkGravity(empty_p, 0.0, 2)(r) == zeros(3)
    duplicate_p = sampling_params((SM.InverseSquaredGravityModel(), SM.InverseSquaredGravityModel()), plan)
    @test SE._ArmLinkGravity(duplicate_p, 0.0, 2)(r) ≈ 2point_mass

    # J2 requests a planet frame; the shared helper must still sample it at the
    # link's position after the arm and articulated implementations are combined.
    j2 = SM.InverseSquaredJ2GravityModel()
    j2_p = sampling_params((j2,), plan)
    sample = SM.EffectorSampling.StateSample(r, zero(r), 1.0;
        spacecraft=j2_p.args.dynamics_model.spacecraft[2])
    frame = SE.sample_planet_frame(sample, j2_p, 2, 0.0)
    env = SM.EffectorSampling.EnvironmentSample(planet; planet_frame=frame)
    expected_j2 = SM.gravity_backbone_acceleration_ii(j2, sample, env, 0.0)
    @test SE._ArmLinkGravity(j2_p, 0.0, 2)(r) ≈ expected_j2
    @test norm(expected_j2 - point_mass) > 1.0e-5
end

end # module
