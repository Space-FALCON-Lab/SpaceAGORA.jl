module CompliantMultibodyNamingTests

using Test, SpaceAGORA, StaticArrays, ComponentArrays

const SM = SpaceAGORA.SimulationModel
const HS = SpaceAGORA.HYPRServices

# Keep this inventory independent of the exported API documentation: deleting an
# old export and its documentation together must not silently pass this contract.
const KINEMATIC_ALIASES = (
    (:RobotArmBasePose, :ClothArmBasePose),
    (:RobotArmLink, :ClothArmLink),
    (:RobotArmJoint, :ClothArmJoint),
    (:RobotArmModel, :ClothArmModel),
    (:RobotArmPose, :ClothArmPose),
    (:RobotArmState, :ClothArmState),
    (:default_robot_arm_model, :default_cloth_arm_model),
    (:robot_arm_fk, :cloth_fk),
    (:robot_arm_fk_state, :cloth_fk_state),
    (:robot_arm_ik, :cloth_ik),
    (:robot_arm_total_reach, :cloth_total_reach),
    (:robot_arm_end_effector_pose, :cloth_end_effector_pose),
)
const DYNAMICS_ALIASES = (
    (:CompliantRobotArmSimulation, :ClothRobotArmSimulation),
    (:CompliantRobotArmReferenceState, :ClothRobotArmReferenceState),
    (:compliant_robot_arm_reference_state, :cloth_reference_state),
    (:compliant_robot_arm_multibody, :cloth_robot_arm_multibody),
    (:compliant_robot_arm_initial_state, :cloth_robot_arm_initial_state),
    (:compliant_robot_arm_rest_quaternions, :cloth_robot_arm_rest_quaternions),
    (:compliant_robot_arm_end_effector, :cloth_robot_arm_end_effector),
    (:compliant_robot_arm_actuators, :cloth_robot_arm_actuators),
    (:simulate_compliant_robot_arm_plan, :simulate_cloth_robot_arm_plan),
    (:coupled_compliant_robot_arm_state_shape, :coupled_cloth_robot_arm_state_shape),
    (:initialize_coupled_compliant_robot_arm_state!, :initialize_coupled_cloth_robot_arm_state!),
    (:assign_coupled_compliant_robot_arm_rhs!, :assign_coupled_cloth_robot_arm_rhs!),
)

@testset "Canonical names retain the historical objects" begin
    for (canonical, historical) in KINEMATIC_ALIASES
        @test canonical in names(SpaceAGORA)
        @test historical in names(SpaceAGORA)
        @test getproperty(SpaceAGORA, canonical) === getproperty(SpaceAGORA, historical)
        @test getproperty(SM, canonical) === getproperty(SM.Robotics, historical)
    end
    for (canonical, historical) in DYNAMICS_ALIASES
        @test canonical in names(SpaceAGORA)
        @test getproperty(SpaceAGORA, canonical) === getproperty(SM.ClothRobotArmDynamics, historical)
        @test getproperty(SM, canonical) === getproperty(SM.ClothRobotArmDynamics, historical)
        # Reference-state helpers were not historical root exports.
        if historical in (:ClothRobotArmReferenceState, :cloth_reference_state)
            continue
        end
        @test historical in names(SpaceAGORA)
        @test getproperty(SpaceAGORA, historical) === getproperty(SpaceAGORA, canonical)
    end
    @test SpaceAGORA.CompliantMultibody === SM.CompliantMultibody === SM.ClothMultibody
    @test SpaceAGORA.CompliantRobotArmDynamics === SM.CompliantRobotArmDynamics === SM.ClothRobotArmDynamics
    @test :CompliantMultibody in names(SpaceAGORA)
    @test :CompliantRobotArmDynamics in names(SpaceAGORA)
    @test !Base.isexported(SpaceAGORA, :CompliantMultibody)
    @test !Base.isexported(SpaceAGORA, :CompliantRobotArmDynamics)
    # Defining owners remain intact for existing dispatch and serialized type paths.
    @test parentmodule(CompliantMultibodyModel) === SM.ClothMultibody
    @test parentmodule(CompliantRobotArmSimulation) === SM.ClothRobotArmDynamics
    @test parentmodule(RobotArmModel) === SM.Robotics
    @test SM.CompliantRobotArmDynamics.CompliantRobotArmDynamicsMode === SM.ClothRobotArmDynamics.ClothRobotArmDynamicsMode === Symbol
end

@testset "Canonical entrypoints build and evaluate a small arm" begin
    arm = default_robot_arm_model(link_lengths_m=(0.2,), link_radii_m=(0.01,),
        link_masses_kg=(1.0,), joint_axes=((0.0, 0.0, 1.0),))
    base = RobotArmBasePose([0.0, 0.0, 0.0])
    @test arm isa RobotArmModel
    @test only(arm.links) isa RobotArmLink
    @test only(arm.joints) isa RobotArmJoint
    @test robot_arm_total_reach(arm) == 0.2
    pose = robot_arm_fk(arm, base, [pi / 2])
    @test pose isa RobotArmPose
    @test pose.end_effector_position ≈ SVector(0.0, 0.2, 0.0) atol=1e-15
    state = robot_arm_fk_state(arm, base, [pi / 2]; dq=[0.0])
    @test state isa RobotArmState
    @test robot_arm_end_effector_pose(state) == robot_arm_end_effector_pose(pose)
    @test robot_arm_ik(arm, base, pose.end_effector_position; q_seed=[pi / 2]) == [pi / 2]
    @test_throws ArgumentError default_robot_arm_model(link_lengths_m=(0.2,), link_radii_m=())

    tip = robot_arm_fk(arm, base, [0.0]).end_effector_position
    plan = plan_robot_arm_motion(arm, base, [0.0], tip;
        config=RobotArmPlannerConfig(dt_s=0.01, duration_s=0.01))
    @test plan.planner === :cloth_quintic
    reference = compliant_robot_arm_reference_state(plan, 0.0)
    @test reference isa CompliantRobotArmReferenceState
    @test reference.end_effector_position == tip
    x = compliant_robot_arm_initial_state(plan)
    @test length(x) == 13
    @test compliant_robot_arm_end_effector(plan, x) == tip
    @test length(compliant_robot_arm_rest_quaternions(plan, 0.0)) == 1
    multibody = compliant_robot_arm_multibody(plan)
    @test only(multibody.joints).name == :cloth_arm_joint_1
    @test only(compliant_robot_arm_actuators(plan)).name == :cloth_arm_actuator_1

    sc = ComponentVector(merge((pos=zeros(3), vel=zeros(3),
        q=[0.0, 0.0, 0.0, 1.0], ω=zeros(3)), coupled_compliant_robot_arm_state_shape(plan)))
    initialize_coupled_compliant_robot_arm_state!(sc, plan)
    @test sc.arm_r[:, 1] == x[1:3]
    sc.arm_r[1, 1] += 0.01
    du = zero(sc)
    forces, torques = zeros(3), zeros(3)
    assign_coupled_compliant_robot_arm_rhs!(du, sc, plan, 0.0, forces, torques;
        k_translation_n_m=100.0, c_translation_n_s_m=0.0,
        k_rotation_n_m_rad=0.0, c_rotation_n_m_s_rad=0.0)
    @test forces[1] > 0.0
    @test du.arm_v[1, 1] < 0.0

    simulation = simulate_compliant_robot_arm_plan(plan; dt_s=0.01, duration_s=0.01,
        integrator=:rk4, k_translation_n_m=0.0, c_translation_n_s_m=0.0,
        k_rotation_n_m_rad=0.0, c_rotation_n_m_s_rad=0.0)
    @test simulation isa CompliantRobotArmSimulation
    @test simulation.trajectory.t_s == [0.0, 0.01]
    @test simulation.trajectory.states[1] == simulation.trajectory.states[2]
end

@testset "HYPR contract 1 retains its historical service names" begin
    @test HS.CONTRACT_VERSION == v"1.0.0"
    @test HS.RobotArm.ClothArmBasePose === RobotArmBasePose
    @test HS.RobotArm.ClothArmModel === RobotArmModel
    @test HS.RobotArm.cloth_fk === robot_arm_fk
    @test HS.RobotArm.cloth_ik === robot_arm_ik
    @test HS.Cloth.ClothRobotArmDynamics.simulate_cloth_robot_arm_plan === simulate_compliant_robot_arm_plan
    @test HS.Cloth.ClothRobotArmDynamics.assign_coupled_cloth_robot_arm_rhs! === assign_coupled_compliant_robot_arm_rhs!
    @test HS.Cloth.ClothMultibody.compliant_state_parts === SpaceAGORA.CompliantMultibody.compliant_state_parts
    @test HS.CLOTH_SERVICES == ((:ClothRobotArmDynamics, :simulate_cloth_robot_arm_plan),
        (:ClothRobotArmDynamics, :assign_coupled_cloth_robot_arm_rhs!), (:ClothMultibody, :compliant_state_parts))
    config = RobotArmHYPRConfig(retime_cloth_physics_enable=false, retime_cloth_dt_s=0.05)
    @test config.retime_cloth_physics_enable === false
    @test config.retime_cloth_dt_s == 0.05
end

end # module CompliantMultibodyNamingTests
