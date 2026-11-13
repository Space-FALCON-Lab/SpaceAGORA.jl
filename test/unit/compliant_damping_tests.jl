module CompliantDampingTests

using Test, LinearAlgebra, StaticArrays, Random, ComponentArrays, SpaceAGORA
const CM = SpaceAGORA.SimulationModel.ClothMultibody
const CA = SpaceAGORA.SimulationModel.ClothRobotArmDynamics
const RB = SpaceAGORA.SimulationModel.Robotics
const RP = SpaceAGORA.SimulationModel.RobotArmPlanning
const Z = SVector(0.0, 0.0, 0.0)
const Q = SVector(0.0, 0.0, 0.0, 1.0)

function fixture(; c=0.2, spin=Z, boost=Z)
    # Separated anchors and transverse motion reproduce the old internal moment.
    r = [SVector(0.3, -0.2, 0.1), SVector(0.49, -0.2, 0.1)]
    nodes = [CM.CompliantTopologyNode(Symbol(:b, i); mass_kg=i,
        inertia_body_kg_m2=0.1i, position=r[i],
        velocity=(i == 1 ? SVector(0.0, 0.1, 0.0) : Z) + boost + cross(spin, r[i]),
        angular_velocity=spin) for i in 1:2]
    edge = CM.CompliantTopologyEdge(:joint, 1, 2;
        parent_point_body=SVector(0.1, 0.0, 0.0),
        child_point_body=SVector(-0.1, 0.0, 0.0), c_translation_n_s_m=c)
    return CM.build_compliant_topology(nodes, [edge])
end

function ledger(model, x)
    j = only(model.joints)
    load = only(CM.compliant_joint_loads(model, x))
    a = CM._parent_kinematics(model, x, j.parent, j.parent_point_body)
    b = CM._parent_kinematics(model, x, j.child, j.child_point_body)
    midpoint = (a.point + b.point) / 2
    va = a.point_velocity + cross(a.ω_world, midpoint - a.point)
    vb = b.point_velocity + cross(b.ω_world, midpoint - b.point)
    fp, fc = load.translation_force_parent_world, load.translation_force_child_world
    tp, tc = load.compliance_torque_parent_world, load.compliance_torque_child_world
    moment = cross(a.point, fp) + cross(b.point, fc) + tp + tc
    power = dot(fp, a.point_velocity) + dot(fc, b.point_velocity) +
        dot(tp, a.ω_world) + dot(tc, b.ω_world)
    return (; fp, fc, moment, power, slip=va-vb, load)
end

@testset "compliant damping balances moment and dissipates midpoint slip" begin
    base = fixture()
    l = ledger(base.model, base.initial_state)
    @test l.fp ≈ SVector(0.0, -0.02, 0.0)
    @test norm(l.moment) < 1e-15
    @test l.power ≈ -0.002 atol=1e-15
    # A common rigid rotation/translation must not change the damper's slip/load.
    boosted = fixture(spin=SVector(0.2, -0.3, 2.0), boost=SVector(3.0, 2.0, -1.0))
    lb = ledger(boosted.model, boosted.initial_state)
    @test lb.fp ≈ l.fp atol=1e-14
    @test norm(lb.moment) < 1e-14
    @test lb.power ≈ l.power atol=1e-14
    rigid=copy(boosted.initial_state);rigid[9]-=.1
    @test norm(ledger(boosted.model,rigid).fp) < 1e-14

    rng = MersenneTwister(4026)
    for trial in 1:32
        C=iseven(trial) ? .2Matrix{Float64}(I,3,3) : [.2 .03 0; .03 .3 0; 0 0 .4]
        base=fixture(c=C)
        x = copy(base.initial_state)
        for i in 1:2
            o = 13(i-1)
            x[o+1:o+3] .= randn(rng, 3)
            x[o+4:o+7] .= normalize(randn(rng, 4))
            x[o+8:o+13] .= randn(rng, 6)
        end
        l = ledger(base.model, x)
        @test norm(l.fp+l.fc) < 1e-14
        @test norm(l.moment) < 1e-13
        @test l.power ≈ -dot(l.slip, C*l.slip) atol=1e-13
        @test l.power <= 1e-14
        # The compiled attachment route must apply the same COM wrenches.
        F, T = zeros(3,2), zeros(3,2)
        pos, quat, vel, omega = zeros(3,2), zeros(4,2), zeros(3,2), zeros(3,2)
        for i in 1:2
            s = CM.compliant_state_parts(x, i)
            pos[:,i] .= s.r; quat[:,i] .= s.q; vel[:,i] .= s.v; omega[:,i] .= s.ω
        end
        mount = CM.compliant_mount_kinematics(Z, Z, Q, Z, Z, Q, Z, Q)
        CM.compliant_joint_loads_in_place!(F,T,CM.compile_compliant_model(base.model),0,
            pos,quat,vel,omega,mount,[Q])
        dx = CM.compliant_multibody_dynamics(base.model,x)
        for i in 1:2
            b = base.model.bodies[i]; o = 13(i-1)
            @test F[:,i]/b.mass_kg ≈ dx[o+8:o+10] atol=1e-13
            @test b.inertia_body_kg_m2 \ (T[:,i]-cross(omega[:,i],b.inertia_body_kg_m2*omega[:,i])) ≈ dx[o+11:o+13] atol=1e-13
        end
    end
end

function invariants(model, x)
    P, H, E = Z, Z, 0.0
    for (i,b) in enumerate(model.bodies)
        s = CM.compliant_state_parts(x,i)
        P += b.mass_kg*s.v
        H += cross(s.r,b.mass_kg*s.v) + CM._rot(s.q)*(b.inertia_body_kg_m2*s.ω)
        E += (b.mass_kg*dot(s.v,s.v) + dot(s.ω,b.inertia_body_kg_m2*s.ω))/2
    end
    return P,H,E
end

@testset "compiled damping transfers the complete wrench to a moving mount" begin
    r=SVector(.49,-.2,.1); v=SVector(.1,-.2,.05); omega=SVector(.3,.1,-.2)
    node=CM.CompliantTopologyNode(:child;mass_kg=1.0,inertia_body_kg_m2=.1,
        position=r,velocity=v,angular_velocity=omega)
    edge=CM.CompliantTopologyEdge(:mount,0,1;parent_point_body=SVector(.1,0.,0.),
        child_point_body=SVector(-.1,0.,0.),c_translation_n_s_m=.2)
    b=CM.build_compliant_topology([node],[edge])
    rm=SVector(.3,-.2,.1); vm=SVector(0.,.1,0.); wm=SVector(.2,-.1,.3)
    mount=CM.compliant_mount_kinematics(rm,vm,Q,wm,Z,Q,Z,Q)
    F,T=zeros(3,1),zeros(3,1)
    fm,tm=CM.compliant_joint_loads_in_place!(F,T,CM.compile_compliant_model(b.model),0,
        reshape(collect(r),3,1),reshape(collect(Q),4,1),reshape(collect(v),3,1),
        reshape(collect(omega),3,1),mount,[Q])
    @test norm(fm+F[:,1]) < 1e-14
    @test norm(cross(rm,fm)+tm+cross(r,F[:,1])+T[:,1]) < 1e-14
    pa=rm+edge.parent_point_body; pb=r+edge.child_point_body; mid=(pa+pb)/2
    slip=vm+cross(wm,mid-rm)-v-cross(omega,mid-r)
    power=dot(fm,vm)+dot(tm,wm)+dot(F[:,1],v)+dot(T[:,1],omega)
    @test power ≈ -.2dot(slip,slip) atol=1e-14
end

@testset "free damped bodies conserve momenta and lose mechanical energy" begin
    b = fixture(spin=SVector(0.1,-0.2,0.3))
    x = copy(b.initial_state)
    p0,h0,e0 = invariants(b.model,x)
    maxp,maxh = 0.0,0.0
    previous = e0
    slip=ledger(b.model,x).slip
    previous_power=-.2dot(slip,slip); work=0.0
    for k in 1:1000
        x = CM.step_compliant_multibody_rk4(b.model,x,(k-1)*0.001,0.001)
        p,h,e = invariants(b.model,x)
        maxp=max(maxp,norm(p-p0)); maxh=max(maxh,norm(h-h0))
        @test e <= previous + 1e-13
        previous=e
        slip=ledger(b.model,x).slip
        power=-.2dot(slip,slip)
        work += .0005*(previous_power+power)
        previous_power=power
    end
    @test maxp < 1e-12
    @test maxh < 1e-11
    @test previous < e0-1e-4
    @test abs(previous-e0-work) < 1e-9
    @info "Damping regression (1 s, RK4 dt=0.001 s)" maximum_linear_momentum_error=maxp maximum_angular_momentum_error=maxh energy_work_residual=abs(previous-e0-work)
end

@testset "robot-arm damping includes the moving mount reaction" begin
    model = RB.default_cloth_arm_model(link_lengths_m=(0.2,),link_radii_m=(0.01,),
        link_masses_kg=(1.0,),joint_axes=((0.0,0.0,1.0),))
    tip=SVector(0.2,0.0,0.0)
    plan=RP.RobotArmPlan(model,RB.ClothArmBasePose(zeros(3)),[0.0,1.0],
        zeros(1,2),zeros(1,2),zeros(1,2),repeat(reshape(collect(tip),3,1),1,2),
        [0.0],[0.0],tip,0.0,:cloth_quintic)
    sc=ComponentVector(merge((pos=zeros(3),vel=[0.0,0.1,0.0],q=collect(Q),ω=[0.2,-0.1,0.3]),
        CA.coupled_cloth_robot_arm_state_shape(plan)))
    CA.initialize_coupled_cloth_robot_arm_state!(sc,plan)
    sc.arm_r[1,1] += 0.01
    sc.arm_v[:,1] .= [0.0,0.0,0.0]
    sc.arm_ω[:,1] .= [-0.1,0.3,0.2]
    du=zero(sc); fp=zeros(3); tp=zeros(3)
    CA.assign_coupled_cloth_robot_arm_rhs!(du,sc,plan,0.0,fp,tp;
        k_translation_n_m=0.0,c_translation_n_s_m=0.2,
        k_rotation_n_m_rad=0.0,c_rotation_n_m_s_rad=0.0)
    link=only(model.links); inertia=CA._link_inertia(link)
    fc=link.mass_kg*du.arm_v[:,1]
    tc=inertia*du.arm_ω[:,1]+cross(sc.arm_ω[:,1],inertia*sc.arm_ω[:,1])
    @test norm(fp+fc) < 1e-14
    @test norm(cross(sc.pos,fp)+cross(sc.arm_r[:,1],fc)+tp+tc) < 1e-13
    a=CA._coupled_parent_kinematics(sc,0,model.mount_offset_body)
    b=CA._coupled_parent_kinematics(sc,1,-link.com_offset_parent)
    mid=(a.point+b.point)/2
    slip=a.point_velocity+cross(a.ω_world,mid-a.point)-b.point_velocity-cross(b.ω_world,mid-b.point)
    power=dot(fp,sc.vel)+dot(tp,sc.ω)+dot(fc,sc.arm_v[:,1])+dot(tc,sc.arm_ω[:,1])
    @test power ≈ -0.2dot(slip,slip) atol=1e-13
end

end
