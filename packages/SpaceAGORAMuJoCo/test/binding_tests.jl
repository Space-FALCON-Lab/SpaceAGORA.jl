# Offsets in binding.jl come from the 3.11.0 headers. Each one is checked here against a value that MuJoCo
# reports through an mj_* call or that the MJCF below fixes by construction.
@testset "binding: header offsets against mj_* readbacks" begin
    B = SpaceAGORAMuJoCo.Binding
    xml = """
    <mujoco>
      <option timestep="0.005" gravity="0 0 -9.81" integrator="implicitfast"/>
      <worldbody>
        <body name="a" pos="1 2 3" euler="0 0 90">
          <freejoint name="ja"/>
          <inertial pos="0.1 0.2 0.3" mass="2.5" diaginertia="0.1 0.2 0.3"/>
          <body name="arm" pos="1 0 0">
            <joint name="h" type="hinge" axis="0 0 1" pos="0 0 0"/>
            <inertial pos="0 0 0" mass="0.5" diaginertia="0.01 0.02 0.03"/>
          </body>
        </body>
        <body name="b" pos="-4 5 6">
          <freejoint name="jb"/>
          <inertial pos="0 0 0" mass="7" diaginertia="1 1 1"/>
        </body>
      </worldbody>
      <actuator><motor name="m" joint="h"/></actuator>
      <equality>
        <weld name="w1" body1="a" body2="b" active="true"/>
        <weld name="w2" body1="a" body2="b" active="false"/>
      </equality>
    </mujoco>"""
    m = B.load_xml_string(xml)
    d = B.MjData(m)
    SIG(bit) = 1 << bit
    getstate(sig) = B.get_state!(Vector{Float64}(undef, B.state_size(m, sig)), m, d, sig)

    @test B.version() == 3_011_000
    # sizes: counts agree with mj_stateSize of the matching state element
    @test B.nq(m) == B.state_size(m, SIG(1)) == 14 + 1
    @test B.nv(m) == B.state_size(m, SIG(2)) == 6 + 1 + 6
    @test B.nu(m) == B.state_size(m, SIG(6)) == 1
    @test B.neq(m) == B.state_size(m, SIG(9)) == 2
    @test B.na(m) == B.state_size(m, SIG(3)) == 0
    @test B.njnt(m) == 3
    @test B.nbody(m) == 4
    @test B.name2id(m, B.OBJ_BODY, "b") == 3
    @test B.id2name(m, B.OBJ_BODY, 2) == "arm"
    # model arrays
    @test B.body_mass(m) ≈ [0.0, 2.5, 0.5, 7.0]
    @test sum(B.body_mass(m)) ≈ B.total_mass(m)
    @test collect(B.body_rootid(m)) == [0, 1, 1, 3]
    @test collect(B.body_dofnum(m)) == [0, 6, 1, 6]
    @test collect(B.body_jntnum(m)) == [0, 1, 1, 1]
    @test collect(B.jnt_type(m)) == [0, 3, 0]            # free, hinge, free
    @test collect(B.jnt_qposadr(m)) == [0, 7, 8]
    @test collect(B.jnt_dofadr(m)) == [0, 6, 7]
    @test B.body_inertia(m)[4:6] ≈ [0.1, 0.2, 0.3]       # body 1 (principal axes of its inertial frame)
    # options
    @test B.timestep(m) == 0.005
    @test B.gravity(m) == (0.0, 0.0, -9.81)
    @test B.integrator(m) == Int(B.INT_IMPLICITFAST)
    B.set_timestep!(m, 0.02); B.set_integrator!(m, B.INT_EULER); B.set_gravity!(m, (1.0, 2.0, 3.0))
    @test (B.timestep(m), B.integrator(m), B.gravity(m)) == (0.02, 0, (1.0, 2.0, 3.0))
    B.set_gravity!(m, (0.0, 0.0, 0.0))
    # data: reads and writes go to the memory MuJoCo itself uses
    q = B.qpos(m, d); v = B.qvel(m, d); c = B.ctrl(m, d); x = B.xfrc_applied(m, d)
    @test q[1:3] == [1.0, 2.0, 3.0] && q[4] == cosd(45)   # qpos0 of body a
    q[8] = 0.25; v[7] = -0.5; c[1] = 0.125; x[:, 4] .= 1:6
    @test getstate(SIG(1))[8] == 0.25 && getstate(SIG(2))[7] == -0.5
    @test getstate(SIG(6)) == [0.125]
    @test getstate(SIG(8)) == vec(copy(x))
    @test B.eq_active(m, d) == UInt8[1, 0] && getstate(SIG(9)) == [1.0, 0.0]
    B.forward!(m, d)
    xi = B.xipos(m, d); xq = B.xquat(m, d)
    # body a is yawed +90 deg; its COM sits at pos + R * (0.1, 0.2, 0.3) = (1 - 0.2, 2 + 0.1, 3.3)
    @test xi[:, 2] ≈ [0.8, 2.1, 3.3] atol = 1e-12
    @test xq[:, 2] ≈ [cosd(45), 0, 0, sind(45)] atol = 1e-12
    @test xi[:, 4] ≈ [-4.0, 5.0, 6.0]
    # time advances by the timestep through step1/step2 (offset of data.time, and the opt.timestep write)
    B.set_timestep!(m, 0.01)
    B.step1!(m, d); B.step2!(m, d)
    @test B.data_time(d) ≈ 0.01
    # a second model from mj_copyModel is independent
    m2 = B.copy_model(m)
    B.set_timestep!(m2, 0.5)
    @test B.timestep(m) == 0.01 && B.timestep(m2) == 0.5
    # loading from a file agrees with loading from a string
    path = tempname() * ".xml"; write(path, xml)
    m3 = B.load_xml_file(path); rm(path)
    @test B.body_mass(m3) == B.body_mass(m)
    @test_throws ErrorException B.load_xml_string("<mujoco><oops/></mujoco>")
end

# mjData.warning[i].number at offsetof(mjData, warning) + 8 i + 4, checked against what MuJoCo itself records.
@testset "binding: warning counters" begin
    B = SpaceAGORAMuJoCo.Binding
    m = B.load_xml_string("<mujoco><worldbody><body><freejoint/><inertial pos='0 0 0' mass='1' diaginertia='1 1 1'/></body></worldbody></mujoco>")
    d = B.MjData(m)
    @test all(B.warning_count(d, w) == 0 for w in 0:6)
    B.xfrc_applied(m, d)[1, 2] = Inf
    B.step1!(m, d); B.step2!(m, d)          # MuJoCo logs, resets the data and carries on
    @test B.warning_count(d, B.WARN_BADQACC) == 1
    @test B.warning_count(d, B.WARN_BADQPOS) == 0 && B.warning_count(d, B.WARN_BADQVEL) == 0
    @test sum(B.warning_count(d, w) for w in 0:6) == 1    # no other slot aliases the BADQACC counter
end
