# Shared setup for the mounted-actuator demo and its solver comparison.
# Include after examples/common.jl; each call makes fresh vehicle/controller state.
using StaticArrays, LinearAlgebra

const WHEEL_TORQUE_NM = 2.0e-3
const WHEEL_MOMENTUM_NMS = 0.02
const ROD_DIPOLE_AM2 = 0.2
const LASER_POWER_W = 1.0           # optical power transmitted
const LASER_EFFICIENCY = 0.2        # so it draws 5 W

# The controller. Its levels are, in mounting order: each wheel's torque on the
# body about its spin axis (N m), each rod's dipole (A m^2), the terminal's
# optical power (W). The wheels damp the body rate and cancel the rods' torque,
# so the rods' cross-product law (dipole = gain * h x B / |B|^2) takes momentum
# out of the wheels without turning the body. The field comes from an ideal
# magnetometer: the field the run's rods act in.
mutable struct DetumbleAndUnload <: AbstractControlEffectorModel
    rate_gain::Float64          # N m s: wheel torque -rate_gain * ω
    unload_gain::Float64        # 1/s
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::DetumbleAndUnload) = :mounted
SpaceAGORA.actuator_levels(c::DetumbleAndUnload, i::Int) = i == 1 ? c.levels : nothing

function SpaceAGORA.update_control!(c::DetumbleAndUnload, u, p, t::Float64, i::Int64)
    i == 1 || return nothing
    state = u.sc[1]
    ω = SVector{3, Float64}(state.ω)
    h = SVector{3, Float64}(state.h_wheels)     # the wheels' integrated momenta; the wheels lie along x, y, z
    B = SM.rod_field_body(MagneticTorqueRodModel(), SVector{3, Float64}(state.pos), SVector{4, Float64}(state.q), p)
    dipole = clamp.(c.unload_gain .* cross(h, B) ./ dot(B, B), -ROD_DIPOLE_AM2, ROD_DIPOLE_AM2)
    c.levels[1:3] .= -c.rate_gain .* ω .- cross(dipole, B)
    c.levels[4:6] .= dipole
    c.levels[7] = LASER_POWER_W
    return nothing
end


function mounted_actuator_configuration(; duration_s=120.0, sample_s=nothing,
    solver=nothing, tolerances=nothing)
    planet = Earth()

    # A 10 kg 6U-class bus in a 500 km sun-synchronous orbit, tumbling at about
    # 3 degrees per second about each axis.
    ic = InitialCondition(planet.Rp_e + 500e3, 1.0e-4, 97.4, 0.0, 30.0, 0.0,
        SVector(0.0, 0.0, 0.0, 1.0), SVector(0.05, -0.05, 0.05))
    bus = Link(root=true, m=10.0, dims=MVector(0.1, 0.2, 0.3))
    spacecraft = SpacecraftModel(links=[bus], root=bus, prop_mass=0.0, inertia_tensor=bus.inertia,
        initial_condition=ic, id=1)

    # Mounting order: wheels first, then rods, then the terminal (thrusters,
    # wheels, rods, laser terminals, hinges, thruster gimbals, CMGs, then user-defined kinds).
    for axis in (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        add_actuator!(spacecraft, bus, ReactionWheel(spin_axis=axis,
            max_torque_nm=WHEEL_TORQUE_NM, max_momentum_nms=WHEEL_MOMENTUM_NMS))
    end
    for axis in (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        add_actuator!(spacecraft, bus, Magnetorquer(axis=axis, max_dipole_am2=ROD_DIPOLE_AM2))
    end
    add_actuator!(spacecraft, bus, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0),
        max_power_w=LASER_POWER_W, efficiency=LASER_EFFICIENCY))

    controller = DetumbleAndUnload(0.01, 0.005, zeros(7))
    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=duration_s,
        initial_time=InitialTime(year=2024, month=3, day=20, hour=0, minute=0, second=0.0),
        dynamic_effectors=(InverseSquaredGravityModel(),),
        ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=true,
        keplerian=false,
        EI_km=100.0,
        verbose=false,
        results=false,
    )
    args = SM.SimConfig._with_configuration(base;
        control_model=ControlModel(control_effectors=(controller,), control_rates=[1.0]),
        solver_config=solver,
        integration_tolerances=isnothing(tolerances) ? base.integration_tolerances : tolerances,
        mission_configuration=isnothing(sample_s) ? base.mission_configuration :
            SM.SimConfig._with_fields(base.mission_configuration; data_rate=sample_s))
    return args
end
