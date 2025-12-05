include("../simulation/Run.jl")
include("../config.jl") #TODO:Figure out how to run multiple times without having to comment this line out
include("../utils/maneuver_plans.jl")
include("../utils/attitude_control_plans.jl")
include("../utils/quaternion_utils.jl")
# include("SpacecraftModel.jl")

import .config
import .ref_sys
using Profile
using Interpolations
using Arrow
using Plots
plotly()

# --- Type Alias for Scalar-Last Quaternion ---
# We define a type alias for clarity: a 4-element vector representing [qx, qy, qz, qw]
const QuatSL = Vector{Float64}

# --- Helper Functions for Quaternion Math (Scalar-Last) ---

"""
    quat_mult(q1::QuatSL, q2::QuatSL)

Performs quaternion multiplication q_result = q1 * q2 for scalar-last quaternions [x, y, z, w].
"""
function quat_mult(q1::QuatSL, q2::QuatSL)::QuatSL
    v1, w1 = q1[1:3], q1[4]
    v2, w2 = q2[1:3], q2[4]

    # Scalar part: w_res = w1*w2 - v1 . v2
    w_res = w1 * w2 - dot(v1, v2)

    # Vector part: v_res = w1*v2 + w2*v1 + cross(v1, v2)
    v_res = w1 * v2 + w2 * v1 + cross(v1, v2)

    return [v_res[1], v_res[2], v_res[3], w_res]
end

"""
    quat_inv(q::QuatSL)

Calculates the inverse of a unit quaternion q.
For a unit quaternion, the inverse is the conjugate: q_inv = [-qx, -qy, -qz, qw].
"""
function quat_inv(q::QuatSL)::QuatSL
    return [-q[1], -q[2], -q[3], q[4]]
end

# --- DCM to Quaternion Conversion (Scalar-Last, Robust Trace Method) ---

"""
    dcm_to_quat_scalar_last(C::AbstractMatrix{Float64})::QuatSL

Converts a 3x3 Direction Cosine Matrix (DCM) to a scalar-last quaternion [x, y, z, w].
Uses a robust method based on the maximum diagonal element.
"""
function dcm_to_quat_scalar_last(C::AbstractMatrix{Float64})::QuatSL
    T = C[1,1] + C[2,2] + C[3,3]

    local x, y, z, w, S

    if T > 0.0
        # Case 1: Trace is largest (most common)
        S = sqrt(T + 1.0) * 2.0  # S = 4*w
        w = 0.25 * S
        x = (C[3,2] - C[2,3]) / S
        y = (C[1,3] - C[3,1]) / S
        z = (C[2,1] - C[1,2]) / S
    elseif (C[1,1] > C[2,2]) && (C[1,1] > C[3,3])
        # Case 2: C[1,1] is largest
        S = sqrt(1.0 + C[1,1] - C[2,2] - C[3,3]) * 2.0  # S = 4*x
        x = 0.25 * S
        w = (C[3,2] - C[2,3]) / S
        y = (C[1,2] + C[2,1]) / S
        z = (C[1,3] + C[3,1]) / S
    elseif C[2,2] > C[3,3]
        # Case 3: C[2,2] is largest
        S = sqrt(1.0 + C[2,2] - C[1,1] - C[3,3]) * 2.0  # S = 4*y
        y = 0.25 * S
        w = (C[1,3] - C[3,1]) / S
        x = (C[1,2] + C[2,1]) / S
        z = (C[3,2] + C[2,3]) / S
    else
        # Case 4: C[3,3] is largest
        S = sqrt(1.0 + C[3,3] - C[1,1] - C[2,2]) * 2.0  # S = 4*z
        z = 0.25 * S
        w = (C[2,1] - C[1,2]) / S
        x = (C[1,3] + C[3,1]) / S
        y = (C[2,3] + C[3,2]) / S
    end

    # Normalize to ensure unit quaternion (optional, but good practice)
    q = [x, y, z, w]
    return q / norm(q)
end


"""
    get_lvlh_to_eci_dcm(r_eci::Vector{Float64}, v_eci::Vector{Float64})

Calculates the Direction Cosine Matrix C_LE that transforms vectors from
LVLH (Local Vertical/Local Horizontal) to ECI. The columns of C_LE are the
LVLH basis vectors expressed in ECI coordinates.

NEW LVLH Frame Axes Convention (Right-Handed System):
- Z_LVLH: Nadir pointing (-R, opposite of position vector r)
- Y_LVLH: Opposite Orbit Normal (-N, opposite of r x v)
- X_LVLH: Transverse (T, completes the system: Y x Z)
"""
function get_lvlh_to_eci_dcm(r_eci::Vector{Float64}, v_eci::Vector{Float64})::Matrix{Float64}
    # 1. Calculate the standard R, N, T basis vectors (normalized)
    
    # R-bar (Radial): Unit vector along position r
    R_bar = normalize(r_eci)
    
    # H (Orbit Normal): Unnormalized normal vector h = r x v
    H_vec = cross(r_eci, v_eci)

    # N-bar (Normal): Unit vector along orbit normal h
    N_bar = normalize(H_vec)

    # T-bar (Transverse): Completes the standard R-T-N system: T = N x R
    T_bar = cross(N_bar, R_bar)
    
    # 2. Define the NEW LVLH Axes based on the user's request
    
    # X_LVLH: Transverse (T_bar)
    X_LVLH = T_bar

    # Y_LVLH: Opposite Orbit Normal (-N_bar)
    Y_LVLH = -N_bar

    # Z_LVLH: Nadir pointing (-R_bar)
    Z_LVLH = -R_bar
    
    # Verification of right-hand rule for the new frame: X = Y x Z
    # T_bar = (-N_bar) x (-R_bar) = N_bar x R_bar = T_bar. (Correct)

    # 3. Form the DCM C_LE = [X_LVLH | Y_LVLH | Z_LVLH] (columns are LVLH axes in ECI)
    C_LE = [
        X_LVLH[1] Y_LVLH[1] Z_LVLH[1];
        X_LVLH[2] Y_LVLH[2] Z_LVLH[2];
        X_LVLH[3] Y_LVLH[3] Z_LVLH[3]
    ]

    return C_LE
end

# ==============================================================================
#                      MAIN CONVERSION FUNCTION
# ==============================================================================

"""
    eci_to_lvlh_attitude(q_eci_body::QuatSL, r_eci::Vector{Float64}, v_eci::Vector{Float64})

Converts a spacecraft attitude quaternion from ECI-to-Body (q_EB) to LVLH-to-Body (q_LB).

The conversion follows the quaternion chain rule:
q_LB = q_EB * q_LE
where:
- q_EB: Spacecraft attitude (ECI -> Body), provided as input.
- q_LE: LVLH -> ECI attitude, which is the inverse of the ECI -> LVLH quaternion (q_EL).

# Arguments
- `q_eci_body`: The input ECI-to-Body attitude quaternion [x, y, z, w].
- `r_eci`: Spacecraft position vector in ECI frame [x, y, z] (m).
- `v_eci`: Spacecraft velocity vector in ECI frame [vx, vy, vz] (m/s).

# Returns
- A scalar-last quaternion [x, y, z, w] representing the LVLH-to-Body attitude (q_LB).
"""
function eci_to_lvlh_attitude(q_eci_body::QuatSL, r_eci::Vector{Float64}, v_eci::Vector{Float64})::QuatSL
    # 1. Calculate the DCM C_LE (LVLH -> ECI) from the state vectors using the new T, -N, -R convention
    C_LE = get_lvlh_to_eci_dcm(r_eci, v_eci)

    # 2. Convert the DCM C_LE to the LVLH -> ECI quaternion (q_LE)
    # The DCM for ECI -> LVLH (C_EL) is the transpose of C_LE: C_EL = C_LE'
    C_EL = C_LE'

    # Convert C_EL to q_EL (ECI -> LVLH)
    q_EL = dcm_to_quat_scalar_last(C_EL)

    # 3. Calculate q_LE (LVLH -> ECI) as the inverse (conjugate) of q_EL
    q_LE = quat_inv(q_EL)

    # 4. Apply the quaternion chain rule: q_LB = q_EB * q_LE
    q_lvlh_body = quat_mult(q_eci_body, q_LE)

    return q_lvlh_body
end
"""
    get_lvlh_angular_velocity(r_eci::Vector{Float64}, v_eci::Vector{Float64})::Vector{Float64}

Calculates the angular velocity vector of the LVLH frame (relative to ECI)
expressed in ECI coordinates (omega_L/I).
"""
function get_lvlh_angular_velocity(r_eci::Vector{Float64}, v_eci::Vector{Float64})::Vector{Float64}
    # Standard R-bar, N-bar, T-bar basis vectors
    r_norm = norm(r_eci)
    R_bar = r_eci / r_norm
    H_vec = cross(r_eci, v_eci)
    N_bar = normalize(H_vec)
    T_bar = cross(N_bar, R_bar) # T = N x R

    # Components of omega_L/I (angular velocity of LVLH frame in R-T-N axes)
    # Omega_R = 0 (Rotation about Radial axis)
    omega_T = dot(r_eci, v_eci) / (r_norm^2)
    omega_N = norm(H_vec) / (r_norm^2)

    # omega_L/I expressed in ECI frame is: omega_T * T_bar + omega_N * N_bar
    omega_L_I_eci = omega_T * T_bar + omega_N * N_bar

    return omega_L_I_eci
end

# import .SpacecraftModel
# Define spacecraft model
spacecraft = config.SpacecraftModel()
# Add bodies to the spacecraft model
# p = SVector{3, Float64}([0.1, 0.2, -0.3])
# q = 1/(1+norm(p)^2)*SVector{4, Float64}([2*p; 1-norm(p)^2])
# skew = (ω) -> SMatrix{3, 3, Float64}([0 -ω[3] ω[2];
#                                    ω[3] 0 -ω[1];
#                                    -ω[2] ω[1] 0])
# dcm = (q[4]^2 - norm(q[1:3])^2)*I(3) - 2*q[4]*skew(q[1:3]) + 2*q[1:3]*q[1:3]' # DCM from quaternion
# ω_body = SVector{3, Float64}([0.001, -0.01, 0.03]) # Reference angular velocity
# ω_ref = dcm'*ω_body
# h = sqrt(30.0/7.0)
# w = sqrt(6.0)
# d = sqrt(66.0/7.0)

rw_torques_data = DataFrame(Arrow.Table("cygnss_rw_momentum_derivatives_no_filter.feather"))

println("RW torque data loaded from file, time range: $(minimum(rw_torques_data[!, 1])) to $(maximum(rw_torques_data[!, 1])) seconds.")
rw_1_itp = cubic_spline_interpolation(range(rw_torques_data[1, 1], rw_torques_data[end, 1], length(rw_torques_data[!, 2])), rw_torques_data[!, 2])
rw_2_itp = cubic_spline_interpolation(range(rw_torques_data[1, 1], rw_torques_data[end, 1], length(rw_torques_data[!, 3])), rw_torques_data[!, 3])
rw_3_itp = cubic_spline_interpolation(range(rw_torques_data[1, 1], rw_torques_data[end, 1], length(rw_torques_data[!, 4])), rw_torques_data[!, 4])
rw_torque_itp = (t) -> SVector{3, Float64}([rw_1_itp(t), rw_2_itp(t), rw_3_itp(t)])
println(rw_torque_itp(2))

rw_torques_cloth = DataFrame(Arrow.Table("slew_maneuver_torques_unfiltered.feather"))
cloth_times = rw_torques_cloth[1, 1]:rw_torques_cloth[2, 1] - rw_torques_cloth[1, 1]:rw_torques_cloth[end, 1]
println(length(cloth_times))
println(length(rw_torques_cloth[!, :rho_dot_1]))
rw_1_itp_cloth = cubic_spline_interpolation(range(rw_torques_cloth[1, 1], stop=rw_torques_cloth[end, 1], length=length(rw_torques_cloth[!, :rho_dot_1])), rw_torques_cloth[!, :rho_dot_1])
rw_2_itp_cloth = cubic_spline_interpolation(range(rw_torques_cloth[1, 1], stop=rw_torques_cloth[end, 1], length=length(rw_torques_cloth[!, :rho_dot_1])), rw_torques_cloth[!, :rho_dot_2])
rw_3_itp_cloth = cubic_spline_interpolation(range(rw_torques_cloth[1, 1], stop=rw_torques_cloth[end, 1], length=length(rw_torques_cloth[!, :rho_dot_1])), rw_torques_cloth[!, :rho_dot_3])
rw_torque_itp_cloth = (t) -> SVector{3, Float64}([rw_1_itp_cloth(t), rw_2_itp_cloth(t), rw_3_itp_cloth(t)])
Plots.plot(cloth_times, rw_1_itp_cloth(cloth_times))
Plots.plot!(cloth_times, rw_2_itp_cloth(cloth_times))
Plots.display(Plots.plot!(cloth_times, rw_3_itp_cloth(cloth_times)))

torque_rods_data = DataFrame(Arrow.Table("cygnss_dipole_commands.feather"))
println("Torque rod data loaded from file, time range: $(minimum(torque_rods_data[!, 1])) to $(maximum(torque_rods_data[!, 1])) seconds.")
tr_x_itp = cubic_spline_interpolation(range(torque_rods_data[1, 1], torque_rods_data[end, 1], length(torque_rods_data[!, 2])), torque_rods_data[!, 2])
tr_y_itp = cubic_spline_interpolation(range(torque_rods_data[1, 1], torque_rods_data[end, 1], length(torque_rods_data[!, 3])), torque_rods_data[!, 3])
tr_z_itp = cubic_spline_interpolation(range(torque_rods_data[1, 1], torque_rods_data[end, 1], length(torque_rods_data[!, 4])), torque_rods_data[!, 4])
tr_command_itp = (t) -> MVector{3, Float64}([tr_x_itp(t), tr_y_itp(t), tr_z_itp(t)])
println("dipole command at 2.0 seconds: ", tr_command_itp(2.0))

function CYGNSS_attitude_function(state, m, b::config.Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    # Simple attitude control function for CYGNSS
    axes = [MVector{3, Float64}(-1.0, 0.0, 0.0),
            MVector{3, Float64}(0.0, 1.0, 0.0),
            MVector{3, Float64}(0.0, 0.0, -1.0)]
    # offset = 890.0017919540405 - 0.5
    offset = 0.0
    # b.ω_wheel_derivatives .= rw_torque_itp_cloth(t+offset)
    # if !isempty(b.magnets)
    #     tr_cmd_mags = tr_command_itp(t+offset)
    #     for i in 1:3
    #         b.magnets[i].magnitude = tr_cmd_mags[i]
    #     end
    # end
end

integral_error = MVector{3, Float64}(0.0, 0.0, 0.0)
function CYGNSS_attitude_pid_function(state, m, b::config.Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64)
    r_eci = state[1:3]
    v_eci = state[4:6]
    quat_idx = 8 + length(m.body.links)
    q_eci_body = state[quat_idx:quat_idx+3]
    omega_eci_body = state[quat_idx+4:quat_idx+6]
    O_BE = rot(q_eci_body) # Orientation matrix of body relative to ECI
    # Determine current LVLH quaternion
    q_lvlh = eci_to_lvlh_attitude(q_eci_body, r_eci, v_eci)
    # Determine target LVLH quaternion
    q_target_lvlh = t < 901.75 ? SVector{4, Float64}(0.0, 0.0, 0.0, 1.0) : SVector{4, Float64}(-0.087156, 0.0, 0.0, 0.99619)
    # Calculate the error quaternion
    q_inv = SVector{4, Float64}(-q_lvlh[1:3]..., q_lvlh[4])
    qe = SVector{4, Float64}(hcat(Ψ(q_target_lvlh), q_target_lvlh) * q_inv)
    # Calculate LVLH-relative angular velocity
    omega_lvlh = get_lvlh_angular_velocity(r_eci, v_eci) # Angular velocity of the LVLH frame
    omega_e = omega_eci_body - O_BE*omega_lvlh # Angular velocity error in resolved in body frame
    # PID control gains
    Kq = 0.0 # Quaternion error gain
    Kw = 0.0 # Angular velocity error gain
    Ki = 0.0 # Integral gain
    q_error = qe[1:3] #* qe[4] # Error quaternion compensation calculation
    angular_error_threshold = deg2rad(5.0) # 0.1 degree threshold
    global integral_error
    integral_error += q_error * b.attitude_control_rate
    tau_cmd = MVector{3, Float64}(0.0, 0.0, 0.0)
    euler_gyro = cross(omega_eci_body, (b.inertia * omega_eci_body + b.J_rw * b.rw)) # eulerian gyroscopic term in torque calculation
    omega_error_term = Kw * omega_e # angular velocity error term in torque calculation
    tau_cmd = omega_error_term + euler_gyro
    angle_error = acos(clamp(qe[4], -1.0, 1.0)) * 2.0
    if abs(angle_error) < angular_error_threshold
        println("Integral control active at time $(t) seconds.")
        tau_cmd .+= Kq * q_error + Ki * integral_error
    else
        integral_error .= MVector{3, Float64}(0.0, 0.0, 0.0)
    end
    # Apply torque command to reaction wheels
    b.ω_wheel_derivatives .= pinv(b.J_rw) * tau_cmd
end

# q = SVector{4, Float64}([0.0, 0.0, sin(pi/4), cos(pi/4)]) # Quaternion for the main bus
main_bus = config.Link(root=true, 
                        r=SVector{3, Float64}(0.0, 0.0, 0.0), # Body z-axis points down, origin is at bottom, CoM from engineering drawing 
                        # q=SVector{4, Float64}(q),
                        # q=SVector{4, Float64}([0.28047528 -0.17599893  0.9414761  -0.06309311]),
                        # q=SVector{4, Float64}([ -0.000178090669, 0.000196625584, -0.000787386924,0.999990655]), # Initial quaternion, from slew data, LVLH
                        # q=SVector{4, Float64}([-1.78090669e-04, 1.96625584e-04, -7.87386924e-04, 9.99999655e-01]), # Stating from ~900s, from slew data, assumes scalar first, LVLH
                        q=SVector{4, Float64}([-0.769326211835314, -0.0287409968395995, 0.368405744863056, 0.521141383921568]), # Initial quaternion, from slew data, ECI
                        # q=SVector{4, Float64}([-0.49694828, -0.27337817, 0.69593993, 0.44042526]), # Starting from ~900s, from slew data, ECI
                        ṙ=SVector{3, Float64}([0.0, 0.0, 0.0]), 
                        # ω=SVector{3, Float64}([0.00029976 -0.00091251  0.00051997]), # Initial angular velocity rad/s, from CYGNSS documentation
                        # ω=SVector{3, Float64}([-9.789142632171234e-5, -8.82827330140926e-5, 0.00012436837964648057]), # Initial angular velocity rad/s, from slew data, LVLH
                        # ω=SVector{3, Float64}([4.94494778455986e-07, -1.896896037665177e-05, -2.7264315784859147e-06]), # Starting from ~900s, from slew data, LVLH
                        ω=SVector{3, Float64}([-9.167855927546074e-05, -0.0011950697001943617, 0.0001295462530991909]), # Initial angular velocity rad/s, from slew data, ECI
                        # ω=SVector{3, Float64}([-1.249549492745143e-06, -0.0011257840337701133, -2.340685254945601e-06]), # Starting from ~900s, from slew data, ECI
                        dims=SVector{3, Float64}([20.222e-2, 52.12e-2, 64.09e-2]), 
                        ref_area=0.1129753, # m^2
                        # ref_area=0.0,
                        m=28.94,
                        gyro=3,
                        attitude_control_rate=0.1, # seconds
                        rw=SVector{3, Float64}(25.790585298785203/6000.0*18.0e-3, -11.830065546928386/6000.0*18.0e-3, -532.0841089980587/6000.0*18.0e-3), # Initial RW angular momentum Nms, from CYGNSS documentation
                        # rw=SVector{3, Float64}(-85.75349176655784/6000.0*18.0e-3, 571.8118011223758/6000.0*18.0e-3, -262.8650773012001/6000.0*18.0e-3), # Initial RW angular momentum starting from 900s, from slew data
                        max_torque=0.65536e-3, # max torque from RW datasheet, Nm
                        max_h=18.0e-3, # max angular momentum from RW datasheet, Nms
                        J_rw=SMatrix{3, 3, Float64}([ 0.8164  0.4083  -0.4083;
                                                     -0.5774  0.5773  -0.5773;
                                                      0.0000 -0.7071  -0.7071]),
                        # J_rw=SMatrix{3, 3, Float64}([ 0.81171104  0.41303245 -0.41295202;
                        #                               0.58395316 -0.56044371  0.5872832 ;
                        #                              -0.01113066  0.7178489   0.69610995]),
                        # attitude_control_function=(m, b::config.Link, root_index::Int, vel_pp_rw::SVector{3, Float64}, h_pp_hat::SVector{3, Float64}, aerobraking_phase::Int, t::Float64) -> (b.ω_wheel_derivatives .= pinv(b.J_rw) * rw_torque_itp_cloth(t+890.0017919540405))) # cloth attitude control
                        attitude_control_function=CYGNSS_attitude_pid_function) # CYGNSS data attitude control

L_panel = config.Link(r=SVector{3, Float64}(0.0, 56.9e-2, -(20.222 - 13.1)*1.0e-2),
                      q=SVector{4, Float64}(0.0, sqrt(2.0)/2.0, 0.0, sqrt(2.0)/2.0),
                      dims=SVector{3, Float64}([49.71e-2, 0.05, 52.12e-2]),
                      ref_area=49.71e-4,
                      m=0.01)
                      
R_panel = config.Link(r=SVector{3, Float64}(0.0, -56.9e-2, -(20.222 - 13.1)*1.0e-2),
                      q=SVector{4, Float64}(0.0, sqrt(2.0)/2.0, 0.0, sqrt(2.0)/2.0),
                      dims=SVector{3, Float64}([49.71e-2, 0.05, 52.12e-2]),
                      ref_area=49.71e-4,
                      m=0.01)
config.add_body!(spacecraft, main_bus, prop_mass=0.0)
config.add_body!(spacecraft, L_panel)
config.add_body!(spacecraft, R_panel)

L_panel_joint = config.Joint(main_bus, L_panel)
R_panel_joint = config.Joint(R_panel, main_bus)
config.add_joint!(spacecraft, L_panel_joint)
config.add_joint!(spacecraft, R_panel_joint)
inertia_tensor = [1.4e6 -1.71e4 8.08e3;
                  -1.71e4 8.19e5 -5.35e3;
                  8.08e3 -5.35e3 1.95e6] * 1e-6
config.set_inertia_tensor!(spacecraft, main_bus, 
                        SMatrix{3, 3, Float64}(inertia_tensor))

# Add torque rods
# Dipole commands in Am^2, locations in m
# TRX = config.Magnet(axis=MVector{3, Float64}(-1.0, 0.0, 0.0), magnitude=0.0)
# TRY = config.Magnet(axis=MVector{3, Float64}(0.0, 1.0, 0.0), magnitude=0.0)
# TRZ = config.Magnet(axis=MVector{3, Float64}(0.0, 0.0, -1.0), magnitude=0.0)
# config.add_magnet!(main_bus, TRX)
# config.add_magnet!(main_bus, TRY)
# config.add_magnet!(main_bus, TRZ)
# Residual dipole moment
# sc_dipole = config.Magnet(axis=normalize(MVector{3, Float64}(-0.0101, -0.0417, 0.1697)), magnitude=norm(MVector{3, Float64}(-0.0101, -0.0417, 0.1697)))
# config.add_magnet!(main_bus, sc_dipole)

println("Spacecraft model initialized with $(length(spacecraft.links)) bodies.")
# println("Spacecraft roots: $spacecraft.roots")
println("Spacecraft COM: $(config.get_COM(spacecraft, main_bus))")
println("Spacecraft MOI: $(config.get_inertia_tensor(spacecraft, main_bus))")
lenXHub = 52.12
lenYHub = 64.09
lenZHub = 29.21
lenYPanel = 49.71
faceArea = 1129.753 # Area of house-shaped face, cm^2
bus_facet_area_list = [faceArea, # Calculated from weirdly shaped face
                       faceArea,
                       25.55*lenXHub,
                       25.55*lenXHub,
                       9.388*lenXHub,
                       9.388*lenXHub,
                       lenYHub*lenXHub,
                       20.275*lenXHub] * 1e-4
panel_facet_area_list = [lenYPanel*lenXHub,
                         lenYPanel*lenXHub] * 1e-4

bus_facet_attitude_list = [q_from_phi(deg2rad(0.0) * [1.0, 0.0, 0.0]),
                           q_from_phi(deg2rad(180.0) * [0.0, 0.0, 1.0]),
                           q_from_phi(deg2rad(63.435) * [1.0, 0.0, 0.0]),
                           q_from_phi(deg2rad(-63.435) * [1.0, 0.0, 0.0]),
                           q_from_phi(deg2rad(0.0) * [1.0, 0.0, 0.0]),
                           q_from_phi(deg2rad(180.0) * [0.0, 0.0, 1.0]),
                           q_from_phi(deg2rad(180.0) * [1.0, 0.0, 0.0]),
                           q_from_phi(deg2rad(0.0) * [1.0, 0.0, 0.0])]
panel_facet_attitude_list = [q_from_phi(deg2rad(180.0) * [1.0, 0.0, 0.0]),
                             q_from_phi(deg2rad(0.0) * [1.0, 0.0, 0.0])]

bus_facet_normal_vectors = [SVector{3, Float64}([1.0, 0.0, 0.0]),
                            SVector{3, Float64}([1.0, 0.0, 0.0]),
                            SVector{3, Float64}([0.0, 1.0, 0.0]),
                            SVector{3, Float64}([0.0, 1.0, 0.0]),
                            SVector{3, Float64}([0.0, 1.0, 0.0]),
                            SVector{3, Float64}([0.0, 1.0, 0.0]),
                            SVector{3, Float64}([0.0, 0.0, 1.0]),
                            SVector{3, Float64}([0.0, 0.0, 1.0])]
panel_facet_normal_vectors = [SVector{3, Float64}([0.0, 0.0, 1.0]),
                              SVector{3, Float64}([0.0, 0.0, 1.0])]

bus_facet_locs = [SVector{3, Float64}([lenXHub*1e-2 * 0.5, 0.0, 0.0]),
                  SVector{3, Float64}([-lenXHub*1e-2 * 0.5, 0.0, 0.0]),
                  SVector{3, Float64}([0.0, 21.0915e-2, -(6.2592 - 13.1)*1.0e-2]),
                  SVector{3, Float64}([0.0, -21.0915e-2, -(6.2592 - 13.1)*1.0e-2]),
                  SVector{3, Float64}([0.0, 32.045e-2, -(12.519 + 9.388/2.0 - 13.1)*1.0e-2]),
                  SVector{3, Float64}([0.0, -32.045e-2, -(12.519 + 9.388/2.0 - 13.1)*1.0e-2]),
                  SVector{3, Float64}([0.0, 0.0, -(12.519 + 9.388 - 13.1)*1.0e-2]),
                  SVector{3, Float64}([0.0, 0.0, 13.1e-2])]
panel_facet_locs = [SVector{3, Float64}([0.0, 0.0, 0.0]),
                    SVector{3, Float64}([0.0, 0.0, 0.0])]

bus_specular_coeffs = [0.336, 0.336, 0.336, 0.336, 0.336, 0.336, 0.336, 0.336]
panel_specular_coeffs = [0.16, 0.0]
bus_diffuse_coeffs = [0.139, 0.139, 0.139, 0.139, 0.139, 0.139, 0.139, 0.139]
panel_diffuse_coeffs = [0.16, 0.56]
bus_facet_names = ["front_hub", "back_hub", "right_slant", "left_slant", "right_vert", "left_vert", "top_hub", "bot_hub"]
panel_facet_names_left = ["left_top_panel", "left_bot_panel"]
panel_facet_names_right = ["right_top_panel", "right_bot_panel"]
bus_facets = config.create_facet_list(bus_facet_area_list,
                                      bus_facet_attitude_list,
                                      bus_facet_normal_vectors,
                                      bus_facet_locs,
                                      bus_diffuse_coeffs,
                                      bus_specular_coeffs,
                                      bus_facet_names)
panel_facets_R = config.create_facet_list(panel_facet_area_list,
                                        panel_facet_attitude_list,
                                        panel_facet_normal_vectors,
                                        panel_facet_locs,
                                        panel_diffuse_coeffs,
                                        panel_specular_coeffs,
                                        panel_facet_names_right)
panel_facets_L = config.create_facet_list(panel_facet_area_list,
                                        panel_facet_attitude_list,
                                        panel_facet_normal_vectors,
                                        panel_facet_locs,
                                        panel_diffuse_coeffs,
                                        panel_specular_coeffs,
                                        panel_facet_names_left)
config.add_facet!(main_bus, bus_facets)
config.add_facet!(L_panel, panel_facets_L)
config.add_facet!(R_panel, panel_facets_R)
args = Dict(# Misc Simulation
            :results => 1,                                                                                      # Generate csv file for results True=1, False=0
            :passresults => false,                                                                                  # Pass results as output True=1, False=0
            :print_res => true,                                                                                    # Print some lines True=1, False=0
            :directory_results => "output/cygnss_comparison_slew_eci",                # Directory where to save the results
            :directory_Gram => "GRAMpy",                                                    # Directory where Gram is
            :directory_Gram_data => "GRAM_Data",                                            # Directory where Gram data is
            :directory_Spice => "GRAM_Data/SPICE",                                          # Directory where SPICE files are located
            :Gram_version => 0,                                                                                 # MarsGram x file to use
            :montecarlo_analysis => 0,                                                                          # Generate csv file for Montecarlo results True=1, False=0
            :plot => 0,                                                                                         # Generate pdf plots of results True=1, False=0
            :filename => 0,                                         # Filename with specifics of simulation, True =1, False=0
            :machine => "",                                         # choices=['Laptop' , 'Cluster' , 'Aero' , 'Desktop_Home','Karnap_Laptop']
            :integrator => "Julia",                                 # choices=['Costumed', 'Julia'] Costumed customed integrator, Julia DifferentialEquations.jl library integrator, only for drag passage, others phases use RK4
            :normalize => 1,                                       # Normalize the integration True=1, False=0
            :closed_form => 0,                                     # Closed form solution True=1, False=0
            :save_csv => false,
            # Type of Mission
            :type_of_mission => "Time",                           # choices=['Drag Passage' , 'Orbits' , 'Aerobraking Campaign']
            :keplerian => true,                                        # Do not include drag passage: True=1, False=0
            :number_of_orbits => 10,                                 # Number of aerobraking passage
            :mission_time => 3598.0,                                  # Mission time in seconds, used only for Time mission type
            # :mission_time => 1000.0,                                  # Mission time in seconds, used only for Time mission type
            :orientation_sim => true,                                  # Orientation simulation True=1, False=0, if false, will only propagate position
            :num_steps_to_save => 10000,                            # Number of timesteps between saves

            # Physical Model
            :planet => 0,                                           # Earth = 0, Mars = 1, Venus = 2
            :planettime => 0.0,                                     # Initial time of the mission, sec. Important for J2 effect and rotation of the planet
            :gravity_model => "Inverse Squared and J2 effect",      # choices=['Constant' , 'Inverse Squared' , 'Inverse Squared and J2 effect', 'GRAM']
            :density_model => "None",                               # choices=['Constant' , 'Exponential' , 'Gram']
            :topography_model => "None",                             # choices=['None' , 'Spherical Harmonics']
            :topography_harmonics_file => "Topography_harmonics_data/Earth2012.csv", # File with the topography harmonics coefficients
            :topo_degree => 90,                                     # Maximum degree of the topography harmonics (Defined in the file)
            :topo_order => 90,                                      # Maximum order of the topography harmonics (Defined in the file)
            :wind => 1,                                             # Wind calculation only if density model is Gram True=1, False=0
            :aerodynamic_model => "Mach-dependent",                 # choices=['Cd and Cl Constant' , 'Mach-dependent' , 'No-Ballistic flight with axial coefficient']: "Mach-dependent" specific for spacecraft shape, "No-Ballistic flight" specific for blunted-cone shape
            :thermal_model => "Maxwellian Heat Transfer",           # choices=['Maxwellian Heat Transfer' , 'Convective and Radiative']: "Maxwellian Heat Transfer" specific for spacecraft shape, "Convective and Radiative" specific for blunted-cone shape
            
            # Perturbations
            :n_bodies => ["Sun", "Moon"],                                        # Add names of bodies you want to simulate the gravity of to a list. Keep list empty if not required to simulate extra body gravity.
            :srp => true,                                             # Solar Radiation Pressure true/false
            :eclipse => true,                                         # Whether to include eclipse conditions in SRP calculation
            :gravity_gradient => true,                                   # Gravity Gradient true/false
            :gravity_harmonics => 1,                                            # Gravity Spherical harmonics True=1, False=0
            :gravity_harmonics_file => "Gravity_harmonics_data/EarthGGM05C.csv", # File with the gravity harmonics coefficients
            :L => 50,                                              # Maximum degree of the gravity harmonics (Defined in the file)
            :M => 50,                                              # Maximum order of the gravity harmonics (Defined in the file)
            :magnetic_field => false,                                    # Magnetic field True=1, False=0

            # Rates
            :trajectory_rate => 100.0,                              # Rate at which the trajectory in drag passage integrate using RK4
            :flash1_rate => 3.0,                                    # Rate at which Control Mode-1 is called
            :save_rate => 1.0,                                      # Rate at which the data trajectory are saved
            
            # Body
            :body_shape => "Spacecraft",                            # choices=['Spacecraft' , 'Blunted Cone']
            :max_heat_rate => 0.15,                                 # Max heat rate the heat rate control will start to react to
            :max_heat_load => 30.0,                                 # Max heat load the heat load control will not be overcomed
            # :dry_mass => 411.0,                                     # Initial dry mass of body in kg
            # :prop_mass => 50.0,                                     # Initial propellant mass of body in kg
            :reflection_coefficient => 0.9,                         # Diffuse reflection sigma =0, for specular reflection sigma = 1
            :thermal_accomodation_factor => 1.0,                    # Thermal accomodation factor, Shaaf and Chambre
            :α => 90.0,                                             # Max angle of attack of solar panels

            # # Fill for Spacecraft body shape only
            # :length_sat => 2.2,                                     # Length of the satellite in m
            # :height_sat => 1.7,                                     # Height of the satellite in m
            # :width_sat => 2.6,                                      # Width of the satellite in m
            # :length_sp => 3.76,                                     # Length of the solar panels in m
            # :height_sp => 1.93,                                     # Height of the solar panels in m

            # # Fill for Blunted Cone body shape only
            # :cone_angle => 70.0,                                    # Cone angle of the blunted cone in deg
            # :base_radius => 2.65/2,                                 # Base radius of the blunted cone in m
            # :nose_radius => 0.6638,                                 # Nose radius of the blunted cone in m
            :spacecraft_model => spacecraft,                            # Spacecraft model with bodies and joints
            
            # Engine
            :thrust => 4.0,                                         # Maximum magnitude thrust in N
            
            # Control Mode
            :control_mode => 0,                                     # Use Rotative Solar Panels Control:  False=0, Only heat rate=1, Only heat load=2, Heat rate and Heat load = 3
            :security_mode => 1,                                    # Security mode that set the angle of attack to 0 deg if predicted heat load exceed heat load limit
            :second_switch_reevaluation => 1,                       # Reevaluation of the second switch time when the time is closer to it
            :control_in_loop => 1,                                  # Control in loop, control called during integration of trajectory, full state knowledge
            :flash2_through_integration => 0,                       # Integration of the equations of motion and lambda to define time switches and revaluation second time switch
            
            # Initial Conditions
            :initial_condition_type => 2,                           # Initial Condition ra,hp = 0, Initial Condition v, gamma = 1
            :ra_initial_a => 15000.0e3,                # Initial Apoapsis Radius for for-loop in m
            :ra_initial_b => 50000e3,                               # Final Apoapsis Radius for for-loop in m
            :ra_step => 5e10,                                       # Step Apoapsis Radius for for-loop in m
            :hp_initial_a => 145.0e3,                                 # Initial Periapsis Altitude for for-loop in m
            :hp_initial_b => 1590000.0e3,                              # Final Periapsis Altitude for for-loop in m
            :hp_step => 1e12,                                 # Step Periapsis Radius for for-loop in m
            :v_initial_a => 4500.0,                                 # Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_initial_a => 4500.0,                                 # Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_initial_b => 5000.0,                                 # Final Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_step => 1000.0,                                       # Step Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :a_initial_a => 6815519.860683523,                # Initial Semi-major axis for for-loop in m, starting from 0s
            # :a_initial_a => 6813679.740083234, # starting from 900s
            # :a_initial_a => 6179921.801554798,
            # :a_initial_a => 10000.0e3,
            :a_initial_b => 6920.0e3,                               # Final Semi-major axis for for-loop in m
            # :a_initial_b => 11000.0e3,
            :a_step => 5e10,                                       # Step Semi-major axis for for-loop in m
            :e_initial_a => 0.001160182015198709,                                   # Initial Eccentricity for for-loop in m, starting from 0s
            # :e_initial_a => 0.0011523634328820804, # Starting from 900s
            # :e_initial_a => 0.10332925554563428,
            # :e_initial_a => 0.01,
            :e_initial_b => 0.1,                                   # Final Eccentricity for for-loop in m
            :e_step => 0.1,                                       # Step Eccentricity for for-loop in m
            
            :orientation_type => 2,                                   # Initial Condition orientation = 0, Initial Condition orientation and velocity = 1
            :γ_initial_a => -2.5,                                    # Initial Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_initial_b => 7.0,                                    # Final Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_step => 100,                                         # Step Gamma (deg) for for-loop if initial conditions are in v and gamma
            :inclination => 35.006305312328244,                                   # Inclination Orbit, deg, starting from 0s
            # :inclination => 35.00691817898257, # Starting from 900s
            # :inclination => 35.60900229798397,
            :ω => 175.8048936813804,                                              # AOP, deg, starting from 0s
            # :ω => 227.97599785070602, # Starting from 900s
            # :ω =>235.59080763543037,
            :Ω => 143.51935099761045,                                              # RAAN, deg, starting from 0s
            # :Ω => 143.50901715025634, # Starting from 900s
            # :Ω => 179.05921043028263,
            :ν => 345.19661998242003,                                               # True Anomaly, deg, starting from 0s
            # :ν => 350.3995039488173,# Starting from 900s
            # :ν => 180.24980892300408,
            # :ν => 40.0,                                               # True Anomaly, deg
            :EI => 160.0,                                           # Entry Interface, km
            :AE => 160.0,                                           # Atmospheric Exit, km
            :year => 2025,                                          # Mission year
            :month => 10,                                           # Mission month
            :day => 4,                                             # Mission day
            :hours => 0,                                           # Mission hour
            :minutes => 56,                                         # Mission minute
            :secs => 59.0,                                          # Mission second
            
            # Final Conditions
            :final_apoapsis => 3390.0e3+503e3, # 5088116.837416616, # 4905.974818462152e3                  # Final apoapsis radius if aerobraking campaign

            # Do not change
            :heat_load_sol => 0,                                    # Heat load solution leave it to 0 and change it only for control mode = 2:  Max energy depletaion=0, Min energy depletion=1, One switch max-min=2, One switch min-max = 3
            :thrust_control => "None",                              # choices=['None' , 'Aerobraking Maneuver' , 'Drag Passage Firing']
            :phi => 180.0,                                          # Thrust Angle, deg
            :delta_v => 0,                                          # Delta-v of Aerobraking Manuver,m/s
            :apoapsis_targeting => 0,                               # Apoapsis Targeting Enabled
            :ra_fin_orbit => 25000e3,                               # Target final apoapsis for the orbit, m
            :maneuver_plan => Odyssey_firing_plan,                # Maneuver plan function
            
            # Monte Carlo Simulations
            :montecarlo => 0,                                       # Run Monte Carlo simulation True=1, False=0
            :monte_carlo_run => 0,
            :initial_montecarlo_number => 1,                        # Initial Monte Carlo sample number
            :montecarlo_size => 1000,                               # number of Monte Carlo samples
            
            # Monte Carlo Perturbations
            :CD_dispersion => 10.0,                                 # Max dispersion of CD for Uniform Distribution, %
            :CL_dispersion => 10.0,                                 # Max dispersion of CL for Uniform Distribution, %
            :rp_dispersion => 87.0*0.05/3,                                  # Max dispersion for initial vacuum periapsis radius following uniform distribution, km
            :ra_dispersion => 28559.0*0.05/3,                                  # Max dispersion for initial apoapsis radius following uniform distribution, km
            :i_dispersion => 0.25,                                  # Max dispersion for initial inclination following uniform distribution, deg
            :Ω_dispersion => 0.25,                                  # Max dispersion for initial right ascension of the ascending node following uniform distribution, deg
            :ω_dispersion => 0.25,                                  # Max dispersion for initial argument of periapsis following uniform distribution, deg
            :vi_dispersion => 0.025,                                # Max dispersion for initial true anomaly following uniform distribution, deg
            
            # MonteCarlo Perturbation Guidance - Closed Form Solution (only for online)
            :ra_dispersion_gnc => 0.25,                             # Max dispersion for initial apoapsis radius used by gnc following uniform distribution, km
            :rp_dispersion_gnc => 0.25,                             # Max dispersion for initial periapsis radius used by gnc following uniform distribution, km
            :i_dispersion_gnc => 0.025,                             # Max dispersion for initial inclination used by gnc following uniform distribution, deg
            :Ω_dispersion_gnc => 0.025,                             # Max dispersion for initial right ascension of the ascending node used by gnc following uniform distribution, deg
            :ω_dispersion_gnc => 0.0,                               # Max dispersion for initial argument of periapsis used by gnc following uniform distribution, deg
            :vi_dispersion_gnc => 0.0,                              # Max dispersion for initial true anomaly used by gnc following uniform distribution, deg
            
            # Online trajectory control (heat rate)
            :ρ_mudispersion_gnc => 0.0,                             # Mean dispersion of rho for Gaussian Distribution, %
            :ρ_sigmadispersion_gnc => 1.0,                          # Std dispersion of rho for Gaussian Distribution, %
            :T_mudispersion_gnc => 0.0,                             # Mean dispersion of T for Gaussian Distribution, %
            :T_sigmadispersion_gnc => 1.0,                          # Std dispersion of T for Gaussian Distribution, %
            :S_mudispersion_gnc => 0.0,                             # Mean dispersion of S for Gaussian Distribution, %
            :S_sigmadispersion_gnc => 1.0,                          # Std dispersion of S for Gaussian Distribution, %
            :multiplicative_factor_heatload => 1.0,                 # Multiplicative factor for heat rate prediction when calculated heat load

            :a_tol => 1e-5,                                         # Absolute tolerance for integration
            :r_tol => 1e-3,                                         # Relative tolerance for integration
            :a_tol_orbit => 1e-10,                                    # Absolute tolerance for orbit integration (outside atmosphere, i.e., step 1 and step 3)
            :r_tol_orbit => 1e-8,                                    # Relative tolerance for orbit integration (outside atmosphere, i.e., step 1 and step 3)
            :a_tol_drag => 1e-10,                                       # Absolute tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :r_tol_drag => 1e-8,                                       # Relative tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :a_tol_quaternion => 1e-11,                                  # Absolute tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :r_tol_quaternion => 1e-9,                                  # Relative tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :dt_max => 1.0,                                         # Maximum time step for integration, s
            :dt_max_orbit => 1.0,                                   # Maximum time step for orbit integration (outside atmosphere, i.e., step 1 and step 3), s
            :dt_max_drag => 1.0,                                    # Maximum time step for drag passage

            :Odyssey_sim => 0                                      # Simulate Odyssey Mission
            )

# # Calculating time of simulation
# @profview run_analysis(args)
mc_runs = 1000
mc_dispersion_pct = 0.001 # Should be roughly on the order of kms for position values
mc_dispersion_deg = 0.1

nominal_a = args[:a_initial_a]
nominal_e = args[:e_initial_a]
nominal_i = args[:inclination]
nominal_ω = args[:ω]
nominal_Ω = args[:Ω]
nominal_ν = args[:ν]

# for i in 92:mc_runs
#     println("Starting MC $i/$mc_runs")
    t = @elapsed begin          
        # Run the simulation
        # args[:monte_carlo_run] = i
        # args[:a_initial_a] = nominal_a * (1 + rand((-1, 1))*rand()*mc_dispersion_pct)
        # args[:e_initial_a] = nominal_e * (1 + rand((-1, 1))*rand()*mc_dispersion_pct)
        # args[:i] = nominal_i + rand((-1, 1))*rand()*mc_dispersion_deg
        # args[:ω] = nominal_ω + rand((-1, 1))*rand()*mc_dispersion_deg
        # args[:Ω] = nominal_Ω + rand((-1, 1))*rand()*mc_dispersion_deg
        # args[:ν] = nominal_ν * (1 + rand((-1, 1))*rand()*mc_dispersion)
        sol = run_analysis(args)
        if Bool(args[:passresults])
            println("Ra initial = " * string((sol.orientation.oe[1][1] * (1 + sol.orientation.oe[2][1]))* 1e-3) * " km, Ra new = " * string((sol.orientation.oe[1][end] * (1 + sol.orientation.oe[2][end]))* 1e-3) * " km - Actual periapsis altitude = " * string(minimum(sol.orientation.alt) * 1e-3) * " km - Target Ra = " * string(args[:final_apoapsis] * 1e-3) * " km")
        end
    end

    println("COMPUTATIONAL TIME = " * string(t) * " s")
# end