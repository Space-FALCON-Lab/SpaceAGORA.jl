using LinearAlgebra
using StaticArrays
using SparseArrays

### Quaternion Function ###
# Quaternion multiplication
function quat_mult(q::AbstractVector, p::AbstractVector)
    q0, qv = q[4], q[1:3]
    p0, pv = p[4], p[1:3]

    s = q0*p0 - dot(qv, pv)
    v = q0*pv + p0*qv + cross(qv, pv)
    return [v; s]
end

## skew matrix
hatIs = [1; 1; 2; 2; 3; 3]
hatJs = [2; 3; 1; 3; 1; 2]
function hat(ω::AbstractVector{<:Float64})
    """
    skew matrix, The matrix equivalent of cross product
    params: 
    ω -> vector ∈ R³
    return:
    [ω]x 
    """

    # Vs = [-ω[3]; ω[2]; ω[3]; -ω[1]; -ω[2]; ω[1]]
    return SMatrix{3, 3, Float64}([0 -ω[3] ω[2];
                                   ω[3] 0 -ω[1];
                                   -ω[2] ω[1] 0])
    # return [0 -ω[3] ω[2]
    #     ω[3] 0 -ω[1]
    #     -ω[2] ω[1] 0]
end

function L!(Lmat, Q)
    Lmat .= [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I+hat(Q[2:4])]
end

## Left
lIs = [1; 1; 1; 1; 2; 2; 2; 2; 3; 3; 3; 3; 4; 4; 4; 4]
lJs = [1; 2; 3; 4; 1; 2; 3; 4; 1; 2; 3; 4; 1; 2; 3; 4]
function L(Q)
    """
    Left Multiply quaternion
    q2 ⨂ q1 -> L(q2)q1
    params:
    Q -> unit quaternion
    returns: 
    L(Q) -> 4x4 matrix
    """
    # Lmat = zeros(eltype(Q), 4, 4)
    # L!(Lmat, Q)
    tmp = hat(Q[2:4])
    Vs = [Q[1]; -Q[2]; -Q[3]; -Q[4]; Q[2]; Q[1]; tmp[1, 2]; tmp[1, 3]; Q[3]; tmp[2, 1]; Q[1]; tmp[2, 3]; Q[4]; tmp[3, 1]; tmp[3,2]; Q[1]]
    return sparse(lIs, lJs, Vs)
    # return [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I+hat(Q[2:4])]
end

function R!(Rmat, Q)
    Rmat .= [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I-hat(Q[2:4])]
end
## Right
function R(Q)
    """
    Right Multiply quaternion
    q2 ⨂ q1 -> R(q1)q2
    params:
    Q -> unit quaternion
    returns: 
    R(Q) -> 4x4 matrix
    """
    # Rmat = zeros(eltype(Q), 4, 4)
    # R!(Rmat, Q)
    # s = @view Q[1]
    # v = @view Q[2:4]
    tmp = -hat(Q[1:3])
    Vs = [Q[1]; -Q[2]; -Q[3]; -Q[4]; Q[2]; Q[1]; tmp[1, 2]; tmp[1, 3]; Q[3]; tmp[2, 1]; Q[1]; tmp[2, 3]; Q[4]; tmp[3, 1]; tmp[3,2]; Q[1]]
    return sparse(lIs, lJs, Vs)
    # return [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I-hat(Q[2:4])]
end

function rot(q)
    """
    Convert quaternion to rotation matrix
    """
    q1, q2, q3, q4 = q
    return SMatrix{3, 3, Float64}([q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4);
                                    2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4);
                                    2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2])
end

function error_quaternion(current::SVector{4, Float64}, target::SVector{4, Float64})

    q = current[4]
    q_bar = current[1:3]'
    q_bar_cross = hat(SVector{3, Float64}(q_bar))
    q_matrix = SMatrix{4, 4, Float64}(Float64[q*I(3)-q_bar_cross q_bar'; -q_bar q])

    qd = SVector{4, Float64}(Float64[-target[1:3]; target[4]])
    return normalize(q_matrix * qd)

    # q = -current
    # q_bar = current[1:3]'
    # q_bar_cross = [0 -q_bar[3] q_bar[2]; q_bar[3] 0 -q_bar[1]; -q_bar[2] q_bar[1] 0]
    # q_matrix = [q[4]*I(3)-q_bar_cross q_bar'; -q_bar q[4]]

    # qd = [-target[1]; -target[2]; -target[3]; target[4]]
    # qe_neg = q_matrix * qd
    # qe = qe_pos[4] > qe_neg[4] ? qe_pos : qe_neg
    # return SVector{4, Float64}(normalize(qe))
end

## take angle components (4, 3)
H = SMatrix{4,3}([zeros(1, 3); I]);
S = SMatrix{4,3}([[1 0 0]; zeros(3,3)])
## Inverse
T = SMatrix{4,4}(Diagonal([1.0; -1; -1; -1]))

## Attitude Jacobian
function G(Q)
    return L(Q) * H
end

function Q̄(q1, q2)
    """
    Average quaternion
    params:
    q1 -> unit quaternion
    q2 -> unit quaternion
    returns: 
    q̄ -> unit quaternion
    """
    q̂ = L(q1)' * q2
    normv = norm(q̂[1:3])
    if normv == 0.0
        return q1
    else
        r = q̂[1:3] / normv
        θ = (2 * atan(normv, q̂[4]))
    end
    # r = q̂[2:4] / (norm(q̂[2:4]) + 1e-10)
    # θ = 2 * atan(norm(q̂[2:4]), q̂[1])
    q̄ = L(q1) * [cos(θ / 4); sin(θ / 4) * r]
    # values, vectors = eigen(0.5.*([q1 q2] * [q1 q2]'))
    q̄ ./= norm(q̄)

    return q̄

end


function phi_from_q(q)
    """
    axis angle from quaternion (scalar last)
    params: 
    q -> unit quaternion
    returns: 
    ϕ -> axis (r) * angle (θ) 
    """
    # cayley_map(q)
    v = @view q[1:3]
    # v = SVector(q[2], q[3], q[4])
    s = q[4]
    normv = norm(v)

    # if isapprox(normv, 0, atol=eps())
    r = v / (normv + eps())
    θ = (2 * atan(normv, s))
    return r * θ
    # end

    # return zeros(3)
end

function q_from_phi(ϕ)
    """
    axis angle from quaternion (scalar last)
    params: 
    ϕ -> axis (r) * angle (θ) 
    returns: 
    q -> unit quaternion
    """

    θ = norm(ϕ)

    if θ == 0.0
        return [1, 0, 0, 0]
    else
        r = ϕ / θ
        return [cos(θ/2); r*sin(θ/2)]
    end
end

function dphi_from_q(q)
    v = @views q[1:3]
    # v = SVector(q[2], q[3], q[4])
    s = q[4]
    normv = norm(v)
    return [-2*v/(s^2*(1.0 + normv^2/s^2)) (2 * atan(normv, s))/(normv+1e-10)*Matrix(I, 3,3) .- v*v'*(2 * atan(normv, s))/(normv+1e-10)^3 .+ 2*v*v'/((normv+1e-10)*s*(normv+1e-10)*(1.0 + normv^2/s^2))] 
end

function cayley_map(q) 
    return q[1:3]/q[4]
end

function qToEulerAngles(q)
    # // this implementation assumes normalized quaternion
    # // converts to Euler angles in 3-2-1 sequence
    qx, qy, qz, qw = q
    # roll (x-axis rotation)
    sinr_cosp = 2 * (qw * qx + qy * qz);
    cosr_cosp = 1 - 2 * (qx * qx + qy * qy);
    roll = atan(cosr_cosp,sinr_cosp);# - π/2;

    # pitch (y-axis rotation)
    sinp = sqrt(1 + 2 * (qw * qy - qx * qz));
    cosp = sqrt(1 - 2 * (qw * qy - qx * qz));
    pitch = 2 * atan(cosp,sinp) - π/2;

    # yaw (z-axis rotation)
    siny_cosp = 2 * (qw * qz + qx * qy);
    cosy_cosp = 1 - 2 * (qy * qy + qz * qz);
    yaw = atan(cosy_cosp,siny_cosp);# - π/2;

    return [roll, pitch, yaw]
end


function EulerAnglesToq(θ)
    roll, pitch, yaw = θ[1],θ[2],θ[3]
    cr = cos(roll * 0.5)
    sr = sin(roll * 0.5)
    cp = cos(pitch * 0.5)
    sp = sin(pitch * 0.5)
    cy = cos(yaw * 0.5)
    sy = sin(yaw * 0.5)

    qw = cr * cp * cy + sr * sp * sy
    qx = sr * cp * cy - cr * sp * sy
    qy = cr * sp * cy + sr * cp * sy
    qz = cr * cp * sy - sr * sp * cy

    # check for norm close to 1 (allowing a small tolerance)
    norm = sqrt(qw^2 + qx^2 + qy^2 + qz^2)
    if abs(norm - 1) > 1e-6
        println("Warning: Quaternion is not normalized. Norm is: ", norm)
    end
    return [qx, qy, qz, qw]
end

function rotation_matrix(link)
    q = link.q  # Assuming link.q stores the quaternion
    # Quaternion components
    qx, qy, qz, qw = q[1], q[2], q[3], q[4]
    
    # Rotation matrix from quaternion
    R = [1 - 2*qy^2 - 2*qz^2    2*qx*qy - 2*qz*qw    2*qx*qz + 2*qy*qw;
         2*qx*qy + 2*qz*qw    1 - 2*qx^2 - 2*qz^2    2*qy*qz - 2*qx*qw;
         2*qx*qz - 2*qy*qw    2*qy*qz + 2*qx*qw    1 - 2*qx^2 - 2*qy^2]
    return R
end

function dcm_to_quaternion(dcm::SMatrix{3, 3, Float64})
    """
    Convert a Direction Cosine Matrix (DCM) to a quaternion.
    
    # Arguments
    - `dcm`: A 3x3 Direction Cosine Matrix (DCM) represented as a StaticVector.
    
    # Returns
    - A quaternion represented as a StaticVector{4, Float64}.
    """
    tr = dcm[1, 1] + dcm[2, 2] + dcm[3, 3]
    i = argmax([dcm[1, 1], dcm[2, 2], dcm[3, 3], tr])
    if i == 4
        q = SVector{4, Float64}([dcm[2, 3] - dcm[3, 2], dcm[3, 1] - dcm[1, 3], dcm[1, 2] - dcm[2, 1], 1+tr])
        return q / norm(q)
    else
        q = zeros(4)
        q[i] = 1+2*dcm[i, i] - tr
        for j in 1:4
            if j != i
                if j != 4
                    q[j] = dcm[i, j] + dcm[j, i]
                else
                    idx = i+1 > 3 ? (i+1)%3 : i+1
                    second_idx = i+2 > 3 ? (i+2)%3 : i+2
                    
                    q[j] = dcm[idx, second_idx] - dcm[second_idx, idx]
                end
            end
        end
        return q / norm(q)
    end
end

function Ψ(q::SVector{4, Float64})
    """
    Convert quaternion to Euler angles in 3-2-1 sequence.

    # Arguments
    - `q`: A quaternion represented as a StaticVector{4, Float64}.

    # Returns
    - A StaticVector{3, Float64} containing the Euler angles (roll, pitch, yaw).
    """
    return SMatrix{4, 3, Float64}([q[4]*I(3) - hat(q[1:3]); -q[1:3]'])
end

function Ξ(q::SVector{4, Float64})
    """
    Convert quaternion to kinematics matrix Ξ.

    # Arguments
    - `q`: A quaternion represented as a StaticVector{4, Float64}.

    # Returns
    - A StaticMatrix{4, 3, Float64} containing the quaternion kinematics matrix.
    """
    return SMatrix{4, 3, Float64}([q[4]*I(3) + hat(q[1:3]); -q[1:3]'])
end
