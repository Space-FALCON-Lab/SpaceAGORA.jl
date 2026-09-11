#!/usr/bin/env julia

using LinearAlgebra
using Printf
using StaticArrays

using PlotlyJS
using SpaceAGORA

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT_DIR = joinpath(REPO_ROOT, "output", "solar_panel_cloth_deployment_demo")

const PANEL_COUNT = 4
const PANEL_TILES_X = 3
const PANEL_TILES_Y = 4
const PANEL_LENGTH_M = 0.42
const PANEL_WIDTH_M = 0.28
const PANEL_THICKNESS_M = 0.008
const PANEL_MASS_KG = 0.18
const TILE_LENGTH_M = PANEL_LENGTH_M / PANEL_TILES_X
const TILE_WIDTH_M = PANEL_WIDTH_M / PANEL_TILES_Y
const TILE_MASS_KG = PANEL_MASS_KG / (PANEL_TILES_X * PANEL_TILES_Y)
const TILE_INERTIA_FLOOR_KG_M2 = 2.0e-3
const DEPLOY_DURATION_S = 6.0
const SIM_DURATION_S = 60.0
const ROTATION_DURATION_S = 8.0
const ROTATION_SIM_DURATION_S = 30.0
const ROTATION_TARGET_RAD = 0.5 * pi
const DT_S = 0.005

@inline function _unit_quat(q)
    qv = SVector{4, Float64}(q)
    nq = norm(qv)
    return nq > eps(Float64) ? qv / nq : SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
end

@inline function _quat_conj(q)
    qv = _unit_quat(q)
    return SVector{4, Float64}(-qv[1], -qv[2], -qv[3], qv[4])
end

@inline function _quat_mul(a, b)
    aq = _unit_quat(a)
    bq = _unit_quat(b)
    ax, ay, az, aw = aq
    bx, by, bz, bw = bq
    return _unit_quat(SVector{4, Float64}(
        aw * bx + bw * ax + ay * bz - az * by,
        aw * by + bw * ay + az * bx - ax * bz,
        aw * bz + bw * az + ax * by - ay * bx,
        aw * bw - ax * bx - ay * by - az * bz,
    ))
end

@inline function _rot(q_raw)
    q = _unit_quat(q_raw)
    x, y, z, w = q
    return SMatrix{3, 3, Float64}(
        1 - 2y^2 - 2z^2, 2x * y + 2z * w, 2x * z - 2y * w,
        2x * y - 2z * w, 1 - 2x^2 - 2z^2, 2y * z + 2x * w,
        2x * z + 2y * w, 2y * z - 2x * w, 1 - 2x^2 - 2y^2,
    )
end

@inline _quat_y(theta_rad) =
    SVector{4, Float64}(0.0, sin(0.5 * theta_rad), 0.0, cos(0.5 * theta_rad))

@inline _quat_x(theta_rad) =
    SVector{4, Float64}(sin(0.5 * theta_rad), 0.0, 0.0, cos(0.5 * theta_rad))

@inline function _smoothstep(u)
    s = clamp(Float64(u), 0.0, 1.0)
    return s * s * (3.0 - 2.0 * s)
end

@inline function _node_index(panel_idx::Int, tile_x::Int, tile_y::Int)
    return ((panel_idx - 1) * PANEL_TILES_X * PANEL_TILES_Y) +
        ((tile_y - 1) * PANEL_TILES_X) +
        tile_x
end

function _panel_state_from_relative_angles(relative_angles)
    positions = SVector{3, Float64}[]
    quaternions = SVector{4, Float64}[]
    q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    r = SVector{3, Float64}(0.5 * PANEL_LENGTH_M, 0.0, 0.0)
    push!(positions, r)
    push!(quaternions, q)

    for theta in relative_angles
        parent_r = positions[end]
        parent_q = quaternions[end]
        hinge_world = parent_r + _rot(parent_q) * SVector{3, Float64}(0.5 * PANEL_LENGTH_M, 0.0, 0.0)
        q = _quat_mul(parent_q, _quat_y(theta))
        r = hinge_world - _rot(q) * SVector{3, Float64}(-0.5 * PANEL_LENGTH_M, 0.0, 0.0)
        push!(positions, r)
        push!(quaternions, q)
    end
    return positions, quaternions
end

function _tile_center(panel_center, panel_q, tile_x::Int, tile_y::Int)
    local_x = -0.5 * PANEL_LENGTH_M + (tile_x - 0.5) * TILE_LENGTH_M
    local_y = -0.5 * PANEL_WIDTH_M + (tile_y - 0.5) * TILE_WIDTH_M
    return panel_center + _rot(panel_q) * SVector{3, Float64}(local_x, local_y, 0.0)
end

function _edge_rest_quats(model, hinge_joint_groups, initial_relative_angles, t_s)
    deploy_fraction = _smoothstep(t_s / DEPLOY_DURATION_S)
    rests = [joint.rest_child_parent_quat for joint in model.joints]
    for hinge_idx in 1:(PANEL_COUNT - 1)
        rest = _quat_y((1.0 - deploy_fraction) * initial_relative_angles[hinge_idx])
        for joint_idx in hinge_joint_groups[hinge_idx]
            rests[joint_idx] = rest
        end
    end
    return rests
end

function _build_folded_panel_model()
    initial_relative_angles = deg2rad.([105.0, -105.0, 105.0])
    panel_positions, panel_quats = _panel_state_from_relative_angles(initial_relative_angles)
    physical_tile_inertia = SpaceAGORA.rectangular_prism_inertia(
        TILE_MASS_KG,
        (TILE_LENGTH_M, TILE_WIDTH_M, PANEL_THICKNESS_M),
    )
    tile_inertia = SMatrix{3, 3, Float64}(
        max(physical_tile_inertia[1, 1], TILE_INERTIA_FLOOR_KG_M2), 0.0, 0.0,
        0.0, max(physical_tile_inertia[2, 2], TILE_INERTIA_FLOOR_KG_M2), 0.0,
        0.0, 0.0, max(physical_tile_inertia[3, 3], TILE_INERTIA_FLOOR_KG_M2),
    )

    nodes = SpaceAGORA.CompliantTopologyNode[]
    for panel_idx in 1:PANEL_COUNT, tile_y in 1:PANEL_TILES_Y, tile_x in 1:PANEL_TILES_X
        push!(nodes, SpaceAGORA.CompliantTopologyNode(
            Symbol("solar_panel_$(panel_idx)_tile_$(tile_x)_$(tile_y)");
            mass_kg=TILE_MASS_KG,
            inertia_body_kg_m2=tile_inertia,
            position=_tile_center(panel_positions[panel_idx], panel_quats[panel_idx], tile_x, tile_y),
            quaternion=panel_quats[panel_idx],
        ))
    end

    edges = SpaceAGORA.CompliantTopologyEdge[]
    hinge_joint_groups = [Int[] for _ in 1:(PANEL_COUNT - 1)]
    root_joint_indices = Int[]

    internal_Kx = SMatrix{3, 3, Float64}(4.0I)
    internal_Cx = SMatrix{3, 3, Float64}(0.16I)
    internal_Kr = SMatrix{3, 3, Float64}(Diagonal([0.08, 0.12, 0.08]))
    internal_Cr = SMatrix{3, 3, Float64}(Diagonal([0.018, 0.026, 0.018]))
    hinge_Kx = SMatrix{3, 3, Float64}(6.0I)
    hinge_Cx = SMatrix{3, 3, Float64}(0.32I)
    hinge_Kr = SMatrix{3, 3, Float64}(Diagonal([0.05, 0.20, 0.05]))
    hinge_Cr = SMatrix{3, 3, Float64}(Diagonal([0.020, 0.060, 0.020]))

    for tile_y in 1:PANEL_TILES_Y
        push!(edges, SpaceAGORA.CompliantTopologyEdge(
            Symbol("root_hinge_row_$tile_y"),
            0,
            _node_index(1, 1, tile_y);
            parent_point_body=SVector{3, Float64}(0.0, -0.5 * PANEL_WIDTH_M + (tile_y - 0.5) * TILE_WIDTH_M, 0.0),
            child_point_body=SVector{3, Float64}(-0.5 * TILE_LENGTH_M, 0.0, 0.0),
            k_translation_n_m=SMatrix{3, 3, Float64}(8.0I),
            c_translation_n_s_m=SMatrix{3, 3, Float64}(0.38I),
            k_rotation_n_m_rad=SMatrix{3, 3, Float64}(0.25I),
            c_rotation_n_m_s_rad=SMatrix{3, 3, Float64}(0.070I),
        ))
        push!(root_joint_indices, length(edges))
    end

    for panel_idx in 1:PANEL_COUNT
        for tile_y in 1:PANEL_TILES_Y
            for tile_x in 1:(PANEL_TILES_X - 1)
                push!(edges, SpaceAGORA.CompliantTopologyEdge(
                    Symbol("panel_$(panel_idx)_x_edge_$(tile_x)_$(tile_y)"),
                    _node_index(panel_idx, tile_x, tile_y),
                    _node_index(panel_idx, tile_x + 1, tile_y);
                    parent_point_body=SVector{3, Float64}(0.5 * TILE_LENGTH_M, 0.0, 0.0),
                    child_point_body=SVector{3, Float64}(-0.5 * TILE_LENGTH_M, 0.0, 0.0),
                    k_translation_n_m=internal_Kx,
                    c_translation_n_s_m=internal_Cx,
                    k_rotation_n_m_rad=internal_Kr,
                    c_rotation_n_m_s_rad=internal_Cr,
                ))
            end
        end
        for tile_x in 1:PANEL_TILES_X
            for tile_y in 1:(PANEL_TILES_Y - 1)
                push!(edges, SpaceAGORA.CompliantTopologyEdge(
                    Symbol("panel_$(panel_idx)_y_edge_$(tile_x)_$(tile_y)"),
                    _node_index(panel_idx, tile_x, tile_y),
                    _node_index(panel_idx, tile_x, tile_y + 1);
                    parent_point_body=SVector{3, Float64}(0.0, 0.5 * TILE_WIDTH_M, 0.0),
                    child_point_body=SVector{3, Float64}(0.0, -0.5 * TILE_WIDTH_M, 0.0),
                    k_translation_n_m=internal_Kx,
                    c_translation_n_s_m=internal_Cx,
                    k_rotation_n_m_rad=internal_Kr,
                    c_rotation_n_m_s_rad=internal_Cr,
                ))
            end
        end
    end

    for hinge_idx in 1:(PANEL_COUNT - 1)
        parent_panel = hinge_idx
        child_panel = hinge_idx + 1
        for tile_y in 1:PANEL_TILES_Y
            push!(edges, SpaceAGORA.CompliantTopologyEdge(
                Symbol("deploy_hinge_$(parent_panel)_$(child_panel)_row_$tile_y"),
                _node_index(parent_panel, PANEL_TILES_X, tile_y),
                _node_index(child_panel, 1, tile_y);
                parent_point_body=SVector{3, Float64}(0.5 * TILE_LENGTH_M, 0.0, 0.0),
                child_point_body=SVector{3, Float64}(-0.5 * TILE_LENGTH_M, 0.0, 0.0),
                k_translation_n_m=hinge_Kx,
                c_translation_n_s_m=hinge_Cx,
                k_rotation_n_m_rad=hinge_Kr,
                c_rotation_n_m_s_rad=hinge_Cr,
            ))
            push!(hinge_joint_groups[hinge_idx], length(edges))
        end
    end

    build = SpaceAGORA.build_compliant_topology(nodes, edges)
    actuators = SpaceAGORA.CompliantJointActuator[]
    for hinge_idx in 1:(PANEL_COUNT - 1), joint_idx in hinge_joint_groups[hinge_idx]
        push!(actuators, SpaceAGORA.CompliantJointActuator(
            Symbol("hinge_motor_$(hinge_idx)_joint_$joint_idx"),
            joint_idx;
            torque_limit_n_m=(0.001, 0.003, 0.001),
            kp_n_m_rad=SMatrix{3, 3, Float64}(Diagonal([0.0, 0.35, 0.0])),
            kd_n_m_s_rad=SMatrix{3, 3, Float64}(Diagonal([0.0, 0.075, 0.0])),
        ))
    end
    return build, initial_relative_angles, hinge_joint_groups, root_joint_indices, actuators
end

function _tile_inertia_matrix()
    physical_tile_inertia = SpaceAGORA.rectangular_prism_inertia(
        TILE_MASS_KG,
        (TILE_LENGTH_M, TILE_WIDTH_M, PANEL_THICKNESS_M),
    )
    return SMatrix{3, 3, Float64}(
        max(physical_tile_inertia[1, 1], TILE_INERTIA_FLOOR_KG_M2), 0.0, 0.0,
        0.0, max(physical_tile_inertia[2, 2], TILE_INERTIA_FLOOR_KG_M2), 0.0,
        0.0, 0.0, max(physical_tile_inertia[3, 3], TILE_INERTIA_FLOOR_KG_M2),
    )
end

function _build_end_point_x_axis_rotation_model()
    tile_inertia = _tile_inertia_matrix()
    panel_centers = [
        SVector{3, Float64}((panel_idx - 0.5 * (PANEL_COUNT + 1)) * PANEL_LENGTH_M, 0.0, 0.0)
        for panel_idx in 1:PANEL_COUNT
    ]
    panel_q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

    nodes = SpaceAGORA.CompliantTopologyNode[]
    for panel_idx in 1:PANEL_COUNT, tile_y in 1:PANEL_TILES_Y, tile_x in 1:PANEL_TILES_X
        push!(nodes, SpaceAGORA.CompliantTopologyNode(
            Symbol("deployed_panel_$(panel_idx)_tile_$(tile_x)_$(tile_y)");
            mass_kg=TILE_MASS_KG,
            inertia_body_kg_m2=tile_inertia,
            position=_tile_center(panel_centers[panel_idx], panel_q, tile_x, tile_y),
            quaternion=panel_q,
        ))
    end

    internal_Kx = SMatrix{3, 3, Float64}(5.0I)
    internal_Cx = SMatrix{3, 3, Float64}(0.20I)
    internal_Kr = SMatrix{3, 3, Float64}(Diagonal([0.10, 0.16, 0.10]))
    internal_Cr = SMatrix{3, 3, Float64}(Diagonal([0.022, 0.035, 0.022]))
    axis_Kx = SMatrix{3, 3, Float64}(12.0I)
    axis_Cx = SMatrix{3, 3, Float64}(0.55I)
    axis_Kr = SMatrix{3, 3, Float64}(Diagonal([0.48, 0.08, 0.08]))
    axis_Cr = SMatrix{3, 3, Float64}(Diagonal([0.140, 0.025, 0.025]))

    edges = SpaceAGORA.CompliantTopologyEdge[]
    for panel_idx in 1:PANEL_COUNT
        for tile_y in 1:PANEL_TILES_Y
            for tile_x in 1:(PANEL_TILES_X - 1)
                push!(edges, SpaceAGORA.CompliantTopologyEdge(
                    Symbol("rotation_panel_$(panel_idx)_x_edge_$(tile_x)_$(tile_y)"),
                    _node_index(panel_idx, tile_x, tile_y),
                    _node_index(panel_idx, tile_x + 1, tile_y);
                    parent_point_body=SVector{3, Float64}(0.5 * TILE_LENGTH_M, 0.0, 0.0),
                    child_point_body=SVector{3, Float64}(-0.5 * TILE_LENGTH_M, 0.0, 0.0),
                    k_translation_n_m=internal_Kx,
                    c_translation_n_s_m=internal_Cx,
                    k_rotation_n_m_rad=internal_Kr,
                    c_rotation_n_m_s_rad=internal_Cr,
                ))
            end
        end
        for tile_x in 1:PANEL_TILES_X
            for tile_y in 1:(PANEL_TILES_Y - 1)
                push!(edges, SpaceAGORA.CompliantTopologyEdge(
                    Symbol("rotation_panel_$(panel_idx)_y_edge_$(tile_x)_$(tile_y)"),
                    _node_index(panel_idx, tile_x, tile_y),
                    _node_index(panel_idx, tile_x, tile_y + 1);
                    parent_point_body=SVector{3, Float64}(0.0, 0.5 * TILE_WIDTH_M, 0.0),
                    child_point_body=SVector{3, Float64}(0.0, -0.5 * TILE_WIDTH_M, 0.0),
                    k_translation_n_m=internal_Kx,
                    c_translation_n_s_m=internal_Cx,
                    k_rotation_n_m_rad=internal_Kr,
                    c_rotation_n_m_s_rad=internal_Cr,
                ))
            end
        end
    end

    for panel_idx in 1:(PANEL_COUNT - 1), tile_y in 1:PANEL_TILES_Y
        push!(edges, SpaceAGORA.CompliantTopologyEdge(
            Symbol("rotation_panel_join_$(panel_idx)_$(panel_idx + 1)_row_$tile_y"),
            _node_index(panel_idx, PANEL_TILES_X, tile_y),
            _node_index(panel_idx + 1, 1, tile_y);
            parent_point_body=SVector{3, Float64}(0.5 * TILE_LENGTH_M, 0.0, 0.0),
            child_point_body=SVector{3, Float64}(-0.5 * TILE_LENGTH_M, 0.0, 0.0),
            k_translation_n_m=internal_Kx,
            c_translation_n_s_m=internal_Cx,
            k_rotation_n_m_rad=internal_Kr,
            c_rotation_n_m_s_rad=internal_Cr,
        ))
    end

    drive_joint_indices = Int[]
    total_panel_length = PANEL_COUNT * PANEL_LENGTH_M
    push!(edges, SpaceAGORA.CompliantTopologyEdge(
        :end_point_x_axis_drive,
        0,
        _node_index(1, 1, 2);
        parent_point_body=SVector{3, Float64}(-0.5 * total_panel_length, 0.0, 0.0),
        child_point_body=SVector{3, Float64}(-0.5 * TILE_LENGTH_M, 0.5 * TILE_WIDTH_M, 0.0),
        k_translation_n_m=axis_Kx,
        c_translation_n_s_m=axis_Cx,
        k_rotation_n_m_rad=axis_Kr,
        c_rotation_n_m_s_rad=axis_Cr,
    ))
    push!(drive_joint_indices, length(edges))

    build = SpaceAGORA.build_compliant_topology(nodes, edges)
    actuators = SpaceAGORA.CompliantJointActuator[
        SpaceAGORA.CompliantJointActuator(
            Symbol("end_point_x_axis_motor_$joint_idx"),
            joint_idx;
            torque_limit_n_m=(0.008, 0.0015, 0.0015),
            kp_n_m_rad=SMatrix{3, 3, Float64}(Diagonal([0.65, 0.0, 0.0])),
            kd_n_m_s_rad=SMatrix{3, 3, Float64}(Diagonal([0.16, 0.0, 0.0])),
        ) for joint_idx in drive_joint_indices
    ]
    return build, drive_joint_indices, actuators
end

function _signed_relative_x_angle(parent_q, child_q)
    q_rel = _quat_mul(_quat_conj(parent_q), child_q)
    v_norm = norm(q_rel[1:3])
    angle = 2.0 * atan(v_norm, q_rel[4])
    sign_x = abs(q_rel[1]) > 1.0e-12 ? sign(q_rel[1]) : 1.0
    return sign_x * angle
end

function _signed_relative_y_angle(parent_q, child_q)
    q_rel = _quat_mul(_quat_conj(parent_q), child_q)
    v_norm = norm(q_rel[1:3])
    angle = 2.0 * atan(v_norm, q_rel[4])
    sign_y = abs(q_rel[2]) > 1.0e-12 ? sign(q_rel[2]) : 1.0
    return sign_y * angle
end

function _relative_hinge_angles(x, hinge_joint_groups)
    angles = zeros(PANEL_COUNT - 1)
    for hinge_idx in 1:(PANEL_COUNT - 1)
        row_angles = Float64[]
        for tile_y in 1:PANEL_TILES_Y
            parent = SpaceAGORA.compliant_state_parts(x, _node_index(hinge_idx, PANEL_TILES_X, tile_y))
            child = SpaceAGORA.compliant_state_parts(x, _node_index(hinge_idx + 1, 1, tile_y))
            push!(row_angles, _signed_relative_y_angle(parent.q, child.q))
        end
        angles[hinge_idx] = sum(row_angles) / length(row_angles)
    end
    return angles
end

function _tile_corners(x, node_idx)
    state = SpaceAGORA.compliant_state_parts(x, node_idx)
    R = _rot(state.q)
    local_corners = SVector{3, Float64}[
        SVector{3, Float64}(-0.5 * TILE_LENGTH_M, -0.5 * TILE_WIDTH_M, 0.0),
        SVector{3, Float64}(0.5 * TILE_LENGTH_M, -0.5 * TILE_WIDTH_M, 0.0),
        SVector{3, Float64}(0.5 * TILE_LENGTH_M, 0.5 * TILE_WIDTH_M, 0.0),
        SVector{3, Float64}(-0.5 * TILE_LENGTH_M, 0.5 * TILE_WIDTH_M, 0.0),
    ]
    return [state.r + R * p for p in local_corners]
end

function _panel_tip(x)
    pts = SVector{3, Float64}[]
    for tile_y in 1:PANEL_TILES_Y
        corners = _tile_corners(x, _node_index(PANEL_COUNT, PANEL_TILES_X, tile_y))
        push!(pts, 0.5 .* (corners[2] .+ corners[3]))
    end
    return sum(pts) / length(pts)
end

function _simulate_deployment(build, initial_relative_angles, hinge_joint_groups, actuators)
    times = collect(0.0:DT_S:SIM_DURATION_S)
    states = Vector{Vector{Float64}}(undef, length(times))
    states[1] = copy(build.initial_state)
    actual_angles = zeros(PANEL_COUNT - 1, length(times))
    commanded_angles = zeros(PANEL_COUNT - 1, length(times))
    compliance_torque_y = zeros(PANEL_COUNT - 1, length(times))
    actuator_torque_y = zeros(PANEL_COUNT - 1, length(times))
    saturation = zeros(PANEL_COUNT - 1, length(times))
    tip = zeros(3, length(times))

    for (k, t) in pairs(times)
        rests = _edge_rest_quats(build.model, hinge_joint_groups, initial_relative_angles, t)
        actual_angles[:, k] .= rad2deg.(_relative_hinge_angles(states[k], hinge_joint_groups))
        commanded_angles[:, k] .= rad2deg.((1.0 - _smoothstep(t / DEPLOY_DURATION_S)) .* initial_relative_angles)
        loads = SpaceAGORA.compliant_joint_loads(build.model, states[k]; joint_rest_quaternions=rests, joint_actuators=actuators)
        for hinge_idx in 1:(PANEL_COUNT - 1)
            joint_ids = hinge_joint_groups[hinge_idx]
            compliance_torque_y[hinge_idx, k] = sum(loads[j].compliance_torque_child_body[2] for j in joint_ids)
            actuator_torque_y[hinge_idx, k] = sum(loads[j].actuator_torque_child_body[2] for j in joint_ids)
            hinge_actuators = filter(a -> a.joint in joint_ids, actuators)
            limit_sum = sum(a.torque_limit_n_m[2] for a in hinge_actuators)
            saturation[hinge_idx, k] = limit_sum > 0.0 ? abs(actuator_torque_y[hinge_idx, k]) / limit_sum : 0.0
        end
        tip[:, k] .= _panel_tip(states[k])

        if k < length(times)
            states[k + 1] = SpaceAGORA.step_compliant_multibody_rk4(
                build.model,
                states[k],
                t,
                times[k + 1] - t;
                joint_rest_quaternions=rests,
                joint_actuators=actuators,
            )
        end
    end

    return (
        times=times,
        states=states,
        actual_angles_deg=actual_angles,
        commanded_angles_deg=commanded_angles,
        compliance_torque_y_nm=compliance_torque_y,
        actuator_torque_y_nm=actuator_torque_y,
        actuator_saturation=saturation,
        tip_m=tip,
    )
end

function _end_point_x_axis_rest_quats(model, drive_joint_indices, t_s)
    slew_fraction = _smoothstep(t_s / ROTATION_DURATION_S)
    rest = _quat_x(slew_fraction * ROTATION_TARGET_RAD)
    rests = [joint.rest_child_parent_quat for joint in model.joints]
    for joint_idx in drive_joint_indices
        rests[joint_idx] = rest
    end
    return rests
end

function _average_panel_rotation_deg(x)
    angles = Float64[]
    for panel_idx in 1:PANEL_COUNT, tile_y in 1:PANEL_TILES_Y, tile_x in 1:PANEL_TILES_X
        state = SpaceAGORA.compliant_state_parts(x, _node_index(panel_idx, tile_x, tile_y))
        push!(angles, _signed_relative_x_angle(SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), state.q))
    end
    return rad2deg(sum(angles) / length(angles))
end

function _free_end_edge_span(x)
    lower_corners = _tile_corners(x, _node_index(PANEL_COUNT, PANEL_TILES_X, 1))
    upper_corners = _tile_corners(x, _node_index(PANEL_COUNT, PANEL_TILES_X, PANEL_TILES_Y))
    return lower_corners[2], upper_corners[3]
end

function _simulate_end_point_x_axis_rotation(build, drive_joint_indices, actuators)
    times = collect(0.0:DT_S:ROTATION_SIM_DURATION_S)
    states = Vector{Vector{Float64}}(undef, length(times))
    states[1] = copy(build.initial_state)
    commanded_angle = zeros(length(times))
    actual_angle = zeros(length(times))
    compliance_torque_x = zeros(length(times))
    actuator_torque_x = zeros(length(times))
    saturation = zeros(length(times))
    lower_free_edge = zeros(3, length(times))
    upper_free_edge = zeros(3, length(times))

    for (k, t) in pairs(times)
        rests = _end_point_x_axis_rest_quats(build.model, drive_joint_indices, t)
        commanded_angle[k] = rad2deg(_smoothstep(t / ROTATION_DURATION_S) * ROTATION_TARGET_RAD)
        actual_angle[k] = _average_panel_rotation_deg(states[k])
        loads = SpaceAGORA.compliant_joint_loads(build.model, states[k]; joint_rest_quaternions=rests, joint_actuators=actuators)
        compliance_torque_x[k] = sum(loads[j].compliance_torque_child_body[1] for j in drive_joint_indices)
        actuator_torque_x[k] = sum(loads[j].actuator_torque_child_body[1] for j in drive_joint_indices)
        limit_sum = sum(a.torque_limit_n_m[1] for a in actuators)
        saturation[k] = limit_sum > 0.0 ? abs(actuator_torque_x[k]) / limit_sum : 0.0
        lower_edge, upper_edge = _free_end_edge_span(states[k])
        lower_free_edge[:, k] .= lower_edge
        upper_free_edge[:, k] .= upper_edge

        if k < length(times)
            states[k + 1] = SpaceAGORA.step_compliant_multibody_rk4(
                build.model,
                states[k],
                t,
                times[k + 1] - t;
                joint_rest_quaternions=rests,
                joint_actuators=actuators,
            )
        end
    end

    return (
        times=times,
        states=states,
        commanded_angle_deg=commanded_angle,
        actual_angle_deg=actual_angle,
        compliance_torque_x_nm=compliance_torque_x,
        actuator_torque_x_nm=actuator_torque_x,
        actuator_saturation=saturation,
        lower_free_edge_m=lower_free_edge,
        upper_free_edge_m=upper_free_edge,
    )
end

function _panel_mesh_geometry(x, panel_idx)
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    is = Int[]
    js = Int[]
    ks = Int[]
    wire_x = Float64[]
    wire_y = Float64[]
    wire_z = Float64[]

    function add_wire(a, b)
        push!(wire_x, a[1], b[1], NaN)
        push!(wire_y, a[2], b[2], NaN)
        push!(wire_z, a[3], b[3], NaN)
    end

    for tile_y in 1:PANEL_TILES_Y, tile_x in 1:PANEL_TILES_X
        corners = _tile_corners(x, _node_index(panel_idx, tile_x, tile_y))
        base = length(xs)
        append!(xs, (p[1] for p in corners))
        append!(ys, (p[2] for p in corners))
        append!(zs, (p[3] for p in corners))
        push!(is, base, base)
        push!(js, base + 1, base + 2)
        push!(ks, base + 2, base + 3)
        add_wire(corners[1], corners[2])
        add_wire(corners[2], corners[3])
        add_wire(corners[3], corners[4])
        add_wire(corners[4], corners[1])
    end
    return (x=xs, y=ys, z=zs, i=is, j=js, k=ks, wire_x=wire_x, wire_y=wire_y, wire_z=wire_z)
end

function _panel_mesh_trace(x, panel_idx, name, color; opacity=0.55, showlegend=false)
    geom = _panel_mesh_geometry(x, panel_idx)
    return PlotlyJS.mesh3d(
        x=geom.x,
        y=geom.y,
        z=geom.z,
        i=geom.i,
        j=geom.j,
        k=geom.k,
        name=name,
        color=color,
        opacity=opacity,
        showlegend=showlegend,
    )
end

function _panel_wire_trace(x, panel_idx, name; opacity=0.75, showlegend=false)
    geom = _panel_mesh_geometry(x, panel_idx)
    return PlotlyJS.scatter3d(
        x=geom.wire_x,
        y=geom.wire_y,
        z=geom.wire_z,
        mode="lines",
        name=name,
        line=PlotlyJS.attr(width=1.6, color=@sprintf("rgba(18,28,36,%.3f)", opacity)),
        showlegend=showlegend,
    )
end

function _deployment_3d_plot(sim)
    traces = PlotlyJS.GenericTrace[]
    idxs = unique(round.(Int, range(1, length(sim.times), length=18)))
    wire_idxs = Set(unique(round.(Int, range(1, length(sim.times), length=10))))
    palette = ["rgb(30,88,145)", "rgb(42,135,112)", "rgb(232,177,61)", "rgb(202,83,73)"]
    for (sidx, idx) in pairs(idxs)
        alpha = 0.10 + 0.50 * sidx / length(idxs)
        for panel_idx in 1:PANEL_COUNT
            push!(traces, _panel_mesh_trace(
                sim.states[idx],
                panel_idx,
                idx == first(idxs) && panel_idx == 1 ? "tile mesh snapshots" : "",
                palette[panel_idx];
                opacity=alpha,
                showlegend=idx == first(idxs) && panel_idx == 1,
            ))
            if idx in wire_idxs
                push!(traces, _panel_wire_trace(
                    sim.states[idx],
                    panel_idx,
                    idx == first(idxs) && panel_idx == 1 ? "tile boundaries" : "";
                    opacity=0.25 + 0.45 * sidx / length(idxs),
                    showlegend=idx == first(idxs) && panel_idx == 1,
                ))
            end
        end
    end
    push!(traces, PlotlyJS.scatter3d(
        x=sim.tip_m[1, :],
        y=sim.tip_m[2, :],
        z=sim.tip_m[3, :],
        mode="lines",
        name="outer tip path",
        line=PlotlyJS.attr(width=6, color="rgb(30,30,30)"),
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=[0.0],
        y=[0.0],
        z=[0.0],
        mode="markers",
        name="root hinge",
        marker=PlotlyJS.attr(size=6, color="rgb(120,35,35)", symbol="diamond"),
    ))

    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="SpaceAGORA Cloth Solar Panel Deployment Mesh",
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis=PlotlyJS.attr(title="x (m)"),
                yaxis=PlotlyJS.attr(title="y (m)"),
                zaxis=PlotlyJS.attr(title="z (m)"),
            ),
            legend=PlotlyJS.attr(x=0.02, y=0.98),
            margin=PlotlyJS.attr(l=0, r=0, t=55, b=0),
        ),
    )
end

function _progress_plot(sim)
    traces = PlotlyJS.GenericTrace[]
    for i in 1:(PANEL_COUNT - 1)
        push!(traces, PlotlyJS.scatter(
            x=sim.times,
            y=sim.commanded_angles_deg[i, :],
            mode="lines",
            name="hinge $i command deg",
            xaxis="x1",
            yaxis="y1",
            line=PlotlyJS.attr(width=2, dash="dash"),
        ))
        push!(traces, PlotlyJS.scatter(
            x=sim.times,
            y=sim.actual_angles_deg[i, :],
            mode="lines",
            name="hinge $i actual deg",
            xaxis="x1",
            yaxis="y1",
            line=PlotlyJS.attr(width=3),
        ))
        push!(traces, PlotlyJS.scatter(
            x=sim.times,
            y=sim.compliance_torque_y_nm[i, :],
            mode="lines",
            name="hinge $i compliance Nm",
            xaxis="x2",
            yaxis="y2",
        ))
        push!(traces, PlotlyJS.scatter(
            x=sim.times,
            y=sim.actuator_torque_y_nm[i, :],
            mode="lines",
            name="hinge $i actuator Nm",
            xaxis="x3",
            yaxis="y3",
        ))
        push!(traces, PlotlyJS.scatter(
            x=sim.times,
            y=sim.actuator_saturation[i, :],
            mode="lines",
            name="hinge $i saturation",
            xaxis="x4",
            yaxis="y4",
        ))
    end
    push!(traces, PlotlyJS.scatter(
        x=sim.times,
        y=sim.tip_m[1, :],
        mode="lines",
        name="tip x m",
        xaxis="x5",
        yaxis="y5",
    ))
    push!(traces, PlotlyJS.scatter(
        x=sim.times,
        y=sim.tip_m[3, :],
        mode="lines",
        name="tip z m",
        xaxis="x5",
        yaxis="y5",
    ))

    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="Solar Panel Mesh Deployment Progress",
            grid=PlotlyJS.attr(rows=3, columns=2, pattern="independent"),
            xaxis=PlotlyJS.attr(title="time (s)"),
            yaxis=PlotlyJS.attr(title="hinge angle (deg)"),
            xaxis2=PlotlyJS.attr(title="time (s)"),
            yaxis2=PlotlyJS.attr(title="aggregate compliance torque about hinge (Nm)"),
            xaxis3=PlotlyJS.attr(title="time (s)"),
            yaxis3=PlotlyJS.attr(title="aggregate actuator torque about hinge (Nm)"),
            xaxis4=PlotlyJS.attr(title="time (s)"),
            yaxis4=PlotlyJS.attr(title="actuator saturation"),
            xaxis5=PlotlyJS.attr(title="time (s)"),
            yaxis5=PlotlyJS.attr(title="outer tip position (m)"),
            legend=PlotlyJS.attr(orientation="h", y=-0.22),
            margin=PlotlyJS.attr(l=75, r=35, t=70, b=95),
        ),
    )
end

function _rotation_3d_plot(sim)
    traces = PlotlyJS.GenericTrace[]
    idxs = unique(round.(Int, range(1, length(sim.times), length=16)))
    wire_idxs = Set(unique(round.(Int, range(1, length(sim.times), length=8))))
    palette = ["rgb(30,88,145)", "rgb(42,135,112)", "rgb(232,177,61)", "rgb(202,83,73)"]
    for (sidx, idx) in pairs(idxs)
        alpha = 0.10 + 0.52 * sidx / length(idxs)
        for panel_idx in 1:PANEL_COUNT
            push!(traces, _panel_mesh_trace(
                sim.states[idx],
                panel_idx,
                idx == first(idxs) && panel_idx == 1 ? "90 deg rotation snapshots" : "",
                palette[panel_idx];
                opacity=alpha,
                showlegend=idx == first(idxs) && panel_idx == 1,
            ))
            if idx in wire_idxs
                push!(traces, _panel_wire_trace(
                    sim.states[idx],
                    panel_idx,
                    idx == first(idxs) && panel_idx == 1 ? "tile boundaries" : "";
                    opacity=0.25 + 0.45 * sidx / length(idxs),
                    showlegend=idx == first(idxs) && panel_idx == 1,
                ))
            end
        end
    end
    push!(traces, PlotlyJS.scatter3d(
        x=sim.lower_free_edge_m[1, :],
        y=sim.lower_free_edge_m[2, :],
        z=sim.lower_free_edge_m[3, :],
        mode="lines",
        name="lower free-edge path",
        line=PlotlyJS.attr(width=5, color="rgb(35,35,35)"),
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=sim.upper_free_edge_m[1, :],
        y=sim.upper_free_edge_m[2, :],
        z=sim.upper_free_edge_m[3, :],
        mode="lines",
        name="upper free-edge path",
        line=PlotlyJS.attr(width=5, color="rgb(95,45,45)"),
    ))
    total_panel_length = PANEL_COUNT * PANEL_LENGTH_M
    push!(traces, PlotlyJS.scatter3d(
        x=[-0.5 * total_panel_length, 0.5 * total_panel_length],
        y=[0.0, 0.0],
        z=[0.0, 0.0],
        mode="lines",
        name="long x rotation axis",
        line=PlotlyJS.attr(width=7, color="rgb(135,35,95)"),
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=[-0.5 * total_panel_length],
        y=[0.0],
        z=[0.0],
        mode="markers",
        name="end-point drive",
        marker=PlotlyJS.attr(size=7, color="rgb(120,35,35)", symbol="diamond"),
    ))

    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="Already-Deployed Solar Panel 90 Degree End-Point X-Axis Rotation",
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis=PlotlyJS.attr(title="x (m)"),
                yaxis=PlotlyJS.attr(title="y (m)"),
                zaxis=PlotlyJS.attr(title="z (m)"),
            ),
            legend=PlotlyJS.attr(x=0.02, y=0.98),
            margin=PlotlyJS.attr(l=0, r=0, t=55, b=0),
        ),
    )
end

function _rotation_progress_plot(sim)
    traces = PlotlyJS.GenericTrace[
        PlotlyJS.scatter(
            x=sim.times,
            y=sim.commanded_angle_deg,
            mode="lines",
            name="command deg",
            xaxis="x1",
            yaxis="y1",
            line=PlotlyJS.attr(width=2, dash="dash"),
        ),
        PlotlyJS.scatter(
            x=sim.times,
            y=sim.actual_angle_deg,
            mode="lines",
            name="mesh average deg",
            xaxis="x1",
            yaxis="y1",
            line=PlotlyJS.attr(width=3),
        ),
        PlotlyJS.scatter(
            x=sim.times,
            y=sim.compliance_torque_x_nm,
            mode="lines",
            name="x-axis compliance torque Nm",
            xaxis="x2",
            yaxis="y2",
        ),
        PlotlyJS.scatter(
            x=sim.times,
            y=sim.actuator_torque_x_nm,
            mode="lines",
            name="x-axis actuator torque Nm",
            xaxis="x3",
            yaxis="y3",
        ),
        PlotlyJS.scatter(
            x=sim.times,
            y=sim.actuator_saturation,
            mode="lines",
            name="x-axis actuator saturation",
            xaxis="x4",
            yaxis="y4",
        ),
        PlotlyJS.scatter(
            x=sim.times,
            y=sim.lower_free_edge_m[3, :],
            mode="lines",
            name="lower free-edge z m",
            xaxis="x5",
            yaxis="y5",
        ),
        PlotlyJS.scatter(
            x=sim.times,
            y=sim.upper_free_edge_m[3, :],
            mode="lines",
            name="upper free-edge z m",
            xaxis="x5",
            yaxis="y5",
        ),
    ]

    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="End-Point X-Axis Rotation Progress",
            grid=PlotlyJS.attr(rows=3, columns=2, pattern="independent"),
            xaxis=PlotlyJS.attr(title="time (s)"),
            yaxis=PlotlyJS.attr(title="rotation angle (deg)"),
            xaxis2=PlotlyJS.attr(title="time (s)"),
            yaxis2=PlotlyJS.attr(title="compliance torque about x axis (Nm)"),
            xaxis3=PlotlyJS.attr(title="time (s)"),
            yaxis3=PlotlyJS.attr(title="actuator torque about x axis (Nm)"),
            xaxis4=PlotlyJS.attr(title="time (s)"),
            yaxis4=PlotlyJS.attr(title="actuator saturation"),
            xaxis5=PlotlyJS.attr(title="time (s)"),
            yaxis5=PlotlyJS.attr(title="free-edge z position (m)"),
            legend=PlotlyJS.attr(orientation="h", y=-0.22),
            margin=PlotlyJS.attr(l=75, r=35, t=70, b=95),
        ),
    )
end

function _save_html(path, plot)
    mkpath(dirname(path))
    open(path, "w") do io
        show(io, MIME("text/html"), plot; include_plotlyjs="cdn", full_html=true)
    end
    return abspath(path)
end

function run_solar_panel_cloth_deployment_demo()
    mkpath(OUT_DIR)
    build, initial_relative_angles, hinge_joint_groups, _, actuators = _build_folded_panel_model()
    sim = _simulate_deployment(build, initial_relative_angles, hinge_joint_groups, actuators)
    deployment_html = _save_html(joinpath(OUT_DIR, "solar_panel_deployment_3d.html"), _deployment_3d_plot(sim))
    progress_html = _save_html(joinpath(OUT_DIR, "solar_panel_deployment_progress.html"), _progress_plot(sim))

    rotation_build, drive_joint_indices, rotation_actuators = _build_end_point_x_axis_rotation_model()
    rotation_sim = _simulate_end_point_x_axis_rotation(rotation_build, drive_joint_indices, rotation_actuators)
    rotation_html = _save_html(joinpath(OUT_DIR, "solar_panel_endpoint_x_axis_rotation_3d.html"), _rotation_3d_plot(rotation_sim))
    rotation_progress_html = _save_html(joinpath(OUT_DIR, "solar_panel_endpoint_x_axis_rotation_progress.html"), _rotation_progress_plot(rotation_sim))

    final_angles = sim.actual_angles_deg[:, end]
    final_rotation = rotation_sim.actual_angle_deg[end]
    println("Solar panel Cloth mesh deployment demo complete.")
    println("  tile bodies             = ", length(build.model.bodies))
    println("  compliant joints        = ", length(build.model.joints))
    println("  final hinge angles deg  = ", join((@sprintf("%.2f", a) for a in final_angles), ", "))
    println("  3D deployment HTML      = ", deployment_html)
    println("  progress HTML           = ", progress_html)
    println("  endpoint x-axis final deg = ", @sprintf("%.2f", final_rotation))
    println("  rotation 3D HTML        = ", rotation_html)
    println("  rotation progress HTML  = ", rotation_progress_html)
    return (
        build=build,
        simulation=sim,
        deployment_html=deployment_html,
        progress_html=progress_html,
        rotation_build=rotation_build,
        rotation_simulation=rotation_sim,
        rotation_html=rotation_html,
        rotation_progress_html=rotation_progress_html,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_solar_panel_cloth_deployment_demo()
end
