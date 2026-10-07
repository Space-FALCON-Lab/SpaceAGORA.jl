"""
Optional MuJoCo proximity scenes for SpaceAGORA: rigid-body contact-capable multibody scenes (MJCF) that
fly in a chief-centered inertial frame while SpaceAGORA's own gravity models drive the chief orbit and the
per-body differential gravity.

Load with `using SpaceAGORAMuJoCo`. The native library is MuJoCo 3.11.0 from a lazy artifact (Linux x86_64);
see `README.md` for the dependency review. Stage 1 is the standalone scene runner: [`ProximityScene`](@ref),
[`scene_reset!`](@ref), [`scene_step!`](@ref), [`scene_state`](@ref), [`scene_set_state!`](@ref). Stage 2 runs a
scene inside `run_simulation` through [`ProximitySceneDynamics`](@ref).
"""
module SpaceAGORAMuJoCo

include("binding.jl")
include("proximity_scene.jl")
include("engine_coupling.jl")

export ProximityScene, SceneBodyState, SceneState, body_state_from_initial_condition
export sa_to_mujoco_quaternion, mujoco_to_sa_quaternion, ProximitySceneDynamics, ProximitySceneRuntime, scene_body_pose_save_field
export scene_reset!, scene_step!, scene_state, scene_set_state!, scene_time, scene_chief, scene_body_names,
    scene_body_state, scene_ctrl

end # module SpaceAGORAMuJoCo
