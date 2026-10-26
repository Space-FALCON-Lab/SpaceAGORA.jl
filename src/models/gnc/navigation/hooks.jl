module NavigationHooks
    using LinearAlgebra
    using StaticArrays

    export update_navigation!
    export RPOStationGeometry, RPOCubeSatGeometry, RPOReferenceGeometry
    export nearest_station_distance_sq, nearest_station_point
    export rpo_clearance_distance_to_station, rpo_clearance_to_station, rpo_path_clearance_stats
    export rpo_surface_normal_from_pointcloud, rpo_goal_standoff_point

    """
        update_navigation!(model, u, p, t, spacecraft_index) -> nothing

    The navigation tick, which every effector in a `NavigationModel` (a sensor
    or a state estimator, usually an `AbstractNavigationModel`) implements for
    its own type. The effector takes its measurements or updates its estimate
    and keeps them in its own state for guidance and control effectors to read;
    the engine ignores the return value. SpaceAGORA has no built-in navigation
    effector, so this function has no methods until one is added.

    - `model`: the navigation effector.
    - `u`: the full integrated state. `u.sc[spacecraft_index]` is the spacecraft, with
      `pos` (m) and `vel` (m/s) relative to the planet's centre in the inertial
      frame and `mass` (kg); with attitude propagation also `q = (x, y, z, w)`,
      scalar last, with `v_body = rot(q) * v_inertial` (body means bus axes), and
      `ω`, the bus angular velocity in rad/s. These are the true values; any sensor
      noise is the effector's own.
      Sensors mounted with `add_sensor!` are sampled before navigation ticks at
      the same time; read them with `sensor_measurement(p.args.dynamics_model.spacecraft[spacecraft_index], k)`.
    - `p`: the run's parameters; `p.args` is the `SimulationConfiguration`.
    - `t`: elapsed time since the run's start epoch, in seconds (`Float64`).
    - `spacecraft_index`: the spacecraft's position in the run (`Int`, from 1), not its ID.

    Effector `k` runs every `navigation_rates[k]` seconds of simulated time, once
    for each spacecraft in index order, starting one period after the solve
    starts. When ticks coincide, navigation effectors run before guidance and
    control effectors.

    Before integration starts, `run_simulation` checks that every navigation
    effector has a method of this function for these arguments and stops with
    an `ArgumentError` naming the effector otherwise.
    """
    function update_navigation! end

    # Reference/collision geometry is planning data, not an estimator. Retain
    # these historical NavigationHooks identities for existing clients.
    include(joinpath(@__DIR__, "..", "planning", "rpo", "geometry", "station_geometry.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "geometry", "cubesat_geometry.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "geometry", "rpo_reference_geometry.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "geometry", "mesh_distance.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "geometry", "clearance.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "geometry", "surface_frames.jl"))
    include(joinpath(@__DIR__, "..", "planning", "rpo", "geometry", "rpo_distance_queries.jl"))
end
