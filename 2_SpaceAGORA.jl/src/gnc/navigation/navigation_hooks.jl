module NavigationHooks
    using LinearAlgebra
    using StaticArrays

    export calcNavigationEffect!
    export RPOStationGeometry, RPOCubeSatGeometry, RPOReferenceGeometry
    export nearest_station_distance_sq, nearest_station_point
    export rpo_clearance_distance_to_station, rpo_clearance_to_station, rpo_path_clearance_stats
    export rpo_surface_normal_from_pointcloud, rpo_goal_standoff_point

    """
        calcNavigationEffect!(model, u, p, t, sat_idx)

    Hook point for navigation/sensor-estimator models executed by typed periodic callbacks.
    Extend this method for your navigation effector type.
    """
    function calcNavigationEffect!(model, u, p, t::Float64, sat_idx::Int)
        throw(MethodError(calcNavigationEffect!, (model, u, p, t, sat_idx)))
    end

    include(joinpath(@__DIR__, "rpo_nav", "reference_geometry", "station_geometry.jl"))
    include(joinpath(@__DIR__, "rpo_nav", "reference_geometry", "cubesat_geometry.jl"))
    include(joinpath(@__DIR__, "rpo_nav", "reference_geometry", "rpo_reference_geometry.jl"))
    include(joinpath(@__DIR__, "rpo_nav", "distances", "mesh_distance.jl"))
    include(joinpath(@__DIR__, "rpo_nav", "distances", "clearance.jl"))
    include(joinpath(@__DIR__, "rpo_nav", "distances", "surface_frames.jl"))
    include(joinpath(@__DIR__, "rpo_nav", "distances", "rpo_distance_queries.jl"))
end
