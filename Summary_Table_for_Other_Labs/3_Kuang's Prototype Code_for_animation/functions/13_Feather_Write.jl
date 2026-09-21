import Arrow
include(joinpath(@__DIR__, "..", "..", "2_SpaceAGORA.jl", "ORACLE", "functions", "FeatherRecording.jl"))
using .FeatherRecording

function make_prototype_recorder(initial_state, params)
    order = vcat(params[:target_ids], params[:helper_ids])
    inverse_order = invperm(order)
    state_adapter(state) = reshape(state[1:6*params[:N]], 6, params[:N])[:, order]
    function active_pairs(states, time)
        forces, _ = laser_forces(vec(states[:, inverse_order]), params)
        return unique([minmax(inverse_order[source], inverse_order[target])
            for ((source, target), force) in forces if source != target && any(!iszero, force)])
    end
    function gates(states, time)
        values = Float64[]
        for first in 1:params[:N]-1, second in first+1:params[:N]
            position_a = SVector{3,Float64}(states[1:3, first])
            position_b = SVector{3,Float64}(states[1:3, second])
            if params[:min_range] > 0
                push!(values, norm(position_b-position_a) - params[:min_range])
            end
            if params[:use_los]
                metrics = los_metrics(position_a, position_b; R_atm=params[:R_atm])
                push!(values, metrics.rclosest_norm - max(params[:R_atm], params[:R_atm]+params[:atm_clearance]))
            end
        end
        return values
    end
    return IntervalRecorder(state_adapter, active_pairs, initial_state, params[:max_range]; gates)
end

"""
    save_timeseries_feather(sol, p; feather_dir)

Write `trajectory.feather` using SpaceAGORA's native Arrow IPC column
names and SI units. Targets precede helpers, so the single ORACLE target is sc1.
Samples follow the runner's 10-second output grid, including the final endpoint.
The native non-attitude schema is preserved; fields not recorded by the
prototype are nullable Float64 columns containing missing, not fabricated zeros.
Laser delta-v uses the prototype's trapezoidal RTN diagnostic on these samples.
"""
function save_timeseries_feather(sol, p;
                                 feather_dir=normpath(joinpath(@__DIR__, "..", "output", "feather")),
                                 scenario=basename(feather_dir))
    length(p[:target_ids]) == 1 ||
        throw(ArgumentError("ORACLE Feather output requires exactly one target"))
    satellite_order = vcat(p[:target_ids], p[:helper_ids])
    sort(satellite_order) == collect(1:p[:N]) ||
        throw(ArgumentError("target_ids and helper_ids must partition all satellites"))
    columns = Pair{Symbol,AbstractVector}[:time => Float64.(sol.t)]
    for (field, offset) in (("pos", 0), ("vel", 3))
        for (output_id, prototype_id) in enumerate(satellite_order)
            for component in 1:3
                name = Symbol("sc$(output_id)_$(field)_$(component)")
                values = Float64[state[idx(prototype_id, offset + component)] for state in sol.u]
                push!(columns, name => values)
            end
        end
    end
    for field in ("altitude", "latitude_deg", "longitude_deg", "mass",
                  "wind", "drag", "lift", "cross", "periapsis_altitude",
                  "heat_rate", "heat_load")
        for (output_id, prototype_id) in enumerate(satellite_order)
            if field == "mass"
                push!(columns, Symbol("sc$(output_id)_mass") =>
                    fill(Float64(p[:masses][prototype_id]), length(sol.t)))
            elseif field == "periapsis_altitude"
                altitudes = map(sol.u) do state
                    position = SVector{3,Float64}(state[idx(prototype_id, 1):idx(prototype_id, 3)])
                    velocity = SVector{3,Float64}(state[idx(prototype_id, 4):idx(prototype_id, 6)])
                    elements = rv2coe(position, velocity, p[:mu])
                    elements.a * (1 - elements.e) - R_EARTH
                end
                push!(columns, Symbol("sc$(output_id)_periapsis_altitude") => altitudes)
            else
                suffixes = field in ("wind", "drag", "lift", "cross") ? ("_1", "_2", "_3") : ("",)
                for suffix in suffixes
                    push!(columns, Symbol("sc$(output_id)_$(field)$(suffix)") =>
                        fill!(Vector{Union{Missing,Float64}}(undef, length(sol.t)), missing))
                end
            end
        end
    end
    target_id = only(p[:target_ids])
    _, delta_v = delta_v_RTN_time_series(sol, p)
    for (component, axis) in enumerate(("r", "t", "n"))
        push!(columns, Symbol("dv_$(axis)_accumulated") => delta_v[target_id][component, :])
    end
    active_helpers = map(sol.u) do state
        forces, _ = laser_forces(state, p)
        helper_index = findfirst(helper -> haskey(forces, (helper, target_id)), p[:helper_ids])
        helper_index === nothing ? 0 : helper_index + 1
    end
    push!(columns, :laser_active_helper => active_helpers)
    mkpath(feather_dir)
    feather_path = joinpath(feather_dir, "trajectory.feather")
    Arrow.write(feather_path, (; columns...); metadata=Dict(
        "source" => "Kuang prototype",
        "scenario" => scenario,
        "output_interval_s" => string(get(p, :output_interval_s, 10.0)),
        "maximum_range_m" => string(p[:max_range]),
        "prototype_satellite_order" => join(satellite_order, ","),
        "unavailable_fields" => "altitude,latitude_deg,longitude_deg,wind,drag,lift,cross,heat_rate,heat_load",
        "delta_v_method" => "trapezoidal RTN integration on solver samples",
    ))
    if haskey(p, :feather_recorder)
        write_intervals(p[:feather_recorder], feather_dir; source="Kuang prototype", scenario,
            laser_timing="range/LOS/min-range roots and accepted-step force selection; GVE selection checked at accepted steps")
    end
    println("Saved Feather to: ", feather_path)
    return feather_path
end