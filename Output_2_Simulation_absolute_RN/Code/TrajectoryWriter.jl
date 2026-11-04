mutable struct TrajectoryWriter{W}
    writer::W
    columns::Vector{Symbol}
    buffer::Matrix{Float64}
    target::Int
    used::Int
    rows::Int
    last_time::Float64
    closed::Bool
end

function TrajectoryWriter(path, spacecraft_count, run_id; chunk_rows=512)
    columns = trajectory_columns(spacecraft_count)
    metadata = Dict("run_id" => run_id, "schema_version" => string(SCHEMA_VERSION),
        "target_id" => string(spacecraft_count), "frame" => "J2000 ECI",
        "units" => "time:s; pos:m; vel,laser_dv,RTN:m/s; work:J",
        "rtn_definition" => "integral of applied laser acceleration in instantaneous target RTN",
        "rtn_branch_definition" => RTN_BRANCH_DEFINITION)
    writer = open(Arrow.Writer, path; file=true, compress=:zstd, ntasks=0, metadata)
    TrajectoryWriter(writer, columns, Matrix{Float64}(undef, chunk_rows, length(columns)),
        spacecraft_count, 0, 0, -Inf, false)
end

function flush_trajectory!(output::TrajectoryWriter)
    output.used == 0 && return
    # Arrow's writer may still own a batch while the next samples arrive.
    table = DataFrame([name => copy(view(output.buffer, 1:output.used, i))
        for (i, name) in enumerate(output.columns)]; copycols=false)
    Arrow.write(output.writer, table)
    output.used = 0
    nothing
end

function save_sample!(output::TrajectoryWriter, time, state)
    time > output.last_time || error("Trajectory times must be strictly increasing")
    all(isfinite, state) || error("Nonfinite simulation state at t=$time")
    output.used += 1
    row, column = output.used, 1
    output.buffer[row, column] = time
    for sc in state.sc
        for field in (:pos, :vel, :laser_dv), axis in 1:3
            column += 1
            output.buffer[row, column] = getproperty(sc, field)[axis]
        end
        column += 1
        output.buffer[row, column] = sc.laser_work
    end
    target = state.sc[output.target]
    for value in (target.laser_rtn..., target.laser_dv_T_positive, target.laser_dv_T_negative,
        target.laser_dv_R_positive, target.laser_dv_R_negative,
        target.laser_dv_N_positive, target.laser_dv_N_negative)
        column += 1
        output.buffer[row, column] = value
    end
    column == length(output.columns) || error("Trajectory column count mismatch")
    output.rows += 1
    output.last_time = time
    output.used == size(output.buffer, 1) && flush_trajectory!(output)
    nothing
end

function close_trajectory!(output::TrajectoryWriter)
    output.closed && return
    try
        flush_trajectory!(output)
    finally
        output.closed = true
        close(output.writer)
    end
    nothing
end
