# Helper: save timeseries CSV
function save_timeseries_csv(sol, p, helper_oe, target_orbit;
                             csv_dir=normpath(joinpath(@__DIR__, "..", "output", "CSV")),
                             include_B=false, include_Pin=false, include_Range=false,
                             include_J2=false, include_toy=false)
    N = p[:N]
    header = ["t"]
    for i in 1:N
        append!(header, ["x$i","y$i","z$i","vx$i","vy$i","vz$i"])
    end
    rows = Vector{Float64}[]
    for (k, t) in enumerate(sol.t)
        u = sol.u[k]
        row = Float64[t]
        for i in 1:N
            push!(row, u[idx(i,1)])
            push!(row, u[idx(i,2)])
            push!(row, u[idx(i,3)])
            push!(row, u[idx(i,4)])
            push!(row, u[idx(i,5)])
            push!(row, u[idx(i,6)])
        end
        push!(rows, row)
    end
    mkpath(csv_dir)
    sim_time_s = sol.t[end] - sol.t[1]
    helper_alt_m = helper_oe[1].a_m - R_EARTH
    target_alt_m = target_orbit.a_m - R_EARTH
    helper_i_deg = helper_oe[1].i_deg
    target_i_deg = target_orbit.i_deg
    fname = @sprintf("timeseries_N%d_T%.0fs_h%.0fkm_t%.0fkm_ih%.1fdeg_it%.1fdeg", N, sim_time_s, helper_alt_m/1000, target_alt_m/1000, helper_i_deg, target_i_deg)
    if include_B || include_Pin
        first_cav = first(values(p[:cavity]))
        if include_B
            fname *= @sprintf("_B%.4g", first_cav[:B])
        end
        if include_Pin
            fname *= @sprintf("_Pin%.4g", first_cav[:Pin])
        end
    end
    if include_Range
        fname *= @sprintf("_rmin%.0fm_rmax%.4g", p[:min_range], p[:max_range])
    end
    if include_J2
        fname *= get(p, :use_J2, false) ? "_J2T" : "_J2F"
    end
    if include_toy
        fname *= "_toyT"
    end
    fname *= ".csv"
    csv_path = joinpath(csv_dir, fname)
    open(csv_path, "w") do io
        # Save metadata
        println(io, "# satellites=", N)
        println(io, "# sim_time_s=", sim_time_s)
        println(io, "# helper_altitude_m=", helper_alt_m)
        println(io, "# target_altitude_m=", target_alt_m)
        # Save additional parameters from p
        println(io, "# mu=", p[:mu])
        println(io, "# c=", p[:c])
        println(io, "# use_los=", p[:use_los])
        println(io, "# R_atm=", p[:R_atm])
        println(io, "# atm_clearance=", p[:atm_clearance])
        println(io, "# min_range=", p[:min_range])
        println(io, "# max_range=", p[:max_range])
        println(io, "# use_J2=", get(p, :use_J2, false))
        println(io, "# useDrag=", get(p, :useDrag, false))
        # Save masses as comma-separated
        println(io, "# masses=", join(p[:masses], ","))
        println(io, "# cavity=", p[:cavity])
        println(io, "# Pmatrix=", p[:Pmatrix])
        println(io, "# helper_ids=", join(p[:helper_ids], ","))
        println(io, "# target_ids=", join(p[:target_ids], ","))
        writedlm(io, permutedims(header), ',')
        for row in rows
            writedlm(io, row', ',')
        end
    end
    println("Saved CSV to: ", csv_path)
    return csv_path
end

# Helper: load timeseries CSV and reconstruct solution
function load_timeseries_csv(csv_path)
    metadata = Dict{String, Any}()
    
    # Read metadata and data
    lines = readlines(csv_path)
    data_start_idx = 1
    
    # Parse metadata from comment lines
    for (i, line) in enumerate(lines)
        if startswith(line, "# satellites=")
            metadata["N"] = parse(Int, split(line, "=")[2])
        elseif startswith(line, "# sim_time_s=")
            metadata["sim_time_s"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# helper_altitude_m=")
            metadata["helper_altitude_m"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# target_altitude_m=")
            metadata["target_altitude_m"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# mu=")
            metadata["mu"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# c=")
            metadata["c"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# use_los=")
            metadata["use_los"] = parse(Bool, split(line, "=")[2])
        elseif startswith(line, "# R_atm=")
            metadata["R_atm"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# atm_clearance=")
            metadata["atm_clearance"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# min_range=")
            metadata["min_range"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# max_range=")
            metadata["max_range"] = parse(Float64, split(line, "=")[2])
        elseif startswith(line, "# use_J2=")
            metadata["use_J2"] = parse(Bool, split(line, "=")[2])
        elseif startswith(line, "# useDrag=")
            metadata["useDrag"] = parse(Bool, split(line, "=")[2])
        elseif startswith(line, "# masses=")
            mass_str = split(line, "=")[2]
            metadata["masses"] = parse.(Float64, split(mass_str, ","))
        elseif !startswith(line, "#")
            data_start_idx = i
            break
        end
    end
    
    # Read CSV data (skip metadata lines and header)
    data = readdlm(csv_path, ',', Float64; skipstart=data_start_idx)
    
    # Extract time and state vectors
    N = metadata["N"]
    t = data[:, 1]
    u = Vector{Vector{Float64}}()
    
    for row_idx in 1:size(data, 1)
        state = zeros(6 * N)
        for i in 1:N
            state[idx(i,1)] = data[row_idx, 1 + (i-1)*6 + 1]  # x
            state[idx(i,2)] = data[row_idx, 1 + (i-1)*6 + 2]  # y
            state[idx(i,3)] = data[row_idx, 1 + (i-1)*6 + 3]  # z
            state[idx(i,4)] = data[row_idx, 1 + (i-1)*6 + 4]  # vx
            state[idx(i,5)] = data[row_idx, 1 + (i-1)*6 + 5]  # vy
            state[idx(i,6)] = data[row_idx, 1 + (i-1)*6 + 6]  # vz
        end
        push!(u, state)
    end
    
    # Create a simple named tuple mimicking ODE solution structure
    sol = (t=t, u=u)
    
    cav = Dict{Tuple{Int, Int}, Dict{Symbol, Any}}() # Initialize the cav dictionary
    for i in 1:(N-1) # Construct cav for helper satellites (1 to helper_num) linking to the target satellite (index N)
        cav[(i, N)] = Dict(:B => 100.0, :Pin => 1e4)
    end

    # Reconstruct parameter dictionary p
    p = Dict{Symbol, Any}(
        :N => N,
        :mu => get(metadata, "mu", MU),
        :c => get(metadata, "c", C),
        :masses => get(metadata, "masses", fill(227.0, N)),
        :use_los => get(metadata, "use_los", true),
        :R_atm => get(metadata, "R_atm", R_EARTH + 100_000.0),
        :atm_clearance => get(metadata, "atm_clearance", 5_000.0),
        :min_range => get(metadata, "min_range", 0.0),
        :max_range => get(metadata, "max_range", Inf),
        :use_J2 => get(metadata, "use_J2", false),
        :useDrag => get(metadata, "useDrag", false),
        :cavity => cav,
        :Pmatrix => zeros(Float64, N, N),
        :helper_ids => collect(1:N-1),
        :target_ids => [N]
    )
    
    println("Loaded CSV from: ", csv_path)
    println("  Satellites: ", N)
    println("  Time points: ", length(t))

    return sol, metadata, p
end