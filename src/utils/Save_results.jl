using .SimulationModel: ODEParams, IntermediateSolution, Solution

if !isdefined(@__MODULE__, :_legacy_get_save_results_runtime_state)
    @inline function _legacy_get_save_results_runtime_state(args=nothing; cnf=nothing, solution=nothing, model=nothing)
        module_runtime = isdefined(@__MODULE__, :config) ? getfield(@__MODULE__, :config) : nothing

        cnf_state = if cnf !== nothing
            cnf
        elseif args isa AbstractDict && haskey(args, :cnf)
            args[:cnf]
        elseif module_runtime !== nothing && hasproperty(module_runtime, :cnf)
            getproperty(module_runtime, :cnf)
        else
            nothing
        end

        solution_state = if solution !== nothing
            solution
        elseif args isa AbstractDict && haskey(args, :solution)
            args[:solution]
        elseif module_runtime !== nothing && hasproperty(module_runtime, :solution)
            getproperty(module_runtime, :solution)
        else
            nothing
        end

        model_state = if model !== nothing
            model
        elseif args isa AbstractDict && haskey(args, :model)
            args[:model]
        elseif module_runtime !== nothing && hasproperty(module_runtime, :model)
            getproperty(module_runtime, :model)
        else
            nothing
        end

        if cnf_state === nothing || solution_state === nothing || model_state === nothing
            throw(ArgumentError("Legacy runtime state missing; provide `cnf`, `solution`, and `model`."))
        end

        return (cnf=cnf_state, solution=solution_state, model=model_state)
    end
end

function save_results(time::Vector{Float64}, ratio::Float64, params::ODEParams)
    initial_time = 0
    solution = params.solution
    cnf = params.cnf
    runtime = _legacy_get_save_results_runtime_state(nothing; cnf=cnf, solution=solution, model=params.m)

    if length(solution.orientation.time) == 0
        cnf.prev_step_integrator = 0.0
        initial_time = cnf.prev_step_integrator
        cnf.initial_time_saved = 0.0
    else
        if cnf.counter_integrator == 1
            cnf.prev_step_integrator =  0.0
            intial_time = cnf.prev_step_integrator
            cnf.initial_time_saved = solution.orientation.time[end]
        end
    end

    n_variable_to_save = length(cnf.solution_intermediate[1]) - 1
    range_time = [item[1] for item in cnf.solution_intermediate]
    unique!(time) # Remove duplicate time entries
    results = Matrix{Float64}(zeros(n_variable_to_save, Int(ceil(length(time)/ratio))))


    t = zeros(Float64, Int(ceil(length(time)/ratio)))

    index_prev = 1

    for (i, true_time) in enumerate(time)
        # if isapprox(true_time, range_time[index_prev], atol = 1e-8) # Skip duplicate time entries
        #     continue
        # end
        if isapprox((i-1) % ratio, 0, atol = 0.1) || i == length(time) # Save every 'ratio' time steps or the last time step
            index = findfirst(x -> x == true_time, range_time[index_prev:end])
            t[Int(floor(i/ratio))] = i == length(time) ? true_time : range_time[index+index_prev] + initial_time
            # range_solution = SVector{n_variable_to_save + 1, Float64}([runtime.cnf.solution_intermediate[index+index_prev], 0])

            if length(t) == 1
                if index + index_prev <= length(runtime.cnf.solution_intermediate)
                    results[:, 1] .= runtime.cnf.solution_intermediate[index+index_prev][2:end]
                else
                    results[:, 1] .= runtime.cnf.solution_intermediate[end][2:end]
                end
            elseif true_time != time[end]
                results[:, Int(floor(i/ratio))] .= runtime.cnf.solution_intermediate[index+index_prev][2:end]
            else
                results[:, Int(floor(i/ratio))] .= runtime.cnf.solution_intermediate[end][2:end]
            end
            index_prev = index
        end

        # i += 1
    end

    time_0 = time[end]
    runtime.cnf.prev_step_integrator = time_0
    runtime.cnf.solution_intermediate = Vector{Number}[]
    runtime.cnf.solution_intermediate = Vector{Number}[]

    t = [i + runtime.cnf.initial_time_saved for i in t]

    ## SAVE RESULTS GLOBAL
    append!(runtime.solution.orientation.time, t)
    append!(runtime.solution.orientation.year, results[1,:])
    append!(runtime.solution.orientation.month, results[2,:])
    append!(runtime.solution.orientation.day, results[3,:])
    append!(runtime.solution.orientation.hour, results[4,:])
    append!(runtime.solution.orientation.minute, results[5,:])
    append!(runtime.solution.orientation.second, results[6,:])
    append!(runtime.solution.orientation.number_of_passage, results[7,:])
    append!(runtime.solution.orientation.pos_ii[1], results[8,:])
    append!(runtime.solution.orientation.pos_ii[2], results[9,:])
    append!(runtime.solution.orientation.pos_ii[3], results[10,:])
    append!(runtime.solution.orientation.vel_ii[1], results[11,:])
    append!(runtime.solution.orientation.vel_ii[2], results[12,:])
    append!(runtime.solution.orientation.vel_ii[3], results[13,:])
    append!(runtime.solution.orientation.pos_ii_mag, results[14,:])
    append!(runtime.solution.orientation.vel_ii_mag, results[15,:])
    append!(runtime.solution.orientation.quaternion[1], results[91,:])
    append!(runtime.solution.orientation.quaternion[2], results[92,:])
    append!(runtime.solution.orientation.quaternion[3], results[93,:])
    append!(runtime.solution.orientation.quaternion[4], results[94,:])
    append!(runtime.solution.orientation.ω[1], results[95,:])
    append!(runtime.solution.orientation.ω[2], results[96,:])
    append!(runtime.solution.orientation.ω[3], results[97,:])

    append!(runtime.solution.orientation.pos_pp[1], results[16,:])
    append!(runtime.solution.orientation.pos_pp[2], results[17,:])
    append!(runtime.solution.orientation.pos_pp[3], results[18,:])
    append!(runtime.solution.orientation.pos_pp_mag, results[19,:])
    append!(runtime.solution.orientation.vel_pp[1], results[20,:])
    append!(runtime.solution.orientation.vel_pp[2], results[21,:])
    append!(runtime.solution.orientation.vel_pp[3], results[22,:])
    append!(runtime.solution.orientation.vel_pp_mag, results[23,:])

    append!(runtime.solution.orientation.oe[1], results[24,:])
    append!(runtime.solution.orientation.oe[2], results[25,:])
    append!(runtime.solution.orientation.oe[3], results[26,:])
    append!(runtime.solution.orientation.oe[4], results[27,:])
    append!(runtime.solution.orientation.oe[5], results[28,:])
    append!(runtime.solution.orientation.oe[6], results[29,:])

    append!(runtime.solution.orientation.lat, results[30,:])
    append!(runtime.solution.orientation.lon, results[31,:])
    append!(runtime.solution.orientation.alt, results[32,:])
    append!(runtime.solution.orientation.γ_ii, results[33,:])
    append!(runtime.solution.orientation.γ_pp, results[34,:])

    append!(runtime.solution.orientation.h_ii[1], results[35,:])
    append!(runtime.solution.orientation.h_ii[2], results[36,:])
    append!(runtime.solution.orientation.h_ii[3], results[37,:])
    append!(runtime.solution.orientation.h_pp[1], results[38,:])
    append!(runtime.solution.orientation.h_pp[2], results[39,:])
    append!(runtime.solution.orientation.h_pp[3], results[40,:])
    append!(runtime.solution.orientation.h_ii_mag, results[41,:])
    append!(runtime.solution.orientation.h_pp_mag, results[42,:])

    append!(runtime.solution.orientation.uD[1], results[43,:])
    append!(runtime.solution.orientation.uD[2], results[44,:])
    append!(runtime.solution.orientation.uD[3], results[45,:])
    append!(runtime.solution.orientation.uE[1], results[46,:])
    append!(runtime.solution.orientation.uE[2], results[47,:])
    append!(runtime.solution.orientation.uE[3], results[48,:])
    append!(runtime.solution.orientation.uN[1], results[49,:])
    append!(runtime.solution.orientation.uN[2], results[50,:])
    append!(runtime.solution.orientation.uN[3], results[51,:])
    append!(runtime.solution.orientation.vN, results[52,:])
    append!(runtime.solution.orientation.vE, results[53,:])
    append!(runtime.solution.orientation.azi_pp, results[54,:])

    # Physical properties
    append!(runtime.solution.physical_properties.ρ, results[55,:])
    append!(runtime.solution.physical_properties.T, results[56,:])
    append!(runtime.solution.physical_properties.p, results[57,:])
    append!(runtime.solution.physical_properties.wind[1], results[58,:])
    append!(runtime.solution.physical_properties.wind[2], results[59,:])
    append!(runtime.solution.physical_properties.wind[3], results[60,:])
    append!(runtime.solution.physical_properties.cL, results[61,:])
    append!(runtime.solution.physical_properties.cD, results[62,:])
    append!(runtime.solution.physical_properties.S, results[63,:])

    append!(runtime.solution.physical_properties.α_control, results[98,:])
    append!(runtime.solution.physical_properties.inertia_tensor[1], results[99,:]) # inertia tensor components
    append!(runtime.solution.physical_properties.inertia_tensor[2], results[100,:])
    append!(runtime.solution.physical_properties.inertia_tensor[3], results[101,:])
    append!(runtime.solution.physical_properties.inertia_tensor[4], results[102,:])
    append!(runtime.solution.physical_properties.inertia_tensor[5], results[103,:])  
    append!(runtime.solution.physical_properties.inertia_tensor[6], results[104,:])
    append!(runtime.solution.physical_properties.inertia_tensor[7], results[105,:])
    append!(runtime.solution.physical_properties.inertia_tensor[8], results[106,:])
    append!(runtime.solution.physical_properties.inertia_tensor[9], results[107,:])
    append!(runtime.solution.physical_properties.τ_rw[1], results[108,:]) # total reaction wheel torque τ_rw_x
    append!(runtime.solution.physical_properties.τ_rw[2], results[109,:]) # total reaction wheel torque τ_rw_y
    append!(runtime.solution.physical_properties.τ_rw[3], results[110,:]) # total reaction wheel torque τ_rw_z

    # Initialize α and β if they are not already initialized
    n_bodies = length(runtime.model.body.links)
    if isempty(runtime.solution.physical_properties.α)
        for i in 1:n_bodies
            append!(runtime.solution.physical_properties.α, [[]])
            append!(runtime.solution.physical_properties.β, [[]])
            append!(runtime.solution.performance.heat_rate, [[]])
            append!(runtime.solution.performance.heat_load, [[]])
        end
    end

    # Append α and β for each link
    for i in 1:n_bodies
        append!(runtime.solution.physical_properties.α[i], results[110 + i,:]) # α
        append!(runtime.solution.physical_properties.β[i], results[110 + n_bodies + i,:]) # β
        append!(runtime.solution.performance.heat_rate[i], results[110 + 2*n_bodies + i,:]) # heat rate
        append!(runtime.solution.performance.heat_load[i], results[110 + 3*n_bodies + i,:]) # heat load
    end


    n_reaction_wheels = runtime.model.body.n_reaction_wheels
    n_thrusters = runtime.model.body.n_thrusters
    # Initialize the reaction wheel properties if they are not already initialized
    if isempty(runtime.solution.physical_properties.rw_h)
        for i in 1:n_reaction_wheels
            append!(runtime.solution.physical_properties.rw_h, [[]])
            append!(runtime.solution.physical_properties.rw_τ, [[]])
        end
    end

    # Append reaction wheel angular momentum and torque for each reaction wheel
    for i in 1:n_reaction_wheels
        append!(runtime.solution.physical_properties.rw_h[i], results[110 + 4*n_bodies + i,:])
        append!(runtime.solution.physical_properties.rw_τ[i], results[110 + 4*n_bodies + n_reaction_wheels + i,:]) # rw_τ
    end

    # Initialize thruster forces if they are not already initialized
    if isempty(runtime.solution.physical_properties.thruster_forces)
        for i in 1:n_thrusters
            append!(runtime.solution.physical_properties.thruster_forces, [[]])
        end
    end
    
    # Append thruster forces for each thruster
    for i in 1:n_thrusters
        append!(runtime.solution.physical_properties.thruster_forces[i], results[110 + 4*n_bodies + 2*n_reaction_wheels + i,:]) # thruster forces
    end
    

    # Performance
    append!(runtime.solution.performance.mass, results[64,:])
    # append!(runtime.solution.performance.heat_rate, results[65,:])
    # append!(runtime.solution.performance.heat_load, results[66,:])
    append!(runtime.solution.performance.T_r, results[65,:])
    append!(runtime.solution.performance.q, results[66,:])

    # Forces
    append!(runtime.solution.forces.gravity_ii[1], results[67,:])
    append!(runtime.solution.forces.gravity_ii[2], results[68,:])
    append!(runtime.solution.forces.gravity_ii[3], results[69,:])
    append!(runtime.solution.forces.drag_pp[1], results[70,:])
    append!(runtime.solution.forces.drag_pp[2], results[71,:])
    append!(runtime.solution.forces.drag_pp[3], results[72,:])
    append!(runtime.solution.forces.drag_ii[1], results[73,:])
    append!(runtime.solution.forces.drag_ii[2], results[74,:])
    append!(runtime.solution.forces.drag_ii[3], results[75,:])
    append!(runtime.solution.forces.lift_pp[1], results[76,:])
    append!(runtime.solution.forces.lift_pp[2], results[77,:])
    append!(runtime.solution.forces.lift_pp[3], results[78,:])
    append!(runtime.solution.forces.lift_ii[1], results[79,:])
    append!(runtime.solution.forces.lift_ii[2], results[80,:])
    append!(runtime.solution.forces.lift_ii[3], results[81,:])
    append!(runtime.solution.forces.force_ii[1], results[82,:])
    append!(runtime.solution.forces.force_ii[2], results[83,:])
    append!(runtime.solution.forces.force_ii[3], results[84,:])
    append!(runtime.solution.forces.τ_ii[1], results[85,:])
    append!(runtime.solution.forces.τ_ii[2], results[86,:])
    append!(runtime.solution.forces.τ_ii[3], results[87,:])
    append!(runtime.solution.forces.energy, results[88,:])

    # Simulation
    append!(runtime.solution.simulation.MC_seed, results[89,:])
    append!(runtime.solution.simulation.drag_passage, results[90,:])

    return time_0
end

function save_results(params::ODEParams)
    # println("Saveing results at time: ", params.intermediate_solution.time)
    solution = params.solution
    intermediate_solution = params.intermediate_solution
    model = params.m
    cnf = params.cnf
    push!(solution.orientation.time, intermediate_solution.time)
    push!(solution.orientation.year, intermediate_solution.year)
    push!(solution.orientation.month, intermediate_solution.month)
    push!(solution.orientation.day, intermediate_solution.day)
    push!(solution.orientation.hour, intermediate_solution.hour)
    push!(solution.orientation.minute, intermediate_solution.minute)
    push!(solution.orientation.second, intermediate_solution.second)
    push!(solution.orientation.number_of_passage, intermediate_solution.number_of_passage)
    push!(solution.orientation.pos_ii[1], intermediate_solution.pos_ii[1])
    push!(solution.orientation.pos_ii[2], intermediate_solution.pos_ii[2])
    push!(solution.orientation.pos_ii[3], intermediate_solution.pos_ii[3])
    push!(solution.orientation.vel_ii[1], intermediate_solution.vel_ii[1])
    push!(solution.orientation.vel_ii[2], intermediate_solution.vel_ii[2])
    push!(solution.orientation.vel_ii[3], intermediate_solution.vel_ii[3])
    push!(solution.orientation.pos_ii_mag, intermediate_solution.pos_ii_mag)
    push!(solution.orientation.vel_ii_mag, intermediate_solution.vel_ii_mag)
    push!(solution.orientation.quaternion[1], intermediate_solution.quaternion[1])
    push!(solution.orientation.quaternion[2], intermediate_solution.quaternion[2])
    push!(solution.orientation.quaternion[3], intermediate_solution.quaternion[3])
    push!(solution.orientation.quaternion[4], intermediate_solution.quaternion[4])
    push!(solution.orientation.ω[1], intermediate_solution.ω[1])
    push!(solution.orientation.ω[2], intermediate_solution.ω[2])
    push!(solution.orientation.ω[3], intermediate_solution.ω[3])

    push!(solution.orientation.pos_pp[1], intermediate_solution.pos_pp[1])
    push!(solution.orientation.pos_pp[2], intermediate_solution.pos_pp[2])
    push!(solution.orientation.pos_pp[3], intermediate_solution.pos_pp[3])
    push!(solution.orientation.pos_pp_mag, intermediate_solution.pos_pp_mag)
    push!(solution.orientation.vel_pp[1], intermediate_solution.vel_pp[1])
    push!(solution.orientation.vel_pp[2], intermediate_solution.vel_pp[2])
    push!(solution.orientation.vel_pp[3], intermediate_solution.vel_pp[3])
    push!(solution.orientation.vel_pp_mag, intermediate_solution.vel_pp_mag)

    push!(solution.orientation.oe[1], intermediate_solution.oe[1])
    push!(solution.orientation.oe[2], intermediate_solution.oe[2])
    push!(solution.orientation.oe[3], intermediate_solution.oe[3])
    push!(solution.orientation.oe[4], intermediate_solution.oe[4])
    push!(solution.orientation.oe[5], intermediate_solution.oe[5])
    push!(solution.orientation.oe[6], intermediate_solution.oe[6])

    push!(solution.orientation.lat, intermediate_solution.lat)
    push!(solution.orientation.lon, intermediate_solution.lon)
    push!(solution.orientation.alt, intermediate_solution.alt)
    push!(solution.orientation.γ_ii, intermediate_solution.γ_ii)
    push!(solution.orientation.γ_pp, intermediate_solution.γ_pp)

    push!(solution.orientation.h_ii[1], intermediate_solution.h_ii[1])
    push!(solution.orientation.h_ii[2], intermediate_solution.h_ii[2])
    push!(solution.orientation.h_ii[3], intermediate_solution.h_ii[3])
    push!(solution.orientation.h_pp[1], intermediate_solution.h_pp[1])
    push!(solution.orientation.h_pp[2], intermediate_solution.h_pp[2])
    push!(solution.orientation.h_pp[3], intermediate_solution.h_pp[3])
    push!(solution.orientation.h_ii_mag, intermediate_solution.h_ii_mag)
    push!(solution.orientation.h_pp_mag, intermediate_solution.h_pp_mag)

    push!(solution.orientation.uD[1], intermediate_solution.uD[1])
    push!(solution.orientation.uD[2], intermediate_solution.uD[2])
    push!(solution.orientation.uD[3], intermediate_solution.uD[3])
    push!(solution.orientation.uE[1], intermediate_solution.uE[1])
    push!(solution.orientation.uE[2], intermediate_solution.uE[2])
    push!(solution.orientation.uE[3], intermediate_solution.uE[3])
    push!(solution.orientation.uN[1], intermediate_solution.uN[1])
    push!(solution.orientation.uN[2], intermediate_solution.uN[2])
    push!(solution.orientation.uN[3], intermediate_solution.uN[3])
    push!(solution.orientation.vN, intermediate_solution.vN)
    push!(solution.orientation.vE, intermediate_solution.vE)
    push!(solution.orientation.azi_pp, intermediate_solution.azi_pp)

    # Physical properties
    push!(solution.physical_properties.ρ, intermediate_solution.ρ)
    push!(solution.physical_properties.T, intermediate_solution.T)
    push!(solution.physical_properties.p, intermediate_solution.p)
    push!(solution.physical_properties.wind[1], intermediate_solution.wind[1])
    push!(solution.physical_properties.wind[2], intermediate_solution.wind[2])
    push!(solution.physical_properties.wind[3], intermediate_solution.wind[3])
    push!(solution.physical_properties.cL, intermediate_solution.cL)
    push!(solution.physical_properties.cD, intermediate_solution.cD)
    push!(solution.physical_properties.S, intermediate_solution.S)

    push!(solution.physical_properties.α_control, intermediate_solution.α_control)
    push!(solution.physical_properties.inertia_tensor[1], intermediate_solution.inertia_tensor[1]) # inertia tensor components
    push!(solution.physical_properties.inertia_tensor[2], intermediate_solution.inertia_tensor[2])
    push!(solution.physical_properties.inertia_tensor[3], intermediate_solution.inertia_tensor[3])
    push!(solution.physical_properties.inertia_tensor[4], intermediate_solution.inertia_tensor[4])
    push!(solution.physical_properties.inertia_tensor[5], intermediate_solution.inertia_tensor[5])  
    push!(solution.physical_properties.inertia_tensor[6], intermediate_solution.inertia_tensor[6])
    push!(solution.physical_properties.inertia_tensor[7], intermediate_solution.inertia_tensor[7])
    push!(solution.physical_properties.inertia_tensor[8], intermediate_solution.inertia_tensor[8])
    push!(solution.physical_properties.inertia_tensor[9], intermediate_solution.inertia_tensor[9])
    push!(solution.physical_properties.τ_rw[1], intermediate_solution.τ_rw[1]) # total reaction wheel torque τ_rw_x
    push!(solution.physical_properties.τ_rw[2], intermediate_solution.τ_rw[2]) # total reaction wheel torque τ_rw_y
    push!(solution.physical_properties.τ_rw[3], intermediate_solution.τ_rw[3]) # total reaction wheel torque τ_rw_z

    # Initialize α and β if they are not already initialized
    n_bodies = length(model.body.links)
    if isempty(solution.physical_properties.α)
        for i in 1:n_bodies
            append!(solution.physical_properties.α, [Float64[]])
            append!(solution.physical_properties.β, [Float64[]])
            append!(solution.performance.heat_rate, [Float64[]])
            append!(solution.performance.heat_load, [Float64[]])
        end
    end
    # Append α and β for each link
    for i in 1:n_bodies
        push!(solution.physical_properties.α[i], intermediate_solution.α[i]) # α
        push!(solution.physical_properties.β[i], intermediate_solution.β[i]) # β
        push!(solution.performance.heat_rate[i], intermediate_solution.heat_rate[i]) # heat rate
        push!(solution.performance.heat_load[i], intermediate_solution.heat_load[i]) # heat load
    end

    n_reaction_wheels = model.body.n_reaction_wheels
    n_thrusters = model.body.n_thrusters
    # Initialize the reaction wheel properties if they are not already initialized
    if isempty(solution.physical_properties.rw_h)
        for i in 1:n_reaction_wheels
            append!(solution.physical_properties.rw_h, [Float64[]])
            append!(solution.physical_properties.rw_τ, [Float64[]])
        end
    end
    # Append reaction wheel properties for each reaction wheel
    for i in 1:n_reaction_wheels
        push!(solution.physical_properties.rw_h[i], intermediate_solution.rw_h[i])
        push!(solution.physical_properties.rw_τ[i], intermediate_solution.rw_τ[i])
    end


    # Initialize thruster forces if they are not already initialized
    if isempty(solution.physical_properties.thruster_forces)
        for i in 1:n_thrusters
            append!(solution.physical_properties.thruster_forces, [Float64[]])
        end
    end
    # Append thruster forces for each thruster
    for i in 1:n_thrusters
        push!(solution.physical_properties.thruster_forces[i], intermediate_solution.thruster_forces[i])
    end

    # Performance
    push!(solution.performance.mass, intermediate_solution.mass)
    push!(solution.performance.T_r, intermediate_solution.T_r)
    push!(solution.performance.q, intermediate_solution.dynamic_pressure)

    # Forces
    push!(solution.forces.gravity_ii[1], intermediate_solution.gravity_ii[1])
    push!(solution.forces.gravity_ii[2], intermediate_solution.gravity_ii[2])
    push!(solution.forces.gravity_ii[3], intermediate_solution.gravity_ii[3])
    push!(solution.forces.drag_pp[1], intermediate_solution.drag_pp[1])
    push!(solution.forces.drag_pp[2], intermediate_solution.drag_pp[2])
    push!(solution.forces.drag_pp[3], intermediate_solution.drag_pp[3])
    push!(solution.forces.drag_ii[1], intermediate_solution.drag_ii[1])
    push!(solution.forces.drag_ii[2], intermediate_solution.drag_ii[2])
    push!(solution.forces.drag_ii[3], intermediate_solution.drag_ii[3])
    push!(solution.forces.lift_pp[1], intermediate_solution.lift_pp[1])
    push!(solution.forces.lift_pp[2], intermediate_solution.lift_pp[2])
    push!(solution.forces.lift_pp[3], intermediate_solution.lift_pp[3])
    push!(solution.forces.lift_ii[1], intermediate_solution.lift_ii[1])
    push!(solution.forces.lift_ii[2], intermediate_solution.lift_ii[2])
    push!(solution.forces.lift_ii[3], intermediate_solution.lift_ii[3])
    push!(solution.forces.force_ii[1], intermediate_solution.force_ii[1])
    push!(solution.forces.force_ii[2], intermediate_solution.force_ii[2])
    push!(solution.forces.force_ii[3], intermediate_solution.force_ii[3])
    push!(solution.forces.τ_body[1], intermediate_solution.τ_body[1])
    push!(solution.forces.τ_body[2], intermediate_solution.τ_body[2])
    push!(solution.forces.τ_body[3], intermediate_solution.τ_body[3])
    push!(solution.forces.energy, intermediate_solution.energy)

    # Simulation
    push!(solution.simulation.MC_seed, intermediate_solution.MC_index)
    push!(solution.simulation.drag_passage, intermediate_solution.drag_state)

    # return intermediate_solution.time[end]
end

function clean_results()
    return Solution()
end
