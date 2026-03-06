using CSV
using DataFrames
using Arrow
using .SimulationModel: Model, Cnf, Solution

function save_csv(filename::String, args::Dict, arrow_filename::String, params::Tuple)
    cnf = params[1]
    model = params[2]
    solution = params[3]
    if args[:save_csv]
        touch(filename)

        writer = open(filename, "a")
    end
    # periapsis_alt = vcat(cnf.altitude_periapsis, zeros(length(solution.orientation.time) - length(cnf.altitude_periapsis)))
    println("Preparing data to save to CSV and Arrow files...")
    println("Max time: ", maximum(solution.orientation.time))
    data_push = DataFrame(time = solution.orientation.time,
                          year = solution.orientation.year,
                          month = solution.orientation.month,
                          day = solution.orientation.day,
                          hour = solution.orientation.hour,
                          min = solution.orientation.minute,
                          sec = solution.orientation.second,
                          number_of_passage = solution.orientation.number_of_passage,
                          pos_ii_1 = solution.orientation.pos_ii[1],
                          pos_ii_2 = solution.orientation.pos_ii[2],
                          pos_ii_3 = solution.orientation.pos_ii[3],
                          vel_ii_1 = solution.orientation.vel_ii[1],
                          vel_ii_2 = solution.orientation.vel_ii[2],
                          vel_ii_3 = solution.orientation.vel_ii[3],
                          pos_ii_mag = solution.orientation.pos_ii_mag,
                          vel_ii_mag = solution.orientation.vel_ii_mag,
                          q_1 = solution.orientation.quaternion[1],
                          q_2 = solution.orientation.quaternion[2],
                          q_3 = solution.orientation.quaternion[3],
                          q_4 = solution.orientation.quaternion[4],
                          omega_1 = solution.orientation.ω[1],
                          omega_2 = solution.orientation.ω[2],
                          omega_3 = solution.orientation.ω[3],
                          pos_pp_1 = solution.orientation.pos_pp[1],
                          pos_pp_2 = solution.orientation.pos_pp[2],
                          pos_pp_3 = solution.orientation.pos_pp[3],
                          pos_pp_mag = solution.orientation.pos_pp_mag,
                          vel_pp_1 = solution.orientation.vel_pp[1],
                          vel_pp_2 = solution.orientation.vel_pp[2],
                          vel_pp_3 = solution.orientation.vel_pp[3],
                          vel_pp_mag = solution.orientation.vel_pp_mag,
                          a = solution.orientation.oe[1],
                          e = solution.orientation.oe[2],
                          i = solution.orientation.oe[3],
                          OMEGA = solution.orientation.oe[4],
                          omega = solution.orientation.oe[5],
                          vi = solution.orientation.oe[6],
                          lat = solution.orientation.lat,
                          lon = solution.orientation.lon,
                          alt = solution.orientation.alt,
                          gamma_ii = solution.orientation.γ_ii,
                          gamma_pp = solution.orientation.γ_pp,
                          h_ii_1 = solution.orientation.h_ii[1],
                          h_ii_2 = solution.orientation.h_ii[2],
                          h_ii_3 = solution.orientation.h_ii[3],
                          h_pp_1 = solution.orientation.h_pp[1],
                          h_pp_2 = solution.orientation.h_pp[2],
                          h_pp_3 = solution.orientation.h_pp[3],
                          h_ii_mag = solution.orientation.h_ii_mag,
                          h_pp_mag = solution.orientation.h_pp_mag,
                          uD_1 = solution.orientation.uD[1],
                          uD_2 = solution.orientation.uD[2],
                          uD_3 = solution.orientation.uD[3],
                          uE_1 = solution.orientation.uE[1],
                          uE_2 = solution.orientation.uE[2],
                          uE_3 = solution.orientation.uE[3],
                          uN_1 = solution.orientation.uN[1],
                          uN_2 = solution.orientation.uN[2],
                          uN_3 = solution.orientation.uN[3],
                          vN = solution.orientation.vN,
                          vE = solution.orientation.vE,
                          azi_pp = solution.orientation.azi_pp,
                          rho = solution.physical_properties.ρ,
                          T = solution.physical_properties.T,
                          p = solution.physical_properties.p,
                          wind_1 = solution.physical_properties.wind[1],
                          wind_2 = solution.physical_properties.wind[2],
                          wind_3 = solution.physical_properties.wind[3],
                          cL = solution.physical_properties.cL,
                          cD = solution.physical_properties.cD,
                          aoa_control = solution.physical_properties.α_control,
                          S = solution.physical_properties.S,
                          mass = solution.performance.mass,
                        #   heat_rate = solution.performance.heat_rate,
                        #   heat_load = solution.performance.heat_load,
                          T_r = solution.performance.T_r,
                          q = solution.performance.q,
                          gravity_ii_1 = solution.forces.gravity_ii[1],
                          gravity_ii_2 = solution.forces.gravity_ii[2],
                          gravity_ii_3 = solution.forces.gravity_ii[3],
                          drag_pp_1 = solution.forces.drag_pp[1],
                          drag_pp_2 = solution.forces.drag_pp[2],
                          drag_pp_3 = solution.forces.drag_pp[3],
                          drag_ii_1 = solution.forces.drag_ii[1],
                          drag_ii_2 = solution.forces.drag_ii[2],
                          drag_ii_3 = solution.forces.drag_ii[3],
                          lift_pp_1 = solution.forces.lift_pp[1],
                          lift_pp_2 = solution.forces.lift_pp[2],
                          lift_pp_3 = solution.forces.lift_pp[3],
                          lift_ii_1 = solution.forces.lift_ii[1],
                          lift_ii_2 = solution.forces.lift_ii[2],
                          lift_ii_3 = solution.forces.lift_ii[3],
                          force_ii_1 = solution.forces.force_ii[1],
                          force_ii_2 = solution.forces.force_ii[2],
                          force_ii_3 = solution.forces.force_ii[3],
                          tau_ii_1 = solution.forces.τ_body[1],
                          tau_ii_2 = solution.forces.τ_body[2],
                          tau_ii_3 = solution.forces.τ_body[3],
                          J_ii_1 = solution.physical_properties.inertia_tensor[1],
                          J_ii_2 = solution.physical_properties.inertia_tensor[2],
                          J_ii_3 = solution.physical_properties.inertia_tensor[3],
                          J_ii_4 = solution.physical_properties.inertia_tensor[4],
                          J_ii_5 = solution.physical_properties.inertia_tensor[5],
                          J_ii_6 = solution.physical_properties.inertia_tensor[6],
                          J_ii_7 = solution.physical_properties.inertia_tensor[7],
                          J_ii_8 = solution.physical_properties.inertia_tensor[8],
                          J_ii_9 = solution.physical_properties.inertia_tensor[9],
                          energy = solution.forces.energy,
                          rw_torque_ii_1 = solution.physical_properties.τ_rw[1],
                          rw_torque_ii_2 = solution.physical_properties.τ_rw[2],
                          rw_torque_ii_3 = solution.physical_properties.τ_rw[3],
                          t_cf = zeros(length(solution.orientation.time)),
                          h_cf = zeros(length(solution.orientation.time)),
                          gamma_cf = zeros(length(solution.orientation.time)),
                          v_cf = zeros(length(solution.orientation.time)))
    
    # Save physical properties with varying length based on the number of bodies
    n_bodies = length(model.body.links)
    for i in 1:n_bodies
        data_push[!, Symbol("link_$(i)_aoa")] = solution.physical_properties.α[i]
        data_push[!, Symbol("link_$(i)_sideslip")] = solution.physical_properties.β[i]
        data_push[!, Symbol("link_$(i)_heat_rate")] = solution.performance.heat_rate[i]
        data_push[!, Symbol("link_$(i)_heat_load")] = solution.performance.heat_load[i]
    end

    # Save properties based on the number of reaction wheels
    for i in 1:model.body.n_reaction_wheels
        data_push[!, Symbol("rw_h_$(i)")] = solution.physical_properties.rw_h[i]
        data_push[!, Symbol("rw_tau_$(i)")] = solution.physical_properties.rw_τ[i]
    end

    for i in 1:model.body.n_thrusters
        data_push[!, Symbol("thruster_force_$(i)")] = solution.physical_properties.thruster_forces[i]
    end

    println("data_push: ", size(data_push[!, :t_cf]), " cf: ", length(solution.closed_form.t_cf))
    
    if args[:closed_form] == 1
        data_push[!, :t_cf] = solution.closed_form.t_cf
        data_push[!, :h_cf] = solution.closed_form.h_cf
        data_push[!, :gamma_cf] = solution.closed_form.γ_cf
        data_push[!, :v_cf] = solution.closed_form.v_cf
    end
    if filesize(filename) == 0.0 && args[:save_csv]
        # Write header if file is empty
        CSV.write(filename, data_push, writeheader=true)
    elseif args[:save_csv]
        CSV.write(filename, data_push, append=true)
    end
    # Write to Arrow file for plotting
    if args[:print_res]
        println("Writing data to Arrow file...")
    end
    # println(solution.orientation.number_of_passage[1])
    # temp_file = joinpath(temp_name, "data$(uuid4()).arrow")
    Arrow.append(arrow_filename, data_push)
    if args[:print_res]
        println("Data written to Arrow file.")
    end
end