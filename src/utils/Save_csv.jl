using CSV
using DataFrames

function save_csv(filename, args)

    touch(filename)

    writer = open(filename, "a")
    # periapsis_alt = vcat(config.cnf.altitude_periapsis, zeros(length(config.solution.orientation.time) - length(config.cnf.altitude_periapsis)))

    data_push = DataFrame(time = config.solution.orientation.time,
                          year = config.solution.orientation.year,
                          month = config.solution.orientation.month,
                          day = config.solution.orientation.day,
                          hour = config.solution.orientation.hour,
                          min = config.solution.orientation.minute,
                          sec = config.solution.orientation.second,
                          number_of_passage = config.solution.orientation.number_of_passage,
                          pos_ii_1 = config.solution.orientation.pos_ii[1],
                          pos_ii_2 = config.solution.orientation.pos_ii[2],
                          pos_ii_3 = config.solution.orientation.pos_ii[3],
                          vel_ii_1 = config.solution.orientation.vel_ii[1],
                          vel_ii_2 = config.solution.orientation.vel_ii[2],
                          vel_ii_3 = config.solution.orientation.vel_ii[3],
                          pos_ii_mag = config.solution.orientation.pos_ii_mag,
                          vel_ii_mag = config.solution.orientation.vel_ii_mag,
                          q_1 = config.solution.orientation.quaternion[1],
                          q_2 = config.solution.orientation.quaternion[2],
                          q_3 = config.solution.orientation.quaternion[3],
                          q_4 = config.solution.orientation.quaternion[4],
                          omega_1 = config.solution.orientation.ω[1],
                          omega_2 = config.solution.orientation.ω[2],
                          omega_3 = config.solution.orientation.ω[3],
                          pos_pp_1 = config.solution.orientation.pos_pp[1],
                          pos_pp_2 = config.solution.orientation.pos_pp[2],
                          pos_pp_3 = config.solution.orientation.pos_pp[3],
                          pos_pp_mag = config.solution.orientation.pos_pp_mag,
                          vel_pp_1 = config.solution.orientation.vel_pp[1],
                          vel_pp_2 = config.solution.orientation.vel_pp[2],
                          vel_pp_3 = config.solution.orientation.vel_pp[3],
                          vel_pp_mag = config.solution.orientation.vel_pp_mag,
                          a = config.solution.orientation.oe[1],
                          e = config.solution.orientation.oe[2],
                          i = config.solution.orientation.oe[3],
                          OMEGA = config.solution.orientation.oe[4],
                          omega = config.solution.orientation.oe[5],
                          vi = config.solution.orientation.oe[6],
                          lat = config.solution.orientation.lat,
                          lon = config.solution.orientation.lon,
                          alt = config.solution.orientation.alt,
                          gamma_ii = config.solution.orientation.γ_ii,
                          gamma_pp = config.solution.orientation.γ_pp,
                          h_ii_1 = config.solution.orientation.h_ii[1],
                          h_ii_2 = config.solution.orientation.h_ii[2],
                          h_ii_3 = config.solution.orientation.h_ii[3],
                          h_pp_1 = config.solution.orientation.h_pp[1],
                          h_pp_2 = config.solution.orientation.h_pp[2],
                          h_pp_3 = config.solution.orientation.h_pp[3],
                          h_ii_mag = config.solution.orientation.h_ii_mag,
                          h_pp_mag = config.solution.orientation.h_pp_mag,
                          uD_1 = config.solution.orientation.uD[1],
                          uD_2 = config.solution.orientation.uD[2],
                          uD_3 = config.solution.orientation.uD[3],
                          uE_1 = config.solution.orientation.uE[1],
                          uE_2 = config.solution.orientation.uE[2],
                          uE_3 = config.solution.orientation.uE[3],
                          uN_1 = config.solution.orientation.uN[1],
                          uN_2 = config.solution.orientation.uN[2],
                          uN_3 = config.solution.orientation.uN[3],
                          vN = config.solution.orientation.vN,
                          vE = config.solution.orientation.vE,
                          azi_pp = config.solution.orientation.azi_pp,
                          rho = config.solution.physical_properties.ρ,
                          T = config.solution.physical_properties.T,
                          p = config.solution.physical_properties.p,
                          wind_1 = config.solution.physical_properties.wind[1],
                          wind_2 = config.solution.physical_properties.wind[2],
                          wind_3 = config.solution.physical_properties.wind[3],
                          cL = config.solution.physical_properties.cL,
                          cD = config.solution.physical_properties.cD,
                          aoa_control = config.solution.physical_properties.α_control,
                          S = config.solution.physical_properties.S,
                          mass = config.solution.performance.mass,
                        #   heat_rate = config.solution.performance.heat_rate,
                        #   heat_load = config.solution.performance.heat_load,
                          T_r = config.solution.performance.T_r,
                          q = config.solution.performance.q,
                          gravity_ii_1 = config.solution.forces.gravity_ii[1],
                          gravity_ii_2 = config.solution.forces.gravity_ii[2],
                          gravity_ii_3 = config.solution.forces.gravity_ii[3],
                          drag_pp_1 = config.solution.forces.drag_pp[1],
                          drag_pp_2 = config.solution.forces.drag_pp[2],
                          drag_pp_3 = config.solution.forces.drag_pp[3],
                          drag_ii_1 = config.solution.forces.drag_ii[1],
                          drag_ii_2 = config.solution.forces.drag_ii[2],
                          drag_ii_3 = config.solution.forces.drag_ii[3],
                          lift_pp_1 = config.solution.forces.lift_pp[1],
                          lift_pp_2 = config.solution.forces.lift_pp[2],
                          lift_pp_3 = config.solution.forces.lift_pp[3],
                          lift_ii_1 = config.solution.forces.lift_ii[1],
                          lift_ii_2 = config.solution.forces.lift_ii[2],
                          lift_ii_3 = config.solution.forces.lift_ii[3],
                          force_ii_1 = config.solution.forces.force_ii[1],
                          force_ii_2 = config.solution.forces.force_ii[2],
                          force_ii_3 = config.solution.forces.force_ii[3],
                          energy = config.solution.forces.energy,
                          rw_torque_ii_1 = config.solution.physical_properties.τ_rw[1],
                          rw_torque_ii_2 = config.solution.physical_properties.τ_rw[2],
                          rw_torque_ii_3 = config.solution.physical_properties.τ_rw[3],
                          t_cf = zeros(length(config.solution.orientation.time)),
                          h_cf = zeros(length(config.solution.orientation.time)),
                          gamma_cf = zeros(length(config.solution.orientation.time)),
                          v_cf = zeros(length(config.solution.orientation.time)))
    
    # Save physical properties with varying length based on the number of bodies
    n_bodies = length(config.model.body.links)
    for i in 1:n_bodies
        data_push[!, Symbol("link_$(i)_aoa")] = config.solution.physical_properties.α[i]
        data_push[!, Symbol("link_$(i)_sideslip")] = config.solution.physical_properties.β[i]
        data_push[!, Symbol("link_$(i)_heat_rate")] = config.solution.performance.heat_rate[i]
        data_push[!, Symbol("link_$(i)_heat_load")] = config.solution.performance.heat_load[i]
    end

    # Save properties based on the number of reaction wheels
    for i in 1:config.model.body.n_reaction_wheels
        data_push[!, Symbol("rw_h_$(i)")] = config.solution.physical_properties.rw_h[i]
        data_push[!, Symbol("rw_tau_$(i)")] = config.solution.physical_properties.rw_τ[i]
    end
    if args[:closed_form] == 1
        data_push[!, :t_cf] = config.solution.closed_form.t_cf
        data_push[!, :h_cf] = config.solution.closed_form.h_cf
        data_push[!, :gamma_cf] = config.solution.closed_form.γ_cf
        data_push[!, :v_cf] = config.solution.closed_form.v_cf
    end
    if filesize(filename) == 0
        # Write header if file is empty
        CSV.write(filename, data_push, writeheader=true)
    else
        CSV.write(filename, data_push, append=true)
    end

end