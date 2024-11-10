module config

export model, cnf, solution, Body, Planet, Initial_condition, Aerodynamics, Engines, Model

    mutable struct Body
        mass::Float64
        length_SA::Float64
        height_SA::Float64
        area_SA::Float64
        length_SC::Float64
        height_SC::Float64
        area_SC::Float64
        area_tot::Float64
        δ::Float64
        nose_radius::Float64
        base_radius::Float64
    end

    mutable struct Planet
        Rp_e::Float64
        Rp_p::Float64
        Rp_m::Float64
        mass::Float64
        p::Float64
        k::Float64
        ω::Vector{Float64}
        g_ref::Float64
        ρ_ref::Float64
        h_ref::Float64
        H::Float64
        R::Float64
        γ::Float64
        T::Float64
        J2::Float64
        μ::Float64
        μ_fluid::Float64
        Lz::Float64
    end

    mutable struct Aerodynamics
        δ::Float64
        α::Float64
        thermal_accomodation_factor::Float64
        reflection_coefficient::Float64
        thermal_contact::Float64
        heat_rate_limit::Float64
        heat_load_limit::Float64
    end

    mutable struct Engines
        ϕ::Float64
        g_e::Float64
        T::Float64
        Isp::Float64
    end

    mutable struct Initial_condition
        a::Float64
        e::Float64
        i::Float64
        Ω::Float64
        ω::Float64
        vi::Float64
        m::Float64
        year::Int64
        month::Int64
        day::Int64
        hour::Int64
        minute::Int64
        second::Float64
        time_rot::Float64
    end

    mutable struct Model 
        body::Body
        planet::Planet
        aerodynamics::Aerodynamics
        engines::Engines
        initial_condition::Initial_condition
    end

    planet = Planet(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, [0.0, 0.0, 0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    body = Body(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    aerodynamics = Aerodynamics(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    engines = Engines(0.0, 0.0, 0.0, 0.0)
    initial_condition = Initial_condition(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, 0.0)
    model = Model(body, planet, aerodynamics, engines, initial_condition)

    mutable struct Cnf
        impact::Bool
        altitude_periapsis::Vector{Float64}
        max_heatrate::Vector{Float64}
        solution_intermediate::Vector{Any}
        atmospheric_data::Dict{String,Float64}
        previous_atmospheric_data::Dict{String,Float64}
        drag_state::Bool
        ascending_phase::Bool
        evaluate_switch_heat_load::Bool
        security_mode::Bool
        time_IEI::Float64
        time_OEI::Float64
        time_switch_1::Float64
        time_switch_2::Float64
        state_inner_boundary_atmosphere::Vector{Float64}
        count_aerobraking::Int64
        count_dori::Int64
        count_phase::Int64
        count_numberofpassage::Int64
        count_overcome_hr::Int64
        counter_random::Int64
        save_index_heat::Int64
        index_warning_alt::Int64
        index_warning_flow::Int64
        index_Mars_Gram_call::Int64
        index_MonteCarlo::Int64
        index_propellant_mass::Int64
        T_w::Float64
        Δv_man::Float64
        closed_form_solution_off::Int64
        α::Float64
        α_past::Float64

        # Results to delete
        periapsis_list::Vector{Float64}
        Δv_list::Vector{Float64}
        orbit_number_list::Vector{Int64}
        heat_load_past::Float64
        heat_load_ppast::Float64
        state_flesh1::Vector{Float64}
        α_list::Vector{Float64}
        initial_position_closed_form::Vector{Float64}
        continue_simulation::Bool
        timer_revaluation::Float64
        MarsGram_recall::Int64
        heat_rate_prev::Float64
        sensible_loads::Bool
        counter_integrator::Int64
        prev_step_integrator::Float64
        initial_time_saved::Float64

        # Extra variables missing in python version
        counter::Int64
        heat_rate_limit::Float64
        time_OP::Float64
        time_IP::Float64
        MarsGram_justrecalled::Int64
        heat_rate_list::Vector{Float64}
        stop_simulation::Bool
        results_save::Int64
        count_impact::Int64
        count_apoapsisgreaterperiapsis::Int64
    end

    cnf = Cnf(false, [], [], [], Dict(), Dict(), false, false, false, false, 0.0, 0.0, 0.0, 0.0, [], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 4.0, 0.0, 1, pi/2, pi/2, [], [], [], 0.0, 0.0, [], [], [], true, 0.0, 0, 0.0, false, 0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, [], false, 0, 0, 0)

    mutable struct Controller
        guidance_t_eval::Vector{Float64}
        count_controller::Int64
        count_prev_controller::Int64
        stored_state::Int64
        prev_time::Float64
        t::Float64
    end

    controller = Controller([], 0, 0, 0, 0.0, 0.0)

    mutable struct Orientation
        time::Vector{Float64}
        year::Vector{Int64}
        month::Vector{Int64}
        day::Vector{Int64}
        hour::Vector{Int64}
        minute::Vector{Int64}
        second::Vector{Float64}
        number_of_passage::Vector{Int64}
        pos_ii::Vector{Vector{Float64}}
        vel_ii::Vector{Vector{Float64}}
        pos_ii_mag::Vector{Float64}
        vel_ii_mag::Vector{Float64}
        pos_pp::Vector{Vector{Float64}}
        pos_pp_mag::Vector{Float64}
        vel_pp::Vector{Vector{Float64}}
        vel_pp_mag::Vector{Float64}
        oe::Vector{Vector{Float64}}
        lat::Vector{Float64}
        lon::Vector{Float64}
        alt::Vector{Float64}
        γ_ii::Vector{Float64}
        γ_pp::Vector{Float64}
        h_ii::Vector{Vector{Float64}}
        h_pp::Vector{Vector{Float64}}
        h_ii_mag::Vector{Float64}
        h_pp_mag::Vector{Float64}
        uD::Vector{Vector{Float64}}
        uE::Vector{Vector{Float64}}
        uN::Vector{Vector{Float64}}
        vN::Vector{Float64}
        vE::Vector{Float64}
        azi_pp::Vector{Float64}
    end

    mutable struct Physical_properties
        ρ::Vector{Float64}
        T::Vector{Float64}
        p::Vector{Float64}
        wind::Vector{Vector{Float64}}
        cL::Vector{Float64}
        cD::Vector{Float64}
        α::Vector{Float64}
        S::Vector{Float64}
    end

    mutable struct Performance
        mass::Vector{Float64}
        heat_rate::Vector{Float64}
        heat_load::Vector{Float64}
        T_r::Vector{Float64}
        q::Vector{Float64}
    end

    mutable struct Forces
        gravity_ii::Vector{Vector{Float64}}
        drag_pp::Vector{Vector{Float64}}
        drag_ii::Vector{Vector{Float64}}
        lift_pp::Vector{Vector{Float64}}
        lift_ii::Vector{Vector{Float64}}
        force_ii::Vector{Vector{Float64}}
        energy::Vector{Float64}
    end

    mutable struct Simulation
        MC_seed::Vector{Int64}
        drag_passage::Vector{Int64}
    end

    mutable struct Closed_form
        t_cf::Vector{Float64}
        h_cf::Vector{Float64}
        γ_cf::Vector{Float64}
        v_cf::Vector{Float64}
    end

    mutable struct Solution
        orientation::Orientation
        physical_properties::Physical_properties
        performance::Performance
        forces::Forces
        simulation::Simulation
        closed_form::Closed_form
    end

    orientation = Orientation([], [], [], [], [], [], [], [], [[],[],[]], [[],[],[]], [], [], [[],[],[]], [], [[],[],[]], [], [[],[],[],[],[],[]], [], [], [], [], [], [[],[],[]], [[],[],[]], [], [], [[],[],[]], [[],[],[]], [[],[],[]], [], [], [])
    physical_properties = Physical_properties([], [], [], [[],[],[]], [], [], [], [])
    performance = Performance([], [], [], [], [])
    forces = Forces([[],[],[]], [[],[],[]], [[],[],[]], [[],[],[]], [[],[],[]], [[],[],[]], [])
    simulation = Simulation([], [])
    closed_form = Closed_form([], [], [], []) 
    solution = Solution(orientation, physical_properties, performance, forces, simulation, closed_form)

end
