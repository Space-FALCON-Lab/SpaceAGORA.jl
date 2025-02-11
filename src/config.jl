module config

export model, cnf, solution, Body, Planet, Initial_condition, Aerodynamics, Engines, Model

    @kwdef mutable struct Body
        mass::Float64 = 0.0
        length_SA::Float64 = 0.0
        height_SA::Float64 = 0.0
        area_SA::Float64 = 0.0
        length_SC::Float64 = 0.0
        height_SC::Float64 = 0.0
        area_SC::Float64 = 0.0
        area_tot::Float64 = 0.0
        δ::Float64 = 0.0
        nose_radius::Float64 = 0.0
        base_radius::Float64 = 0.0
    end

    @kwdef mutable struct Planet
        Rp_e::Float64 = 0.0
        Rp_p::Float64 = 0.0
        Rp_m::Float64 = 0.0
        mass::Float64  = 0.0
        p::Float64 = 0.0
        k::Float64 = 0.0
        ω::Vector{Float64} = [0.0, 0.0, 0.0]
        g_ref::Float64 = 0.0
        ρ_ref::Float64 = 0.0
        h_ref::Float64 = 0.0
        H::Float64 = 0.0
        R::Float64 = 0.0
        γ::Float64 = 0.0
        T::Float64 = 0.0
        J2::Float64 = 0.0
        μ::Float64 = 0.0
        μ_fluid::Float64 = 0.0
        Lz::Float64 = 0.0
        name::String = ""
    end

    @kwdef mutable struct Aerodynamics
        δ::Float64 = 0.0
        α::Float64 = 0.0
        thermal_accomodation_factor::Float64 = 0.0
        reflection_coefficient::Float64 = 0.0
        thermal_contact::Float64 = 0.0
        heat_rate_limit::Float64 = 0.0
        heat_load_limit::Float64 = 0.0
    end

    @kwdef mutable struct Engines
        ϕ::Float64 = 0.0
        g_e::Float64 = 0.0
        T::Float64 = 0.0
        Isp::Float64 = 0.0
    end

    @kwdef mutable struct Initial_condition
        a::Float64 = 0.0
        e::Float64 = 0.0
        i::Float64 = 0.0
        Ω::Float64 = 0.0
        ω::Float64 = 0.0
        vi::Float64 = 0.0
        m::Float64 = 0.0
        year::Int64 = 0
        month::Int64 = 0
        day::Int64 = 0
        hour::Int64 = 0
        minute::Int64 = 0
        second::Float64 = 0.0
        time_rot::Float64 = 0.0
    end

    @kwdef mutable struct Model 
        body::Body = Body()
        planet::Planet = Planet()
        aerodynamics::Aerodynamics = Aerodynamics()
        engines::Engines = Engines()
        initial_condition::Initial_condition = Initial_condition()
    end

    model = Model()

    @kwdef mutable struct Cnf
        impact::Bool = false
        altitude_periapsis::Vector{Float64} = []
        max_heatrate::Vector{Float64} = []
        solution_intermediate::Vector{Any} = []
        atmospheric_data::Dict{String,Float64} = Dict()
        previous_atmospheric_data::Dict{String,Float64} = Dict()
        drag_state::Bool = false
        ascending_phase::Bool = false
        evaluate_switch_heat_load::Bool = false
        security_mode::Bool = false
        time_IEI::Float64 = 0.0
        time_OEI::Float64 = 0.0
        time_switch_1::Float64 = 0.0
        time_switch_2::Float64 = 0.0
        state_inner_boundary_atmosphere::Vector{Float64} = []
        count_aerobraking::Int64 = 0
        count_dori::Int64 = 0
        count_phase::Int64 = 0
        count_numberofpassage::Int64 = 0
        count_overcome_hr::Int64 = 0
        counter_random::Int64 = 0
        save_index_heat::Int64 = 0
        index_warning_alt::Int64 = 0
        index_warning_flow::Int64 = 0
        index_Mars_Gram_call::Int64 = 0
        index_MonteCarlo::Int64 = 1
        index_propellant_mass::Int64 = 1
        T_w::Float64 = 4.0
        Δv_man::Float64 = 0.0
        closed_form_solution_off::Int64 = 1
        α::Float64 = pi/2
        α_past::Float64 = pi/2

        # Results to delete
        periapsis_list::Vector{Float64} = []
        Δv_list::Vector{Float64} = []
        orbit_number_list::Vector{Int64} = []
        heat_load_past::Float64 = 0.0
        heat_load_ppast::Float64 = 0.0
        state_flesh1::Vector{Vector{Float64}} = [[]]
        α_list::Vector{Float64} = []
        initial_position_closed_form::Vector{Float64} = []
        continue_simulation::Bool = true
        timer_revaluation::Float64 = 0.0
        MarsGram_recall::Int64 = 0
        heat_rate_prev::Float64 = 0.0
        sensible_loads::Bool = false
        counter_integrator::Int64 = 0
        prev_step_integrator::Float64 = 0.0
        initial_time_saved::Float64 = 0.0

        # Extra variables missing in python version
        counter::Int64 = 0
        heat_rate_limit::Float64 = 0.0
        time_OP::Float64 = 0.0
        time_IP::Float64 = 0.0
        Gram_justrecalled::Int64 = 0
        Gram_directory::String = ""
        heat_rate_list::Vector{Float64} = []
        stop_simulation::Bool = false
        results_save::Int64 = 0
        count_eventfirststep::Int64 = 0
        eventfirststep_periapsis::Int64 = 0
        count_eventsecondstep::Int64 = 0
        count_reached_EI::Int64 = 0
        count_reached_AE::Int64 = 0
        count_out_drag_passage::Int64 = 0
        count_in_drag_passage::Int64 = 0
        count_in_drag_passage_nt::Int64 = 0
        count_apoapsispoint::Int64 = 0
        count_periapsispoint::Int64 = 0
        count_impact::Int64 = 0
        count_apoapsisgreaterperiapsis::Int64 = 0
        count_stop_firing::Int64 = 0
        count_guidance::Int64 = 0
        count_heat_rate_check::Int64 = 0
        count_heat_load_check_exit::Int64 = 0
        t_out_drag_passage::Float64 = 0.0

        n_bodies_list::Vector{Planet} = []
        DU::Float64 = 0.0
        TU::Float64 = 0.0
        MU::Float64 = 0.0   
    end

    cnf = Cnf()

    @kwdef mutable struct Controller
        guidance_t_eval::Vector{Float64} = []
        count_controller::Int64 = 0
        count_prev_controller::Int64 = 0
        stored_state::Int64 = 0
        prev_time::Float64 = 0.0
        t::Float64 = 0.0
    end

    controller = Controller()

    @kwdef mutable struct Orientation
        time::Vector{Float64} = []
        year::Vector{Float64} = []
        month::Vector{Float64} = []
        day::Vector{Float64} = []
        hour::Vector{Float64} = []
        minute::Vector{Float64} = []
        second::Vector{Float64} = []
        number_of_passage::Vector{Int64} = []
        pos_ii::Vector{Vector{Float64}} = [[], [], []]
        vel_ii::Vector{Vector{Float64}} = [[], [], []]
        pos_ii_mag::Vector{Float64} = []
        vel_ii_mag::Vector{Float64} = []
        pos_pp::Vector{Vector{Float64}} = [[], [], []]
        pos_pp_mag::Vector{Float64} = []
        vel_pp::Vector{Vector{Float64}} = [[], [], []]
        vel_pp_mag::Vector{Float64} = []
        oe::Vector{Vector{Float64}} = [[], [], [], [], [], []]
        lat::Vector{Float64} = []
        lon::Vector{Float64} = []
        alt::Vector{Float64} = []
        γ_ii::Vector{Float64} = []
        γ_pp::Vector{Float64} = []
        h_ii::Vector{Vector{Float64}} = [[], [], []]
        h_pp::Vector{Vector{Float64}} = [[], [], []]
        h_ii_mag::Vector{Float64} = []
        h_pp_mag::Vector{Float64} = []
        uD::Vector{Vector{Float64}} = [[], [], []]
        uE::Vector{Vector{Float64}} = [[], [], []]
        uN::Vector{Vector{Float64}} = [[], [], []]
        vN::Vector{Float64} = []
        vE::Vector{Float64} = []
        azi_pp::Vector{Float64} = []
    end

    @kwdef mutable struct Physical_properties
        ρ::Vector{Float64} = []
        T::Vector{Float64} = []
        p::Vector{Float64} = []
        wind::Vector{Vector{Float64}} = [[], [], []] 
        cL::Vector{Float64} = []
        cD::Vector{Float64} = []
        α::Vector{Float64} = []
        S::Vector{Float64} = []
    end

    @kwdef mutable struct Performance
        mass::Vector{Float64} = []
        heat_rate::Vector{Float64} = []
        heat_load::Vector{Float64} = []
        T_r::Vector{Float64} = []
        q::Vector{Float64} = []
    end

    @kwdef mutable struct Forces
        gravity_ii::Vector{Vector{Float64}} = [[], [], []]
        drag_pp::Vector{Vector{Float64}} = [[], [], []]
        drag_ii::Vector{Vector{Float64}} = [[], [], []]
        lift_pp::Vector{Vector{Float64}} = [[], [], []]
        lift_ii::Vector{Vector{Float64}} = [[], [], []]
        force_ii::Vector{Vector{Float64}} = [[], [], []]
        energy::Vector{Float64} = []
    end

    @kwdef mutable struct Simulation
        MC_seed::Vector{Int64} = []
        drag_passage::Vector{Int64} = []
    end

    @kwdef mutable struct Closed_form
        t_cf::Vector{Float64} = []
        h_cf::Vector{Float64} = []
        γ_cf::Vector{Float64} = []
        v_cf::Vector{Float64} = []
    end

    @kwdef mutable struct Solution
        orientation::Orientation = Orientation()
        physical_properties::Physical_properties = Physical_properties()
        performance::Performance = Performance()
        forces::Forces = Forces()
        simulation::Simulation = Simulation()
        closed_form::Closed_form = Closed_form()
    end

    solution = Solution()

end