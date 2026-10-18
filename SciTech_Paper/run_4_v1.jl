using ComponentArrays
using DiffEqBase
using LinearAlgebra
using OrdinaryDiffEq
using StaticArrays
using TOML

const PAPER_DIR = @__DIR__
const V1_ROOT = normpath(joinpath(PAPER_DIR, "..", "4_SpaceAGORA.jl-main_v1"))
length(ARGS) in (0, 2, 4) ||
    error("Usage: run_4_v1.jl [helper_count orbit_count [schedule target_inclination_deg]]")
const N_HELPERS = isempty(ARGS) ? 20 : parse(Int, ARGS[1])
const ORBITS = isempty(ARGS) ? 5 : parse(Int, ARGS[2])
const SCHEDULE = length(ARGS) == 4 ? Symbol(ARGS[3]) : :gve_sma
const TARGET_INC_DEG = length(ARGS) == 4 ? parse(Float64, ARGS[4]) : 0.0
N_HELPERS > 0 && ORBITS > 0 || error("Helper and orbit counts must be positive.")
SCHEDULE in (:gve_sma, :gve_inc) || error("Supported comparison schedules: gve_sma, gve_inc.")
isfinite(TARGET_INC_DEG) && 0 <= TARGET_INC_DEG <= 180 ||
    error("Target inclination must be finite and between 0 and 180 degrees.")
const CASE_DIR = isempty(ARGS) ? "" : "$(N_HELPERS)_helpers_$(ORBITS)_orbits" *
    (length(ARGS) == 4 ? "_$(SCHEDULE)_target_inc_$(TARGET_INC_DEG)_deg" : "")
const MASS_KG = 227.0
const LASER_POWER_W = 10_000.0
const LASER_MAGNIFICATION = 100.0
const LASER_RANGE_M = 200e3
const DT_MAX_S = 10.0
const MU = 3.98600436233e14
const EARTH_RADIUS_M = 6_378_136.6

module V1Bridge
module SimulationModel
module AbstractTypes
abstract type AbstractForceTorqueModel end
end
module DynamicEffectors
function calcForceTorque end
function solver_partition end
include(joinpath(@__DIR__, "..", "4_SpaceAGORA.jl-main_v1", "src", "dynamics",
    "coupled", "force_torque_models", "laser_link_effectors.jl"))
end
end
end

const Laser = V1Bridge.SimulationModel.DynamicEffectors.LaserLinkEffectors

function circular_state(radius, true_anomaly_deg, inclination_deg=0.0)
    anomaly = deg2rad(true_anomaly_deg)
    inclination = deg2rad(inclination_deg)
    position = SVector(radius * cos(anomaly),
        radius * sin(anomaly) * cos(inclination), radius * sin(anomaly) * sin(inclination))
    speed = sqrt(MU / radius)
    velocity = SVector(-speed * sin(anomaly),
        speed * cos(anomaly) * cos(inclination), speed * cos(anomaly) * sin(inclination))
    return position, velocity
end

function inclination_deg_from_state(state)
    h = cross(state.pos, state.vel)
    return rad2deg(atan(hypot(h[1], h[2]), h[3]))
end

function sma_from_state(state)
    return inv(2 / norm(state.pos) - dot(state.vel, state.vel) / MU)
end

function rhs!(du, u, model, t)
    fill!(du, 0.0)
    positions, _ = Laser._state_vectors(u)
    force_totals = zeros(Float64, 6, length(positions))
    Laser.accumulate_laser_link_forces!(force_totals, model, positions, trues(length(positions)))
    for satellite in eachindex(positions)
        position = SVector{3, Float64}(u.sc[satellite].pos)
        du.sc[satellite].pos .= u.sc[satellite].vel
        du.sc[satellite].vel .= -MU * position / norm(position)^3 +
            SVector{3, Float64}(force_totals[1:3, satellite]) / u.sc[satellite].mass
    end
    return nothing
end

function run_case(; output_root=PAPER_DIR)
    helper_radius = EARTH_RADIUS_M + 1000e3
    target_radius = EARTH_RADIUS_M + 1050e3
    target_period = 2pi * sqrt(target_radius^3 / MU)
    duration = ORBITS * target_period
    positions_and_velocities = [circular_state(
        index == 1 ? target_radius : helper_radius,
        index == 1 ? 0.0 : 360.0 * (index - 2) / N_HELPERS,
        index == 1 ? TARGET_INC_DEG : 0.0)
        for index in 1:(N_HELPERS + 1)]
    initial_state = ComponentVector(sc=[
        (pos=collect(position), vel=collect(velocity), mass=MASS_KG)
        for (position, velocity) in positions_and_velocities])
    model = Laser.OpenCavityLaserLinkModel(target_idx=1,
        helper_indices=collect(2:(N_HELPERS + 1)), range_m=LASER_RANGE_M,
        power_w=LASER_POWER_W, magnification=LASER_MAGNIFICATION,
        beta=1.0, eta=1.0, schedule=SCHEDULE)

    schedule = Tuple{Float64, Int}[]
    function update_and_record!(integrator)
        helper = Laser.update_laser_link_schedule!(model, integrator.u)
        (isempty(schedule) || schedule[end][1] != integrator.t) &&
            push!(schedule, (Float64(integrator.t), helper))
        return nothing
    end
    scheduler = DiffEqBase.DiscreteCallback(
        (u, t, integrator) -> true,
        update_and_record!;
        initialize=(callback, u, t, integrator) -> update_and_record!(integrator))
    problem = ODEProblem(rhs!, initial_state, (0.0, duration), model)
    output_times = range(0.0, duration; length=1001)
    elapsed = @elapsed solution = solve(problem, Tsit5(); callback=scheduler,
        reltol=1e-9, abstol=1e-11, dtmax=DT_MAX_S, saveat=output_times)
    string(solution.retcode) == "Success" || error("Simulation failed: $(solution.retcode)")
    solution.t[end] == duration || error("Simulation ended before the requested duration.")

    final_state = solution.u[end].sc[1]
    initial_inc = inclination_deg_from_state(initial_state.sc[1])
    final_inc = inclination_deg_from_state(final_state)
    isapprox(initial_inc, TARGET_INC_DEG; atol=1e-10) ||
        error("Initial target inclination does not match the requested value.")
    activations = model.link_activation_count
    initial_helper = first(schedule)[2]
    initial_selection_end = findfirst(entry -> entry[2] != initial_helper, schedule)
    startup_activations = initial_helper > 0 && initial_selection_end !== nothing &&
        schedule[initial_selection_end][1] < 1.0 ? 1 : 0
    physical_activations = activations - startup_activations
    result = (
        code="4_SpaceAGORA.jl-main_v1",
        helper_count=N_HELPERS,
        schedule=string(SCHEDULE),
        initial_target_altitude_m=1050e3,
        initial_helper_altitude_m=1000e3,
        initial_helper_inclination_deg=0.0,
        initial_inclination_deg=initial_inc,
        final_inclination_deg=final_inc,
        delta_inclination_deg=final_inc - initial_inc,
        requested_target_orbits=ORBITS,
        duration_s=duration,
        target_period_s=target_period,
        initial_sma_m=target_radius,
        final_sma_m=sma_from_state(final_state),
        delta_sma_m=sma_from_state(final_state) - target_radius,
        final_eccentricity=norm(cross(final_state.vel,
            cross(final_state.pos, final_state.vel)) / MU -
            final_state.pos / norm(final_state.pos)),
        scheduler_updates=length(schedule),
        link_activations=activations,
        physical_link_activations=physical_activations,
        startup_link_activations=startup_activations,
        active_updates=count(entry -> entry[2] > 0, schedule),
        final_helper=model.active_helper_idx,
        runtime_s=elapsed,
        retcode=string(solution.retcode),
    )

    result_dir = joinpath(output_root, "results", CASE_DIR, "4_v1")
    mkpath(result_dir)
    open(joinpath(result_dir, "summary.toml"), "w") do io
        TOML.print(io, Dict{String, Any}(string(k) => v for (k, v) in pairs(result)))
    end
    open(joinpath(result_dir, "schedule.csv"), "w") do io
        println(io, "time_s,target_orbits,active_helper")
        for (time, helper) in schedule
            println(io, time, ",", time / target_period, ",", helper)
        end
    end
    println(result)
    return result
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_case()
