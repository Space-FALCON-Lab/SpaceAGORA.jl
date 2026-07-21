struct PaperObservation
    drag_passage_time_s::Float64
    apoapsis_radius_m::Float64
    periapsis_altitude_m::Float64
    argument_of_periapsis_rad::Float64
    raan_rad::Float64
    inclination_rad::Float64
    date_ordinal::Float64
    max_density_kg_m3::Float64
    max_heat_rate_w_cm2::Float64
end

paper_observation_names() = [
    :drag_passage_time_s,
    :apoapsis_radius_m,
    :periapsis_altitude_m,
    :argument_of_periapsis_rad,
    :raan_rad,
    :inclination_rad,
    :date_ordinal,
    :max_density_kg_m3,
    :max_heat_rate_w_cm2,
]

raw_observation_vector(obs::PaperObservation) = [
    obs.drag_passage_time_s,
    obs.apoapsis_radius_m,
    obs.periapsis_altitude_m,
    obs.argument_of_periapsis_rad,
    obs.raan_rad,
    obs.inclination_rad,
    obs.date_ordinal,
    obs.max_density_kg_m3,
    obs.max_heat_rate_w_cm2,
]

python_ordinal(date::Date) = Dates.value(date - Date(1, 1, 1)) + 1
python_ordinal(dt::DateTime) = python_ordinal(Date(dt))
