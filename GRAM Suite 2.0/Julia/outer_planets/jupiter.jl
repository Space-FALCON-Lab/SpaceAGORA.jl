_define_outer_planet_bindings(:jupiter, :J, BODY_JUPITER, "Jupiter")

function jupiter_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_JUPITER, "Jupiter")
    return (state = get_gases_state(atmos),)
end
