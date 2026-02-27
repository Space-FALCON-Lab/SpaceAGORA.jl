_define_outer_planet_bindings(:uranus, :U, BODY_URANUS, "Uranus")

function uranus_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_URANUS, "Uranus")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    dihydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    methane = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    helium = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall((:getGasesState_U, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GasesStateC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}), atmos.ptr, state, dihydrogen, methane, helium)
    return (state = state[], dihydrogen = dihydrogen[], methane = methane[], helium = helium[])
end
