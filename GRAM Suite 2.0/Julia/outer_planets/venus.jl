_define_outer_planet_bindings(:venus, :V, BODY_VENUS, "Venus")

function venus_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_VENUS, "Venus")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    hydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    helium = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    oxygen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    nitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dinitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_monoxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_dioxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall(
        (:getGasesState_V, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Ref{GasesStateC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}),
        atmos.ptr, state, hydrogen, helium, oxygen, nitrogen, dinitrogen, carbon_monoxide, carbon_dioxide
    )
    return (
        state = state[],
        hydrogen = hydrogen[],
        helium = helium[],
        oxygen = oxygen[],
        nitrogen = nitrogen[],
        dinitrogen = dinitrogen[],
        carbon_monoxide = carbon_monoxide[],
        carbon_dioxide = carbon_dioxide[]
    )
end
