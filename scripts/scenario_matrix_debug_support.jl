# Maintained definition loader for the diagnostic entrypoints. No matrix runs,
# plots, exports, or generated source writes occur when this module is included.
module ScenarioMatrixDebugSupport

const SCENARIO_MATRIX_DEFINITIONS_ONLY = true
include(joinpath(@__DIR__, "..", "test", "gmat_scenario_matrix.jl"))

function require_matrix_inputs(; stk::Bool=false)
    _basilisk_reference_available() || throw(ArgumentError(
        "Basilisk parity references missing from $(_BASILISK_REFERENCE_DIR); $(_FETCH_REFERENCES_HINT)"))
    scenarios = _active_basilisk_expected_scenario_names()
    isempty(scenarios) && throw(ArgumentError(
        "No matrix reference scenarios found in $(_BASILISK_REFERENCE_DIR); $(_FETCH_REFERENCES_HINT)"))
    if stk && !_stk_reference_available()
        throw(ArgumentError("STK references missing from $(_STK_RESULTS_DIR); supply the STK comparison CSVs before running this diagnostic."))
    end
    for scenario in sort!(collect(scenarios))
        paths = stk ? (_scenario_basilisk_path(scenario), _scenario_stk_path(scenario)) :
                      (_scenario_basilisk_path(scenario),)
        for path in paths
            isfile(path) || throw(ArgumentError("Missing matrix reference: $path"))
        end
    end
    for relative_path in (_GMAT_HARMONICS_EARTH_FILE, _GMAT_HARMONICS_MARS_FILE,
                          _GMAT_HARMONICS_VENUS_FILE, _GMAT_HARMONICS_MOON_FILE)
        path = joinpath(_GMAT_REPO_ROOT, relative_path)
        isfile(path) || throw(ArgumentError("Missing gravity coefficients: $path"))
    end
    _gmat_planetary_kernel_relpath()
    return nothing
end

end # module ScenarioMatrixDebugSupport
