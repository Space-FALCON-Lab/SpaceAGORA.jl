const _VENUS_MATRIX_SCENARIOS = join([
    "venus_j0_tbfalse",
    "venus_j0_tbtrue",
    "venus_j2_tbfalse",
    "venus_j2_tbtrue",
    "venus_j50_tbfalse",
    "venus_j50_tbtrue",
], ",")

ENV["SPACEAGORA_GMAT_SCENARIOS"] = _VENUS_MATRIX_SCENARIOS

include(joinpath(@__DIR__, "gmat_scenario_matrix.jl"))
