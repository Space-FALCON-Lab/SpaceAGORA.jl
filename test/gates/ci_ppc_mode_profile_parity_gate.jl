# Compare the nine base benchmark modes with the actual shipped profiles.
# Attribution/route-exploration arms and the intentional R3 backend override
# are outside this callback-mode and scheduler contract.

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODES_PATH = joinpath(REPO_ROOT, "benchmarks", "studies",
                            "parallelization_performance", "modes.jl")
const PROFILE_DEFINITIONS_PATH = joinpath(REPO_ROOT, "src", "orchestration", "parallel", "routing",
                                          "profile_definitions.jl")
isfile(MODES_PATH) || error("modes.jl not found at $(MODES_PATH)")
isfile(PROFILE_DEFINITIONS_PATH) ||
    error("profile_definitions.jl not found at $(PROFILE_DEFINITIONS_PATH)")

# These trusted repository definitions construct their values with Base alone.
# Fresh modules avoid package/native loading and check constructed values,
# including the separately constructed serial mode and default scheduler fields.
const PROFILE_SOURCE = Module(gensym(:PPCParityProfiles))
Base.include(PROFILE_SOURCE, PROFILE_DEFINITIONS_PATH)
const MODE_SOURCE = Module(gensym(:PPCParityModes))
# PPCConfig occurs only in signatures of environment/reporting helpers. Those
# helpers are not called by ppc_mode_specs(); no benchmark launcher is loaded.
Core.eval(MODE_SOURCE, :(struct PPCConfig end))
Base.include(MODE_SOURCE, MODES_PATH)
const MODE_SPECS = Base.invokelatest(getfield(MODE_SOURCE, :ppc_mode_specs))

const CHECKED = [
    ("serial", "R0"), ("outer_threads", "R1_a"), ("outer_process", "R1_b"),
    ("inner_only", "R2"), ("outer_inner_static", "R3"),
    ("outer_inner_adaptive", "R4"), ("full_smart", "R5"),
    ("policy_v2", "R6"), ("predictive", "R7"),
]
const CHECKED_FIELDS = (
    :thermal => :thermal_mode,
    :density => :density_mode,
    :control => :control_mode,
    :multibody => :multibody_mode,
    :effector => :effector_mode,
    :scheduler => :inner_scheduler,
)

violations = String[]
for (mode, profile) in CHECKED
    # Preserve the missing-mode diagnostic used by the existing PR297 controls.
    haskey(MODE_SPECS, mode) || error(
        "Benchmark mode table has no PPCModeSpec for checked mode $(repr(mode)) in modes.jl.")
    spec = MODE_SPECS[mode]
    cfg = Base.invokelatest(getfield(PROFILE_SOURCE, :profile_config), profile)
    if !hasproperty(spec, :profile)
        push!(violations, "modes.jl mode $(repr(mode)) is missing profile")
    else
        actual_profile = try
            Base.invokelatest(getfield(PROFILE_SOURCE, :parse_parallel_profile), spec.profile)
        catch err
            push!(violations, "modes.jl mode $(repr(mode)) has invalid profile " *
                              "$(repr(spec.profile)): $(sprint(showerror, err))")
            nothing
        end
        if actual_profile !== nothing && actual_profile != cfg.profile
            push!(violations, "modes.jl mode $(repr(mode)) names profile " *
                              "$(repr(spec.profile)) but requires $(profile)")
        end
    end
    for (field, shipped_field) in CHECKED_FIELDS
        if !hasproperty(spec, field)
            push!(violations, "modes.jl mode $(repr(mode)) is missing $(field)")
        elseif !hasproperty(cfg, shipped_field)
            push!(violations, "profile $(profile) is missing $(shipped_field)")
        else
            harness = getproperty(spec, field)
            shipped = getproperty(cfg, shipped_field)
            if harness != shipped
                push!(violations,
                      "modes.jl mode $(repr(mode)) sets $(field)=$(repr(harness)) but " *
                      "profile $(profile) ships $(field)=$(repr(shipped))")
            end
        end
    end
end

if !isempty(violations)
    error("Benchmark mode table disagrees with the shipped profiles:\n" *
          join(violations, "\n") *
          "\n\nFix modes.jl (or the profile), do not silence this: the harness " *
          "is what backs the paper's routing numbers.")
end

println("ppc_mode_profile_parity_gate_ok")
