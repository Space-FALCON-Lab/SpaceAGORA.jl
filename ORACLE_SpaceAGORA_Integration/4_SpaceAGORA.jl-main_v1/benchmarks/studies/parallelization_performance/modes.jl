Base.@kwdef struct PPCModeSpec
    name::String
    profile::String
    backend::String
    outer_active::Bool
    policy_adaptive::Bool
    rhs_batch::String
    density::String
    control::String
    thermal::String
    multibody::String
    effector::String
    scheduler::String = "static"
    persistent::Bool = false
    allow_inner_with_outer::Bool = false
end

function ppc_mode_specs()::Dict{String, PPCModeSpec}
    serial = PPCModeSpec(
        name="serial",
        profile="R0",
        backend="none",
        outer_active=false,
        policy_adaptive=false,
        rhs_batch="off",
        density="off",
        control="off",
        thermal="off",
        multibody="off",
        effector="off"
    )
    return Dict(
        "serial" => serial,
        "outer_threads" => PPCModeSpec(
            name="outer_threads",
            profile="R1a",
            backend="threads",
            outer_active=true,
            policy_adaptive=false,
            rhs_batch="auto",
            density="off",
            control="off",
            thermal="off",
            multibody="off",
            effector="off"
        ),
        "outer_process" => PPCModeSpec(
            name="outer_process",
            profile="R1b",
            backend="process",
            outer_active=true,
            policy_adaptive=false,
            rhs_batch="auto",
            density="off",
            control="off",
            thermal="off",
            multibody="off",
            effector="off"
        ),
        "inner_only" => PPCModeSpec(
            name="inner_only",
            profile="R2",
            backend="none",
            outer_active=false,
            policy_adaptive=false,
            rhs_batch="auto",
            density="auto",
            control="auto",
            thermal="auto",
            multibody="auto",
            effector="auto"
        ),
        "outer_inner_static" => PPCModeSpec(
            name="outer_inner_static",
            profile="R3",
            backend="threads",
            outer_active=true,
            policy_adaptive=false,
            rhs_batch="auto",
            density="auto",
            control="auto",
            thermal="auto",
            multibody="auto",
            effector="auto",
            allow_inner_with_outer=true
        ),
        "outer_inner_adaptive" => PPCModeSpec(
            name="outer_inner_adaptive",
            profile="R4",
            backend="threads",
            outer_active=true,
            policy_adaptive=true,
            rhs_batch="auto",
            density="auto",
            control="auto",
            thermal="auto",
            multibody="auto",
            effector="auto",
            scheduler="dynamic",
            allow_inner_with_outer=true
        ),
        "full_smart" => PPCModeSpec(
            name="full_smart",
            profile="R5",
            backend="threads",
            outer_active=true,
            policy_adaptive=true,
            rhs_batch="auto",
            density="auto",
            control="auto",
            thermal="on",
            multibody="auto",
            effector="auto",
            scheduler="dynamic",
            persistent=true,
            allow_inner_with_outer=true
        )
    )
end

function ppc_mode_env_pairs(mode::PPCModeSpec, cfg::PPCConfig)::Vector{Pair{String, Union{Nothing, String}}}
    persist = mode.persistent ? "1" : "0"
    return Pair{String, Union{Nothing, String}}[
        "SPACEAGORA_PARALLEL_PROFILE" => mode.profile,
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => mode.backend,
        "SPACEAGORA_PERF_PROCS" => string(cfg.process_workers),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => (mode.outer_active ? "1" : "0"),
        "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE" => (mode.policy_adaptive ? "1" : "0"),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => (mode.policy_adaptive ? "1" : "0"),
        "SPACEAGORA_RHS_BATCH_PARALLEL" => mode.rhs_batch,
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => mode.density,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => mode.control,
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => mode.thermal,
        "SPACEAGORA_MULTIBODY_PARALLEL" => mode.multibody,
        "SPACEAGORA_EFFECTOR_PARALLEL" => mode.effector,
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => (mode.allow_inner_with_outer ? "1" : "0"),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => (mode.allow_inner_with_outer ? "1" : "0"),
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => (mode.allow_inner_with_outer ? "1" : "0"),
        "SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER" => (mode.allow_inner_with_outer ? "1" : "0"),
        "SPACEAGORA_EFFECTOR_PARALLEL_ALLOW_WITH_OUTER" => (mode.allow_inner_with_outer ? "1" : "0"),
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
        "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => mode.scheduler,
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => persist,
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => persist,
        "SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD" => (mode.policy_adaptive ? "1" : "0"),
        "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => (mode.name == "full_smart" ? "1" : "0"),
        "SPACEAGORA_SOLVER_MODE" => cfg.solver_mode,
        "SPACEAGORA_SAVE_BUNDLE" => "0",
        "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
        "OPENBLAS_NUM_THREADS" => "1",
        "GKSwstype" => "100"
    ]
end

function ppc_effective_env_string(mode::PPCModeSpec, cfg::PPCConfig)::String
    pairs = ppc_mode_env_pairs(mode, cfg)
    return join(["$(p.first)=$(p.second)" for p in pairs], ";")
end
