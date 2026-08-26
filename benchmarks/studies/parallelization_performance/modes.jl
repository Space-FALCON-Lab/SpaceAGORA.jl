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
    # Overrides the SPACEAGORA_RHS_CALIBRATE value that policy_adaptive would
    # otherwise imply. Exists so pre-solve plan calibration can be separated from
    # the rest of the adaptive machinery in an attribution run; `nothing` keeps
    # the derived behaviour and is what every shipped mode uses.
    calibrate::Union{Nothing, String} = nothing
    # Overrides all five inner callback/effector parallel-mode knobs at once
    # (density/control/thermal/multibody/effector). Exists so an attribution run
    # can switch off the *inner* policy surface -- and with it the per-callback
    # decision and observation bookkeeping -- without touching the outer route
    # selection, which is a separate axis. `nothing` keeps each mode's own
    # per-knob values and is what every shipped mode uses.
    inner_modes::Union{Nothing, String} = nothing
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
            # Stays pinned even though R3 also declares :auto -- this mode is the
            # fixed hybrid baseline the adaptive profiles are scored against, so
            # it must not move with the router.
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
            # "auto" hands the outer-route choice to select_outer_route!, which
            # is what R4 declares (profile_definitions.jl: outer_backend=:auto).
            # It was pinned to "threads" here, which silently removed the
            # process route from R4/R5's reachable set. See
            # ppc_resolve_outer_backend.
            backend="auto",
            outer_active=true,
            policy_adaptive=true,
            rhs_batch="auto",
            density="auto",
            control="auto",
            thermal="auto",
            multibody="auto",
            effector="auto",
            # Mirrors profile_definitions.jl, which moved R4/R5 to the static
            # scheduler after the dynamic default measured 7.0% slower on
            # gravity-only constellations (paired probe, 21 pairs, p = 0.0072).
            # The harness must not diverge from the shipped profile or it
            # measures a configuration nobody actually runs. Dynamic is still
            # reachable through calibration, which crosses both schedulers with
            # the allotment ladder and pins whichever measurably wins.
            scheduler="static",
            allow_inner_with_outer=true
        ),
        "full_smart" => PPCModeSpec(
            name="full_smart",
            profile="R5",
            # See outer_inner_adaptive above: R5 declares outer_backend=:auto.
            backend="auto",
            outer_active=true,
            policy_adaptive=true,
            rhs_batch="auto",
            density="auto",
            control="auto",
            thermal="on",
            multibody="auto",
            effector="auto",
            scheduler="static",
            persistent=true,
            allow_inner_with_outer=true
        ),

        # --- Attribution variants of full_smart -------------------------------
        # Not part of any reported ladder. R5 on a process-routed Monte Carlo
        # workload is ~60% slower per sample than the pinned process route, and
        # three mechanisms are switched on together in R5 and off together in
        # outer_process, so the raw pair cannot say which one costs. These two
        # peel them apart one at a time against an otherwise identical R5.
        "full_smart_nocalib" => PPCModeSpec(
            name="full_smart_nocalib",
            profile="R5", backend="auto", outer_active=true, policy_adaptive=true,
            rhs_batch="auto", density="auto", control="auto", thermal="on",
            multibody="auto", effector="auto", scheduler="static", persistent=true,
            allow_inner_with_outer=true,
            # Only difference from full_smart: no pre-solve plan sweep.
            calibrate="off"
        ),
        "full_smart_noinner" => PPCModeSpec(
            name="full_smart_noinner",
            profile="R5", backend="auto", outer_active=true, policy_adaptive=true,
            rhs_batch="auto", density="auto", control="auto", thermal="on",
            multibody="auto", effector="auto", scheduler="static", persistent=true,
            # Only difference from full_smart: the inner policy is not consulted
            # while an outer split is active, which is what outer_process does.
            allow_inner_with_outer=false
        ),
        # Third attribution arm. nocalib and noinner between them recovered only
        # 2-3% of the 63.5% gap to the pinned process route, so neither plan
        # calibration nor inner-policy consultation is the cost; both of those
        # variants left policy_adaptive on. This one turns it off and changes
        # nothing else, which isolates the adaptive decision machinery itself
        # (SPACEAGORA_PARALLEL_POLICY_ADAPTIVE and MEASURED_REWARD, the latter
        # timing candidate configurations rather than merely branching).
        # Deliberately incoherent as a shipped profile -- an adaptive profile
        # with adaptation disabled -- which is exactly what makes it a probe.
        "full_smart_nopolicy" => PPCModeSpec(
            name="full_smart_nopolicy",
            profile="R5", backend="auto", outer_active=true,
            policy_adaptive=false,
            rhs_batch="auto", density="auto", control="auto", thermal="on",
            multibody="auto", effector="auto", scheduler="static", persistent=true,
            allow_inner_with_outer=true,
            # policy_adaptive=false would otherwise drag calibration off with it;
            # pinning it to auto keeps this arm's only difference the policy itself.
            calibrate="auto"
        ),
        # Fourth and fifth attribution arms, added for the post-fix re-run on
        # montecarlo_heavy_aerobraking.
        #
        # The three arms above cannot separate the remaining gap on a
        # *process-routed* workload, because every Distributed worker is launched
        # with --threads=1 and therefore has inner thread budget 1. At budget 1
        # pre-solve calibration returns immediately (rhs_calibration.jl:334) and
        # allow_inner_with_outer cannot change any answer (thread_policy_decision
        # forces use_threads=false), so `nocalib` and `noinner` are expected to be
        # no-ops there by construction rather than by measurement. What is *not*
        # a no-op at budget 1 is the per-callback bookkeeping: the decision path
        # still records telemetry under a lock on every call, and
        # record_policy_observation! -- which has no budget<=1 short-circuit --
        # still runs its adaptive branch whenever the knob is "auto". Both are
        # driven by the five inner mode knobs being non-"off", which is an axis
        # none of the three arms above vary.
        "full_smart_innermodes_off" => PPCModeSpec(
            name="full_smart_innermodes_off",
            profile="R5", backend="auto", outer_active=true, policy_adaptive=true,
            rhs_batch="auto", density="auto", control="auto", thermal="on",
            multibody="auto", effector="auto", scheduler="static", persistent=true,
            allow_inner_with_outer=true,
            # Only difference from full_smart: the five inner parallel-mode knobs
            # are off, so no inner decision or observation bookkeeping is entered
            # at all. Outer route selection is untouched -- backend is still
            # "auto" and SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE is still 1 -- so
            # this arm is still routing adaptively.
            inner_modes="off"
        ),
        "full_smart_nohints" => PPCModeSpec(
            name="full_smart_nohints",
            profile="R5", backend="auto", outer_active=true, policy_adaptive=true,
            rhs_batch="auto", density="auto", control="auto", thermal="on",
            multibody="auto", effector="auto", scheduler="static",
            # Only difference from full_smart: persistent hints (and with them
            # the measured-reward branch, which requires persistent_hints_enabled)
            # are off. Separates the hint/measured-reward machinery from the base
            # per-callback decision+observation bookkeeping that innermodes_off
            # removes wholesale.
            persistent=false,
            allow_inner_with_outer=true
        )
    )
end

# `outer_tasks` is how many outer units of work (MC samples / constellation
# members dispatched by this harness) will run concurrently under this mode.
# Pass it wherever it is known; -1 means "unknown", which keeps the mode's own
# declared value.
#
# It exists because SPACEAGORA_OUTER_PARALLEL_ACTIVE is not a label for "this is
# a parallel profile" -- the RHS router reads it as the factual claim "an
# enclosing outer split already owns the thread pool, so do not start a nested
# one", and responds by clamping the RHS to a single worker (setup.jl's
# outer_serialized branches). Asserting it unconditionally for every
# threads/process-backed mode therefore pinned the RHS at allotment 1 for the
# single-simulation constellation cases that make up most of B1/B2/B3/B5/B6 --
# cases that dispatch exactly one outer task and so have no outer split at all.
# The claim has to track the actual sample count, not the mode name.
function ppc_mode_env_pairs(
    mode::PPCModeSpec,
    cfg::PPCConfig;
    outer_tasks::Int=-1,
)::Vector{Pair{String, Union{Nothing, String}}}
    persist = mode.persistent ? "1" : "0"
    outer_active = mode.outer_active && outer_tasks != 1
    inner_override = mode.inner_modes
    density_mode = inner_override === nothing ? mode.density : inner_override
    control_mode = inner_override === nothing ? mode.control : inner_override
    thermal_mode = inner_override === nothing ? mode.thermal : inner_override
    multibody_mode = inner_override === nothing ? mode.multibody : inner_override
    effector_mode = inner_override === nothing ? mode.effector : inner_override
    return Pair{String, Union{Nothing, String}}[
        "SPACEAGORA_PARALLEL_PROFILE" => mode.profile,
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => mode.backend,
        "SPACEAGORA_PERF_PROCS" => string(cfg.process_workers),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => (outer_active ? "1" : "0"),
        "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE" => (mode.policy_adaptive ? "1" : "0"),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => (mode.policy_adaptive ? "1" : "0"),
        "SPACEAGORA_RHS_BATCH_PARALLEL" => mode.rhs_batch,
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => density_mode,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => control_mode,
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => thermal_mode,
        "SPACEAGORA_MULTIBODY_PARALLEL" => multibody_mode,
        "SPACEAGORA_EFFECTOR_PARALLEL" => effector_mode,
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
        # Pre-solve RHS-plan auto-calibration follows the mode's own adaptive
        # flag rather than being disabled everywhere.
        #
        # It was previously forced off for every mode, on the reasoning that it
        # short-circuits the routing policy via rhs_plan_override and so
        # contaminates a profile-ladder comparison. That reasoning holds for the
        # *static* profiles (R0-R3), which are meant to measure one fixed route.
        # Applying it to R4/R5 as well disabled the very mechanism those profiles
        # exist to evaluate -- calibration IS part of "the full adaptive routing
        # policy SpaceAGORA ships with".
        #
        # The cost of getting this wrong was large. On heavy_1024sat_fullstack_1hr
        # at 12 threads: 29.90 s with calibration off vs 5.71 s with it on (5.2x),
        # i.e. 0.96x vs 4.2x against the serial baseline. Every "multi-effector
        # constellations do not scale" measurement taken with the blanket-off
        # setting understated the shipped configuration by that factor.
        "SPACEAGORA_RHS_CALIBRATE" => (mode.calibrate === nothing ?
            (mode.policy_adaptive ? "auto" : "off") : mode.calibrate),
        "SPACEAGORA_SOLVER_MODE" => cfg.solver_mode,
        "SPACEAGORA_SAVE_BUNDLE" => "0",
        "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
        "OPENBLAS_NUM_THREADS" => "1",
        "GKSwstype" => "100"
    ]
end

function ppc_effective_env_string(mode::PPCModeSpec, cfg::PPCConfig; outer_tasks::Int=-1)::String
    pairs = ppc_mode_env_pairs(mode, cfg; outer_tasks=outer_tasks)
    return join(["$(p.first)=$(p.second)" for p in pairs], ";")
end
