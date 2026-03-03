using Arrow
using Random
using SpaceAGORACalibration
using Tables
using TOML

function stage_rows(table, stage_name::String)
    return [row for row in Tables.rows(table) if String(row.stage) == stage_name]
end

mktempdir() do tmp
    base = default_spec()
    spec = CalibrationSpec(
        schema_version=base.schema_version,
        id="demo_resume_test",
        name=base.name,
        description=base.description,
        output_root=joinpath(tmp, "calibration_out"),
        seed=base.seed,
        objective=base.objective,
        verification_script=base.verification_script,
        manifest_paths=base.manifest_paths,
        scenario_weights=base.scenario_weights,
        parameters=base.parameters,
        budgets=base.budgets
    )

    backend = MockBackend(noise_sigma=0.0, fail_rate=0.0)
    first = run_calibration(spec, backend)
    second = run_calibration(spec, backend)

    @assert !isempty(first.run_id)
    @assert first.run_id == second.run_id
    @assert first.run_dir == second.run_dir
    @assert isfile(first.report_path)
    @assert isfinite(first.best_score)

    run_dir = first.run_dir
    @assert isfile(joinpath(run_dir, "spec.toml"))
    @assert isfile(joinpath(run_dir, "evaluations.arrow"))
    @assert isfile(joinpath(run_dir, "state.toml"))
    @assert isfile(joinpath(run_dir, "best_manifest.toml"))
    @assert isfile(joinpath(run_dir, "report.md"))

    row_count_1 = length(collect(Arrow.Table(joinpath(run_dir, "evaluations.arrow"))))
    row_count_2 = length(collect(Arrow.Table(joinpath(second.run_dir, "evaluations.arrow"))))
    @assert row_count_1 == row_count_2
end

mktempdir() do tmp
    base = default_spec()
    budgets = BudgetSpec(
        initial_samples=2,
        global_iters=0,
        batch_size=1,
        global_acquisition="ei",
        bo_pool_size=32,
        bo_length_scale=0.35,
        bo_noise=1.0e-6,
        bo_kappa=1.96,
        bo_xi=0.01,
        local_refine_strategy="trust_region",
        local_refine_topk=1,
        local_refine_steps=2,
        local_refine_neighbors=3,
        local_refine_init_scale=0.2,
        local_refine_shrink=0.6,
        local_refine_expand=1.2,
        local_refine_min_improvement=1.0e-8,
        robustness_samples=1,
        initial_design="lhs",
        robustness_uncertainty="none",
        uncertainty_atm_scale=0.0,
        uncertainty_ic_scale=0.0,
        robust_ranking="weighted",
        robust_cvar_alpha=0.9,
        robust_p95_weight=0.5,
        robust_fail_weight=5.0
    )
    spec = CalibrationSpec(
        schema_version=base.schema_version,
        id="local_refine_v3_test",
        name=base.name,
        description=base.description,
        output_root=joinpath(tmp, "calibration_out"),
        seed=base.seed,
        objective=base.objective,
        verification_script=base.verification_script,
        manifest_paths=base.manifest_paths,
        scenario_weights=base.scenario_weights,
        parameters=base.parameters,
        budgets=budgets
    )

    backend = MockBackend(noise_sigma=0.0, fail_rate=0.0)
    result = run_calibration(spec, backend)
    table = Arrow.Table(joinpath(result.run_dir, "evaluations.arrow"))
    local_rows = stage_rows(table, "local_refine_full")
    @assert length(local_rows) == budgets.local_refine_topk * budgets.local_refine_steps * budgets.local_refine_neighbors
end

mktempdir() do tmp
    base = default_spec()
    budgets = BudgetSpec(
        initial_samples=2,
        global_iters=0,
        batch_size=1,
        global_acquisition="ei",
        bo_pool_size=32,
        bo_length_scale=0.35,
        bo_noise=1.0e-6,
        bo_kappa=1.96,
        bo_xi=0.01,
        local_refine_strategy="trust_region",
        local_refine_topk=1,
        local_refine_steps=0,
        local_refine_neighbors=2,
        local_refine_init_scale=0.2,
        local_refine_shrink=0.6,
        local_refine_expand=1.2,
        local_refine_min_improvement=1.0e-8,
        robustness_samples=5,
        initial_design="lhs",
        robustness_uncertainty="normal",
        uncertainty_atm_scale=0.05,
        uncertainty_ic_scale=0.002,
        robust_ranking="cvar",
        robust_cvar_alpha=0.8,
        robust_p95_weight=0.5,
        robust_fail_weight=5.0
    )
    spec = CalibrationSpec(
        schema_version=base.schema_version,
        id="robustness_v4_test",
        name=base.name,
        description=base.description,
        output_root=joinpath(tmp, "calibration_out"),
        seed=base.seed,
        objective=base.objective,
        verification_script=base.verification_script,
        manifest_paths=base.manifest_paths,
        scenario_weights=base.scenario_weights,
        parameters=base.parameters,
        budgets=budgets
    )

    backend = MockBackend(noise_sigma=0.0, fail_rate=0.0)
    result = run_calibration(spec, backend)
    table = Arrow.Table(joinpath(result.run_dir, "evaluations.arrow"))
    robust_rows = stage_rows(table, "robustness_validation")
    @assert length(robust_rows) == budgets.local_refine_topk * max(1, budgets.robustness_samples)

    payloads = [
        TOML.parse(String(row.candidate_values))["payload"]
        for row in robust_rows
    ]
    for p in payloads
        @assert haskey(p, "robustness_replica")
    end
    @assert isfinite(result.best_score)
end

for acq in ("ei", "lcb")
    mktempdir() do tmp
        base = default_spec()
        budgets = BudgetSpec(
            initial_samples=3,
            global_iters=2,
            batch_size=1,
            global_acquisition=acq,
            bo_pool_size=48,
            bo_length_scale=0.35,
            bo_noise=1.0e-6,
            bo_kappa=1.96,
            bo_xi=0.01,
            local_refine_topk=1,
            local_refine_steps=0,
            robustness_samples=1,
            initial_design="random",
            robust_p95_weight=0.5,
            robust_fail_weight=5.0
        )
        spec = CalibrationSpec(
            schema_version=base.schema_version,
            id="global_bo_" * acq,
            name=base.name,
            description=base.description,
            output_root=joinpath(tmp, "calibration_out"),
            seed=base.seed,
            objective=base.objective,
            verification_script=base.verification_script,
            manifest_paths=base.manifest_paths,
            scenario_weights=base.scenario_weights,
            parameters=base.parameters,
            budgets=budgets
        )

        backend = MockBackend(noise_sigma=0.0, fail_rate=0.0)
        result = run_calibration(spec, backend)
        @assert isfile(joinpath(result.run_dir, "evaluations.arrow"))

        table = Arrow.Table(joinpath(result.run_dir, "evaluations.arrow"))
        n_global = count(row -> String(row.stage) == "global_search_quick", Tables.rows(table))
        @assert n_global == budgets.initial_samples + budgets.global_iters * budgets.batch_size
    end
end

mktempdir() do tmp
    base = default_spec()
    budgets = BudgetSpec(
        initial_samples=4,
        global_iters=0,
        batch_size=1,
        parallel_evaluations=3,
        global_acquisition="lcb",
        bo_pool_size=32,
        bo_length_scale=0.35,
        bo_noise=1.0e-6,
        bo_kappa=1.96,
        bo_xi=0.01,
        local_refine_topk=1,
        local_refine_steps=0,
        robustness_samples=0,
        initial_design="lhs",
        robust_p95_weight=0.5,
        robust_fail_weight=5.0
    )
    spec = CalibrationSpec(
        schema_version=base.schema_version,
        id="parallel_eval_engine_test",
        name=base.name,
        description=base.description,
        output_root=joinpath(tmp, "calibration_out"),
        seed=base.seed,
        objective=base.objective,
        verification_script=base.verification_script,
        manifest_paths=base.manifest_paths,
        scenario_weights=base.scenario_weights,
        parameters=base.parameters,
        budgets=budgets
    )

    backend = MockBackend(noise_sigma=0.0, fail_rate=0.0)
    result = run_calibration(spec, backend)
    table = Arrow.Table(joinpath(result.run_dir, "evaluations.arrow"))
    n_global = count(row -> String(row.stage) == "global_search_quick", Tables.rows(table))
    @assert n_global == budgets.initial_samples

    spec_doc = SpaceAGORACalibration.Spec.spec_to_dict(spec)
    @assert Int(spec_doc["budgets"]["parallel_evaluations"]) == budgets.parallel_evaluations
end

mktempdir() do tmp
    manifest_path = joinpath(tmp, "base_manifest.toml")
    open(manifest_path, "w") do io
        TOML.print(io, Dict(
            "version" => 1,
            "scenarios" => Any[
                Dict(
                    "name" => "test",
                    "value" => 3.0
                )
            ]
        ))
    end

    verifier_path = joinpath(tmp, "fake_verify.jl")
    open(verifier_path, "w") do io
        write(io, """
        using TOML

        function argval(flag::String)
            prefix = string(flag, "=")
            for a in ARGS
                startswith(a, prefix) && return split(a, "=", limit=2)[2]
            end
            error("missing required argument: " * flag)
        end

        manifest = argval("--manifest")
        out_summary = argval("--out-summary")
        out_errors = argval("--out-errors")

        doc = TOML.parsefile(manifest)
        value = Float64(doc["scenarios"][1]["value"])

        open(out_summary, "w") do io
            println(io, "scenario,event,nmae")
            println(io, "test,metric,\$(value)")
        end
        open(out_errors, "w") do io
            println(io, "scenario,event,error_km")
        end
        """)
    end

    cmd_spec = CalibrationSpec(
        id="command_backend_mapping_test",
        output_root=joinpath(tmp, "out"),
        verification_script=verifier_path,
        manifest_paths=[manifest_path],
        parameters=[
            ParameterSpec(
                name="value_scale",
                lower=2.0,
                upper=2.0,
                manifest_targets=["scenarios[name=test].value"],
                transform="mul"
            )
        ],
        budgets=BudgetSpec(
            initial_samples=1,
            global_iters=0,
            batch_size=1,
            local_refine_topk=1,
            local_refine_steps=0,
            robustness_samples=1
        )
    )

    candidate = Candidate(id=1, values=Dict{String, Any}("value_scale" => 2.0), stage="global_search_quick")
    cmd_backend = CommandBackend(
        julia_cmd=Base.julia_cmd(),
        project_path=joinpath(@__DIR__, ".."),
        verification_script=verifier_path,
        manifest_path=manifest_path,
        profile="quick",
        enforce=false
    )
    ev = evaluate_candidate(cmd_backend, cmd_spec, candidate; stage="global_search_quick", run_dir=tmp)
    @assert ev.success
    @assert isfinite(ev.score)
    @assert ev.score == 6.0

    tuned_manifest = TOML.parsefile(ev.artifacts["manifest_tuned"])
    tuned_value = Float64(tuned_manifest["scenarios"][1]["value"])
    @assert tuned_value == 6.0

    replica_candidate = Candidate(
        id=1,
        values=Dict{String, Any}("value_scale" => 2.0, "robustness_replica" => 2),
        stage="robustness_validation"
    )
    ev_rep = evaluate_candidate(cmd_backend, cmd_spec, replica_candidate; stage="robustness_validation", run_dir=tmp)
    @assert ev_rep.success
    @assert occursin("_r002", ev_rep.artifacts["summary"])
end

mktempdir() do tmp
    manifest_path = joinpath(tmp, "base_manifest.toml")
    open(manifest_path, "w") do io
        TOML.print(io, Dict(
            "version" => 1,
            "scenarios" => Any[
                Dict(
                    "name" => "test",
                    "value" => 5.0
                )
            ]
        ))
    end

    spec = CalibrationSpec(
        id="inprocess_backend_mapping_test",
        output_root=joinpath(tmp, "out"),
        manifest_paths=[manifest_path],
        parameters=[
            ParameterSpec(
                name="value_scale",
                lower=2.0,
                upper=2.0,
                manifest_targets=["scenarios[name=test].value"],
                transform="mul"
            )
        ],
        budgets=BudgetSpec(
            initial_samples=1,
            global_iters=0,
            batch_size=1,
            local_refine_topk=1,
            local_refine_steps=0,
            robustness_samples=1
        )
    )

    seen_requests = Any[]
    fake_runner = function (request)
        push!(seen_requests, request)
        doc = TOML.parsefile(request.manifest_path)
        value = Float64(doc["scenarios"][1]["value"])

        open(request.out_summary, "w") do io
            println(io, "scenario,event,nmae")
            println(io, "test,metric,$(value)")
        end
        open(request.out_errors, "w") do io
            println(io, "scenario,event,error_km")
        end
        return nothing
    end

    backend = InProcessBackend(
        run_verification=fake_runner,
        manifest_path=manifest_path,
        profile="quick",
        enforce=true,
        plots=false
    )
    candidate = Candidate(
        id=2,
        values=Dict{String, Any}("value_scale" => 2.0, "robustness_replica" => 1),
        stage="robustness_validation"
    )

    ev = evaluate_candidate(backend, spec, candidate; stage="robustness_validation", run_dir=tmp)
    @assert ev.success
    @assert ev.score == 10.0
    @assert occursin("_r001", ev.artifacts["summary"])

    @assert length(seen_requests) == 1
    request = seen_requests[1]
    @assert request.profile == :full
    @assert request.enforce == true
    @assert request.generate_plots == false
end

mktempdir() do tmp
    manifest_path = joinpath(tmp, "base_manifest.toml")
    open(manifest_path, "w") do io
        TOML.print(io, Dict(
            "version" => 1,
            "scenarios" => Any[
                Dict(
                    "name" => "test",
                    "value" => 1.0
                )
            ]
        ))
    end

    spec = CalibrationSpec(
        id="full_auto_runtime_controller_smoke",
        output_root=joinpath(tmp, "out"),
        manifest_paths=[manifest_path],
        parameters=[
            ParameterSpec(
                name="value_scale",
                lower=0.8,
                upper=1.2,
                transform="mul",
                manifest_targets=["scenarios[name=test].value"]
            )
        ],
        budgets=BudgetSpec(
            initial_samples=4,
            global_iters=2,
            batch_size=1,
            parallel_evaluations=3,
            local_refine_topk=1,
            local_refine_steps=1,
            local_refine_neighbors=2,
            robustness_samples=1
        )
    )

    fake_runner = function (request)
        doc = TOML.parsefile(request.manifest_path)
        value = Float64(doc["scenarios"][1]["value"])
        open(request.out_summary, "w") do io
            println(io, "scenario,event,nmae,pass,total_runtime_s")
            println(io, "test,metric,$(value),true,0.05")
        end
        open(request.out_errors, "w") do io
            println(io, "scenario,event,error_km")
        end
        return nothing
    end

    backend = InProcessBackend(
        run_verification=fake_runner,
        manifest_path=manifest_path,
        profile="quick",
        parallel_profile="R4_full_auto",
        enforce=false,
        plots=false
    )

    withenv(
        "SPACEAGORA_CALIBRATION_MACHINE_LABEL" => "ci_box",
        "SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS" => "3",
        "SPACEAGORA_CALIBRATION_AUTO_WARMUP_BATCHES" => "2",
        "SPACEAGORA_CALIBRATION_AUTO_REBALANCE_WINDOW" => "1"
    ) do
        result = run_calibration(spec, backend)
        @assert isfinite(result.best_score)
    end

    cache_dir = joinpath(spec.output_root, "runtime_policy_cache")
    @assert isdir(cache_dir)
    cache_files = filter(f -> endswith(f, ".toml"), readdir(cache_dir))
    @assert !isempty(cache_files)
    cache_doc = TOML.parsefile(joinpath(cache_dir, cache_files[1]))
    @assert String(get(cache_doc, "profile_name", "")) == "R4_full_auto"
    @assert !isempty(get(cache_doc, "stats", Any[]))
end

mktempdir() do tmp
    manifest_path = joinpath(tmp, "base_manifest.toml")
    open(manifest_path, "w") do io
        TOML.print(io, Dict(
            "version" => 1,
            "scenarios" => Any[
                Dict(
                    "name" => "test",
                    "value" => 2.0
                )
            ]
        ))
    end

    verifier_path = joinpath(tmp, "fake_verify_robust.jl")
    open(verifier_path, "w") do io
        write(io, """
        function argval(flag::String)
            prefix = string(flag, \"=\")
            for a in ARGS
                startswith(a, prefix) && return split(a, \"=\", limit=2)[2]
            end
            error(\"missing required argument: \" * flag)
        end

        out_summary = argval(\"--out-summary\")
        out_errors = argval(\"--out-errors\")

        open(out_summary, \"w\") do io
            println(io, \"scenario,event,rmse_km,max_abs_km,limit_max_rmse_km,limit_max_abs_km,limit_nmae,nmae,pass,total_runtime_s\")
            println(io, \"a,e1,2.0,1.0,1.0,1.0,1.0,2.0,true,50.0\")
            println(io, \"b,e2,0.5,0.5,1.0,1.0,1.0,0.5,true,50.0\")
        end
        open(out_errors, \"w\") do io
            println(io, \"scenario,event,error_km\")
        end
        """)
    end

    budgets = BudgetSpec(
        initial_samples=1,
        global_iters=0,
        batch_size=1,
        local_refine_topk=1,
        local_refine_steps=0,
        robustness_samples=1,
        objective_huber_delta=1.0,
        objective_lambda_fail=25.0,
        objective_lambda_time=2.0,
        objective_runtime_budget_s=200.0,
        objective_telemetry_noise_frac=0.0
    )
    spec = CalibrationSpec(
        id="command_backend_objective_contract_test",
        output_root=joinpath(tmp, "out"),
        verification_script=verifier_path,
        manifest_paths=[manifest_path],
        scenario_weights=Dict("a" => 2.0, "b" => 1.0),
        parameters=[
            ParameterSpec(
                name="value_scale",
                lower=1.0,
                upper=1.0,
                manifest_targets=["scenarios[name=test].value"],
                transform="mul"
            )
        ],
        budgets=budgets
    )

    candidate = Candidate(id=1, values=Dict{String, Any}("value_scale" => 1.0), stage="global_search_quick")
    cmd_backend = CommandBackend(
        julia_cmd=Base.julia_cmd(),
        project_path=joinpath(@__DIR__, ".."),
        verification_script=verifier_path,
        manifest_path=manifest_path,
        profile="quick",
        enforce=false
    )
    ev = evaluate_candidate(cmd_backend, spec, candidate; stage="global_search_quick", run_dir=tmp)
    @assert ev.success
    @assert abs(ev.score - 3.6875) < 1.0e-6
end

let
    PS = SpaceAGORACalibration.ParamSpace
    K = SpaceAGORACalibration.Spec
    p_cont = ParameterSpec(name="pc", kind=K.continuous, lower=0.0, upper=1.0, manifest_targets=["x"])
    p_int = ParameterSpec(name="pi", kind=K.integer, lower=1.0, upper=3.0, manifest_targets=["x"])
    p_cat = ParameterSpec(name="pk", kind=K.categorical, choices=["a", "b"], manifest_targets=["x"])
    rng = MersenneTwister(5)
    @assert 0.0 <= PS._sample_value(rng, p_cont) <= 1.0
    @assert PS._sample_value(rng, p_int) isa Int
    @assert PS._sample_value(rng, p_cat) in ("a", "b")
    p_int_empty = ParameterSpec(name="pi_bad", kind=K.integer, lower=2.7, upper=2.2, manifest_targets=["x"])
    @assert try
        PS._sample_value(rng, p_int_empty)
        false
    catch e
        e isa ArgumentError
    end

    @assert PS._value_from_unit(p_int, 0.0) == 1
    @assert PS._value_from_unit(p_int, 1.0) == 3
    p_cat_empty = ParameterSpec(name="pk_empty", kind=K.categorical, choices=String[], manifest_targets=["x"])
    @assert try
        PS._value_from_unit(p_cat_empty, 0.5)
        false
    catch e
        e isa ArgumentError
    end

    spec_dup = CalibrationSpec(
        id="paramspace_duplicates",
        parameters=[ParameterSpec(name="locked", kind=K.continuous, lower=1.0, upper=1.0, manifest_targets=["x"])],
        budgets=BudgetSpec(initial_samples=1, global_iters=0, batch_size=1, local_refine_topk=1, local_refine_steps=0, robustness_samples=0)
    )
    pop_dup = PS.sample_initial_population(MersenneTwister(9), spec_dup, 4; design="random")
    @assert length(pop_dup) == 4
    @assert all(c.values["locked"] == 1.0 for c in pop_dup)

    spec_mix = CalibrationSpec(
        id="paramspace_perturb",
        parameters=[
            ParameterSpec(name="ci", kind=K.continuous, lower=0.0, upper=10.0, manifest_targets=["x"]),
            ParameterSpec(name="ii", kind=K.integer, lower=1.0, upper=5.0, manifest_targets=["x"])
        ],
        budgets=BudgetSpec(initial_samples=1, global_iters=0, batch_size=1, local_refine_topk=1, local_refine_steps=0, robustness_samples=0)
    )
    base_mix = Candidate(id=1, values=Dict{String, Any}("ci" => 5.0, "ii" => 3))
    pert_mix = PS.perturb_candidate(MersenneTwister(11), spec_mix, base_mix, 2; scale=0.5, perturb_discrete=true)
    @assert 0.0 <= pert_mix.values["ci"] <= 10.0
    @assert 1 <= pert_mix.values["ii"] <= 5

    spec_cat = CalibrationSpec(
        id="paramspace_cat",
        parameters=[ParameterSpec(name="cc", kind=K.categorical, choices=["a", "b"], manifest_targets=["x"])],
        budgets=BudgetSpec(initial_samples=1, global_iters=0, batch_size=1, local_refine_topk=1, local_refine_steps=0, robustness_samples=0)
    )
    base_cat = Candidate(id=1, values=Dict{String, Any}("cc" => "a"))
    saw_changed = false
    saw_kept = false
    for seed in 1:128
        cand = PS.perturb_candidate(MersenneTwister(seed), spec_cat, base_cat, seed; perturb_discrete=true)
        if cand.values["cc"] == "a"
            saw_kept = true
        else
            saw_changed = true
        end
        saw_changed && saw_kept && break
    end
    @assert saw_changed
    @assert saw_kept
end

mktempdir() do tmp
    B = SpaceAGORACalibration.Backend
    K = SpaceAGORACalibration.Spec

    p_int = ParameterSpec(name="i", kind=K.integer, lower=1.0, upper=9.0, manifest_targets=["x"], env_targets=["CAL_INT_ENV"])
    p_cat = ParameterSpec(name="mode", kind=K.categorical, choices=["on", "off"], manifest_targets=["x"])
    spec = CalibrationSpec(
        id="backend_helper_cov",
        seed=77,
        parameters=[p_int, p_cat],
        budgets=BudgetSpec(
            initial_samples=1,
            global_iters=0,
            batch_size=1,
            local_refine_topk=1,
            local_refine_steps=0,
            robustness_samples=1,
            robustness_uncertainty="normal",
            uncertainty_atm_scale=0.05,
            uncertainty_ic_scale=0.01,
            objective_telemetry_noise_frac=0.1
        )
    )

    @assert B._replica_id(Candidate(id=1, values=Dict{String, Any}("robustness_replica" => 2))) == 2
    @assert B._replica_id(Candidate(id=1, values=Dict{String, Any}("robustness_replica" => "oops"))) === nothing

    @assert B._lookup_parameter(spec, "i").name == "i"
    @assert try
        B._lookup_parameter(spec, "missing")
        false
    catch e
        e isa ArgumentError
    end

    @assert B._coerce_parameter_value(p_int, 4.0) == 4
    @assert try
        B._coerce_parameter_value(p_int, 4.8)
        false
    catch e
        e isa InexactError
    end
    @assert B._coerce_parameter_value(p_cat, "on") == "on"
    @assert try
        B._coerce_parameter_value(p_cat, "invalid_choice")
        false
    catch e
        e isa ArgumentError
    end

    @assert B._parse_segment("node[*]") == ("node", (:all, nothing))
    @assert B._parse_segment("node[2]") == ("node", (:index, 2))
    @assert B._parse_segment("node[name=test]") == ("node", (:field_eq, ("name", "test")))
    @assert try
        B._parse_segment("node[bad")
        false
    catch e
        e isa ArgumentError
    end
    @assert try
        B._parse_segment("node[0]")
        false
    catch e
        e isa ArgumentError
    end

    arr = Any[Dict("name" => "a"), Dict("name" => "b")]
    @assert length(B._select_items(arr, (:all, nothing), "seg[*]")) == 2
    @assert B._select_items(arr, (:index, 1), "seg[1]")[1]["name"] == "a"
    @assert length(B._select_items(arr, (:field_eq, ("name", "b")), "seg[name=b]")) == 1
    @assert try
        B._select_items(arr, (:index, 4), "seg[4]")
        false
    catch e
        e isa ArgumentError
    end
    @assert try
        B._select_items(arr, (:unsupported, nothing), "seg")
        false
    catch e
        e isa ArgumentError
    end

    @assert B._apply_transform(1.0, 3.0, "set") == 3.0
    @assert B._apply_transform([1.0, 2.0], 2.0, "add") == [3.0, 4.0]
    @assert B._apply_transform([1.0, 2.0], 3.0, "mul") == [3.0, 6.0]
    @assert try
        B._apply_transform(1.0, 2.0, "bad")
        false
    catch e
        e isa ArgumentError
    end

    fail_eval = evaluate_candidate(MockBackend(noise_sigma=0.0, fail_rate=1.0), spec, Candidate(id=2, values=Dict{String, Any}("i" => 2, "mode" => "on")); stage="global_search_quick", run_dir=tmp)
    @assert !fail_eval.success
    @assert !isfinite(fail_eval.score)

    @assert isnan(B._safe_float(missing))
    @assert B._safe_float("3.25") == 3.25
    @assert isnan(B._safe_float("not_a_number"))
    @assert B._safe_bool(1) == true
    @assert B._safe_bool("off"; default=true) == false
    @assert B._safe_bool("not_bool"; default=true) == true

    rows_runtime_secondary = [(simulation_runtime_s=12.0,)]
    @assert B._runtime_from_rows(rows_runtime_secondary, 99.0) == 12.0
    @assert B._runtime_from_rows(Any[], 5.5) == 5.5

    rows_fallback = [(nmae="2.0", pass="false"), (nmae=1.0, pass=true)]
    obj_fallback = B._objective_from_rows(rows_fallback, spec; run_failed=false, runtime_s=10.0)
    @assert isfinite(obj_fallback.score)
    @assert obj_fallback.failed_rows == 1.0
    @assert obj_fallback.all_pass == false

    rows_robust = [(scenario="unknown", rmse_km=2.0, max_abs_km=1.0, limit_nmae=1.0, nmae=1.5, limit_max_abs_km=2.0, pass=false)]
    obj_robust = B._objective_from_rows(rows_robust, spec; run_failed=false, runtime_s=10.0, noise_rng=MersenneTwister(12))
    @assert isfinite(obj_robust.score)
    @assert obj_robust.failed_rows == 1.0

    # Display-normalized columns should take precedence when present.
    rows_display = [(
        scenario="unknown",
        rmse_km=10.0,
        max_abs_km=10.0,
        rmse_display=1.0,
        max_abs_display=1.0,
        limit_max_rmse_km=1.0,
        limit_max_abs_km=1.0,
        limit_max_rmse_display=1.0,
        limit_max_abs_display=1.0,
        limit_nmae=1.0,
        nmae=1.0,
        pass=true
    )]
    obj_display = B._objective_from_rows(rows_display, spec; run_failed=false, runtime_s=10.0)
    @assert isapprox(obj_display.score, 0.75; atol=1e-12, rtol=0.0)

    scenario = Dict{String, Any}(
        "name" => "s1",
        "ra_m" => 7000e3,
        "rp_altitude_m" => 120e3,
        "i_deg" => 30.0,
        "aop_deg" => 355.0,
        "raan_deg" => 359.0,
        "ta_deg" => 358.0,
        "atmosphere_truth" => Dict("gram_seed" => 1234, "gram_perturbation_scales" => [0.2, 0.0, 0.0, 0.5])
    )
    B._apply_uncertainty_to_scenario!(scenario, MersenneTwister(1), spec)
    @assert scenario["ra_m"] > 0.0
    @assert 0.0 <= scenario["i_deg"] <= 180.0
    @assert 0.0 <= scenario["aop_deg"] < 360.0
    @assert 0.0 <= scenario["raan_deg"] < 360.0
    @assert 0.0 <= scenario["ta_deg"] < 360.0

    manifest_doc = Dict{String, Any}(
        "scenarios" => Any[
            scenario,
            "skip_non_dict"
        ]
    )
    B._apply_uncertainty_to_manifest!(manifest_doc, spec, 7, 2)
    @assert manifest_doc["scenarios"][1] isa AbstractDict

    @assert B._format_env_value(true) == "1"
    @assert B._format_env_value(4) == "4"
    @assert occursin("2.5", B._format_env_value(2.5))
    @assert B._env_with_default("SPACEAGORA_CAL_TEST_UNSET", "x") == "x"
    @assert length(B._candidate_runtime_policy_env_pairs()) >= 5

    env_pairs = B._candidate_env_pairs(spec, Candidate(id=3, values=Dict{String, Any}("i" => 4, "mode" => "off")))
    env_map = Dict(env_pairs)
    @assert env_map["CAL_INT_ENV"] == "4"
    @assert haskey(env_map, "SPACEAGORA_INNER_THREAD_BUDGET")

    missing_csv_path = joinpath(tmp, "does_not_exist.csv")
    @assert isempty(B._read_summary_rows(missing_csv_path))

    manifest_path = joinpath(tmp, "manifest.toml")
    open(manifest_path, "w") do io
        TOML.print(io, Dict("value" => 1.0))
    end
    failing_spec = CalibrationSpec(
        id="backend_command_failure_cov",
        verification_script="missing_script.jl",
        manifest_paths=[manifest_path],
        parameters=[ParameterSpec(name="i", kind=K.integer, lower=1.0, upper=1.0, transform="set", manifest_targets=["value"])],
        budgets=BudgetSpec(initial_samples=1, global_iters=0, batch_size=1, local_refine_topk=1, local_refine_steps=0, robustness_samples=1)
    )
    cmd_backend_fail = CommandBackend(
        julia_cmd=Base.julia_cmd(),
        project_path=joinpath(@__DIR__, ".."),
        verification_script=joinpath(tmp, "missing_verify.jl"),
        manifest_path=manifest_path,
        profile="quick",
        enforce=false
    )
    fail_command_eval = evaluate_candidate(
        cmd_backend_fail,
        failing_spec,
        Candidate(id=9, values=Dict{String, Any}("i" => 1));
        stage="global_search_quick",
        run_dir=tmp
    )
    @assert !fail_command_eval.success
    @assert !isempty(fail_command_eval.error_message)
end

println("SpaceAGORACalibration smoke tests passed")
