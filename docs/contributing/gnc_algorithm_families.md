# Adding a GNC algorithm family

This page is for SpaceAGORA maintainers. It lists what a new family of the
general GNC algorithm contract needs in the core repository. A user who writes
an algorithm *for* an existing family follows
[Writing a GNC Algorithm](../src/user/gnc_algorithm.md) instead; a family is
trusted, core-owned code, so it never lives in a student's package.

## What a family is

A family is the trusted definition of one planning problem: the request an
algorithm receives, the product it proposes, how a proposal is validated, where
an accepted product is installed, when an installed product stops being current,
and what happens on rejection. The request driver (`run_algorithm_request!`)
runs the same steps for every family; the family supplies the content of each
step. The generic layer owns only the order of the steps, the copy rules, the
stream, the records and the error type. Everything else is family-specific:
the request's content and its constructor checks, the product's representation,
the validation's metrics and tolerances, the currency rule, the failure policy
and the run-report fields.

The families today:

| Family | Module and file | Request, product | On rejection | Binding |
|---|---|---|---|---|
| `RPOReferenceFamily` | `RPOPlannerInterfaces`, `src/models/gnc/planning/rpo/contract/rpo_planner.jl` | `RPOPlanningRequest`, `RPOReference` | `:fail`, or `:install_fallback` when the binding configures it | `PlannerGuidance`, the RPO pilot |
| `BurnPlanFamily` | `BurnPlanning`, `src/models/gnc/control/maneuvers/burn_planning.jl` | `BurnPlanRequest`, `PropulsiveBurnPlan` | `:keep_previous` | `AlgorithmGuidance` |
| `JointReferenceFamily` | `JointReferencePlanning`, `src/models/gnc/planning/robot_arm/joint_reference_family.jl` | `JointReferenceRequest`, `RobotArmPlan` | `:fail` | `AlgorithmGuidance`, with the arm controller's `RobotArmPlanBuffer` as its link |
| `EDGGuidanceFamily` | `EDGGuidanceFamilies`, `src/models/gnc/guidance/aerobraking/edg_guidance_family.jl` | `EDGGuidanceRequest`, `EDGGuidanceDecision` | `:fail` | `AlgorithmGuidance`, with `EDGGuidanceSettings`; its algorithm comes from the SpaceAGORAEDG companion through `EDGServices`, whose contract version is the family's |

The burn-plan family sits under `control/` because it writes the burn
scheduler's command buffer, which the typed-command boundary gate
(`test/gates/ci_gnc_typed_command_boundary_gate.jl`) allows only there.

## The family's obligations

A family is a subtype of `AbstractAlgorithmFamily` with methods of these
functions (all in `GNCAlgorithms`, `src/models/gnc/contracts/gnc_algorithms.jl`):

| Function | Required | Notes |
|---|---|---|
| `validate_proposal(family, request, result, trusted)` | yes | against the lifecycle's own copy of the request; never installs, never draws from a stream, never reads diagnostics as geometry |
| `install!(link, family, product, request, t; validation, diagnostics)` | yes | the mutated link first; atomic, with no yield or callback |
| `is_current(family, request, product, t)` | yes | whether the consumer may still use the product, including any preview it needs |
| `on_rejection(family)` | default `:fail` | `:keep_previous` or `:install_fallback` |
| `fallback_product(family, request)` | with `:install_fallback` | such as a hold |
| `validate_capabilities(family, algorithm, request; checkpoint)` | default accepts | for a family with capability rules |
| `default_capabilities(family)` | default empty | the family's capability type, supporting nothing |
| `algorithm_result_type(family)` | default `AbstractAlgorithmResult` | `AlgorithmResult` unless the family needs its own |
| `is_revision_request(family, request)` | default `false` | requests sent to `propose_revision!` |
| `proposed_product(result)` | for a result type of its own | must return, not throw, for a result validation refuses |
| `installed_record`, `rejected_record` | defaults | the run-report fields |
| `family_contract_version(family)` | for a family a companion serves | so the companion can require the version it was built for |

A family that runs through the generic binding, `AlgorithmGuidance`, also
implements `family_request(family, settings, u, p, t, spacecraft_index, reason,
request_id)` (the owned request, built by trusted code), `family_link` (the
install target, when the family finds its own in the run) and
`preflight_family` (the configuration check before the run, such as a missing
consumer). A family with a binding of its own implements the binding functions
instead: `algorithm_runtime`, `bound_algorithm`, `bound_family`, `bound_link`,
`build_algorithm_request`, `trusted_queries` and `fail_algorithm_request!`.

Give the family an abstract algorithm type that declares `algorithm_family`
(`AbstractBurnPlanner` does), so that its algorithms inherit the declaration.

## The requests

- **Owned snapshots.** The request is an immutable copy in SI units with a
  positive `request_id`, the time and its validity. It never holds live
  configuration, the integrator or a model.
- **Constructor checks.** Refuse a malformed request where it is built, not in
  validation.
- **Stale observations** are refused, not planned from.

## Steps for the pull request

1. **The module.** Add the family next to its consumer, with its request and
   product types, the obligations above, docstrings on every public name, and
   at least one algorithm of the family. Include it in `src/SpaceAGORA.jl` and
   re-export its public names from the root.
2. **The public surface.** Add a section for the family to
   `docs/public_api_symbols.jl` and its own `API_REFERENCE_PAGES` entry in
   `docs/make.jl` (as "Burn-plan family" has), because `generated/public_api.md`
   is close to Documenter's size limit. Docstrings `@ref` only names their own
   module can see.
3. **The hook registry.** Add rows to `docs/extension_hooks.jl` for each
   function an algorithm or the family implements that is not already a row,
   with `page` naming the handwritten user guide that explains the hook.
   Add the section title to `HOOK_SECTIONS` when it holds hooks.
4. **The conformance kit.** Add the family's representative requests (a
   `<family>_conformance_requests()` builder, including one request the family
   refuses), its `conformance_link`, and when it runs in a simulation a
   `<family>_conformance_simulation()`, under `src/interfaces/analysis/conformance/`.
   Add the family to `family_instances()` and its algorithms to
   `algorithm_instances()` in `test/fixtures/conformance_instances.jl`, so the
   core runs `check_algorithm_family` and `check_algorithm` on them.
5. **Tests.** Contract tests in `test/unit/gnc/` for validation (every refusal
   reason), installation, currency and the rejection policy, and a run through
   the binding with no family-specific orchestration.
6. **The docs.** A row in the families table of
   [Planning Algorithms](../src/user/gnc_planning.md#The-planning-profile),
   and how its algorithms run.
7. **Review.** A review record under `test/ai_reviews/` (the CI artifact gate).
   For a family a companion serves, qualify the pair as
   [Optional adapter maintenance](adapter_maintenance.md#Cross-repository-qualification)
   describes, and advance the pins only after both repositories have qualified it.

[Contributor kit integrations](contributor_kit_integrations.md) has the
checklist for this and for a new service facade.
