# Independent review of the current-upstream optional KS adapter port

Read-only source review against upstream `a66436549`, checkpoint-base `615a447c0`, and old qualified prototype `9c6e4b8`. The initial inspection encountered the expected in-progress Project/solver-setting cherry-pick conflicts; the final captured port is clean at `2ccc99f4f693124e2cd19b589f90e766c2940d62`. The JSON captures exact inspected file hashes. Implementation, branches and historical evidence were not changed by this reviewer. No numerical run, timing experiment or scientific requalification was performed.

The core trial boundary and current configuration resolution order are retained correctly. At the initial inspected port, two guard gaps required correction before treating the new source as qualified: newly reachable onboard data models bypass KS admission, and central-gravity counting ignores an already-resolved scope mask. The latter is a latent resolved-table contract issue; normal public one-spacecraft scope resolution does not create an empty central scope.

## Actionable guard findings

### P1: newly mounted onboard data bypasses the KS model contract

`src/orchestration/simulation/engine/ks_solver.jl:52–54` rejects mounted actuator/power/sensor fields but omits `Link.data_stores` and `Link.data_sources`, which now exist at `src/models/vehicle/spacecraft/model.jl:308–309`. These fields are absent from the old prototype's Link layout.

This is an actual full-RHS extension path, not inert metadata. `dynamics_rhs.jl:3697–3699` always appends power/data/output work after the translational dispatch. `dynamics_rhs.jl:3726–3728` takes the storage route whenever a store exists, including when `interlink_model === nothing`. That route calls `_data_flows` (`dynamics_rhs.jl:3757–3762`), which evaluates mounted production models through `DataStorage._production_bps`. The public custom `data_production_bps` hook is evaluated for every source at `src/models/vehicle/data/data_storage.jl:74–75`. Its type is an arbitrary `AbstractDataSource`, not one of the dynamic/control effectors checked by the KS purity loop at `ks_solver.jl:63–64`.

A one-link gravity/NoAtmosphere configuration can therefore pass the current KS guard while introducing unqualified data auxiliary equations, storage-boundary callbacks, and a custom source whose implementation might consume query history. The generic data API documents a repeatability requirement, but the optional KS route currently neither checks that owner nor excludes this new capability. A full-state coupled solver's ability to carry an auxiliary component is not evidence that every model producing it is safe.

Required initial-scope repair: reject nonempty stores and sources before optional extension loading/solve, leaving ordinary solver data support intact. Expanding support is possible later, but needs explicit source purity, full/empty-store boundary, same-time perturbed trial, rejection, checkpoint and restart evidence. Do not infer it from the older mass/heat lifecycle qualification.

Focused tests should show (1) an otherwise supported configuration still passes; (2) a built-in store/source is explicitly refused; (3) an unknown custom source, with and without a positive `trial_rhs_safe` declaration elsewhere, cannot bypass the mounted-owner guard; and (4) ordinary DP8/Tsit5 data behavior is unchanged. A counting/history-consuming source is an appropriate negative contract fixture, not a new supported physics model.

### P2: central model must be applicable to spacecraft 1

`ks_solver.jl:60–62` counts inverse-square/J2 types in the global tuple, whereas the actual RHS skips entries whose scope mask is false (`dynamics_rhs.jl:355–356`, `383–389`, and `439–440`; the common predicate is `engine/effector_scopes.jl:13–15`).

Normal public behavior is sound: `run_simulation` resolves dynamic scopes before KS validation (`execution.jl:813–816`, then `860`). `_resolve_effector_scopes` rejects a selection resolving to no spacecraft (`effector_scopes.jl:127–129`) and collapses a full one-spacecraft selection into unscoped bare effectors (`130–141`). Thus a valid `scoped(central, BoundSpacecraft((1,)))` need not be rejected, and this audit does not claim the normal resolver creates empty scopes.

However, an already-resolved `DynamicsModel(..., EffectorScopes)` is accepted after size checks at `effector_scopes.jl:85–94`. `_check_scoped_run` tests `any(gravity) && !all(gravity)` (`171–173`), which permits the all-false case for a one-spacecraft run. A stored/manual internal table with a false mask for its sole central model therefore passes the current global count while the physical RHS omits that central force. The generic engine permits gravity-free models; the KS route has a stricter advertised central-model contract.

Required repair: count the central models actually applicable to spacecraft 1 using the same resolved-scope predicate as the RHS and require exactly one. Check resolved-table consistency, or explicitly reject an empty/nonapplicable selected central. A manually malformed inactive table is a useful negative guard fixture. Preserve ordinary unscoped and valid scoped-to-1 cases. Also retain tests for no central, two applicable centrals, empty/out-of-range public selections, and unsupported coupled/multi-spacecraft configurations. This is a contract/admission finding, not a demonstrated numerical failure in the normal scoped public path.

## Correct ownership already present; preserve it during reconciliation

- **Effective configuration precedes validation.** `execution.jl:535–556` folds per-spacecraft GNC into configuration-level tuples using bound copies and carries the resolved effector-scope table. `run_simulation` calls it before scope resolution and before KS preflight. Models declared on spacecraft cannot simply be checked by inspecting the original caller's `args.control_model`; validate the effective bound model objects. Current order is correct.
- **Shallow configuration rebuilding preserves new fields.** `_with_fields` and `_with_configuration` (`config/simulation_configuration.jl:41–81`) avoid old positional reconstruction dropping scopes/new settings. Keep those current-upstream helpers rather than reinstating an old prototype container layout.
- **Typed/environment solver overrides have one owner.** `_run_solver_config` in `solver_policy.jl` reads explicit `args.solver_config` first and otherwise uses the current environment adapter. `_engine_env_get` and its scoped-override dictionary (`engine/adapters/from_env.jl:238–272`, `314–364`) must remain authoritative. Do not add a direct ENV-only KS selector that bypasses a typed `SimulationEngineConfig`.
- **Current ordinary solver modes remain intact.** Add `:ks` without dropping `:multirate_mri`, its SplitODEProblem route, or ordinary data/power support. KS must still use the ordinary full-state ODEParams route, not a split or externally owned shadow route. Conflict resolution is not permission to restore the old mode list wholesale.
- **Per-spacecraft scopes remain lifecycle state.** GNC member lists are resolved at model startup, checked against later changes and recorded in checkpoints (`engine/gnc_run_scopes.jl:1–31`, `68–77`). Dynamic scope checkpoint identities are separately checked (`effector_scopes.jl:207–240`). Retain that ownership; do not flatten away checkpoint scope identities to make KS work.

## Trial RHS and cache requirements

The copied core boundary is present at `dynamics_rhs.jl:4345–4364`. It invalidates both `density_sample_t` and `in_atmosphere_sample_t` before the normal full RHS. That is necessary because both reusable atmosphere and entry-interface state caches can otherwise identify samples by time alone. The actual density owner is `_density_model_for_sat` (`callbacks/density_callbacks/model_selection.jl:1–17`), which first consults the runtime per-spacecraft model vector.

The initial allowed analytic configurations normally leave that runtime vector empty: `setup.jl:1847–1889` populates it only for GRAM/native/surrogate instances. Therefore this review does not establish an ordinary public custom-density bypass through the vector. Nevertheless solve-start validation should inspect effective runtime instances or enforce the initial no-overrides invariant, rather than treating the configured fallback type as proof of the runtime model. Keep the snapshotted `CallbackEnvConfig` refusals for density history, freeze-per-step, vacuum/track caches; inspect actual runtime settings after typed/environment overrides resolve. Do not silently disable configured frozen semantics.

Current thermal work recomputes its per-link heat-rate buffer (`thermal_callbacks.jl:3–11`, `14–39`) and consumes the effective atmosphere through the normal sampling service. Time-only planet rotations can remain reusable. The current initial force publication also checks the exact saved state (`callbacks/save_fields.jl:66–81`), so a same-epoch perturbed trial cannot overwrite an initialization row. Keep this current-upstream check.

The optional extension correctly re-evaluates the start and accepted candidate via the core boundary (`ext/SpaceAGORAKSPropagatorExt.jl:112–118`, `140–145`) and retains dense interval endpoint states/rates independently from FSAL scratch (`62–104`). Those invariants must survive source reconciliation. Rejected/nonlinear trial evaluations can overwrite derived buffers, but must not advance physical/controller/random history; accepted-candidate refresh alone cannot undo history consumption by a custom model.

Required focused checks on the new port: same epoch/different position and velocity; repeated and permuted trial order; current mass and heat derivatives; accepted refresh after a rejected candidate; active phase/cache settings through typed overrides; correct initial force output; root-shortened dense intervals before affect!; command/impulse boundaries; cold fresh-process restart. These are correctness gates, not timing work.

## Solver and callback contracts

The current upstream callback container helper may use tuples or erased vectors while preserving order (`solver_policy.jl` `_solver_callback_set` and `_callbacks_for_cached_integrator`). Preserve current construction/ordering and the checkpoint PR's table cadence fix. Do not restore an old tuple-only or segment-local output mechanism.

The imported algorithm-type check for cached explicit solvers is important (`solver_policy.jl:503–505` in the inspected port): otherwise a previous solver's integrator can satisfy callback/save-option checks while carrying a different numerical algorithm. KS additionally resets its cached integrator before initializing its per-run mu/diagnostics (`solver_policy.jl:848–865`). Test ordinary→KS→ordinary and changed mu/config runs with the same user cache; do not silently reuse another run's algorithm or physical central parameter.

The initial trait `trial_rhs_safe(::Any)=false` and explicit known-model methods (`ks_solver.jl:3–15`) are a repeatability contract, not scientific validation of arbitrary custom models. A user's deliberate custom opt-in is not itself a hidden type-allowlist bug; every actual RHS extension owner must still be checked or excluded. Accepted-step guidance/navigation and ideal burn scheduling may legitimately update held commands outside trial evaluation. Rejecting all such lifecycle models would narrow away the originally qualified controls/checkpoint scope and is not the requested repair.

SRP facets alone are geometry; the unsupported SRP dynamic force is already rejected unless explicitly admitted by a custom trait. Do not confuse passive facet presence with an active unguarded SRP force. Similarly, keep the existing attitude, mounted hardware, robot/attachment, multi-spacecraft, interlink and external-propagator restrictions without using them to suppress ordinary solver capabilities.

## Verdict and next gate

The port's architecture is consistent with the intended ownership split: SpaceAGORA owns physical Cartesian state, configured physics, runtime/lifecycle and trial sampling; KSPropagator owns the numerical map. It is **not yet approved as a requalified new-source adapter**. Correct the data-owner and effective-central guard gaps, preserve the current-upstream contracts above, then run focused normal-source guard/parity/callback regressions before the frozen scientific and lifecycle replays. Existing old prototype success cannot establish current-upstream correctness automatically. This review makes no performance or full scientific validity claim.


## Follow-up source review: guard corrections at0037984e0

Reviewed `0037984e0ab6ed0aa5fd2b82d3e4ca0193e7304f` read-only after the parent
implemented the two corrections. Relative to the initial captured port2ccc99f4,
the commit changes only `engine/ks_solver.jl`: mounted admission now excludes
both data containers (`52–55`), and central counting uses the same resolved-scope
predicate as the physical RHS (`61–69`). Both changes are approved at source
level. Ordinary solver data behavior, valid unscoped/scoped-to-1 support and the
previously qualified mass/heat/control lifecycle are preserved by this design.
No further mandatory source guard/design change was found for the stated initial
contract. Arbitrarily malformed internal scope-table objects are not promoted to
a supported public input contract by this approval.

The parent independently retained clean before reproductions: the data probe
passes0 and fails4 assertions at the unguarded source; the scope probe passes4 and fails2
because the inactive sole central is admitted. Logs are
`../logs/adapter-data-guard-before-attempt1.log` and
`../logs/adapter-scope-guard-before-attempt3.log`. Earlier scope harness-error
attempts remain preserved and are not substituted for those clean failures.
Probe sources are `../adapter-verification/upstream-data-guard.jl` and
`upstream-scope-guard.jl`.

The after probes and new-source force/heat/dense-history/original-guard checks
were still running when this follow-up source verdict was written. Their success
is **not** inferred from the small diff or from prior prototype evidence. The
parent must record their actual final-source outcomes and perform the separate
scientific/lifecycle qualification gates. No timing work was done here.
