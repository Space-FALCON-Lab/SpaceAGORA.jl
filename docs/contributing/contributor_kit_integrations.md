# Contributor kit integrations

Some work cannot be done from a contributor repository: a new GNC algorithm
family, or a new service facade for a companion package. Both are trusted,
core-owned code and need a SpaceAGORA pull request; this checklist covers that
core half. A feedback controller or estimator uses the public model hooks and
needs neither. That work, and any work that uses an existing family or facade,
stays in the contributor's repository, built from a starter; see
[Building Your Own Research Repository](../src/user/research_repository.md).

Start each item from `main`, and name in the pull request the starter or
companion that needs it.

## A new algorithm family

- [ ] The family module, its request and product types and every family
      obligation, as [Adding a GNC algorithm family](gnc_algorithm_families.md)
      lists them.
- [ ] Its public names re-exported from the root, a section for them in
      `docs/public_api_symbols.jl`, and a page of their own: an entry in
      `API_REFERENCE_PAGES` and the generated page under "API Reference" in the
      `pages` of `docs/make.jl`.
- [ ] Rows in `docs/extension_hooks.jl` for the functions an algorithm or the
      family implements. `test/gates/ci_extension_hook_reference_gate.jl` fails
      when an exported function in a section listed in `HOOK_SECTIONS` has
      neither a row nor an entry in `NOT_HOOKS`.
- [ ] A `<family>_conformance_requests()` builder that includes a request the
      family refuses, the family's `conformance_link`, and a
      `<family>_conformance_simulation()` when the family runs in a simulation.
- [ ] The family in `family_instances()` and its algorithms in
      `algorithm_instances()` (`test/fixtures/conformance_instances.jl`), so
      that core runs `check_algorithm_family` and `check_algorithm` on them.
- [ ] A row in the families table of
      [Planning Algorithms](../src/user/gnc_planning.md#the-planning-profile).
- [ ] If the package starter should show an algorithm of the family, its
      example under `templates/starters/package/`. Core CI generates each
      starter and runs its workflow.

## A new service boundary

- [ ] The facade, following the
      [optional adapter standard](../src/reference/adapter_contract.md) and the
      [extension conventions](../src/reference/extension_conventions.md), with
      contract version 1.0.0 as its first: `CONTRACT_VERSION`,
      `CompatibilityError`, a `require_version` that accepts only the exact
      version, `IMPLEMENTED_FUNCTIONS` and `CONSUMED_SERVICES`.
- [ ] A `public` declaration in `src/SpaceAGORA.jl`, an entry in
      `docs/public_api_symbols.jl` and a row in `VERSIONED_FACADES`
      (`docs/extension_hooks.jl`).
- [ ] A reference page under `docs/src/reference/`, listed under "Extend
      SpaceAGORA" in `docs/make.jl`, and rows in the supported-boundaries table
      and the contract history of `docs/src/reference/adapter_contract.md`.
- [ ] A weak-dependency extension under `ext/` if core calls into the companion.
      Loading SpaceAGORA alone must not load it.
- [ ] The facade in `test/unit/adapter_contract_tests.jl`, which checks its
      shape with `check_service_facade`.
- [ ] A companion starter pointed at the new facade
      (`templates/starters/companion/` uses `ExternalPropagationServices` as
      its example).

## When core CI qualifies the pair

A companion that core CI installs and tests also needs:

- [ ] a gitlink under `packages/<domain>/` at a published companion commit;
- [ ] an entry in `packages/OptionalPackages.toml` (name, path, URL, UUID,
      version);
- [ ] its name in `NAMES` and its domain in `DOMAINS` in
      `scripts/optional_packages.jl`;
- [ ] an `optional-<name>-installation` job in `.github/workflows/julia-ci.yml`,
      added to the `needs` of the `tests` job, which fails when any of those
      jobs fails.

## Always

- [ ] A review record, `test/ai_reviews/PR_<number>.md`, when the pull request
      changes a `.jl` file under `src/`. `test/gates/ci_ai_review_artifact_gate.jl`
      requires it and checks that it names every changed source file.
- [ ] For a facade or family a companion serves: qualify the pair as
      [optional adapter maintenance](adapter_maintenance.md#cross-repository-qualification)
      describes; change the facade's contract version for any change to a
      signature, unit, frame, result field, ownership or failure behavior; and
      move the pins forward only after both repositories have qualified the
      pair.
- [ ] No study results, private data or performance numbers in the core
      repository or the starters.
