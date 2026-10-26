# Contributor kit integrations

This checklist is for the core half of work that a contributor repository
cannot do on its own: a new GNC algorithm family, or a new service boundary for
a companion package. Both are trusted, core-owned code and need a SpaceAGORA
pull request. A feedback controller or estimator uses the public model hooks
and needs neither a new family nor a core pull request. That work, and work
that uses an existing family or facade, stays in the
contributor's repository, built from a starter; see
[Building Your Own Research Repository](../src/user/research_repository.md).

Start each item from `origin/main`, and record the starter or companion that
needs it in the pull request.

## A new algorithm family

- [ ] The family module, its request and product types and every family
      obligation, as [Adding a GNC algorithm family](gnc_algorithm_families.md)
      lists them.
- [ ] Public names re-exported from the root, a section in
      `docs/public_api_symbols.jl` and its own `API_REFERENCE_PAGES` entry in
      `docs/make.jl`.
- [ ] Rows in `docs/extension_hooks.jl` for the functions an algorithm or the
      family implements (`ci_extension_hook_reference_gate.jl` fails on a
      missing row).
- [ ] A `<family>_conformance_requests()` builder that includes a request the
      family refuses, the family's `conformance_link`, and a
      `<family>_conformance_simulation()` when it runs in a simulation.
- [ ] The family in `family_instances()` and its algorithms in
      `algorithm_instances()` (`test/fixtures/conformance_instances.jl`), so the
      core runs `check_algorithm_family` and `check_algorithm` on them.
- [ ] A row in the families table of
      [Planning Algorithms](../src/user/gnc_planning.md#The-planning-profile).
- [ ] If the package starter should show an algorithm of the family, its
      example under `templates/starters/package/` (core CI generates and tests
      each starter).

## A new service boundary

- [ ] The facade, following the
      [Optional Adapter Standard](../src/reference/adapter_contract.md) and the
      [Extension Conventions](../src/reference/extension_conventions.md) from
      contract 1.0.0: `CONTRACT_VERSION`, an exact `require_version`,
      `IMPLEMENTED_FUNCTIONS` and `CONSUMED_SERVICES`.
- [ ] A `public` declaration, an entry in `docs/public_api_symbols.jl` and a row
      in `VERSIONED_FACADES` (`docs/extension_hooks.jl`).
- [ ] A reference page under `docs/src/reference/`, linked from "Extend
      SpaceAGORA" in `docs/make.jl`, and a row in the supported-boundaries table
      and the contract history of `adapter_contract.md`.
- [ ] A weak-dependency extension under `ext/` if the core calls into the
      companion; loading SpaceAGORA alone must not load it.
- [ ] The facade in `test/unit/adapter_contract_tests.jl`, which checks its shape
      with `check_service_facade`.
- [ ] A companion starter pointed at the new facade (`templates/starters/companion/`
      serves `ExternalPropagationServices` as its example).

## When core CI qualifies the pair

A companion that core CI installs and tests also needs:

- [ ] a gitlink under `packages/<domain>/` at a published companion commit;
- [ ] an entry in `packages/OptionalPackages.toml` (name, path, URL, UUID, version);
- [ ] the name in `NAMES` and its domain in `DOMAINS` of
      `scripts/optional_packages.jl`;
- [ ] an `optional-<name>-installation` job in `.github/workflows/julia-ci.yml`.

## Always

- [ ] A review record under `test/ai_reviews/` (the CI artifact gate) when the
      pull request changes `src/`.
- [ ] For a facade or family a companion serves: qualify the pair as
      [Optional adapter maintenance](adapter_maintenance.md#Cross-repository-qualification)
      describes, change the facade's contract version for any change of a
      signature, unit, frame, result field, ownership or failure semantics, and
      advance pins only after both repositories have qualified the pair.
- [ ] No study results, private data or performance numbers in the core
      repository or the starters.
