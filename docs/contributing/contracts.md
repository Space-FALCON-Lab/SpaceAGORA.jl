# Architecture and Quality Contracts

These are the maintainers' contracts: the rules for repository structure,
public API naming, verification gates and runtime analysis. Read them before
you move code ownership, change the topology, decide where behaviour belongs,
or change a CI gate or public name. To run SpaceAGORA, start from the
[User Guide](../src/index.md) instead.

The contracts are repository documents. The published user site does not
render them, and the generated API pages cover only the public names; the
links below go to the files in this repository. For documentation setup and
checks, use the [local build procedure](documentation_policy.md#build-the-docs-locally).

See also the [Contributor Guide](index.md), the
[Documentation Policy](documentation_policy.md) and the
[API Policy](../src/public_api_policy.md).

## Architecture

- [Canonical topology contract](../architecture/canonical_topology_contract.md): the file-placement and ownership model for the repository.
- [GNC aerobraking boundary contract](../architecture/gnc_aerobraking_boundary_contract.md): typed EDG ownership and the typed-command separation rules.
- [Canonical owner audit](../architecture/src_canonical_owner_audit.md): where runnable examples, plotting and report entry points, shared helpers and runtime locks live.
- [Completeness contract](../architecture/src_completeness_contract.md): the source-file categories and the paths forbidden in canonical code.

## Quality

- [API naming contract](../quality/api_naming_contract.md): naming rules and file-level ownership for stable surfaces.
- [Verification contract](../quality/verification_contract.md): required tests, release gates and coverage policy.
- [Performance runtime analysis contract](../quality/performance_runtime_analysis.md): the supported runtime-analysis study and its outputs.
- [Geometry contract](../quality/geometry_contract.md): who owns the shared geometry helpers and which caller-specific policies stay with their callers.
- [Quaternion contract](../quality/quaternion_contract.md): the canonical, robotics and compliant-multibody quaternion helpers and how they differ.
- [GRAM loading contract](../quality/gram_loading_contract.md): the GRAM consumers' loading policies and the discovery they share.
- [Legacy plotting contracts](../quality/legacy_plotting_contracts.md): the saved-data contracts of the retained legacy plotting code.
