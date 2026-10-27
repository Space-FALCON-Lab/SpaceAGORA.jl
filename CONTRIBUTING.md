# Contributing to SpaceAGORA

To change SpaceAGORA itself, start with the
[contributor guide and checklist](docs/contributing/index.md). The checklist
asks you to name the source owners and the API boundary a change touches, keep
configuration and numerical conventions, say who owns mutable state and random
streams between runs, and keep the evidence for the change with it.

The pages it points to most often:

- [Architecture and source ownership](docs/contributing/architecture.md)
- [Configuration ownership](docs/contributing/configuration_ownership.md)
- [Documentation policy and release checks](docs/contributing/documentation_policy.md)
- [Architecture and quality contracts](docs/contributing/contracts.md)
- [Optional adapter conformance and release evidence](docs/contributing/adapter_maintenance.md)
- [GNC and planning ownership](docs/contributing/gnc_ownership.md)

If you are adding a model in a project of your own rather than changing
SpaceAGORA, use the [extension guide](docs/src/extensibility.md) and the
[API compatibility reference](docs/src/public_api_policy.md). The published
documentation site covers using SpaceAGORA; contributor procedures stay here in
the repository.
