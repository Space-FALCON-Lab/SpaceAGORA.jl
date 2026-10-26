## Summary

- scope:
- key files:
- risks:

## Docs Drift Checklist

- [ ] If root exports changed, I updated `docs/public_api_symbols.jl`
- [ ] If assets changed, I updated `data/assets_manifest.toml` and `docs/src/assets.md`
- [ ] If new extension points were added, I updated `docs/src/extensibility.md`
- [ ] If a hook was added, renamed or deprecated, I updated `docs/extension_hooks.jl`
- [ ] If new CLI commands or flags were added, I updated `docs/src/cli.md`

## Release Docs Gates

- [ ] `julia --startup-file=no docs/setup.jl`, then `julia --project=docs docs/make.jl`
- [ ] `julia --project=. test/gates/ci_public_api_surface_gate.jl`
- [ ] `julia --project=. test/gates/ci_hpc_extensibility_docs_gate.jl`
