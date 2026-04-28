## Summary

- scope:
- key files:
- risks:

## Docs Drift Checklist

- [ ] If root exports changed, I updated `docs/public_api_symbols.jl`
- [ ] If assets changed, I updated `data/assets_manifest.toml` and `docs/src/assets.md`
- [ ] If new extension points were added, I updated `docs/src/extensibility.md`
- [ ] If new CLI commands or flags were added, I updated `docs/src/cli.md`

## Release Docs Gates

- [ ] `julia --project=docs docs/make.jl`
- [ ] `julia --project=. test/ci_public_api_surface_gate.jl`
- [ ] `julia --project=. test/ci_hpc_extensibility_docs_gate.jl`
