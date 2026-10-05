This version is just a cleaned version of ver3. Ver3 is already a fully functioning version.

Run scripts from their original repository locations, not VS Code attachment
snapshots. Historical scripts under `historic_codes/` load the shared helpers
from `../functions/`; deeper post-processing scripts use the corresponding
relative path anchored at `@__DIR__`. The helper folder is not duplicated inside
the archive.

From the outer integration workspace, validate Julia syntax and static include
paths across all project copies with:

```sh
julia --startup-file=no test/include_paths.jl
```