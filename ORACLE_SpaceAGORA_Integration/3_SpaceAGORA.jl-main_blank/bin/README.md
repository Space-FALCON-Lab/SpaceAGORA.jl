# Bin folder

This folder contains the project launcher.

- The script [spaceagora] finds the repo root and runs Julia with the correct project environment.
- It then starts the CLI entrypoint at [src/cli/main.jl].
- So this folder is just a convenience shortcut for running the app from the terminal.

## Example

### With the launcher

```bash
./bin/spaceagora --help
```

### Without the launcher

```bash
cd /path/to/repo
julia --project=. src/cli/main.jl --help
```

The bin folder saves you from manually typing the repo path and CLI entrypoint each time.
