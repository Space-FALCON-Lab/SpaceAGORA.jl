# Julia Hello World: Application Environment vs Package Project

A Julia project does not have to be a package. The same hello-world behavior can
be organized as a script inside an application environment or as a reusable package.

## 1. Application/Environment Project

This project provides an environment for running a script. It does not define
a package of its own.

### Folder Layout

```text
HelloApp/
    Project.toml
    hello.jl
```

### Project.toml

```toml
[deps]
```

The dependency section is empty because printing text needs no external packages.
No package name, UUID, version, or author fields are required for this environment.

### hello.jl

```julia
println("Hello, world!")
```

### Run It

From inside the `HelloApp` folder:

```bash
julia --project=. hello.jl
```

`--project=.` selects the current folder's environment; `hello.jl` selects the
script to execute. Activating an environment alone does not execute the script.

Output:

```text
Hello, world!
```

There is no `module HelloApp` or `using HelloApp`. Dependencies can be added later
without turning the application itself into a package. For this dependency-free
script, running `julia hello.jl` without a project environment also works.

## 2. Package Project

This project defines a named package with a function that callers can reuse.

### Folder Layout

```text
HelloWorld/
    Project.toml
    src/
        HelloWorld.jl
```

### Project.toml

```toml
name = "HelloWorld"
uuid = "e7a35852-5297-4ed1-9fd7-d79e1f213927"
version = "0.1.0"
```

`name` declares the package name, and `uuid` identifies the package uniquely.
`version` records its version; `authors` is optional descriptive metadata.
For a different new package, generate a fresh UUID, for example with
`using UUIDs; uuid4()` in Julia.

### src/HelloWorld.jl

```julia
module HelloWorld
    export greet

    greet() = println("Hello, world!")
end
```

This is the conventional source entry point matching the package name.

- `module HelloWorld ... end` defines the module, which is a namespace.
- `greet() = ...` defines a function; its body runs when called.
- `export greet` allows callers using `using HelloWorld` to call `greet()` directly.

### Run It

From inside the `HelloWorld` folder:

```bash
julia --project=. -e 'using HelloWorld; greet()'
```

Output:

```text
Hello, world!
```

The active project identifies this package, so Julia can locate its entry point.
`using HelloWorld` loads the package; `greet()` prints the message. Loading alone
does not call `greet()`. No manual `include` is needed in this command.

To use the package from a different project, add the local package to that
environment with `Pkg.develop(path="/path/to/HelloWorld")` after `import Pkg`.
The package does not need to be published or registered online.

## What Makes the Difference?

| Question | Application/Environment Project | Package Project |
| --- | --- | --- |
| Main purpose | Run scripts or an application in a chosen environment. | Provide code loadable as a named package. |
| Defines a package of its own? | No, not in this example. | Yes: `HelloWorld`. |
| Package name and UUID? | Not required. | Identify the package in project metadata. |
| Matching source entry point? | Not required. | Conventionally `src/HelloWorld.jl`. |
| Module declaration? | Not required for the script. | Defines the package's main module. |
| How is hello world executed? | Run the script. | Load the package, then call `greet()`. |
| Can it have dependencies? | Yes. | Yes. |

The five lines containing `module HelloWorld ... end` alone define a module,
not a complete package project. Likewise, metadata alone does not supply the
module's implementation. The package example combines both with a loadable layout.

A module can also be defined directly in a REPL or through `include`, without
being a separate package. `using .HelloWorld` accesses an already-defined local
module; it does not search arbitrary folders for its source file.

One project file can identify one package of its own or no package at all, while
listing many dependency packages. One repository can contain multiple such projects.

These are organizational choices, not a measure of program complexity: a package
can print only hello world, and a large application can use a script-based environment.

For the related loading and namespace mechanisms, see
[Julia Modules: Loading, Imports, Exports, and Aliases](julia_module_mechanisms.md).