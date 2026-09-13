# Julia Modules: Loading, Imports, Exports, and Aliases

These mechanisms answer different questions: where code is evaluated, which
names are available, and which names a module exposes to its users.

## Reference Table

| Mechanism | Meaning | Important Detail |
| --- | --- | --- |
| `module SomeModule ... end` | Defines a module: a namespace containing functions, types, constants, and other modules. | A module is not necessarily one file; it can include several files. |
| `include("file.jl")` | Evaluates the file's code in the including module. | Executes top-level statements. Function bodies run only when called. Repeating an include evaluates the file again. |
| `using SomeModule` | Loads/accesses the module and makes its module name and exported names available in the current scope. | Does not expose every internal name or automatically run the functions. |
| `import SomeModule` | Loads/accesses the module and brings only its module name into scope. | Access members using qualified names such as `SomeModule.helper(...)`. |
| `using SomeModule: helper` | Brings only the explicitly selected name into scope. | The name need not be exported. Does not by itself bring the name `SomeModule` into scope. |
| `import SomeModule: helper` | Brings the selected name into scope and, for a function, permits adding methods using its unqualified name. | For merely calling `helper`, selective `using` also works. |
| `export helper` | Marks a name as part of the names made available by `using` the current module. | Does not define, import, or load `helper`. The module must separately establish its binding. |
| `@reexport using SomeModule` | Imports and re-exports another module's exported API through the current module. | Provided by the external `Reexport.jl` package, not a Julia built-in keyword. Does not recursively expose every internal name. |
| `const SM = SimulationModel` | Creates a constant alias referring to the same module object. | Does not copy the module, reload its code, or start a simulation. |
| `SpaceAGORA.SimulationModel` | Accesses the `SimulationModel` binding inside `SpaceAGORA`. | This is a qualified Julia name, not a filesystem path. |
| `using .SimulationModel` | Accesses `SimulationModel` relative to the current module and makes its exports available. | The module binding must already be available, for example through an alias. The dot is not a folder reference. |
| `isdefined(@__MODULE__, :SM)` | Checks whether the current module has a defined binding named `SM`. | `:SM` is a Symbol representing the name; `@__MODULE__` is the current module. |
| `!isdefined(@__MODULE__, :SM)` | Returns true when that binding is not defined. | Leading `!` is logical NOT; a trailing `!` in a function name is only a naming convention. |

## Loading Is Not the Same as Name Visibility

`using SpaceAGORA` loads the package and the modules its source includes. An
example does not need to load each internal module separately. However, loading
the package does not make every internal name available without a prefix.

After these aliases are established:

```julia
using SpaceAGORA
const SimulationModel = SpaceAGORA.SimulationModel
const SM = SimulationModel
```

`SM`, `SimulationModel`, and `SpaceAGORA.SimulationModel` refer to the same module.
For example, `SM.InitialTime` and `SimulationModel.InitialTime` access the same type.

Adding `using .SimulationModel` makes its exported names available directly, so
an example can write `InitialTime(...)` without the module prefix.

## Export Does Not Create a Binding

This module declares an export but never defines the name:

```julia
module IncompleteAPI
    export helper
end
```

The export alone does not make `IncompleteAPI.helper` a defined function. A module
must define the function itself or import an existing binding before exposing it.

Here is an import-and-export chain:

```julia
module ExamplePackage
    module ModelTools
        export helper
        helper(value) = 2 * value
    end

    using .ModelTools: helper
    export helper
end

using .ExamplePackage
helper(3)
```

The final call returns `6`. `ExamplePackage.helper` and
`ExamplePackage.ModelTools.helper` refer to the same function.

With `Reexport.jl` installed, `using Reexport` followed by
`@reexport using .ModelTools` can replace the explicit selective import and export
when the intention is to expose that module's exported API.

## Method Extension and Re-Exports

Re-exporting a function preserves its identity and method table. It does not
create a copy or change which module originally defined it.

- `import SomeModule: helper` permits an unqualified extension such as
  `helper(value::MyType) = ...`.
- With a module name available, a qualified definition such as
  `SomeModule.helper(value::MyType) = ...` can extend that function.
- The qualifier can be a re-exporting module if its binding refers to the same
  function. It does not have to name the original defining module.
- Ordinary calls do not need extension permission. Only add methods when that
  extension is intentional, generally for a type you own.

An earlier chat explanation incorrectly claimed that extending through a
re-exporting module creates a different method slot. That claim is incorrect.

A separately defined forwarding wrapper is different from a re-export:

```julia
run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)
```

When this defines a new function in the current module, it is a wrapper that calls
the engine function, not an alias to the engine function.

## Why Examples Have Shared Setup

The shared example setup selects the project environment, creates convenient
aliases, and imports helpers commonly used by examples. Individual examples can
import additional packages or names for their particular needs.

```julia
import SpaceAGORA.TelemetryVerification:
    make_example_config, make_three_body_spacecraft, run_and_report
```

This imports three helpers; it does not call them. Some examples use all three,
even though the ORACLE runner does not. That does not make the whole telemetry
module unnecessary: ORACLE separately uses `rvtoorbitalelement` from it.

A file in the examples folder is not executed automatically. An example explicitly
includes the shared setup. Library source modules manage their own imports and
should not depend on an example changing the caller's project environment.

## Definition Order Matters

After activating the intended project, the following order ensures each name is
available before it is used:

```julia
using SpaceAGORA

if !isdefined(@__MODULE__, :SimulationModel)
    const SimulationModel = SpaceAGORA.SimulationModel
end

using .SimulationModel

if !isdefined(@__MODULE__, :SM)
    const SM = SimulationModel
end
```

The guards avoid defining an already-defined name; they do not verify that an
existing binding points to the expected module. Validate loading in a fresh Julia
process, because an existing REPL can hide missing imports or ordering errors.