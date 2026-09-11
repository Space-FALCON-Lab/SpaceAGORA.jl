# Why Do You Need Constructors?

Constructors let you control **how** an object gets built — enforcing rules, filling in defaults, or computing derived values — instead of leaving every call site to build the struct "by hand" and possibly get it wrong.

## Example where it gets tricky without one

```julia
struct Percentage
    value::Float64
end
```

Nothing stops code from doing `Percentage(150.0)` or `Percentage(-20.0)` — invalid values silently sneak in, and the bug only surfaces later when some unrelated calculation produces nonsense. Every single place in your codebase that creates a `Percentage` would need to remember to check `0 <= value <= 100` itself — easy to forget, and easy to get inconsistent across files.

With a custom (inner) constructor, you enforce the rule **once**, in one place:

```julia
struct Percentage
    value::Float64
    function Percentage(value)
        0 <= value <= 100 || error("value must be between 0 and 100")
        new(value)
    end
end
```

Now `Percentage(150.0)` throws immediately at creation time, no matter where in the code it's called from.

Constructors exist so validation, defaults, and derived-field logic live in one place instead of being duplicated (and potentially forgotten) everywhere a struct is created.
