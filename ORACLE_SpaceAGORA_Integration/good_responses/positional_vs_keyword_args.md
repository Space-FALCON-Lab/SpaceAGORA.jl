# Positional vs. Keyword Arguments in Julia

Whether a parameter must be called positionally or by keyword depends on which side of a `;` it's on in the *definition* — unlike Python, Julia won't accept the "wrong" style.

## Positional only (no semicolon)

function f(a, b)
    return a + b
end

f(1, 2)        # ✅ works
f(a=1, b=2)    # ❌ MethodError

## Keyword only (semicolon right after the opening paren)

function g(; a, b)
    return a + b
end

g(a=1, b=2)    # ✅ works
g(1, 2)        # ❌ MethodError

`@kwdef` auto-generates this exact "keyword only" method for a struct, based on its fields:

Base.@kwdef struct g
    a::Int
    b::Int = 2
end

# equivalent to writing:
# function g(; a, b=2)
#     g(a, b)   # forwards to the plain positional constructor
# end

g(a=1, b=2)   # ✅ works
g(1, 2)       # ✅ also works — @kwdef keeps the positional constructor too

## Mixing both

function h(a, b; c, d)
    return a + b + c + d
end

h(1, 2, c=3, d=4)   # ✅ works — a,b positional; c,d keyword

**Rule:** everything before `;` is positional-only, everything after is keyword-only.

