# Minimal SVG writer shared by the documentation diagrams in diagrams.jl.
# Standard library only; the palette matches docs/src/assets/architecture_run_flow.svg.
module DocsSvg

export Svg, rect!, text!, rich!, line!, path!, circle!, save_svg,
       SANS, MONO, INK, SOFT, MUTED, ACCENT, FILL, EDGE, PALE, PALE_EDGE, WARM, WARM_INK

const SANS = "DejaVu Sans, Helvetica, Arial, sans-serif"
const MONO = "DejaVu Sans Mono, Menlo, monospace"
const INK, SOFT, MUTED, ACCENT = "#1f2937", "#4b5563", "#6b7280", "#2f5f7f"
const FILL, EDGE, PALE, PALE_EDGE = "#f2f7fc", "#8cbde7", "#fafafa", "#9aa5b1"
const WARM, WARM_INK = "#b45309", "#7c2d12"

struct Svg
    width::Int
    height::Int
    title::String
    desc::String
    parts::Vector{String}
end
Svg(width, height, title, desc) = Svg(width, height, title, desc, String[])

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
num(x) = isinteger(x) ? string(Int(x)) : string(round(Float64(x); digits=2))

function rect!(s::Svg, x, y, w, h; fill=FILL, stroke=EDGE, dash=nothing, rx=0, sw=1)
    d = dash === nothing ? "" : " stroke-dasharray=\"$dash\""
    push!(s.parts, "<rect x=\"$(num(x))\" y=\"$(num(y))\" width=\"$(num(w))\" height=\"$(num(h))\" " *
                   "rx=\"$(num(rx))\" fill=\"$fill\" stroke=\"$stroke\" stroke-width=\"$(num(sw))\"$d/>")
end

function text!(s::Svg, x, y, str; size=11, color=SOFT, bold=false, mono=false, anchor="start", italic=false)
    attrs = ["x=\"$(num(x))\"", "y=\"$(num(y))\"", "font-size=\"$(num(size))\"", "fill=\"$color\""]
    bold && push!(attrs, "font-weight=\"bold\"")
    italic && push!(attrs, "font-style=\"italic\"")
    mono && push!(attrs, "font-family=\"$MONO\"")
    anchor == "start" || push!(attrs, "text-anchor=\"$anchor\"")
    push!(s.parts, "<text $(join(attrs, ' '))>$(esc(str))</text>")
end

# spans: (text, kind) with kind one of :plain, :mono, :bold, :monobold.
function rich!(s::Svg, x, y, spans; size=11, color=SOFT, anchor="start")
    tag = Dict(:plain => "", :mono => " font-family=\"$MONO\"", :bold => " font-weight=\"bold\"",
               :monobold => " font-family=\"$MONO\" font-weight=\"bold\"")
    body = join("<tspan$(tag[kind])>$(esc(t))</tspan>" for (t, kind) in spans)
    a = anchor == "start" ? "" : " text-anchor=\"$anchor\""
    push!(s.parts, "<text x=\"$(num(x))\" y=\"$(num(y))\" font-size=\"$(num(size))\" fill=\"$color\"$a>$body</text>")
end

function line!(s::Svg, x1, y1, x2, y2; color=ACCENT, sw=1.5, dash=nothing, arrow=false)
    d = dash === nothing ? "" : " stroke-dasharray=\"$dash\""
    m = arrow ? " marker-end=\"url(#arrow)\"" : ""
    push!(s.parts, "<line x1=\"$(num(x1))\" y1=\"$(num(y1))\" x2=\"$(num(x2))\" y2=\"$(num(y2))\" " *
                   "stroke=\"$color\" stroke-width=\"$(num(sw))\"$d$m/>")
end

function path!(s::Svg, d; color=ACCENT, sw=1.5, fill="none", dash=nothing)
    da = dash === nothing ? "" : " stroke-dasharray=\"$dash\""
    push!(s.parts, "<path d=\"$d\" stroke=\"$color\" stroke-width=\"$(num(sw))\" fill=\"$fill\"$da/>")
end

circle!(s::Svg, cx, cy, r; fill=ACCENT, stroke="none", sw=1) =
    push!(s.parts, "<circle cx=\"$(num(cx))\" cy=\"$(num(cy))\" r=\"$(num(r))\" fill=\"$fill\" " *
                   "stroke=\"$stroke\" stroke-width=\"$(num(sw))\"/>")

function save_svg(s::Svg, path::AbstractString)
    head = """<?xml version="1.0" encoding="utf-8"?>
    <svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="title desc" width="$(s.width)" height="$(s.height)" viewBox="0 0 $(s.width) $(s.height)" font-family="$SANS">
    <title id="title">$(esc(s.title))</title><desc id="desc">$(esc(s.desc))</desc>
    <defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="$ACCENT"/></marker></defs>
    <rect width="$(s.width)" height="$(s.height)" fill="#ffffff"/>
    """
    open(path, "w") do io
        print(io, head, join(s.parts, "\n"), "\n</svg>\n")
    end
    return path
end

end # module DocsSvg
