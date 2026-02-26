const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const ALLOWLIST_PATH = joinpath(REPO_ROOT, "test", "p1_findings_allowlist.txt")

const SCAN_ROOTS = [
    joinpath(REPO_ROOT, "src"),
    joinpath(REPO_ROOT, "test"),
]

const EXTRA_FILES = [
    joinpath(REPO_ROOT, "README.md"),
]

const ALLOWED_EXTENSIONS = Set([
    ".jl",
    ".md",
    ".txt",
    ".yml",
    ".yaml",
    ".toml",
])

const EXCLUDED_RELATIVE_FILES = Set([
    "test/ci_p1_findings_gate.jl",
    "test/p1_findings_allowlist.txt",
])

const P1_PATTERNS = [
    r"\[P1\]",
    r"\bP1\b",
    r"(?i)\bpriority\s*=\s*1\b",
    r"(?i)\bseverity\s*[:=]\s*P?1\b",
]

struct P1Finding
    file::String
    line::Int
    snippet::String
    key::String
end

function load_allowlist(path::String)
    entries = Set{String}()
    if !isfile(path)
        return entries
    end

    for raw in eachline(path)
        line = strip(raw)
        if isempty(line) || startswith(line, "#")
            continue
        end
        push!(entries, line)
    end

    return entries
end

function iter_scan_files()
    files = String[]

    for root in SCAN_ROOTS
        isdir(root) || continue
        for (dir, _, names) in walkdir(root)
            for name in names
                ext = splitext(name)[2]
                ext in ALLOWED_EXTENSIONS || continue
                file = joinpath(dir, name)
                rel = relpath(file, REPO_ROOT)
                rel in EXCLUDED_RELATIVE_FILES && continue
                push!(files, file)
            end
        end
    end

    for file in EXTRA_FILES
        if isfile(file)
            rel = relpath(file, REPO_ROOT)
            rel in EXCLUDED_RELATIVE_FILES || push!(files, file)
        end
    end

    return sort(unique(files))
end

function collect_findings(files::Vector{String})
    findings = P1Finding[]
    for file in files
        rel = relpath(file, REPO_ROOT)
        for (idx, raw) in enumerate(eachline(file))
            line = strip(raw)
            isempty(line) && continue
            any(p -> occursin(p, line), P1_PATTERNS) || continue

            key = string(rel, "::", line)
            push!(findings, P1Finding(rel, idx, line, key))
        end
    end
    return findings
end

function main()
    allowlist = load_allowlist(ALLOWLIST_PATH)
    files = iter_scan_files()
    findings = collect_findings(files)

    new_findings = [f for f in findings if !(f.key in allowlist)]
    stale_allowlist = sort([k for k in allowlist if all(f -> f.key != k, findings)])

    println("p1_scan_files=$(length(files)) findings=$(length(findings)) allowlisted=$(length(allowlist))")

    if !isempty(stale_allowlist)
        println("stale_allowlist_entries=$(length(stale_allowlist))")
        for entry in stale_allowlist
            println("  - $entry")
        end
    end

    if !isempty(new_findings)
        println("new_p1_findings=$(length(new_findings))")
        for f in new_findings
            println("  - $(f.file):$(f.line): $(f.snippet)")
        end
        error("P1 findings gate failed: unallowlisted P1 markers found.")
    end

    println("p1_findings_gate_ok")
end

main()
