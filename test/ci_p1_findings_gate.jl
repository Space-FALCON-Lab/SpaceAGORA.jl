const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const ALLOWLIST_PATH = joinpath(REPO_ROOT, "test", "p1_findings_allowlist.txt")

using Dates

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
    r"(?i)\bpriority\s*=\s*1\b",
    r"(?i)\bseverity\s*[:=]\s*P?1\b",
]

struct P1Finding
    file::String
    line::Int
    snippet::String
    key::String
end

struct AllowlistEntry
    key::String
    owner::Union{Nothing, String}
    opened::Union{Nothing, Date}
    expires::Union{Nothing, Date}
    line::Int
end

function parse_allowlist_entry(line::String, line_no::Int)
    parts = split(line, "#"; limit=2)
    key = strip(parts[1])
    key_parts = split(key, "::"; limit=2)
    if length(key_parts) != 2 || isempty(strip(key_parts[1])) || isempty(strip(key_parts[2]))
        error("Invalid allowlist entry format at $ALLOWLIST_PATH:$line_no. Expected `path::exact_line`.")
    end

    owner = nothing
    opened = nothing
    expires = nothing

    if length(parts) == 2
        metadata = strip(parts[2])
        for token in split(metadata)
            occursin("=", token) || continue
            k, v = split(token, "="; limit=2)
            if k == "owner"
                owner = v
            elseif k == "opened"
                parsed = tryparse(Date, v)
                parsed === nothing && error("Invalid opened date at $ALLOWLIST_PATH:$line_no: $v (expected YYYY-MM-DD)")
                opened = parsed
            elseif k == "expires"
                parsed = tryparse(Date, v)
                parsed === nothing && error("Invalid expires date at $ALLOWLIST_PATH:$line_no: $v (expected YYYY-MM-DD)")
                expires = parsed
            else
                error("Unknown allowlist metadata key at $ALLOWLIST_PATH:$line_no: $k")
            end
        end
    end

    return AllowlistEntry(key, owner, opened, expires, line_no)
end

function load_allowlist(path::String)
    entries = Dict{String, AllowlistEntry}()
    if !isfile(path)
        return entries
    end

    for (line_no, raw) in enumerate(eachline(path))
        line = strip(raw)
        if isempty(line) || startswith(line, "#")
            continue
        end
        entry = parse_allowlist_entry(line, line_no)
        if haskey(entries, entry.key)
            error("Duplicate allowlist entry at $ALLOWLIST_PATH:$line_no for key: $(entry.key)")
        end
        entries[entry.key] = entry
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

    allowlist_keys = Set(keys(allowlist))
    new_findings = [f for f in findings if !(f.key in allowlist_keys)]
    stale_allowlist = sort([k for k in allowlist_keys if all(f -> f.key != k, findings)])
    expired_allowlist = sort([
        v for v in values(allowlist) if v.expires !== nothing && v.expires < today()
    ]; by=v -> v.line)

    println("p1_scan_files=$(length(files)) findings=$(length(findings)) allowlisted=$(length(allowlist))")

    if !isempty(stale_allowlist)
        println("stale_allowlist_entries=$(length(stale_allowlist))")
        for entry in stale_allowlist
            println("  - $entry")
        end
    end

    if !isempty(expired_allowlist)
        println("expired_allowlist_entries=$(length(expired_allowlist))")
        for entry in expired_allowlist
            println("  - $(entry.key) (line $(entry.line), expires=$(entry.expires))")
        end
    end

    if !isempty(new_findings)
        println("new_p1_findings=$(length(new_findings))")
        for f in new_findings
            println("  - $(f.file):$(f.line): $(f.snippet)")
        end
        error("P1 findings gate failed: unallowlisted P1 markers found.")
    end

    if !isempty(expired_allowlist)
        error("P1 findings gate failed: expired allowlist entries found.")
    end

    println("p1_findings_gate_ok")
end

main()
