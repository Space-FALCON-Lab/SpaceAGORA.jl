const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

tracked = readchomp(`sh -lc "cd '$REPO_ROOT' && git ls-files '*.cov' '.DS_Store'"`)
!isempty(strip(tracked)) && error("Tracked artifact files must be removed:\n$tracked")

println("no_artifact_files_gate_ok")
