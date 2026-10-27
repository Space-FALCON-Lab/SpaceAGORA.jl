"""Location of the managed RPO example environment, separate from source checkouts."""
module RPOExampleEnvironment
using SHA
const ROOT = realpath(joinpath(@__DIR__, "..", ".."))
function directory(; root=ROOT)
    override = get(ENV, "SPACEAGORA_RPO_ENV", "")
    isempty(override) || return abspath(override)
    identity = bytes2hex(sha256(realpath(root)))[1:12]
    return joinpath(first(DEPOT_PATH), "environments", "spaceagora-rpo-" * identity)
end
end
