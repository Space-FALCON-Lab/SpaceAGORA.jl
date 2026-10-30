"""Results root shared by the periapsis border study scripts.

Precedence: --results-root (add_argument below), then env SPACEAGORA_BORDER_RESULTS, then the
`results_root` key of study_config.toml (default output/periapsis_border_switchability under the repo root).
"""
import os, tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))


def add_results_root_argument(parser):
    parser.add_argument("--results-root", default=None,
                        help="results root (default: $SPACEAGORA_BORDER_RESULTS, else results_root in study_config.toml)")


def results_root(cli_value=None):
    if cli_value:
        return os.path.abspath(cli_value)
    env = os.environ.get("SPACEAGORA_BORDER_RESULTS", "")
    if env:
        return os.path.abspath(env)
    with open(os.path.join(HERE, "study_config.toml"), "rb") as f:
        return os.path.join(REPO, tomllib.load(f)["results_root"])
