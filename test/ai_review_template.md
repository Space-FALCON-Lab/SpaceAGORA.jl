# AI Review Template

Use this template for PR review artifacts saved as `test/ai_reviews/PR_<number>.md`.

## Scope
- Describe the intent of the PR and what was reviewed.

## Changed Files
- List every changed `src/**/*.jl` file in this PR.
- Example: `src/control/Propulsive_maneuvers.jl`

## Findings
- List findings by severity (`P1`, `P2`, `P3`) and file/line context.
- If none, state `No findings`.

## P1 Assessment
- State whether any unresolved `P1` findings remain.
- Required final line: `Unresolved P1: Yes` or `Unresolved P1: No`.

## Tests Added/Updated
- List tests added or modified for this PR.
- Include rationale for why these tests cover the risk.

## Residual Risk
- Describe any remaining risk that is accepted and why.
- Include follow-up actions if needed.
