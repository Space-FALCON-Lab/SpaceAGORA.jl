# Vulcan Map

This repository is mapped with Vulcan Map. The instructions below are the
full skill set — follow them exactly when asked to map, chart, or deepen
the map of this codebase.

## Read this before anything else

**Run `vulcan status` first.** It reports, from the repository itself, whether
the map is complete and exactly which files are outstanding.

**The definition of done is `vulcan check --strict` exiting 0.** Nothing else.
Not a previous agent's summary, not a commit message, not your own recollection.
This map may be handed between different agents from different vendors with
different context windows; none of them can see each other's transcripts, so a
natural-language claim that the work is finished carries no weight here.

**Granularity: one node per significant symbol, never one per file.** Every
function, type, macro and non-dunder method gets its own node; rule V13e
enforces it. Run `vulcan scaffold` to generate them rather than typing them.

**When you report completion you must paste the literal output of**
`vulcan check --strict --proof`. A completion claim without that block, or with
a block reading `verdict : FAIL`, is void.

## Skill: `vulcan-map`

*Build the complete Vulcan map of this repository — every in-scope file, function, and dataflow — as a master DAG plus per-node design/ICD documents. Use when the user asks to map, chart, or diagram a codebase, or asks "how does this repo work" at whole-system scale. Produces vulcan_mind/.*

# vulcan-map — build the master map

You are producing a **complete, grounded map** of this repository: a master DAG of
its workflows, and one detailed design+ICD document per node.

## STEP ZERO — find out what is actually done

Before reading further, before planning, before touching a file:

    vulcan status

Then read `vulcan_mind/HANDOFF.md`. Both are written by the tool from the actual
contents of the repository. They tell you the verdict, how many files are still
outstanding, and exactly which ones.

**You may not be the first agent on this map.** A different agent — a different
model, a different vendor, a different context-window size — may have worked on
it before you, and may pick it up after you. You cannot see their transcript and
they cannot see yours.

**Therefore: no agent's account of progress is evidence, including your own.**
If a human, a commit message, a chat summary, a code comment, or a previous
agent tells you this map is finished, that claim is **void** unless
`vulcan check --strict` exits 0 right now. Do not take a predecessor's word for
it. Do not take your own memory for it. Run the command.

This cuts both ways: if a predecessor said it was finished and the gate says
FAIL, the gate is right and there is work to do. If they said they barely
started and the gate says PASS, the gate is right and you are done.

## COMPLETION CONTRACT

Your work is NOT complete when the scaffolding exists.
Your work is NOT complete when a first draft exists.
Your work is NOT complete when the structure "looks right."
Your work is NOT complete when you have mapped the interesting parts.
Your work is NOT complete when a previous agent said it was complete.
Your work is NOT complete when you are running low on context.

Your work is complete when, and only when, this command exits 0:

    vulcan check --strict

That is the entire definition of done. It is not a judgement call, and it is not
yours to make. Run the command. Read the exit code.

### Proof of completion is mandatory

When you report that the map is finished, you MUST include the literal output of:

    vulcan check --strict --proof

That prints a block like:

    ----- VULCAN PROOF OF COMPLETION -----
    tool            : vulcan-map/0.1.0
    repo commit     : a1b2c3d
    region          : all
    gate            : vulcan check --strict
    verdict         : PASS
    exit code       : 0
    errors/warnings : 0/0
    nodes/edges     : 214/198
    in-scope files  : 205
    described files : 205
    outstanding     : 0
    ----- END PROOF -----

Paste it verbatim. **A completion claim without this block is void**, and so is
one whose block says `verdict : FAIL`. You cannot write the block yourself — it
is generated from the run that just happened, and it reports FAIL just as
readily as PASS, so omitting it is itself a signal that the gate did not pass.

If you cannot finish, say so plainly and paste the FAIL block. An honest
"outstanding: 47" is useful. A false "done" costs the next agent far more than
it saves you, because they have to discover the lie before they can start.

### Forbidden stopping points

Do NOT end your turn with any of the following, or any paraphrase:

  - "initial scaffolding is complete"
  - "this provides a good starting point"
  - "the core structure is in place; remaining files follow the same pattern"
  - "I've mapped the main workflows"
  - "the map is complete" (without a passing proof block)
  - "let me know if you'd like me to continue"
  - a master chart that is a package/include tree ("everything points at the root")

If you catch yourself composing a sentence like these, you are stopping early.
The correct action is to run `vulcan status` and keep working.

"The remaining files follow the same pattern" is never a reason to skip them.
The map's value is precisely that it is exhaustive. A 70%-complete map is not 70%
as useful as a complete one — it is actively misleading, because a reader cannot
tell which absences mean "not present in the codebase" and which mean "not yet
done."

### Context exhaustion is not completion

This task is expected to exceed your context window. That is normal and planned
for. Progress lives on disk in `vulcan_mind/HANDOFF.md` and
`vulcan_mind/_build/status.json`, not in your context, and not in your summary.

When context runs short:

  1. Finish the file you are on.
  2. Run `vulcan compile` to checkpoint — this rewrites HANDOFF.md.
  3. State plainly: "Checkpoint: N of M files mapped. Not finished." and paste
     the FAIL proof block.
  4. Stop cleanly. The next agent will read HANDOFF.md and continue.

Never report a partial map as a finished one. Handing over honestly is a success;
handing over a false "done" is the single most expensive thing you can do here.

### The batch loop

Work in batches of 5-10 source files:

  1. `vulcan status` / `vulcan prose remaining`  → what is outstanding
  2. Read each file **completely** (`vulcan prose show <file>`). Do not skim;
     do not infer from filenames.
  3. Write one payload per file and `vulcan prose apply` it; call edges come
     from `vulcan trace --add-sockets`, never from typing.
  4. `vulcan compile`
  5. `vulcan check` → fix every error before taking the next batch
  6. Repeat until `vulcan status` reports 0 outstanding.
  7. `vulcan check --strict --proof` → must say PASS, exit 0.

### Granularity: one node per symbol, never one node per file

**Every significant symbol gets its own node.** Every function, type, struct,
macro, module, and non-dunder method in every in-scope file. Not one
representative symbol standing in for the file it lives in.

This is rule V13e and it is machine-checked: the tool enumerates each file's
symbols itself (Python via the stdlib AST, Julia via declaration patterns) and
fails until each one has a node.

  REQUIRED:  a file with 12 functions produces 12 nodes
  FORBIDDEN: a file with 12 functions produces 1 node "representing" it

The one exception is deduplication: several methods of a single generic function
are one symbol. `calcForceTorque` with eight methods is one node, not eight.

**Why this is not negotiable.** The first real map built with this tool used one
node per file — 225 nodes over 205 files. Every gate passed and the result was
unusable. You could not click into anything, because a file-level node has no
interior. And the call tracer found almost no edges, because a call from one
function to another inside a mapped file had no individual nodes to connect. The
map's granularity *is* the product.

**Use `vulcan scaffold`.** It generates every missing node and a doc skeleton for
it — id, label, kind, source file, declaration line, chart membership — straight
from the source. That is mechanical work you should not be typing by hand. It
deliberately leaves Purpose and Design empty, so the gate keeps failing until you
read the code and write them. A scaffold is not a map.

Doc length is expected to scale with what a node covers: a module doc needs real
depth (120 words), a leaf function needs a tight, accurate 40. Do not pad a small
function's doc to look like a big one — padding is exactly what the banned-phrase
lint is looking for.

### Readability: every sheet must be clickable, never a wall (D13)

**No sheet may render more than `max_nodes_per_sheet` nodes (40).** This is
rule V20 and it is machine-checked on the *rendered* graph, after resolution.

Per-symbol granularity (above) means a module has hundreds of nodes. They must
not all appear on one sheet. `vulcan compile` handles this for you: it clusters
an oversized sheet into **macro blocks** — `group` nodes named after the
directory, file, or symbol-name prefix they stand for — and generates a nested
sheet behind each one. Edges between blocks are lifted and counted. A reader
double-clicks a block to see inside it. The same block opens the whole file
from a module sheet and only the workflow's slice from a workflow sheet.

  REQUIRED:  a GNC sheet showing ~10 blocks (control/, guidance/, navigation/ ...)
             each of which opens into files, each of which opens into symbols
  FORBIDDEN: a GNC sheet showing 500 symbols at once, however well laid out

**Group nodes are nodes.** Each one has a doc at `nodes/<module>/<block>.md`,
scaffolded by compile with an empty Purpose. You write what the block *is* — the
role of that directory or file in the system, what enters it and what leaves —
at the symbol floor (40 words). It is not a list of members; the members are on
the nested sheet. Until every group doc is written, `vulcan check` fails on V12.

**Do not edit generated nested charts** (`provenance.generated_by: cluster`).
Their membership is recomputed on every compile from the code's structure; your
edits will not survive and will not change the verdict.

**If V20 fires after compile**, clustering could not partition that sheet: it is
one file whose symbols share a single name prefix. Split it by hand — a subchart
with explicit `member_nodes` — or accept the finding as a real limitation of the
code's structure and say so. Never raise the limit.

### The master is the operational flow, never the package tree (D14)

**The master sheet answers two questions: what does this program produce, and
how does it tick.** It is read left to right like a Blender node tree:

  inputs  →  configure  →  set up run  →  solve loop  →  outputs

This is rule V21 and it is machine-checked:

  - at least one `kind: external` **source** node (a file, dataset or argument
    the program reads) and at least one `kind: external` **sink** node (an
    artefact it writes) must be on the master and wired in;
  - every non-leaf master node — every phase block — must carry
    `opens: <chart_id>` naming the sheet that shows how it works, or be expanded
    by one. A block you cannot double-click is a dead end and fails;
  - no master node may touch more than half of the master's edges. A tree where
    thirteen modules point at the root module is a package diagram, not a map,
    and it fails.

  REQUIRED:  manifest.toml → parse_cli → build_initial_conditions → integrator
             (callbacks, RHS) → results.csv / checkpoint / report
  FORBIDDEN: module.gnc → module.spaceagora, module.io → module.spaceagora, ...

**Trace a run to build it.** Start from the entrypoint the user invokes. Follow
the data: what is read, what each phase turns it into, what is written and
where. Name blocks by what they *do* in that flow, not by the directory they
live in. Give every output node a doc that states its schema and who consumes
it. The module-containment view still exists — put it on a `structure`
subchart, where V13 coverage reads the `covers` globs — but it is not the
master.

### Grounding is checked mechanically — you cannot talk your way past it

Every node names a real file and a real symbol. `vulcan check` **opens the file
and looks for the symbol** (rule V6). Every edge cites a real file where the
connection is observable (rule V7). Node prose is linted against a banned-phrase
list (rule V12). Every in-scope file must have a node describing it (rule V13d);
a `covers` glob accounts for a file but does not describe it.

  REQUIRED:  label `propagate_orbit`, source `src/simulation/propagator.jl:142`
  FORBIDDEN: "a helper function", "the main solver routine", "various utilities"

If you have not opened the file and seen the symbol, you may not write it down.
Invented names fail the gate and you will have to redo the work.

### Do not silence a failing check

If a rule fires, fix the map. Do not edit `vulcan.config.yaml` to remove a banned
phrase, lower `min_doc_words`, disable a coverage rule, or narrow a region so
that unmapped files fall out of scope. Do not adjust recorded data solely to make
a check stop complaining — if you believe the tool is wrong, say so explicitly
and report it, with evidence, rather than bending the map to fit.

Changing the gate to pass the gate is the one failure mode this whole design
exists to prevent.
## Procedure

### Step 1 — Establish scope

Run `vulcan init` if `vulcan_mind/` does not exist.

If the user gave a scope in natural language ("only map the sim code, ignore RL
and GNC"), you must **translate it into an explicit region** before mapping
anything:

  1. Inspect the real directory structure — do not guess at paths.
  2. Draft a concrete region block with literal globs.
  3. Show it to the user and get confirmation:

         region "sim":
           include: src/simulation/**, src/dynamics/**, src/environment/**
           exclude: src/gnc/**, src/rl/**
         → 47 files in scope, 212 excluded. Proceed?

  4. Write it to `vulcan_mind/vulcan.config.yaml` under `regions:`
     (or `vulcan region add sim --include ... --exclude ...`).
  5. Run everything afterwards with `--region sim`.

Never map from an unwritten interpretation of scope. The config file is the
record of what you were asked to do, and rule V11 enforces it.

Excluded regions are not deleted from the map. Where in-scope code calls into
excluded code, create a node with `kind: external` so the boundary is visible
rather than silently truncated.

### Step 2 — Build the worklist

    vulcan worklist --build --region <name>

This enumerates every in-scope file. It is your ledger. It outlives your context.

### Step 2b — Scaffold the structure and the symbol nodes

    vulcan scaffold

On an empty map this does three things, in order, and prints each:

  1. **Bootstraps the structure sheet** — one covering node per source
     directory (`module.<dir>` with `covers: [<dir>/**]`, or an explicit file
     list for loose files at a root), grounded in the directory's `module`
     declaration where the language has one, on
     `graph/subcharts/structure.graph.json`. Rule V13 reads coverage from here.
  2. **Creates one subchart per covering node** (`graph/subcharts/<dir>.graph.json`)
     so symbols have a sheet to land on.
  3. **Creates a node, a doc skeleton and a containment edge** for every
     significant symbol the tool can enumerate — potentially thousands.

Structure comes from the source; the prose is left empty on purpose, so the
gate keeps failing until you write it. Do this before writing any prose:
hand-typing node entries for a real codebase is slow and gets ids, line numbers
and chart membership wrong.

### Step 2c — Derive the call edges

    vulcan trace --add-sockets
    vulcan compile

`trace` reads every symbol's body and writes a `call` edge (or `feedback` where
it closes a cycle) for each call to another mapped symbol, citing the call site
as evidence. Never type call edges by hand. `--add-sockets` declares the
standard `callers`/`callees` sockets in the frontmatter of any doc that lacks
them (a hand-authored doc, or one whose sockets were later rewritten), so no
observed call is dropped; without it such calls are counted and skipped.
Workflow and generated sheets lay themselves out again whenever the edges they
show change; `vulcan compile --reflow [CHART]` does the same for an authored
sheet on request.

### Step 3 — Trace a run: build the master as the operational flow

The master chart is **not** a diagram of which module includes which. It is the
program in operation, read left to right:

    inputs  →  configure  →  set up run  →  solve loop  →  outputs

Build it by tracing, not by listing directories:

  1. **Start at the entrypoint** the user actually invokes — the CLI `main`, the
     `run_*` function an example script calls, the campaign runner. Read it
     completely.
  2. **Follow the data forward.** At every call that hands work to another
     subsystem, ask: what enters, what is it turned into, what leaves. Each such
     stage is a candidate **phase block** — a `kind: group` node on the master
     named by what it *does* ("Set up run", "Solve loop", "Write results"), never
     by where it lives ("simulation/").
  3. **Make every input explicit.** Every file, dataset, kernel, asset or
     configuration object the program reads is a `kind: external` node with only
     outgoing edges. If a stage reads it, draw the edge, with evidence at the
     read site.
  4. **Make every output explicit.** Every artefact the program writes — result
     tables, bundles, checkpoints, plots, reports, caches — is a `kind: external`
     node with only incoming edges. Its doc states the format, the schema and who
     consumes it. This is how a reader learns what the program *can produce*.
  5. **Give every block a way in.** Each phase block carries
     `opens: <chart_id>` naming the sheet that shows how it works: a module
     subchart, a workflow view, or a generated block sheet. `vulcan find --charts`
     lists every chart id, generated ones included. A block that opens nothing
     is a dead end and fails V21.
  5b. **Make the entrypoint a workflow view.** The block for the solve/run
     phase should open a workflow seeded at the entrypoint, so "how does it
     tick" is answered by the real call graph, clustered into blocks:

         vulcan workflow add run-pipeline --title "How a run ticks" \
             --seed <entrypoint node id> --direction downstream --depth 2

     and set `opens: run-pipeline` on that block.
  6. **Edges are dataflow between phases**, each with `evidence` at the call site
     where the handoff is observable. Containment ("includes", "using") is not an
     edge on the master.

Rule V21 checks all of this mechanically: a source and a sink must exist, every
block must open something, and no node may be a hub. A package tree fails.

**The package structure still exists** — it is just not the master. Put the
module nodes with their `covers:` globs on a `structure` subchart
(`graph/subcharts/structure.graph.json`, `derives_from: master`). Rule V13 reads
coverage from wherever those nodes live, and every module node must still be
expanded to function granularity by a subchart (V13b). A module node may not
claim files outside its own root (V13a).

### Step 3b — Map the symbols, batch by batch

Follow the batch loop in the completion contract. For each file:

  - Read the whole file.
  - Every significant symbol in it already has a scaffolded node (V13e). Your
    job is the prose for each: what it does, its real inputs and outputs, its
    assumptions and limits. One node per symbol — never one standing in for
    the file.
  - Write it in batches with the prose pipeline, one source file per payload:

        vulcan prose remaining          # outstanding docs, grouped by file
        vulcan prose show <file>        # the source of that file's symbols
        vulcan prose apply batch.json   # {"<doc path>": {"purpose","design","limits","math"?}}

    `apply` refuses anything V12 would refuse — banned phrases, under the
    floor — and touches only the author-owned sections, so a bad batch cannot
    poison the map. Fill **every** section; `min_doc_words` is enforced on
    prose only, so generated tables cannot pad a stub past the floor.
  - Declare real inputs/outputs in frontmatter — actual argument names and types.
    Sockets are authored here and lifted into the chart JSON by the compiler;
    never hand-write sockets into JSON.
  - Add edges for every connection you can *point at*: a call site, an argument
    passed, a struct field written. Each edge needs `evidence` naming the file
    and lines where it is observable.

### Step 4 — Converge

Repeat until `vulcan worklist --remaining` reports 0 **and**
`vulcan check --strict` exits 0.

Common failures and what they mean:

| Rule | Meaning | Fix |
|---|---|---|
| V6 | A symbol you named is not in the file you named | Open the file; use the real name |
| V7 | An edge cites a file that does not exist | Point at real code, or drop the edge |
| V10 | A socket is in JSON but not in frontmatter | Declare it in the doc; frontmatter is canonical |
| V12 | Vague prose, or too few words | Write the real content |
| V13 | Files are unaccounted for | You are not finished — keep mapping |
| V13b | A module has no function-level subchart | Run the `vulcan-subchart` skill for it |
| V13d | A file is claimed but not described | Give it a real node |
| V13e | Symbols in a file have no node of their own | `vulcan scaffold`, then write their prose |
| V12 on `nodes/<module>/src_*.md` | A macro block (group node) generated by compile has no prose | Write what the block *is* — the role of that directory or file |
| V20 | A sheet still renders more than the readability limit | Split it by hand with an explicit subchart; never raise the limit |
| V21 | The master is a package tree, lacks an input/output node, or has an unclickable block | Rebuild it as the operational flow (Step 3); add `external` sources/sinks; set `opens` on every block |

### Step 5 — Report

Report only after `--strict` passes, and **paste the literal output of
`vulcan check --strict --proof`** as the first thing in your report. Without that
block your claim is void; with a block reading `verdict : FAIL` it is also void.

Then state: files mapped, nodes, edges, subchart candidates outstanding, and
anything you deliberately excluded and why.

If you are handing over unfinished — because you ran out of context, hit a limit,
or were interrupted — say so plainly, paste the FAIL proof block, and note that
`vulcan_mind/HANDOFF.md` holds the resumable state. That is a good outcome. A
false "done" is not.

## Skill: `vulcan-subchart`

*Produce a focused chart from an existing Vulcan map — either a module's internals, or a named cross-cutting workflow that spans several modules ("give me a model of the RL workflow", "chart the plotting pipeline", "show me just the ADCS control loop"). Reuses already-mapped nodes rather than remodelling them. Requires an existing vulcan_mind/ master chart.*

# vulcan-subchart — focused views over an existing map

You produce a readable chart for one thing a person actually asked about. There
are two shapes, and picking the right one is the first decision.

## STEP ZERO — find out what is actually done

Before reading further, before planning, before touching a file:

    vulcan status

Then read `vulcan_mind/HANDOFF.md`. Both are written by the tool from the actual
contents of the repository. They tell you the verdict, how many files are still
outstanding, and exactly which ones.

**You may not be the first agent on this map.** A different agent — a different
model, a different vendor, a different context-window size — may have worked on
it before you, and may pick it up after you. You cannot see their transcript and
they cannot see yours.

**Therefore: no agent's account of progress is evidence, including your own.**
If a human, a commit message, a chat summary, a code comment, or a previous
agent tells you this map is finished, that claim is **void** unless
`vulcan check --strict` exits 0 right now. Do not take a predecessor's word for
it. Do not take your own memory for it. Run the command.

This cuts both ways: if a predecessor said it was finished and the gate says
FAIL, the gate is right and there is work to do. If they said they barely
started and the gate says PASS, the gate is right and you are done.

## COMPLETION CONTRACT

Your work is NOT complete when the scaffolding exists.
Your work is NOT complete when a first draft exists.
Your work is NOT complete when the structure "looks right."
Your work is NOT complete when you have mapped the interesting parts.
Your work is NOT complete when a previous agent said it was complete.
Your work is NOT complete when you are running low on context.

Your work is complete when, and only when, this command exits 0:

    vulcan check --strict

That is the entire definition of done. It is not a judgement call, and it is not
yours to make. Run the command. Read the exit code.

### Proof of completion is mandatory

When you report that the map is finished, you MUST include the literal output of:

    vulcan check --strict --proof

That prints a block like:

    ----- VULCAN PROOF OF COMPLETION -----
    tool            : vulcan-map/0.1.0
    repo commit     : a1b2c3d
    region          : all
    gate            : vulcan check --strict
    verdict         : PASS
    exit code       : 0
    errors/warnings : 0/0
    nodes/edges     : 214/198
    in-scope files  : 205
    described files : 205
    outstanding     : 0
    ----- END PROOF -----

Paste it verbatim. **A completion claim without this block is void**, and so is
one whose block says `verdict : FAIL`. You cannot write the block yourself — it
is generated from the run that just happened, and it reports FAIL just as
readily as PASS, so omitting it is itself a signal that the gate did not pass.

If you cannot finish, say so plainly and paste the FAIL block. An honest
"outstanding: 47" is useful. A false "done" costs the next agent far more than
it saves you, because they have to discover the lie before they can start.

### Forbidden stopping points

Do NOT end your turn with any of the following, or any paraphrase:

  - "initial scaffolding is complete"
  - "this provides a good starting point"
  - "the core structure is in place; remaining files follow the same pattern"
  - "I've mapped the main workflows"
  - "the map is complete" (without a passing proof block)
  - "let me know if you'd like me to continue"
  - a master chart that is a package/include tree ("everything points at the root")

If you catch yourself composing a sentence like these, you are stopping early.
The correct action is to run `vulcan status` and keep working.

"The remaining files follow the same pattern" is never a reason to skip them.
The map's value is precisely that it is exhaustive. A 70%-complete map is not 70%
as useful as a complete one — it is actively misleading, because a reader cannot
tell which absences mean "not present in the codebase" and which mean "not yet
done."

### Context exhaustion is not completion

This task is expected to exceed your context window. That is normal and planned
for. Progress lives on disk in `vulcan_mind/HANDOFF.md` and
`vulcan_mind/_build/status.json`, not in your context, and not in your summary.

When context runs short:

  1. Finish the file you are on.
  2. Run `vulcan compile` to checkpoint — this rewrites HANDOFF.md.
  3. State plainly: "Checkpoint: N of M files mapped. Not finished." and paste
     the FAIL proof block.
  4. Stop cleanly. The next agent will read HANDOFF.md and continue.

Never report a partial map as a finished one. Handing over honestly is a success;
handing over a false "done" is the single most expensive thing you can do here.

### The batch loop

Work in batches of 5-10 source files:

  1. `vulcan status` / `vulcan prose remaining`  → what is outstanding
  2. Read each file **completely** (`vulcan prose show <file>`). Do not skim;
     do not infer from filenames.
  3. Write one payload per file and `vulcan prose apply` it; call edges come
     from `vulcan trace --add-sockets`, never from typing.
  4. `vulcan compile`
  5. `vulcan check` → fix every error before taking the next batch
  6. Repeat until `vulcan status` reports 0 outstanding.
  7. `vulcan check --strict --proof` → must say PASS, exit 0.

### Granularity: one node per symbol, never one node per file

**Every significant symbol gets its own node.** Every function, type, struct,
macro, module, and non-dunder method in every in-scope file. Not one
representative symbol standing in for the file it lives in.

This is rule V13e and it is machine-checked: the tool enumerates each file's
symbols itself (Python via the stdlib AST, Julia via declaration patterns) and
fails until each one has a node.

  REQUIRED:  a file with 12 functions produces 12 nodes
  FORBIDDEN: a file with 12 functions produces 1 node "representing" it

The one exception is deduplication: several methods of a single generic function
are one symbol. `calcForceTorque` with eight methods is one node, not eight.

**Why this is not negotiable.** The first real map built with this tool used one
node per file — 225 nodes over 205 files. Every gate passed and the result was
unusable. You could not click into anything, because a file-level node has no
interior. And the call tracer found almost no edges, because a call from one
function to another inside a mapped file had no individual nodes to connect. The
map's granularity *is* the product.

**Use `vulcan scaffold`.** It generates every missing node and a doc skeleton for
it — id, label, kind, source file, declaration line, chart membership — straight
from the source. That is mechanical work you should not be typing by hand. It
deliberately leaves Purpose and Design empty, so the gate keeps failing until you
read the code and write them. A scaffold is not a map.

Doc length is expected to scale with what a node covers: a module doc needs real
depth (120 words), a leaf function needs a tight, accurate 40. Do not pad a small
function's doc to look like a big one — padding is exactly what the banned-phrase
lint is looking for.

### Readability: every sheet must be clickable, never a wall (D13)

**No sheet may render more than `max_nodes_per_sheet` nodes (40).** This is
rule V20 and it is machine-checked on the *rendered* graph, after resolution.

Per-symbol granularity (above) means a module has hundreds of nodes. They must
not all appear on one sheet. `vulcan compile` handles this for you: it clusters
an oversized sheet into **macro blocks** — `group` nodes named after the
directory, file, or symbol-name prefix they stand for — and generates a nested
sheet behind each one. Edges between blocks are lifted and counted. A reader
double-clicks a block to see inside it. The same block opens the whole file
from a module sheet and only the workflow's slice from a workflow sheet.

  REQUIRED:  a GNC sheet showing ~10 blocks (control/, guidance/, navigation/ ...)
             each of which opens into files, each of which opens into symbols
  FORBIDDEN: a GNC sheet showing 500 symbols at once, however well laid out

**Group nodes are nodes.** Each one has a doc at `nodes/<module>/<block>.md`,
scaffolded by compile with an empty Purpose. You write what the block *is* — the
role of that directory or file in the system, what enters it and what leaves —
at the symbol floor (40 words). It is not a list of members; the members are on
the nested sheet. Until every group doc is written, `vulcan check` fails on V12.

**Do not edit generated nested charts** (`provenance.generated_by: cluster`).
Their membership is recomputed on every compile from the code's structure; your
edits will not survive and will not change the verdict.

**If V20 fires after compile**, clustering could not partition that sheet: it is
one file whose symbols share a single name prefix. Split it by hand — a subchart
with explicit `member_nodes` — or accept the finding as a real limitation of the
code's structure and say so. Never raise the limit.

### The master is the operational flow, never the package tree (D14)

**The master sheet answers two questions: what does this program produce, and
how does it tick.** It is read left to right like a Blender node tree:

  inputs  →  configure  →  set up run  →  solve loop  →  outputs

This is rule V21 and it is machine-checked:

  - at least one `kind: external` **source** node (a file, dataset or argument
    the program reads) and at least one `kind: external` **sink** node (an
    artefact it writes) must be on the master and wired in;
  - every non-leaf master node — every phase block — must carry
    `opens: <chart_id>` naming the sheet that shows how it works, or be expanded
    by one. A block you cannot double-click is a dead end and fails;
  - no master node may touch more than half of the master's edges. A tree where
    thirteen modules point at the root module is a package diagram, not a map,
    and it fails.

  REQUIRED:  manifest.toml → parse_cli → build_initial_conditions → integrator
             (callbacks, RHS) → results.csv / checkpoint / report
  FORBIDDEN: module.gnc → module.spaceagora, module.io → module.spaceagora, ...

**Trace a run to build it.** Start from the entrypoint the user invokes. Follow
the data: what is read, what each phase turns it into, what is written and
where. Name blocks by what they *do* in that flow, not by the directory they
live in. Give every output node a doc that states its schema and who consumes
it. The module-containment view still exists — put it on a `structure`
subchart, where V13 coverage reads the `covers` globs — but it is not the
master.

### Grounding is checked mechanically — you cannot talk your way past it

Every node names a real file and a real symbol. `vulcan check` **opens the file
and looks for the symbol** (rule V6). Every edge cites a real file where the
connection is observable (rule V7). Node prose is linted against a banned-phrase
list (rule V12). Every in-scope file must have a node describing it (rule V13d);
a `covers` glob accounts for a file but does not describe it.

  REQUIRED:  label `propagate_orbit`, source `src/simulation/propagator.jl:142`
  FORBIDDEN: "a helper function", "the main solver routine", "various utilities"

If you have not opened the file and seen the symbol, you may not write it down.
Invented names fail the gate and you will have to redo the work.

### Do not silence a failing check

If a rule fires, fix the map. Do not edit `vulcan.config.yaml` to remove a banned
phrase, lower `min_doc_words`, disable a coverage rule, or narrow a region so
that unmapped files fall out of scope. Do not adjust recorded data solely to make
a check stop complaining — if you believe the tool is wrong, say so explicitly
and report it, with evidence, rather than bending the map to fit.

Changing the gate to pass the gate is the one failure mode this whole design
exists to prevent.
## Which shape are you building?

| The request | Shape | Chart kind |
|---|---|---|
| "chart the vehicle module", "expand dynamics" | **Module detail** — one module, finer grain | `sub` |
| "model the RL workflow", "show me the aerobraking pipeline" | **Workflow view** — one behaviour, across modules | `workflow` |

A workflow almost never lives inside one module. An RL workflow plausibly touches
dynamics, GNC, simulation and analysis. If you build it as a module subchart you
will either truncate it at the module boundary or start remodelling other
modules inside it — both wrong.

## Prerequisite

A master chart must exist. If `vulcan_mind/graph/master.graph.json` is absent,
stop and tell the user to run the `vulcan-map` skill first. Do not invent a
master.

---

# Shape A — module detail (`chart_kind: sub`)

Members are the module's own symbol nodes; `expands` points each at the module
node it refines. See the granularity rules in the completion contract: every
significant symbol in the module's files already needs a node (V13e), so this
chart is mostly assembling what exists rather than creating.

---

# Shape B — cross-cutting workflow view (`chart_kind: workflow`)

## The mechanism, and why it is split this way

Deciding *what counts as "the RL workflow"* is a reading of intent. No traversal
can derive it, and pretending otherwise would produce confident nonsense. But
deciding *what that workflow touches* is not a judgement call at all — it follows
from the call graph.

So the work splits:

- **You choose the seeds.** The entry points where the workflow begins. This is
  the semantic step, it is yours, and it gets written into the chart file where a
  human can read and correct it.
- **The compiler computes membership.** From the seeds it walks the traced call
  edges and generates `member_nodes`. You never type that list. It is regenerated
  on every compile, so the view cannot drift away from the code it describes.

This is the same pattern as sockets and the Connections block: a human-authored
input, a generated projection, and a rule that fails if they disagree.

## Procedure

### 1. Find what is already mapped — do not start from the source

    vulcan find rl
    vulcan find policy
    vulcan find reward

`find` searches node ids, labels, tags, file paths and doc prose, and tells you
which of those matched. **Read the results before creating anything.** If a
dynamics pipeline is already mapped with real depth, the RL view must reference
those nodes, not build a second, shallower model of dynamics inside itself.

Creating a duplicate node for a symbol another chart already owns is rejected by
rule V19. The map must not disagree with itself about what a symbol is.

### 2. Choose and justify the seeds

Pick the entry points — usually the top-level functions a person would name if
asked "where does this workflow start?". State them back to the user with your
reasoning before creating the chart:

    workflow "rl": seeds
      gnc.rl_train_policy    — the training entry point
      gnc.rl_rollout         — the rollout loop the trainer drives
    traversal: both, depth 4
    → reaches 31 nodes across gnc, dynamics, simulation, analysis. Proceed?

Seeds must already be mapped nodes. If the workflow's entry point is not on the
map yet, map it first with `vulcan-map` — a workflow view borrows nodes, it does
not invent them.

### 3. Create it

    vulcan workflow add rl \
      --title "Reinforcement-learning training workflow" \
      --seed gnc.rl_train_policy \
      --seed gnc.rl_rollout \
      --why "training entry point and the rollout loop it drives" \
      --direction both --depth 4

`--direction downstream` follows what the seeds call; `upstream` follows what
calls them; `both` gives the whole neighbourhood. `--depth` bounds the reach —
raise it if the view stops short of the real terminus, lower it if it swallows
half the repo.

### 4. Compile and read the result

    vulcan compile

Membership is generated here. Then open it and check it against your
understanding of the workflow:

- **Too small?** The call edges may not be traced yet, or the direction is wrong.
- **Too large?** A large workflow is fine — compile nests it (D13). The top
  sheet shows macro blocks (directories and files of the modules the workflow
  crosses) with lifted, counted edges; each block opens into the workflow's
  members of that block. Only reduce `--depth` or narrow the seeds if the
  workflow genuinely spans more than one concern.

### 5. Only now, create nodes — and only for genuine gaps

If part of the workflow is real but nothing on the map covers it, add it as a
`local_node` on the workflow chart with full grounding and a real doc, exactly as
in the master skill. That node then belongs to this chart and is borrowed by any
other view that reaches it.

Never add a local node for a symbol that already has one elsewhere. Borrow it.

## Completeness

For a **module detail** chart: every significant symbol in the module's files has
a node, and the traced call edges between them are present.

For a **workflow view**: the traced closure covers the behaviour end to end, from
the seeds to its terminal effects. If it stops at a boundary, either follow it by
raising the depth, or terminate deliberately at a `kind: external` stub and say
so. Every dangling thread is followed or explicitly marked — there is no third
option.

`vulcan check --strict --proof` must pass before you report either as done.

## Skill: `vulcan-augment`

*Add fidelity to an existing Vulcan chart — deepen one area of a master or sub chart on request, e.g. "augment the plotting subchart with more detail on animations". Additive by default; never silently removes existing map content.*

# vulcan-augment — deepen an existing chart

You are increasing the resolution of a **bounded region** of an existing chart.
You are not rebuilding it.

## STEP ZERO — find out what is actually done

Before reading further, before planning, before touching a file:

    vulcan status

Then read `vulcan_mind/HANDOFF.md`. Both are written by the tool from the actual
contents of the repository. They tell you the verdict, how many files are still
outstanding, and exactly which ones.

**You may not be the first agent on this map.** A different agent — a different
model, a different vendor, a different context-window size — may have worked on
it before you, and may pick it up after you. You cannot see their transcript and
they cannot see yours.

**Therefore: no agent's account of progress is evidence, including your own.**
If a human, a commit message, a chat summary, a code comment, or a previous
agent tells you this map is finished, that claim is **void** unless
`vulcan check --strict` exits 0 right now. Do not take a predecessor's word for
it. Do not take your own memory for it. Run the command.

This cuts both ways: if a predecessor said it was finished and the gate says
FAIL, the gate is right and there is work to do. If they said they barely
started and the gate says PASS, the gate is right and you are done.

## COMPLETION CONTRACT

Your work is NOT complete when the scaffolding exists.
Your work is NOT complete when a first draft exists.
Your work is NOT complete when the structure "looks right."
Your work is NOT complete when you have mapped the interesting parts.
Your work is NOT complete when a previous agent said it was complete.
Your work is NOT complete when you are running low on context.

Your work is complete when, and only when, this command exits 0:

    vulcan check --strict

That is the entire definition of done. It is not a judgement call, and it is not
yours to make. Run the command. Read the exit code.

### Proof of completion is mandatory

When you report that the map is finished, you MUST include the literal output of:

    vulcan check --strict --proof

That prints a block like:

    ----- VULCAN PROOF OF COMPLETION -----
    tool            : vulcan-map/0.1.0
    repo commit     : a1b2c3d
    region          : all
    gate            : vulcan check --strict
    verdict         : PASS
    exit code       : 0
    errors/warnings : 0/0
    nodes/edges     : 214/198
    in-scope files  : 205
    described files : 205
    outstanding     : 0
    ----- END PROOF -----

Paste it verbatim. **A completion claim without this block is void**, and so is
one whose block says `verdict : FAIL`. You cannot write the block yourself — it
is generated from the run that just happened, and it reports FAIL just as
readily as PASS, so omitting it is itself a signal that the gate did not pass.

If you cannot finish, say so plainly and paste the FAIL block. An honest
"outstanding: 47" is useful. A false "done" costs the next agent far more than
it saves you, because they have to discover the lie before they can start.

### Forbidden stopping points

Do NOT end your turn with any of the following, or any paraphrase:

  - "initial scaffolding is complete"
  - "this provides a good starting point"
  - "the core structure is in place; remaining files follow the same pattern"
  - "I've mapped the main workflows"
  - "the map is complete" (without a passing proof block)
  - "let me know if you'd like me to continue"
  - a master chart that is a package/include tree ("everything points at the root")

If you catch yourself composing a sentence like these, you are stopping early.
The correct action is to run `vulcan status` and keep working.

"The remaining files follow the same pattern" is never a reason to skip them.
The map's value is precisely that it is exhaustive. A 70%-complete map is not 70%
as useful as a complete one — it is actively misleading, because a reader cannot
tell which absences mean "not present in the codebase" and which mean "not yet
done."

### Context exhaustion is not completion

This task is expected to exceed your context window. That is normal and planned
for. Progress lives on disk in `vulcan_mind/HANDOFF.md` and
`vulcan_mind/_build/status.json`, not in your context, and not in your summary.

When context runs short:

  1. Finish the file you are on.
  2. Run `vulcan compile` to checkpoint — this rewrites HANDOFF.md.
  3. State plainly: "Checkpoint: N of M files mapped. Not finished." and paste
     the FAIL proof block.
  4. Stop cleanly. The next agent will read HANDOFF.md and continue.

Never report a partial map as a finished one. Handing over honestly is a success;
handing over a false "done" is the single most expensive thing you can do here.

### The batch loop

Work in batches of 5-10 source files:

  1. `vulcan status` / `vulcan prose remaining`  → what is outstanding
  2. Read each file **completely** (`vulcan prose show <file>`). Do not skim;
     do not infer from filenames.
  3. Write one payload per file and `vulcan prose apply` it; call edges come
     from `vulcan trace --add-sockets`, never from typing.
  4. `vulcan compile`
  5. `vulcan check` → fix every error before taking the next batch
  6. Repeat until `vulcan status` reports 0 outstanding.
  7. `vulcan check --strict --proof` → must say PASS, exit 0.

### Granularity: one node per symbol, never one node per file

**Every significant symbol gets its own node.** Every function, type, struct,
macro, module, and non-dunder method in every in-scope file. Not one
representative symbol standing in for the file it lives in.

This is rule V13e and it is machine-checked: the tool enumerates each file's
symbols itself (Python via the stdlib AST, Julia via declaration patterns) and
fails until each one has a node.

  REQUIRED:  a file with 12 functions produces 12 nodes
  FORBIDDEN: a file with 12 functions produces 1 node "representing" it

The one exception is deduplication: several methods of a single generic function
are one symbol. `calcForceTorque` with eight methods is one node, not eight.

**Why this is not negotiable.** The first real map built with this tool used one
node per file — 225 nodes over 205 files. Every gate passed and the result was
unusable. You could not click into anything, because a file-level node has no
interior. And the call tracer found almost no edges, because a call from one
function to another inside a mapped file had no individual nodes to connect. The
map's granularity *is* the product.

**Use `vulcan scaffold`.** It generates every missing node and a doc skeleton for
it — id, label, kind, source file, declaration line, chart membership — straight
from the source. That is mechanical work you should not be typing by hand. It
deliberately leaves Purpose and Design empty, so the gate keeps failing until you
read the code and write them. A scaffold is not a map.

Doc length is expected to scale with what a node covers: a module doc needs real
depth (120 words), a leaf function needs a tight, accurate 40. Do not pad a small
function's doc to look like a big one — padding is exactly what the banned-phrase
lint is looking for.

### Readability: every sheet must be clickable, never a wall (D13)

**No sheet may render more than `max_nodes_per_sheet` nodes (40).** This is
rule V20 and it is machine-checked on the *rendered* graph, after resolution.

Per-symbol granularity (above) means a module has hundreds of nodes. They must
not all appear on one sheet. `vulcan compile` handles this for you: it clusters
an oversized sheet into **macro blocks** — `group` nodes named after the
directory, file, or symbol-name prefix they stand for — and generates a nested
sheet behind each one. Edges between blocks are lifted and counted. A reader
double-clicks a block to see inside it. The same block opens the whole file
from a module sheet and only the workflow's slice from a workflow sheet.

  REQUIRED:  a GNC sheet showing ~10 blocks (control/, guidance/, navigation/ ...)
             each of which opens into files, each of which opens into symbols
  FORBIDDEN: a GNC sheet showing 500 symbols at once, however well laid out

**Group nodes are nodes.** Each one has a doc at `nodes/<module>/<block>.md`,
scaffolded by compile with an empty Purpose. You write what the block *is* — the
role of that directory or file in the system, what enters it and what leaves —
at the symbol floor (40 words). It is not a list of members; the members are on
the nested sheet. Until every group doc is written, `vulcan check` fails on V12.

**Do not edit generated nested charts** (`provenance.generated_by: cluster`).
Their membership is recomputed on every compile from the code's structure; your
edits will not survive and will not change the verdict.

**If V20 fires after compile**, clustering could not partition that sheet: it is
one file whose symbols share a single name prefix. Split it by hand — a subchart
with explicit `member_nodes` — or accept the finding as a real limitation of the
code's structure and say so. Never raise the limit.

### The master is the operational flow, never the package tree (D14)

**The master sheet answers two questions: what does this program produce, and
how does it tick.** It is read left to right like a Blender node tree:

  inputs  →  configure  →  set up run  →  solve loop  →  outputs

This is rule V21 and it is machine-checked:

  - at least one `kind: external` **source** node (a file, dataset or argument
    the program reads) and at least one `kind: external` **sink** node (an
    artefact it writes) must be on the master and wired in;
  - every non-leaf master node — every phase block — must carry
    `opens: <chart_id>` naming the sheet that shows how it works, or be expanded
    by one. A block you cannot double-click is a dead end and fails;
  - no master node may touch more than half of the master's edges. A tree where
    thirteen modules point at the root module is a package diagram, not a map,
    and it fails.

  REQUIRED:  manifest.toml → parse_cli → build_initial_conditions → integrator
             (callbacks, RHS) → results.csv / checkpoint / report
  FORBIDDEN: module.gnc → module.spaceagora, module.io → module.spaceagora, ...

**Trace a run to build it.** Start from the entrypoint the user invokes. Follow
the data: what is read, what each phase turns it into, what is written and
where. Name blocks by what they *do* in that flow, not by the directory they
live in. Give every output node a doc that states its schema and who consumes
it. The module-containment view still exists — put it on a `structure`
subchart, where V13 coverage reads the `covers` globs — but it is not the
master.

### Grounding is checked mechanically — you cannot talk your way past it

Every node names a real file and a real symbol. `vulcan check` **opens the file
and looks for the symbol** (rule V6). Every edge cites a real file where the
connection is observable (rule V7). Node prose is linted against a banned-phrase
list (rule V12). Every in-scope file must have a node describing it (rule V13d);
a `covers` glob accounts for a file but does not describe it.

  REQUIRED:  label `propagate_orbit`, source `src/simulation/propagator.jl:142`
  FORBIDDEN: "a helper function", "the main solver routine", "various utilities"

If you have not opened the file and seen the symbol, you may not write it down.
Invented names fail the gate and you will have to redo the work.

### Do not silence a failing check

If a rule fires, fix the map. Do not edit `vulcan.config.yaml` to remove a banned
phrase, lower `min_doc_words`, disable a coverage rule, or narrow a region so
that unmapped files fall out of scope. Do not adjust recorded data solely to make
a check stop complaining — if you believe the tool is wrong, say so explicitly
and report it, with evidence, rather than bending the map to fit.

Changing the gate to pass the gate is the one failure mode this whole design
exists to prevent.
## Step 1 — Bound the augmentation, explicitly

From the request ("augment the plotting subchart with additional fidelity
regarding animations"), identify and state back:

  - target chart: `plotting`
  - target subgraph: nodes reachable from `viz.animate_sequence`
  - fidelity axis: animation frame lifecycle, timing, and state
  - out of scope: static plotting path (already mapped)

Get agreement before editing. Ambiguous augmentation requests produce sprawl.

## Step 2 — Augment

Fidelity is added in four ways. Use whichever the request calls for:

  a. **Decomposition** — expand a coarse node into finer nodes (`expands` set).
  b. **Socket refinement** — split an over-general socket into the real distinct
     inputs/outputs, with real parameter names and types.
  c. **Edge refinement** — replace one vague edge with the several specific
     dataflows it was standing in for; each needs its own `evidence`.
  d. **Doc deepening** — fill in the math, assumptions, and limitations sections
     with real detail: actual equations, actual numerical tolerances, actual
     failure modes.

## Step 3 — Additive discipline (D6)

Augmentation is **additive by default**.

  - Never delete a node or edge whose `origin` is `"human"` — those are the
    user's own edits. Removing them requires `--allow-removal` and an explicit
    instruction.
  - Never reposition existing nodes. The layout engine leaves `ui.pos` alone and
    so do you.
  - If existing content is genuinely *wrong*, do not silently overwrite it.
    State what is wrong, why, and what you propose — then change it only with
    agreement.

## Step 4 — Converge

`vulcan compile && vulcan check --strict` must exit 0.

Augmentation frequently breaks two rules in particular:

  - **V10 (socket parity)** — you refined sockets in frontmatter but an edge
    still references the old socket id. Update the edge.
  - **V5 (acyclicity)** — a finer decomposition can expose a real cycle. If it is
    genuine recursion, classify the back-edge as `kind: "feedback"` (D5). If it
    is not, the decomposition is wrong. Do not disable the rule.
  - **V12 on a `src_*.md` group doc** — adding symbols pushed a sheet over the
    readability limit, so compile clustered it (D13) and scaffolded a doc for
    the new macro block. Write it. Never edit the generated nested chart.

## Completeness for an augmentation

Complete = the stated fidelity axis is at the requested depth **across the whole
bounded subgraph**, not just the first node you touched. If you deepened
`animate_sequence` but left its three children coarse, you are not done.

