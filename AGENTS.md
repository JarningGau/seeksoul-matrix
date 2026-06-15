# Repository Guidelines

Contributor guide for **seeksoul-matrix** — SeepSpace spatial methylation (DD-MET5, Slide-Tag–like chemistry). This repository is under active development.

Two **reference templates** sit alongside the production codebase — study them when building seeksoul-matrix; do not treat either as the production codebase or modify them unless explicitly updating the template:

- **`dbit-matrix/`** — engineering patterns (thin stage scripts, JSON workflows, `make_cmd.py` drivers, stage contracts, docs layout).
- **`SeekSoulMethyl/`** — SeekSoul's official single-cell methylation + RNA data-processing pipeline; domain reference for DD-MET3/DD-MET5 chemistry, barcode rules, stage logic, and analysis steps.

## Project Structure & Module Organization

```
seeksoul-matrix/
├── README.md              # DD-MET5 chemistry, barcode rules, whitelist path
├── pixi.toml              # Root workspace (linux-64)
├── whitelist/             # Cell-barcode whitelists (e.g. DD-MET5/U3CB_methylation.txt.gz)
├── scripts/               # Stage scripts and workflow driver (make_cmd.py)
├── workflow/              # JSON workflow configs
├── docs/                  # Project documentation and stage contracts
├── dbit-matrix/           # Reference template — engineering patterns
└── SeekSoulMethyl/        # Reference template — SeekSoul official methylation pipeline
```

Barcodes for methylation data must not contain `C` (enzymatic C→T conversion). Fixed sequences (TSO, 17L, ME) retain `C` for conversion-rate QC — see `README.md`.

When adding seeksoul-matrix code, mirror **dbit-matrix** engineering conventions (thin stage scripts, JSON-driven workflows, `make_cmd.py` drivers, explicit stage contracts) and align stage behavior with **SeekSoulMethyl** where chemistry and analysis steps overlap.

## Build, Test, and Development Commands

**Root environment** (primary):

```bash
pixi install
```

**Reference templates** (read-only guidance; not the seeksoul-matrix runtime):

```bash
cd dbit-matrix && pixi install        # engineering patterns
# SeekSoulMethyl uses conda — see SeekSoulMethyl/README.md
```

Useful references:

- **dbit-matrix:** `docs/developers/architecture.md`, `contracts.md`, `AGENTS.md`
- **SeekSoulMethyl:** `README.md` (chemistry, installation), `nf/` (Nextflow workflow), `nf/bin/` (stage scripts), `docs/` (BAM dedup, reference genome)

## Coding Style & Naming Conventions

Follow **dbit-matrix** engineering patterns when implementing seeksoul-matrix; consult **SeekSoulMethyl** for methylation-specific stage logic and chemistry details:

- **Python 3.11**; 4-space indentation; `snake_case` for modules and functions.
- Thin, explicit single-stage scripts; workflow parameters in `workflow/*.json`.
- Every stage must support `--dry-run`.
- New docs in **English** under `docs/`.
- Update `pixi.lock` when root dependencies change.

## Testing Guidelines

No seeksoul-matrix test suite yet. When adding tests, follow the template's regression style in `dbit-matrix/docs/maintenance/`. Before finishing workflow changes, run relevant CLI `--help`, `--version`, and `--dry-run`.

## Lightweight Development Loop

Minimal loop for pipeline work — not a heavyweight engineering process.

### Before editing

Briefly declare the task:

- what will be changed
- why the change is needed
- which files are likely to be touched
- how the result will be checked

### During implementation

- make the smallest change that solves the task
- avoid unrelated refactoring
- preserve existing input/output formats unless explicitly required
- preserve existing command-line behavior when possible
- keep code and documentation consistent with the current project structure

### After editing

1. Run the most relevant lightweight check available — e.g. script on a small input, one representative pipeline step, expected output files, existing tests, or generated reports/tables.
2. Update `docs/developers/logs.md` with a short entry:

   - **date**
   - **task name**
   - **files changed**
   - **summary of changes**
   - **check performed**
   - **status:** `done`, `needs_review`, or `blocked`
   - **notes:** remaining risks, if any

**Validation rule:** if no check was run, status must be `needs_review` (not `done`); record why.

## Commit & Pull Request Guidelines

Use [Conventional Commits](https://www.conventionalcommits.org/): `type: short description` in imperative mood, lowercase, under 72 characters. Common types: `feat`, `fix`, `refactor`, `docs`, `chore`, `data`, `analysis`.

```
feat: add DD-MET5 whitelist validation
fix: reject barcodes containing C in demux
docs: document library structure in README
```

One commit per logical change. PRs should link issues and describe chemistry or workflow impact.

## Agent-Specific Instructions

- Build new functionality at the **repository root**, adapting engineering patterns from `dbit-matrix/` and domain logic from `SeekSoulMethyl/` — do not extend either template in place.
- Do not change stage input/output contracts silently; document contracts in `docs/` as they are defined.
- Prefer the simplest fix that preserves required behavior.
- Do not overwrite user changes without explicit permission.
