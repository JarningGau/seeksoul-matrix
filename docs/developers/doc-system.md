# Documentation System

This page defines how developer documentation should be maintained for seeksoul-matrix.

## Target Audiences

- [`README.md`](../../README.md): chemistry, installation, top-level navigation
- [`docs/developers/`](./): stage contracts, implementation notes, QC definitions, validation log
- Reference templates (`dbit-matrix/`, `SeekSoulMethyl/`): engineering and domain patterns — read-only unless explicitly updating the template

## Source Of Truth

| Topic | Canonical page |
|-------|----------------|
| Stage input/output paths and I/O-level contracts | [`contracts.md`](contracts.md) |
| Read-order vs analysis chunks, gather rules, barcode-selection modes | [`chunk_model.md`](chunk_model.md) |
| QC field definitions, demux stats schema, saturation model, summary columns | [`qc_metrics.md`](qc_metrics.md) |
| Per-stage implementation (defaults, skip logic, toolchain, SeekSoul alignment) | [`stage_notes/`](stage_notes/) |
| Manual validation history | [`logs.md`](logs.md) |

## Update Rules

When behavior changes:

- update [`contracts.md`](contracts.md) if stage inputs, outputs, stage order, or I/O-level contracts change
- update [`chunk_model.md`](chunk_model.md) if chunk/shard boundaries, gather semantics, or barcode-selection routing change
- update [`qc_metrics.md`](qc_metrics.md) if QC column names, metric definitions, or saturation/qc_summary formulas change
- update the matching [`stage_notes/<stage>.md`](stage_notes/) if defaults, skip/re-run logic, CLI flags, or environment setup change without altering I/O paths
- append an entry to [`logs.md`](logs.md) after implementation (date, task, files, check, status)

## Writing Rules

- keep [`contracts.md`](contracts.md) lean: paths, file tables, and normative I/O constraints only
- put algorithm and field semantics in [`qc_metrics.md`](qc_metrics.md) or [`stage_notes/`](stage_notes/), not in contracts
- do not duplicate normative facts across multiple documents — link to the canonical page
- keep stage names, workflow keys, and output paths exact
- write new documentation in English under `docs/`

## Legacy References

Historical log entries may cite `contracts.md` before modularization; treat those as archive context. New work should link to the appropriate layer above.
