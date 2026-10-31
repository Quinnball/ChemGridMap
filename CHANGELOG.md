# Changelog

## 0.5.0 — 2026-10-29

- Apply the offset-aware display rule: an entry touched by a flagged document
  pair is uncertain and pale, with `status_reason=offset_pair`, while all model
  probabilities remain unchanged. Record touched-entry counts in fits and
  manifests, show document-pair details and a distinct bottom-left map marker.
- Rename the user-facing "conflict" flag to "wide record range", retain
  `is_conflicted` for compatibility, and add the equal-valued
  `wide_record_range` alias.
- Add `neighbor_recall_fixed_k` for expected recall under random tie-breaking.
  Keep the existing tie-inclusive `neighbor_recall` unchanged.
- Make UMAP layouts repeatable for a given seed when the kNN graph is
  disconnected: such graphs are rerun from UMAP's seeded uniform
  initialization, because its multi-component spectral initialization was not
  repeatable. Their coordinates therefore differ from earlier runs. Connected
  graphs keep the spectral initialization, and the frozen seed-42 layouts are
  unchanged.
- Add `core_smiles_attached` with `[*]` attachment points for FMCS and selected
  coverage cores, and use it for the browser core drawing.
- Show inspector percentages as `> 99%` above 99.5% and `< 1%` below 0.5%.
  Display σ_doc, degrees of freedom and the number of multi-unit entries.
- Export evidence, offset, recall, group and core functions from
  `chemgridmap.__all__`; add a tested Python API example and document output
  columns.
- Record the package source-tree SHA-256, file count, Git commit and dirty
  state in run manifests.
- Report every offset document pair touching an inspector entry, fall back to a
  Kekulé matched member for aromatic bond-cut cores and export
  `core_smiles_attached_source`, and limit Git provenance to a tracked package
  directory.
- Correct the CA II tutorial: the main map is the activity layer, while class
  support is a separate layer.
- Add document-balanced activity class support, configurable merging of identical
  document values, and systematic document-pair offset diagnostics. Export pair
  statistics and annotate affected entries; show source-disagreement warnings in
  the inspector.
- Make UMAP a core dependency. Use UMAP, neighbor-refined grid assignment,
  coverage-target group cores, and Gaussian activity evidence by default. These
  defaults now differ from 0.4.0 in both the CLI, library and browser interface.
- Add generic activity-record CSV input through compound, pActivity, document,
  SMILES and optional relation columns. Reuse the ChEMBL evidence-unit, variance
  and class-support implementation.
- Add optional FMCS group cores with a configurable timeout and explicit
  machine-dependent timeout diagnostics.
- Write FMCS cores as valid SMILES extracted from matched members, so aromatic
  heterocycle cores no longer cause crashes.
- Expose `gaussian` and `exclude-flagged` activity evidence models. Gaussian
  remains the default after the multi-target model comparison.
- Add deterministic neighbor-guided grid refinement and per-entry tie-aware
  neighbor recall. CLI and browser enable refinement by default.
- Add coverage-target group cores from member scaffolds and scaffold-network ring
  fragments, retaining the most frequent scaffold coverage for comparison.
  Document and validate the CLI coverage target in [0.5, 1] (default 0.8).
- Export activity, class-support, neighbor-recall, and structure-group SVG/PNG
  maps at identical cell positions. Preserve the default activity renderer and
  keep neighbor reliability available without structure grouping.
  Browser builds export additional layers only as SVG and reuse the main map for
  observed activity. Library `map_color_formats` controls additional-layer
  formats; its default preserves CLI SVG/PNG exports.
- Record reliability settings, model diagnostics, and output hashes in run
  manifests. Include diagnostics and map layers in CLI/browser result ZIPs.
- Add regression, integration, frozen-result reproduction, and browser checks.
  The explicit legacy switches retain the 0.4.0 projection, layout and core rules.
- Fix relayout so a missing optional caption source produces a clear skip message
  instead of aborting the public figure build.
