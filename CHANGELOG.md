# Changelog

## 0.5.0 — 2026-10-29

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
