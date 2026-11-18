# Changelog

## 0.6.0 — 2026-11-18

Changes made in response to the simulated review of the v0.6 manuscript:

- Inspector: show class reproducibility (the quantity validated against later
  documents) next to class support, each labelled by the displayed class, and
  the agreement reached by entries of the same class and status in the
  evaluation (manuscript Table 4). The support-layer legend repeats those rates.
- Inspector: list the bound records (`>`, `<`, `≥`, `≤`) that curation excluded
  for the selected entry, with the corresponding pActivity bound
  (`excluded_bound_records`).
- Neighbor evidence: descriptive labels ("Consistent evidence", "Mixed
  evidence", "Insufficient evidence"); a visible statement with the measured
  agreement rates of the time-split evaluation in the legend and the query
  panel; "Same class" / "Different class (possible activity cliff, assay
  difference or data issue)" in place of agrees / disagrees; a note column that
  marks neighbors with the same fingerprint (Tanimoto 1).
- Neighbor-evidence layer: consistent categories get a dark outer ring and a
  coloured inner frame, so that a frame stays visible on a fill of the same
  hue; reference molecules are shown in grey rather than pale class tints; the
  legend states that the fill is the new molecule's own measured class.
- Inserted molecules: `placement_recall`, `placement_recall_shared` and
  `placement_recall_fingerprint_count` in `new_molecules.csv` and the inspector.
- Map: **Emphasize class** fades molecules of the other classes, a cue that does
  not depend on hue. `paper/cvd_separations.py` reports every pairwise colour
  separation under simulated protanopia, deuteranopia and tritanopia.
- Documentation: `class_confidence` is documented as class support and
  distinguished from `class_reproducibility`; guidance for batches and dates as
  evidence units; supported browsers and maintenance.
- Dependencies: the declared minimum versions are now the lowest set that the
  test suite passes with on Python 3.10 (RDKit 2025.3.1, NumPy 1.26, pandas
  2.1, SciPy 1.11.1, scikit-learn 1.3, Matplotlib 3.8, umap-learn 0.5.5). The
  earlier bounds were declared but not tested: older RDKit releases reject the
  distance matrix that structural grouping passes to Butina clustering, and
  umap-learn 0.5.3 and 0.5.4 import `pkg_resources`, which current setuptools no
  longer provides.
- CLI: `--value-unit p --unit-column ...` is rejected as a conflict on every
  supported Python version (Python 3.10 and 3.11 had accepted it because the
  value equalled the default).
- CI: a job with the lowest declared dependency versions on Python 3.10.
- macOS: RDKit 2026.9.2 is excluded, because its macOS wheels raise
  "bad any cast" in `Chem.MurckoDecompose`, which every scaffold and core
  calculation uses (RDKit 2026.9.1 and the Linux and Windows wheels of 2026.9.2
  are not affected).
- CLI: the `--unit-column` help is plain ASCII, so `chemgridmap --help` works
  in Windows consoles that cannot print the Greek letter mu.
- Paper: post hoc analyses for the v0.7 manuscript (`paper/rev07/`) and the
  v0.7 figure sources (Figures 1, 4, 6 and 7, and the Additional file 1
  hERG-pair figure; `make paper-figure1-v07 paper-figures-v07
  paper-figure7-v07`).

Earlier changes since 0.5.0:

- Add the optional `coverage_recap` core rule: ring-containing RECAP fragments
  (dummy atoms removed; members above 80 heavy atoms not decomposed) join the
  scaffold and ring-fragment candidates. Show the generic Murcko scaffold and the
  members sharing it on the group card and as a supplementary last row of the core
  comparison table and CSV. Existing rules, defaults and exports are unchanged.
- Add the analysis scripts of the E28 evaluation (plan hashed before the run): the time-split test of
  neighbor evidence (`paper/u10_neighbor_validation.py`, `paper/u10_e28_report.py`),
  the core-rule comparison on the sampled targets (`paper/u10_core_rules.py`,
  `paper/u10_e28c_summary.py`) and the builder of the Additional file 1 addendum
  (`paper/u10_e28_si.py`).

- Compute independent non-FMCS group cores with shared ordered loky workers in
  browser builds and rule switches; add per-second progress, `CHEMGRIDMAP_JOBS`
  and opt-in library/CLI `jobs` while retaining serial defaults outside the app.
  Keep FMCS serial and retain existing scientific output bytes.
- Save compressed, versioned and hash-verified session caches outside existing
  exports/ZIPs. Open the latest compatible local results without recomputation,
  restoring original settings, core caches and active inserted molecules,
  including manual additions and measured feedback. Restrict deserialization and
  reject incompatible, modified, incomplete or unsafe sessions before applying them.
- Improve browser disagreement diamonds with zoom-aware size, bilingual measured
  versus neighbor-evidence tooltips, and an overlay above selection outlines.
  Keep the legend text and exported SVG unchanged. Add bilingual restore controls
  and help, plus Node VM coverage for progress, lists, restoration and markers.

- Add a bilingual one-molecule form with SMILES preview, automatic IDs and
  append-only placement in the fixed reference grid. Record, edit and remove
  manual measured results; retain imported measurements as read-only. Show
  neighbor-evidence disagreements on every color layer and in a filterable
  list, with manual-entry and feedback CSVs, manifest hashes and updated ZIPs.
  Preserve existing reference and file-only insertion outputs.
- Evaluate coverage cores by descending heavy-atom count with an early failure
  cutoff, completing tied candidates at the first qualifying size. Preserve
  the original exhaustive fallback, candidate generation and selection rules.

- Add retained-source highlighting with distinct mapped-molecule counts,
  source unions and separate inserted-molecule batches on all map layers.
- Add the optional `coverage_brics` core rule, cached map-page rule switches,
  and per-group comparisons with a fixed 10-second FMCS budget. Keep reference
  groups, coordinates and existing exports unchanged; save rule outputs under
  `core_rules/`, with bilingual controls and help.

- Fix reference SVG display after inserting new molecules, bilingual control
  labels and modal help. Show the selected new molecule’s similar references
  and placement link, clarify signed comparisons and export excluded
  measurement records with their original columns.

- Build the reference map first, then insert up to 1,000 new molecules into
  existing empty cells without changing reference outputs. Add browser imports
  for measured records, neighbor-evidence layers and lists, reference-match
  comparisons, bilingual help and independent versioned exports; retain CLI
  query arguments with the new insertion workflow. No activity is predicted.

- Add optional local source lookup from a user-supplied ChEMBL activity export
  for generic records, with parent-structure matches, same-pChEMBL marks and
  bilingual possible-source and audit notes; preserve calculations and exports.
- Resolve unmeasured IDs through existing structure matches to measured sources
  or the first duplicate; add bilingual inspector notes and a source-count table
  in Data audit without changing calculations or exports.
- Add inspector source-document links from the bundled ChEMBL 37 document table,
  optional DOI/PubMed/URL references for generic records, link-coverage audit
  summaries and document links in inspector record downloads; preserve map
  calculations and default exports.
- Add English/Chinese browser switching and question-mark setting explanations;
  keep exports in English, correct Chinese punctuation and include the language
  dictionary and document data in installed packages.
- List the retained records, record CSV and offset-pair notes for generic
  activity records in the browser inspector; show their excluded records and
  the evaluated evidence-unit pairs, flagged pairs first, in Data audit.
- Reload a CSV chosen again, name the wide-range marker in the class-support
  legend, and write group cards as measured members plus unmeasured molecules.
- Add browser uploads for unmeasured CSV/TSV molecules, measured/unmeasured
  filtering, and a query inspector linking to measured Morgan neighbors.
- Expose activity units, evidence-unit columns and explicit assumed σ in the
  browser, with query-check counts and the σ source in Data audit.
- Include all three activity classes and the actual thresholds before
  `Not measured` in static activity legends when queries are present.
- Highlight a selected structure group across browser color layers and views,
  with a core structure and coverage card and an inspector highlight toggle.
- Convert molar concentration units in generic record and molecule tables,
  preserve original values/units, and report exclusion reasons.
- Add explicit assumed-sigma class support for ChEMBL and generic records,
  including one evidence unit per compound when document IDs are unavailable.
- Map unmeasured query molecules with measured references, export input checks
  and tie-inclusive Morgan neighbors, and keep query activity/support blank.
- Preserve reference evidence, groups and cores when queries are added; assign
  queries to eligible reference centers and count query members separately.
- Render unmeasured activity/support cells with a white fill and gray border;
  keep structural diagnostics over all molecules and activity diagnostics over
  measured references. Document the public API, CLI and optional audit outputs.

## 0.5.0 — 2026-10-29

- Make `offset_display="none"` the default in the Python API, CLI and browser.
  Offset document pairs are still reported, exported and marked in the support
  layer, while display status follows support alone; use
  `--offset-display uncertain` to opt into the previous downgrade. In the
  held-out re-evaluation that kept held-out documents out of flag construction,
  the downgrade did not increase the agreement of supported cells with the
  held-out literature.
- Export every evaluated document pair to `document_offsets.csv`, including
  unflagged pairs with `systematic_offset=false`; the manifest and inspector
  continue to use only flagged pairs for offset diagnostics and display.
- Emit an explicit warning when fewer than five multi-unit entries require the
  literature fallback for σ_doc, and omit the inapplicable dominant-entry
  proportion warning in that case.
- Add `chemgridmap --example caii`, which runs the packaged CHEMBL205 dataset
  without a source-tree path and is now the README quick-start command.
- Report the activity evidence entry that contributes the largest share of
  σ_doc degrees of freedom and pooled sum of squares in fits, run manifests and
  CLI summaries, with a warning above 30%. The σ_doc estimate is unchanged.
- Show every CLI option default in help, mention the `inspect` subcommand in
  usage, add the bundled CA II quick-start command, and export
  `prepare_generic_activity_records` as documented public API.
- Compare parent identities by the standard InChIKey after removing explicit
  hydrogens. Keep same-key tautomer, resonance and de-deuterated writing
  variants as one deterministic entry, warn about them, and export
  `structure_variants.csv`; reserve identity errors for different constitutions
  or stereochemistry. Report the package's own version in manifests.
- Make the support level and offset display rule configurable in the Python
  API, CLI and browser. Add configurable minimum shared entries, median
  log-unit difference and sign agreement; validate each value and record all
  five settings in fits, `app_options.json` and run manifests. Record the union
  of entries touched by each flagged document pair and show it in the
  inspector.
- Allow offset display to be disabled without discarding offset diagnostics:
  touched entries remain marked while `class_status` follows support alone.
  Inspector and CLI output distinguish attached markers from display effect.
- Render reliable layers with opaque cells. Full-detail structure cells now use
  an 8% colored frame, a white inset and existing molecule margin; overview
  cells use the full color. Use the same five-step OKLab-checked recall palette,
  corner-marker geometry and Morgan-fingerprint terminology in static and
  browser exports.
- Add rendered-color tests using sRGB-to-OKLab distance and the Machado et al.
  (2009) severity-1.0 deuteranopia simulation.
- Trim `structure_members.csv` to identifiers, grid positions, groups,
  scaffolds and core membership; add a generated complete column dictionary
  and public API reference.
- Check the 10,000-entry exact-grouping limit before projection and assignment.
  Add `identity_conflicts="exclude"` and the matching CLI/browser option to
  remove conflicting parent identities, export the conflicts and report the
  counts.
- Keep FMCS core search non-chiral, use the same chirality-aware containment
  check as the other core rules, and record numba, pynndescent and llvmlite
  versions in run manifests.
- Add `neighbor_kth_similarity` for the kth Morgan fingerprint neighbor and
  `neighbor_kth_distance` for other representations, without changing recall.
- Add bundled CA II UMAP and exact PCA regression tests and document the tested
  dependency/platform matrix and Figure 2 renderer command.
- Keep selected-cell and neighbor outlines at a fixed on-screen width while
  zooming in the browser, and list the offset marker in the legend only when an
  offset document pair touches a displayed entry.
- When `offset_display=uncertain`, an entry touched by a flagged document pair
  is uncertain and pale, with `status_reason=offset_pair`, while all model
  probabilities remain unchanged. Record touched-entry counts in fits and
  manifests, show document-pair details and a distinct bottom-left map marker.
- Rename the user-facing "conflict" flag to "wide record range", retain
  `is_conflicted` for compatibility, and add the equal-valued
  `wide_record_range` alias. Present this note in the inspector with a neutral
  style rather than a warning style.
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
