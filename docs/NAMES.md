# Name map

The same concept uses different names at the CLI, Python API, run-manifest and
output-table boundaries.

| Concept | CLI | Python API | Run manifest | Output column |
|---|---|---|---|---|
| Support level | `--support-level` | `confidence_threshold` | `parameters.confidence_threshold`; `activity_confidence.confidence_threshold` | No direct column; it determines `class_status` |
| Class support | Implied by `--support-level` | `class_confidence` column returned by `annotate_activity_confidence` | `activity_confidence` records the model and threshold | `class_confidence` |
| Class reproducibility | No CLI option | `class_reproducibility` column returned by `annotate_activity_confidence` | `activity_confidence` records the model that generates it | `class_reproducibility` |
| σ_doc | Printed as `sigma_doc` | Fit field `tau`; `sigma_doc` is the persisted alias | `activity_confidence.tau` and `activity_confidence.sigma_doc` | No column; recorded in the manifest |
| Minimum shared entries | `--offset-min-shared` | `offset_min_shared` | `parameters.offset_min_shared`; `activity_confidence.offset_min_shared` and `document_offset_min_shared` | No column; `shared_entries` reports each pair's observed value |
| Median-offset threshold | `--offset-min-difference` | `offset_threshold` | `parameters.offset_threshold`; `activity_confidence.offset_threshold` and `document_offset_threshold` | No column; determines `systematic_offset` |
| Sign-agreement threshold | `--offset-sign-agreement` | `offset_sign_agreement` | `parameters.offset_sign_agreement`; `activity_confidence.offset_sign_agreement` and `document_offset_sign_agreement` | No column; determines `systematic_offset` |
| Offset display | `--offset-display {none,uncertain}` | `offset_display` | `parameters.offset_display`; `activity_confidence.offset_display` | Controls display use of `offset_pair_touched`, `status_reason` and `class_status` |
| Identity-conflict rule | `--identity-rule {inchikey,normalized-structure}` | `identity_rule` | `source.curation_parameters.identity_rule` | No direct column; changes conflict handling and related audit tables |
