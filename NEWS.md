# scProkaR 0.99.2

## Maintenance

* Add the maintainer's verified ORCID to the package author metadata.
* Replace blanket warning suppression with explicit handling of numeric
  labels, undefined correlations and unreachable graph paths. Preserve
  meaningful warnings from model fitting and Harmony backends.
* Extract workflow and clustering stages into internal helpers, retaining
  the public interfaces and processing metadata.
* Wrap vignette prose and align continuation lines consistently while
  preserving the rendered Markdown structure.

# scProkaR 0.99.1

## Fixes

* Resolve the legacy Harmony matrix function only when exported, avoiding
  missing-export warnings with current Harmony versions.
* Use `scrapper::modelGeneVariances()` for optional feature selection,
  ranking genes by residual variance. The updated trend fitter can change
  the selected features relative to the previous `scran` implementation.
* Apply the existing MNN cosine normalization explicitly before `fastMNN()`,
  avoiding its deprecated normalization dependency while preserving the
  normalization, feature subset and caller-supplied backend settings.
* Compute per-cell QC directly with the existing sparse-aware sums,
  avoiding the newly deprecated `scuttle::perCellQCMetrics()` call while
  preserving the metrics for single-cell and multi-cell objects.
* Use exact PCA for tiny matrices or large requested fractions of the
  singular values, avoiding the current `irlba` small-matrix error and
  high-rank warnings. Retain sparse iterative PCA for low-rank large data.
* Use R's default package citation until an associated publication is
  available, removing the custom citation without a publication DOI.

# scProkaR 0.99.0

## New

* `RegisterExistingIntegrationEmbeddings()` registers latent spaces already
  stored in an object, with optional copying and automatic name detection.
* `RegisterIntegrationEmbedding()` protects existing reductions from changes
  unless `overwrite = TRUE` is supplied.
* `BenchmarkIntegration()` discovers unregistered `integrated_*` reductions
  and accepts method and reduction names without regard to case.
* First submission to Bioconductor.
* `CreateBacObject()` builds a `SingleCellExperiment` from a count matrix, a
  10x Genomics directory, or a Seurat object, and `MergeBacObjects()` combines
  several such objects while keeping cell barcodes unique.
* `RunBacQC()` and `FilterBacCells()` compute and apply prokaryote-specific
  quality-control metrics, including ribosomal RNA fraction.
* `IntegrateBacData()`, `RegisterIntegrationEmbedding()`,
  `RunIntegratedUMAP()` and `RunIntegratedClustering()` provide batch
  integration and downstream embedding, with `BenchmarkIntegration()` scoring
  the corrected embeddings using scIB-inspired metrics.
* `aggregate_pseudobulk()`, `normalize_pseudobulk()` and the `run_edger_*()`
  functions provide pseudobulk differential expression built on edgeR,
  covering both pairwise contrasts and spline-based time courses.
* `run_tata()` implements Time Aware Trajectory Analysis (TATA), a
  time-informed trajectory abstraction for perturbation time courses, together
  with branch probability, gene trend, plotting and benchmarking helpers.
* The packaged `sce1` dataset provides a small bacterial single-cell example
  used throughout the vignettes.
