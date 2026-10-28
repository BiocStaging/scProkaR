#' Compute bacterial QC metrics
#'
#' Adds Seurat-style per-cell feature fractions to a `SingleCellExperiment`,
#' including rRNA and ribosomal-protein content.
#'
#' @param sce A `SingleCellExperiment`.
#' @param rrna_pattern Regular expression used to identify rRNA features.
#' @param ribo_pattern Regular expression used to identify ribosomal protein
#'   features.
#' @param gene_class_col Optional `rowData` column describing feature classes.
#'   When present, values matching `"rrna"` are used for rRNA metrics and values
#'   matching `"ribo"` or `"ribosomal_protein"` are used for ribosomal protein
#'   metrics.
#' @param store If `TRUE`, store metrics in `colData(sce)` and return the
#'   modified object. If `FALSE`, return the QC table only.
#'
#' @return A `SingleCellExperiment` when `store = TRUE`, otherwise a data frame
#'   of QC metrics.
#' @export
#'
#' @examples
#' set.seed(1)
#' genes <- c("rrsA", "rrlB", "rrfC", "rplA", "rpsB", paste0("gene", 1:25))
#' counts <- matrix(
#'     rpois(length(genes) * 40, lambda = 4),
#'     nrow = length(genes),
#'     dimnames = list(genes, paste0("cell", seq_len(40)))
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts)
#' )
#'
#' # rRNA and ribosomal protein features are matched on feature names
#' sce <- RunBacQC(sce)
#' head(SummarizedExperiment::colData(sce)[, c(
#'     "total_counts", "detected_features", "rrna_fraction", "pct_ribo"
#' )])
#'
#' # Return the QC table only, leaving the object untouched
#' qc <- RunBacQC(sce, store = FALSE)
#' summary(qc$rrna_fraction)
#'
#' # Feature classes stored in rowData() can be used instead of patterns
#' SummarizedExperiment::rowData(sce)$gene_class <- ifelse(
#'     rownames(sce) %in% c("rrsA", "rrlB", "rrfC", "gene1"),
#'     "rrna",
#'     "other"
#' )
#' qc_class <- RunBacQC(sce, gene_class_col = "gene_class", store = FALSE)
#' summary(qc_class$rrna_fraction)
RunBacQC <- function(
    sce,
    rrna_pattern = "^(rrs|rrl|rrf)",
    ribo_pattern = "^(rpl|rps)",
    gene_class_col = NULL,
    store = TRUE
) {
    counts <- assay(sce, "counts")
    if (is.null(counts)) {
        stop(
            "No `counts` assay found. Run CreateBacObject() first.",
            call. = FALSE
        )
    }

    row_data <- as.data.frame(rowData(sce))
    genes <- rownames(sce)
    rrna_index <- grepl(rrna_pattern, genes, ignore.case = TRUE)
    ribo_index <- grepl(ribo_pattern, genes, ignore.case = TRUE)

    if (!is.null(gene_class_col)) {
        if (!gene_class_col %in% colnames(row_data)) {
            stop(
                "`gene_class_col` was not found in rowData(sce).",
                call. = FALSE
            )
        }
        classes <- tolower(as.character(row_data[[gene_class_col]]))
        rrna_index <- classes %in% c("rrna", "r_rna", "ribosomal_rna")
        ribo_index <- classes %in% c(
            "ribo", "ribosomal_protein", "ribo_protein"
        )
    }

    qc_df <- .scprokar_per_cell_qc(sce, rrna_index, ribo_index)

    if (!isTRUE(store)) {
        return(qc_df)
    }

    col_data <- as.data.frame(colData(sce))
    col_data[colnames(qc_df)] <- qc_df
    colData(sce) <- DataFrame(col_data)

    .scprokar_store_step(
        sce,
        "qc",
        list(
            rrna_pattern = rrna_pattern,
            ribo_pattern = ribo_pattern,
            gene_class_col = gene_class_col,
            metrics = colnames(qc_df)
        )
    )
}

#' Filter cells using bacterial QC thresholds
#'
#' Applies simple threshold-based filtering to QC metrics produced by
#' `RunBacQC()`.
#'
#' @param sce A `SingleCellExperiment` with QC metrics already present.
#' @param min_counts Minimum allowed total counts.
#' @param min_features Minimum allowed number of detected features.
#' @param max_rrna_fraction Maximum allowed rRNA fraction.
#' @param custom_filter Optional custom filter. May be a logical vector, a
#'   column name in `colData(sce)`, or a function returning a logical vector.
#'
#' @return A filtered `SingleCellExperiment`.
#' @export
#'
#' @examples
#' set.seed(1)
#' genes <- c("rrsA", "rrlB", "rplA", "rpsB", paste0("gene", 1:26))
#' counts <- matrix(
#'     rpois(length(genes) * 50, lambda = 3),
#'     nrow = length(genes),
#'     dimnames = list(genes, paste0("cell", seq_len(50)))
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts)
#' )
#' sce <- RunBacQC(sce)
#'
#' filtered <- FilterBacCells(
#'     sce,
#'     min_counts = 80,
#'     min_features = 20,
#'     max_rrna_fraction = 0.2
#' )
#' c(before = ncol(sce), after = ncol(filtered))
#'
#' # A custom filter may be given as a function of colData()
#' custom <- FilterBacCells(
#'     sce,
#'     custom_filter = function(cd) {
#'         cd$detected_features >= stats::median(cd$detected_features)
#'     }
#' )
#' ncol(custom)
FilterBacCells <- function(
    sce,
    min_counts = NULL,
    min_features = NULL,
    max_rrna_fraction = NULL,
    custom_filter = NULL
) {
    cd <- as.data.frame(colData(sce))
    required <- c("total_counts", "detected_features", "rrna_fraction")
    missing_cols <- setdiff(required, colnames(cd))
    if (length(missing_cols) > 0) {
        stop("Run RunBacQC() before filtering cells.", call. = FALSE)
    }

    keep <- rep(TRUE, nrow(cd))
    if (!is.null(min_counts)) {
        keep <- keep & (cd$total_counts >= min_counts)
    }
    if (!is.null(min_features)) {
        keep <- keep & (cd$detected_features >= min_features)
    }
    if (!is.null(max_rrna_fraction)) {
        keep <- keep & (cd$rrna_fraction <= max_rrna_fraction)
    }

    if (!is.null(custom_filter)) {
        if (is.function(custom_filter)) {
            custom_keep <- custom_filter(cd)
        } else if (is.character(custom_filter) && length(custom_filter) == 1L) {
            custom_keep <- cd[[custom_filter]]
        } else {
            custom_keep <- custom_filter
        }
        if (!is.logical(custom_keep) || length(custom_keep) != nrow(cd)) {
            stop(
                "`custom_filter` must resolve to one logical value per cell.",
                call. = FALSE
            )
        }
        keep <- keep & custom_keep
    }

    filtered <- sce[, keep, drop = FALSE]
    .scprokar_store_step(
        filtered,
        "filtering",
        list(
            min_counts = min_counts,
            min_features = min_features,
            max_rrna_fraction = max_rrna_fraction,
            kept_cells = sum(keep),
            removed_cells = sum(!keep)
        )
    )
}

#' Per-cell quality-control metrics
#'
#' Counts library size, detected features, and the share of counts falling on
#' the rRNA and ribosomal-protein feature sets directly from the counts.
#' The same sparse-aware calculation handles any number of cells.
#'
#' @param sce A `SingleCellExperiment` with a `counts` assay.
#' @param rrna_index,ribo_index Logical vectors over `rownames(sce)`.
#'
#' @return A data frame of QC metrics, one row per cell.
#' @keywords internal
#' @noRd
.scprokar_per_cell_qc <- function(sce, rrna_index, ribo_index) {
    counts <- assay(sce, "counts")

    totals <- as.numeric(Matrix::colSums(counts))
    detected <- as.numeric(Matrix::colSums(counts > 0))
    subset_sum <- function(index) {
        if (!any(index)) {
            return(rep(0, ncol(counts)))
        }
        as.numeric(Matrix::colSums(counts[index, , drop = FALSE]))
    }
    rrna_counts <- subset_sum(rrna_index)
    ribo_counts <- subset_sum(ribo_index)

    ## A cell with an empty library has no meaningful share; report 0 rather
    ## than the NaN that dividing by zero would give.
    share <- function(selected) selected / pmax(totals, 1)
    rrna_fraction <- share(rrna_counts)
    ribo_fraction <- share(ribo_counts)

    data.frame(
        total_counts = totals,
        detected_features = detected,
        rrna_counts = rrna_counts,
        rrna_fraction = rrna_fraction,
        pct_rrna = rrna_fraction * 100,
        ribo_counts = ribo_counts,
        ribo_fraction = ribo_fraction,
        pct_ribo = ribo_fraction * 100,
        row.names = colnames(sce)
    )
}
