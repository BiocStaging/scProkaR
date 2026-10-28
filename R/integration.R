#' Integrate microbial single-cell batches
#'
#' Runs one of several integration backends and stores the corrected embedding
#' in `reducedDims(sce)`.
#'
#' @param sce A `SingleCellExperiment`.
#' @param batch_col Column in `colData(sce)` containing batch assignments.
#' @param method Integration backend: `"mnn"` or `"harmony"`.
#' @param feature_set Either `"hvg"` or `"all"`.
#' @param dims Integer vector of dimensions to retain.
#' @param assay_name Assay used for preprocessing.
#' @param integrated_name Optional name for the output reduced dimension.
#' @param ... Additional backend-specific arguments.
#'
#' @return A `SingleCellExperiment` with an integrated reduced dimension.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 300, n_features = 60, n_pcs = 10
#' )
#' # Assign a technical batch label to integrate over.
#' sce$batch <- rep(c("batch1", "batch2"), length.out = ncol(sce))
#'
#' if (requireNamespace("batchelor", quietly = TRUE)) {
#'     sce <- IntegrateBacData(
#'         sce,
#'         batch_col = "batch",
#'         method = "mnn",
#'         dims = 1:8
#'     )
#'     print(SingleCellExperiment::reducedDimNames(sce))
#' }
#'
#' if (requireNamespace("harmony", quietly = TRUE)) {
#'     sce <- IntegrateBacData(
#'         sce,
#'         batch_col = "batch",
#'         method = "harmony",
#'         dims = 1:8
#'     )
#'     harmony_emb <- SingleCellExperiment::reducedDim(
#'         sce, "integrated_harmony"
#'     )
#'     print(dim(harmony_emb))
#' }
IntegrateBacData <- function(
    sce,
    batch_col,
    method = c("mnn", "harmony"),
    feature_set = c("hvg", "all"),
    dims = seq_len(30),
    assay_name = "counts",
    integrated_name = NULL,
    ...
) {
    method <- match.arg(method)
    feature_set <- match.arg(feature_set)
    .scprokar_match_columns(sce, batch_col, label = "batch")

    sce <- .scprokar_prepare_reduction(
        sce,
        assay_name = assay_name,
        feature_set = feature_set,
        dims = dims,
        reduction_name = "PCA"
    )

    if (is.null(integrated_name)) {
        integrated_name <- paste0("integrated_", method)
    }

    batch <- as.factor(colData(sce)[[batch_col]])
    features <- .scprokar_select_features(sce, feature_set = feature_set)

    if (method == "mnn") {
        .scprokar_require("batchelor", "MNN integration")
        split_idx <- split(seq_len(ncol(sce)), batch)
        split_sce <- lapply(split_idx, function(idx) sce[, idx, drop = FALSE])
        mnn_args <- list(...)
        min_batch_cells <- min(vapply(split_idx, length, integer(1)))
        if (!"k" %in% names(mnn_args)) {
            mnn_args$k <- max(1L, min(20L, min_batch_cells - 1L))
        }
        if (!"BSPARAM" %in% names(mnn_args) &&
                requireNamespace("BiocSingular", quietly = TRUE)) {
            mnn_args$BSPARAM <- BiocSingular::ExactParam()
        }
        prepared <- .scprokar_prepare_mnn_batches(
            split_sce, features = features, args = mnn_args
        )
        corrected <- do.call(
            batchelor::fastMNN,
            c(
                prepared$batches,
                list(subset.row = features, d = max(dims)),
                prepared$args
            )
        )
        embedding <- reducedDim(corrected, "corrected")[, dims, drop = FALSE]
        rownames(embedding) <- colnames(corrected)
        embedding <- embedding[colnames(sce), , drop = FALSE]
    } else if (method == "harmony") {
        .scprokar_require("harmony", "Harmony integration")
        embedding <- .scprokar_run_harmony(
            reducedDim(sce, "PCA")[, dims, drop = FALSE],
            meta_data = as.data.frame(colData(sce)),
            batch_col = batch_col,
            ...
        )
    }

    rownames(embedding) <- colnames(sce)
    colnames(embedding) <- paste0(
        toupper(method), "_", seq_len(ncol(embedding))
    )
    reducedDim(sce, integrated_name) <- embedding

    meta <- .scprokar_get_metadata(sce)
    if (is.null(meta$integration)) {
        meta$integration <- list(results = list())
    }
    meta$integration$results[[method]] <- list(
        reduction = integrated_name,
        batch_col = batch_col,
        dims = dims,
        feature_set = feature_set,
        n_features = length(features)
    )
    .scprokar_set_metadata(sce, meta)
}

#' Run UMAP on an existing embedding
#'
#' Computes a two-dimensional UMAP from a stored reduced dimension, which is
#' useful for inspecting integrated embeddings returned by
#' [IntegrateBacData()] or [RegisterIntegrationEmbedding()].
#'
#' @param sce A `SingleCellExperiment`.
#' @param reduction Existing reduced-dimension name to use as UMAP input.
#' @param umap_name Optional name for the output UMAP reduced dimension.
#' @param n_neighbors Number of neighbours passed to `uwot::umap()`.
#' @param min_dist UMAP `min_dist`.
#' @param metric Distance metric passed to `uwot::umap()`.
#' @param ... Additional arguments passed to `uwot::umap()`.
#'
#' @return A `SingleCellExperiment` with the computed UMAP stored in
#'   `reducedDims(sce)`.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#'
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'     sce <- RunIntegratedUMAP(
#'         sce,
#'         reduction = "PCA",
#'         umap_name = "umap_pca",
#'         n_neighbors = 15
#'     )
#'     print(SingleCellExperiment::reducedDimNames(sce))
#'     print(head(SingleCellExperiment::reducedDim(sce, "umap_pca"), 3))
#' }
RunIntegratedUMAP <- function(
    sce,
    reduction,
    umap_name = NULL,
    n_neighbors = 30,
    min_dist = 0.3,
    metric = "cosine",
    ...
) {
    .scprokar_require("uwot", "RunIntegratedUMAP()")

    reduction_name <- .scprokar_resolve_reduction_name(sce, reduction)
    emb <- .scprokar_reduction_embedding(sce, reduction_name, "run UMAP")
    if (nrow(emb) < 3) {
        stop(
            "Need at least 3 cells to run UMAP on reduction '",
            reduction_name, "'.", call. = FALSE
        )
    }

    if (is.null(umap_name)) {
        umap_name <- paste0("umap_", reduction_name)
    }
    n_neighbors <- min(n_neighbors, nrow(emb) - 1L)

    umap <- uwot::umap(
        emb,
        n_neighbors = n_neighbors,
        min_dist = min_dist,
        metric = metric,
        verbose = FALSE,
        ...
    )
    rownames(umap) <- colnames(sce)
    colnames(umap) <- c("UMAP_1", "UMAP_2")
    reducedDim(sce, umap_name) <- umap

    meta <- .scprokar_get_metadata(sce)
    meta$visualization <- c(
        meta$visualization,
        list(
            last_umap = list(
                input_reduction = reduction_name,
                output_reduction = umap_name,
                n_neighbors = n_neighbors,
                min_dist = min_dist,
                metric = metric
            )
        )
    )
    .scprokar_set_metadata(sce, meta)
}

#' Cluster cells on an integrated embedding
#'
#' Builds a k-nearest-neighbour graph on a corrected embedding and detects
#' communities for downstream trajectory and pseudobulk workflows.
#'
#' @param sce A `SingleCellExperiment`.
#' @param reduction Integrated reduced-dimension name. If `NULL`, the function
#'   prefers a stored `integrated_*` embedding.
#' @param cluster_col Column name used to store cluster labels in
#'   `colData(sce)`.
#' @param dims Optional dimensions to retain from the embedding.
#' @param k Number of nearest neighbours used to build the graph.
#' @param algorithm Graph clustering algorithm: `"louvain"` or `"walktrap"`.
#' @param resolution Cluster granularity parameter. This is analogous to
#'   Seurat's `resolution`: higher values typically produce more clusters. It is
#'   used when the selected graph clustering algorithm supports it.
#'
#' @return A `SingleCellExperiment` with cluster labels added to `colData(sce)`.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 300, n_features = 60, n_pcs = 10
#' )
#' sce <- RunIntegratedClustering(
#'     sce,
#'     reduction = "PCA",
#'     cluster_col = "pca_clusters",
#'     dims = 1:8,
#'     k = 15,
#'     algorithm = "louvain",
#'     resolution = 0.8
#' )
#' print(table(sce$pca_clusters))
#' print(table(sce$pca_clusters, sce$condition))
RunIntegratedClustering <- function(
    sce,
    reduction = NULL,
    cluster_col = "integrated_clusters",
    dims = NULL,
    k = 20,
    algorithm = c("louvain", "walktrap", "leiden"),
    resolution = 0.8
) {
    algorithm <- match.arg(algorithm)

    prepared <- .scprokar_clustering_embedding(sce, reduction, dims)
    reduction_name <- prepared$reduction
    emb <- prepared$embedding

    k <- max(1L, min(as.integer(k), nrow(emb) - 1L))
    graph_data <- .scprokar_build_cluster_graph(emb, k = k)
    membership <- .scprokar_cluster_membership(
        graph_data$graph, algorithm, resolution, n_cells = nrow(emb)
    )

    colData(sce)[[cluster_col]] <- factor(membership)

    meta <- .scprokar_get_metadata(sce)
    meta$integration$clustering <- list(
        reduction = reduction_name,
        cluster_col = cluster_col,
        dims = if (is.null(dims)) seq_len(ncol(emb)) else dims,
        k = k,
        algorithm = algorithm,
        resolution = resolution
    )
    .scprokar_set_metadata(sce, meta)
}

#' Plot any stored reduced dimension
#'
#' Creates a quick scatter plot for an existing reduced dimension such as an
#' original Seurat UMAP, an integrated latent space, or a UMAP computed by
#' [RunIntegratedUMAP()].
#'
#' @param sce A `SingleCellExperiment`.
#' @param reduction Reduced-dimension name to plot.
#' @param colour_by Optional `colData(sce)` column used for colouring points.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @param palette Optional colour palette. For discrete variables, this
#'   should be a character vector of colours; for continuous variables, it
#'   should be a gradient vector.
#' @param facet_by Optional `colData(sce)` column used to facet the plot.
#' @param shuffle Whether to shuffle plotting order before drawing points.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 300, n_features = 60, n_pcs = 10
#' )
#' plot <- PlotReduction(
#'     sce,
#'     reduction = "PCA",
#'     colour_by = "condition",
#'     point_size = 0.6,
#'     shuffle = TRUE
#' )
#' print(class(plot))
#' print(plot$labels[c("x", "y")])
#'
#' # Facet the same embedding by sampling time.
#' faceted <- PlotReduction(
#'     sce,
#'     reduction = "PCA",
#'     colour_by = "condition",
#'     facet_by = "timepoint"
#' )
#' print(nrow(faceted$data))
PlotReduction <- function(
    sce,
    reduction,
    colour_by = NULL,
    point_size = 0.4,
    point_alpha = 0.8,
    palette = NULL,
    facet_by = NULL,
    shuffle = FALSE
) {
    reduction_name <- .scprokar_resolve_reduction_name(sce, reduction)
    emb <- .scprokar_reduction_embedding(sce, reduction_name, "plot")

    plot_df <- .scprokar_reduction_plot_data(
        sce,
        emb = emb,
        colour_by = colour_by,
        facet_by = facet_by,
        shuffle = shuffle
    )

    plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dim1, y = dim2)) +
        ggplot2::geom_point(
            mapping = ggplot2::aes(color = .colour),
            size = point_size,
            alpha = point_alpha
        ) +
        ggplot2::labs(
            title = reduction_name,
            x = colnames(emb)[1],
            y = colnames(emb)[2],
            color = if (is.null(colour_by)) NULL else colour_by
        ) +
        ggplot2::theme_classic()

    plot <- plot + .scprokar_reduction_colour_scale(
        colour_values = plot_df$.colour,
        colour_by = colour_by,
        palette = palette
    )

    if (!is.null(facet_by)) {
        plot <- plot + ggplot2::facet_wrap(~.facet)
    }

    plot
}

#' Plot Seurat-style integration result panels
#'
#' Creates a standard set of integration diagnostic plots inspired by common
#' Seurat workflows: unintegrated and integrated views coloured by batch and,
#' optionally, a biological label, plus a split-panel view of the integrated
#' embedding.
#'
#' @param sce A `SingleCellExperiment`.
#' @param unintegrated_reduction Reduction name for the pre-integration view.
#' @param integrated_reduction Reduction name for the integrated view. If
#'   `NULL`, the function prefers a stored `umap_integrated_*` reduction and
#'   otherwise falls back to `integrated_*`.
#' @param batch_col Batch column in `colData(sce)`.
#' @param label_col Optional biological label column in `colData(sce)`.
#' @param split_by Optional metadata column used to facet the integrated plot.
#'   This is the `scProkaR` equivalent of Seurat's `split.by`.
#' @param point_size Point size passed to [PlotReduction()].
#' @param point_alpha Point alpha passed to [PlotReduction()].
#' @param shuffle Whether to shuffle plotting order.
#'
#' @return A named list of `ggplot2` objects.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 300, n_features = 60, n_pcs = 10
#' )
#' sce$batch <- rep(c("batch1", "batch2"), length.out = ncol(sce))
#'
#' # Any corrected embedding can stand in for the integrated view; here a
#' # toy correction is registered so the example needs no extra backend.
#' corrected <- SingleCellExperiment::reducedDim(sce, "PCA")[, 1:2]
#' is_batch2 <- sce$batch == "batch2"
#' corrected[is_batch2, ] <- corrected[is_batch2, ] * 0.9
#' sce <- RegisterIntegrationEmbedding(
#'     sce,
#'     embedding = corrected,
#'     method_name = "toy"
#' )
#'
#' plots <- PlotIntegrationOverview(
#'     sce,
#'     unintegrated_reduction = "PCA",
#'     integrated_reduction = "integrated_toy",
#'     batch_col = "batch",
#'     label_col = "condition",
#'     split_by = "batch"
#' )
#' print(names(plots))
#' print(class(plots$integrated_by_batch))
PlotIntegrationOverview <- function(
    sce,
    unintegrated_reduction = "umap",
    integrated_reduction = NULL,
    batch_col,
    label_col = NULL,
    split_by = NULL,
    point_size = 0.4,
    point_alpha = 0.8,
    shuffle = FALSE
) {
    .scprokar_match_columns(sce, batch_col, label = "batch")
    if (!is.null(label_col)) {
        .scprokar_match_columns(sce, label_col, label = "label")
    }
    if (!is.null(split_by)) {
        .scprokar_match_columns(sce, split_by, label = "split")
    }

    integrated_reduction <- .scprokar_default_integration_reduction(
        sce,
        integrated_reduction = integrated_reduction
    )

    plots <- list(
        unintegrated_by_batch = PlotReduction(
            sce,
            reduction = unintegrated_reduction,
            colour_by = batch_col,
            point_size = point_size,
            point_alpha = point_alpha,
            shuffle = shuffle
        ),
        integrated_by_batch = PlotReduction(
            sce,
            reduction = integrated_reduction,
            colour_by = batch_col,
            point_size = point_size,
            point_alpha = point_alpha,
            shuffle = shuffle
        )
    )

    if (!is.null(label_col)) {
        plots$unintegrated_by_label <- PlotReduction(
            sce,
            reduction = unintegrated_reduction,
            colour_by = label_col,
            point_size = point_size,
            point_alpha = point_alpha,
            shuffle = shuffle
        )
        plots$integrated_by_label <- PlotReduction(
            sce,
            reduction = integrated_reduction,
            colour_by = label_col,
            point_size = point_size,
            point_alpha = point_alpha,
            shuffle = shuffle
        )
    }

    if (!is.null(split_by)) {
        split_colour <- if (is.null(label_col)) batch_col else label_col
        plots$integrated_split <- PlotReduction(
            sce,
            reduction = integrated_reduction,
            colour_by = split_colour,
            facet_by = split_by,
            point_size = point_size,
            point_alpha = point_alpha,
            shuffle = shuffle
        )
    }

    plots
}

#' Register an external embedding for benchmarking and downstream analysis
#'
#' Adds a precomputed embedding to `reducedDims(sce)` and records it in
#' `metadata(sce)$scProkaR$integration`, allowing custom integration outputs to
#' be benchmarked alongside native scProkaR methods.
#'
#' @param sce A `SingleCellExperiment`.
#' @param embedding A matrix-like cell embedding with one row per cell.
#' @param method_name Name used when recording the embedding.
#' @param reduction_name Optional reducedDim name. Defaults to
#'   `paste0("integrated_", method_name)`.
#' @param metadata Optional named list describing the external method.
#' @param overwrite Logical indicating whether to replace an existing reduced
#'   dimension with different values. Defaults to `FALSE`; registering the same
#'   embedding again is allowed.
#'
#' @return A `SingleCellExperiment` with the external embedding registered.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#'
#' # Stand-in for an embedding produced by an external integration tool.
#' external <- SingleCellExperiment::reducedDim(sce, "PCA")[, 1:2]
#' colnames(external) <- c("EXT_1", "EXT_2")
#'
#' sce <- RegisterIntegrationEmbedding(
#'     sce,
#'     embedding = external,
#'     method_name = "my_method",
#'     metadata = list(software = "external_tool", version = "1.0")
#' )
#' print(SingleCellExperiment::reducedDimNames(sce))
#' record <- S4Vectors::metadata(sce)$scProkaR$integration$results
#' print(record$my_method)
RegisterIntegrationEmbedding <- function(
    sce,
    embedding,
    method_name,
    reduction_name = NULL,
    metadata = list(),
    overwrite = FALSE
) {
    if (is.null(reduction_name)) {
        reduction_name <- paste0("integrated_", method_name)
    }

    embedding <- as.matrix(embedding)
    if (nrow(embedding) != ncol(sce)) {
        stop("`embedding` must have one row per cell in `sce`.", call. = FALSE)
    }

    if (!is.null(rownames(embedding))) {
        missing_cells <- setdiff(colnames(sce), rownames(embedding))
        if (length(missing_cells) > 0) {
            stop(
                "`embedding` rownames are missing cells: ",
                paste(missing_cells, collapse = ", "),
                call. = FALSE
            )
        }
        embedding <- embedding[colnames(sce), , drop = FALSE]
    } else {
        rownames(embedding) <- colnames(sce)
    }

    available <- reducedDimNames(sce)
    if (reduction_name %in% available && !isTRUE(overwrite)) {
        existing <- as.matrix(reducedDim(sce, reduction_name))
        same_embedding <- isTRUE(all.equal(
            existing, embedding, check.attributes = FALSE, tolerance = 0
        ))
        if (!same_embedding) {
            stop(
                "Reduction '", reduction_name, "' already exists. ",
                "Use `overwrite = TRUE` or choose a different ",
                "`reduction_name`.",
                call. = FALSE
            )
        }
    }

    reducedDim(sce, reduction_name) <- embedding
    ## Use [[ ]] rather than $: $ partial-matches, so a user key such as
    ## `source_file` would silently be promoted into the source field.
    embedding_source <- if (!is.null(metadata[["source"]])) {
        metadata[["source"]]
    } else {
        "external"
    }
    metadata[["source"]] <- NULL

    meta <- .scprokar_get_metadata(sce)
    if (is.null(meta$integration)) {
        meta$integration <- list(results = list())
    }
    meta$integration$results[[method_name]] <- c(
        list(
            reduction = reduction_name,
            source = embedding_source,
            dims = seq_len(ncol(embedding))
        ),
        metadata
    )
    .scprokar_set_metadata(sce, meta)
}

#' Register existing reduced dimensions as integration embeddings
#'
#' Registers embeddings already stored in a `SingleCellExperiment`, including
#' Seurat CCA/RPCA, Scanorama, scVI, SCALEX, BBKNN, ComBat and LIGER latent
#' spaces. It does not run these external tools or import their native objects.
#' The selected reductions are copied or referenced and recorded in
#' `metadata(sce)$scProkaR$integration`, so [BenchmarkIntegration()] can compare
#' them with scProkaR MNN and Harmony outputs.
#'
#' @param sce A `SingleCellExperiment`.
#' @param reductions Character vector of existing reduced-dimension names.
#'   If `NULL`, likely integration reductions are detected from their names.
#'   Pass names explicitly when automatic detection is unsuitable.
#' @param method_names Optional names to use in benchmark tables. Must be
#'   unique and have the same length as `reductions`. By default, names are
#'   inferred from the reductions and converted to lower case.
#' @param copy Logical. If `TRUE`, reductions not already starting with
#'   `prefix` are copied to `paste0(prefix, method_name)`. If `FALSE`, their
#'   original names are registered directly.
#' @param prefix Prefix used when `copy = TRUE`.
#' @param metadata Optional named list added to each method's metadata. A list
#'   entry named by `method_name`, if present, is used for that method instead.
#'   `source` and `original_reduction` are recorded by this function.
#' @param overwrite Logical indicating whether copied reduced dimensions may
#'   replace existing reductions with different values.
#'
#' @return A `SingleCellExperiment` with registered integration embeddings.
#' @export
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' # Stand-in for a latent space from an external integration tool.
#' SingleCellExperiment::reducedDim(sce, "X_scVI") <-
#'     SingleCellExperiment::reducedDim(sce, "PCA")[, seq_len(3)]
#' sce <- RegisterExistingIntegrationEmbeddings(sce, reductions = "X_scVI")
#' SingleCellExperiment::reducedDimNames(sce)
#' S4Vectors::metadata(sce)$scProkaR$integration$results$scvi
RegisterExistingIntegrationEmbeddings <- function(
    sce,
    reductions = NULL,
    method_names = NULL,
    copy = TRUE,
    prefix = "integrated_",
    metadata = list(),
    overwrite = FALSE
) {
    if (!is(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment.", call. = FALSE)
    }

    available <- reducedDimNames(sce)
    if (length(available) == 0L) {
        stop("No reduced dimensions are available in `sce`.", call. = FALSE)
    }

    if (is.null(reductions)) {
        reductions <- .scprokar_guess_external_integration_reductions(
            available
        )
        ## Never re-register a reduction scProkaR already owns. Doing so would
        ## replace the provenance recorded by IntegrateBacData(), which knows
        ## the batch column and feature set, with a generic external record.
        reductions <- setdiff(
            reductions, .scprokar_registered_reductions(sce)
        )
        if (length(reductions) == 0L) {
            stop(
                "No likely external integration reductions were detected. ",
                "Pass `reductions = ...` explicitly.",
                call. = FALSE
            )
        }
    }

    if (!is.character(reductions) || length(reductions) == 0L ||
            anyNA(reductions) || any(!nzchar(reductions))) {
        stop("`reductions` must contain non-empty names.", call. = FALSE)
    }
    missing <- setdiff(reductions, available)
    if (length(missing) > 0L) {
        stop(
            "The following reductions were not found: ",
            paste(missing, collapse = ", "),
            ". Available reductions: ", paste(available, collapse = ", "),
            call. = FALSE
        )
    }

    if (is.null(method_names)) {
        method_names <- vapply(
            reductions, .scprokar_infer_integration_method_name, character(1)
        )
    }
    if (length(method_names) != length(reductions)) {
        stop(
            "`method_names` must have the same length as `reductions`.",
            call. = FALSE
        )
    }
    if (!is.character(method_names) || anyNA(method_names) ||
            any(!nzchar(method_names)) || anyDuplicated(method_names)) {
        stop("`method_names` must be unique, non-empty names.", call. = FALSE)
    }

    Reduce(function(sce, i) {
        reduction <- reductions[[i]]
        method_name <- method_names[[i]]
        target <- if (isTRUE(copy) && !startsWith(reduction, prefix)) {
            paste0(prefix, method_name)
        } else {
            reduction
        }
        method_metadata <- metadata
        if (is.list(metadata[[method_name]])) {
            method_metadata <- metadata[[method_name]]
        }
        method_metadata$source <- "external_existing_reduction"
        method_metadata$original_reduction <- reduction

        RegisterIntegrationEmbedding(
            sce = sce,
            embedding = reducedDim(sce, reduction),
            method_name = method_name,
            reduction_name = target,
            metadata = method_metadata,
            overwrite = overwrite
        )
    }, seq_along(reductions), sce)
}

#' @keywords internal
.scprokar_guess_external_integration_reductions <- function(reduction_names) {
    lower <- tolower(reduction_names)
    ## Require a recognisable tool name. A bare "^x_" would match every
    ## AnnData reduction, including plain visualisation layouts.
    include <- grepl(
        paste(
            c(
                "^integrated_", "cca", "rpca", "scanorama", "scvi",
                "scgen", "scalex", "bbknn", "combat", "regress", "liger",
                "harmony", "mnn"
            ),
            collapse = "|"
        ),
        lower
    )
    ## Visualisation layouts are never latent spaces, wherever the token sits
    ## in the name: "harmony_umap" is a picture of an embedding, not one.
    exclude <- grepl(
        paste(
            c(
                "umap", "tsne", "t_sne", "fdl", "draw_graph", "phate",
                "densmap", "diffmap", "diffusion", "force"
            ),
            collapse = "|"
        ),
        lower
    ) | lower %in% c("pca", "x_pca", "truth")
    reduction_names[include & !exclude]
}

#' @keywords internal
.scprokar_infer_integration_method_name <- function(reduction_name) {
    out <- sub("^integrated_", "", reduction_name, ignore.case = TRUE)
    out <- sub("^X_", "", out, ignore.case = TRUE)
    out <- sub("^seurat_", "", out, ignore.case = TRUE)
    out <- gsub("[^A-Za-z0-9]+", "_", out)
    out <- gsub("^_+|_+$", "", out)
    tolower(out)
}

#' @keywords internal
.scprokar_prepare_mnn_batches <- function(batches, features, args) {
    cos_norm <- if ("cos.norm" %in% names(args)) args$cos.norm else TRUE
    if (!cos_norm) {
        return(list(batches = batches, args = args))
    }
    assay_type <- if ("assay.type" %in% names(args)) {
        args$assay.type
    } else {
        "logcounts"
    }
    cosine_args <- list(mode = "l2norm", subset.row = features)
    if ("BPPARAM" %in% names(args)) {
        cosine_args$BPPARAM <- args$BPPARAM
    }
    batches <- lapply(batches, function(sce) {
        x <- assay(sce, assay_type)
        l2 <- do.call(batchelor::cosineNorm, c(list(x = x), cosine_args))
        ## Match fastMNN's column scaling, including its zero-norm floor,
        ## without the deprecated scuttle::normalizeCounts() call.
        normalized <- Matrix::t(Matrix::t(x) / pmax(l2, 1e-8))
        dimnames(normalized) <- dimnames(x)
        assay(sce, assay_type) <- normalized
        sce
    })
    ## The input is already cosine-normalized; avoid applying it twice.
    args$cos.norm <- FALSE
    list(batches = batches, args = args)
}

#' @keywords internal
.scprokar_run_harmony <- function(pca_embedding, meta_data, batch_col, ...) {
    extra <- list(...)
    n_cells <- nrow(pca_embedding)
    if (!"nclust" %in% names(extra)) {
        if (n_cells <= 2L) {
            extra$nclust <- 1L
        } else {
            extra$nclust <- max(2L, min(20L, n_cells - 1L))
        }
    }
    if (!"verbose" %in% names(extra)) {
        extra$verbose <- FALSE
    }

    if ("RunHarmony" %in% getNamespaceExports("harmony")) {
        harmony_fun <- getExportedValue("harmony", "RunHarmony")
        args_new <- c(
            list(
                data_mat = pca_embedding,
                meta_data = meta_data,
                vars_use = batch_col
            ),
            extra
        )
        ## Retry version-specific signatures on errors while allowing
        ## convergence, parameter and other backend warnings to reach callers.
        out <- .scprokar_attempt(do.call(harmony_fun, args_new))

        if (is.null(out)) {
            args_legacy <- c(
                list(
                    data_mat = pca_embedding,
                    meta_data = meta_data,
                    vars_use = batch_col,
                    do_pca = FALSE
                ),
                extra
            )
            out <- .scprokar_attempt(do.call(harmony_fun, args_legacy))
        }

        if (!is.null(out)) {
            return(as.matrix(out))
        }
    }

    if ("HarmonyMatrix" %in% getNamespaceExports("harmony")) {
        harmony_fun <- getExportedValue("harmony", "HarmonyMatrix")
        args_matrix <- c(
            list(
                data_mat = pca_embedding,
                meta_data = meta_data,
                vars_use = batch_col,
                do_pca = FALSE
            ),
            extra
        )
        out <- .scprokar_attempt(do.call(harmony_fun, args_matrix))
        if (!is.null(out)) {
            return(as.matrix(out))
        }
    }

    stop(
        "Harmony integration failed for this dataset. Try reducing dimensions ",
        "or using method = \"mnn\".",
        call. = FALSE
    )
}

#' @keywords internal
.scprokar_default_integration_reduction <- function(
    sce,
    integrated_reduction = NULL,
    prefer_umap = TRUE
) {
    if (!is.null(integrated_reduction)) {
        return(.scprokar_resolve_reduction_name(sce, integrated_reduction))
    }

    reduction_names <- reducedDimNames(sce)
    integrated_matches <- grep("^integrated_", reduction_names, value = TRUE)
    umap_matches <- grep("^umap_integrated_", reduction_names, value = TRUE)
    if (isTRUE(prefer_umap)) {
        if (length(umap_matches) > 0L) {
            return(umap_matches[[1L]])
        }
        if (length(integrated_matches) > 0L) {
            return(integrated_matches[[1L]])
        }
    } else {
        if (length(integrated_matches) > 0L) {
            return(integrated_matches[[1L]])
        }
        if (length(umap_matches) > 0L) {
            return(umap_matches[[1L]])
        }
    }

    stop(
        "No integrated reduction was found. Run IntegrateBacData() and ",
        "optionally RunIntegratedUMAP(), ",
        "or pass `integrated_reduction` explicitly.",
        call. = FALSE
    )
}

#' @keywords internal
.scprokar_clustering_embedding <- function(sce, reduction, dims) {
    reduction_name <- .scprokar_default_integration_reduction(
        sce, integrated_reduction = reduction, prefer_umap = FALSE
    )
    emb <- as.matrix(reducedDim(sce, reduction_name))
    if (!is.null(dims)) {
        emb <- emb[, dims, drop = FALSE]
    }
    if (nrow(emb) < 2L) {
        stop("Need at least two cells to run clustering.", call. = FALSE)
    }
    if (ncol(emb) < 2L) {
        stop(
            "Integrated reduction '", reduction_name,
            "' must have at least two dimensions.", call. = FALSE
        )
    }
    list(reduction = reduction_name, embedding = emb)
}

.scprokar_reduction_embedding <- function(sce, reduction_name, purpose) {
    emb <- as.matrix(reducedDim(sce, reduction_name))
    if (ncol(emb) < 2) {
        stop(
            "Reduction '", reduction_name,
            "' must have at least two columns to ", purpose, ".",
            call. = FALSE
        )
    }
    emb
}

#' @keywords internal
.scprokar_cluster_membership <- function(
    graph, algorithm, resolution, n_cells
) {
    if (igraph::gsize(graph) == 0L) {
        return(rep(1L, n_cells))
    }
    if (algorithm == "louvain") {
        cluster_args <- list(graph, weights = igraph::E(graph)$weight)
        if ("resolution" %in% names(formals(igraph::cluster_louvain))) {
            cluster_args$resolution <- resolution
        }
        return(igraph::membership(do.call(
            igraph::cluster_louvain, cluster_args
        )))
    }
    if (algorithm == "leiden") {
        return(.scprokar_leiden_membership(graph, resolution))
    }
    if (!is.null(resolution)) {
        warning(
            "`resolution` is ignored when `algorithm = \"walktrap\"` ",
            "because walktrap does not expose a resolution parameter.",
            call. = FALSE
        )
    }
    igraph::membership(igraph::cluster_walktrap(
        graph, weights = igraph::E(graph)$weight
    ))
}

#' @keywords internal
.scprokar_leiden_membership <- function(graph, resolution) {
    if (!"cluster_leiden" %in% getNamespaceExports("igraph")) {
        stop(
            "The installed igraph does not provide `cluster_leiden()`. ",
            "Use `algorithm = \"louvain\"` or upgrade igraph.",
            call. = FALSE
        )
    }
    ## Modularity keeps resolution on the same scale as Louvain, unlike
    ## Leiden's default CPM objective, which over-partitions a kNN graph.
    cluster_args <- list(
        graph, weights = igraph::E(graph)$weight,
        objective_function = "modularity"
    )
    leiden_formals <- names(formals(get(
        "cluster_leiden", envir = asNamespace("igraph")
    )))
    ## igraph 2.1.0 renamed and deprecated `resolution_parameter`.
    if ("resolution" %in% leiden_formals) {
        cluster_args$resolution <- resolution
    } else if ("resolution_parameter" %in% leiden_formals) {
        cluster_args$resolution_parameter <- resolution
    }
    igraph::membership(do.call(igraph::cluster_leiden, cluster_args))
}

#' @keywords internal
.scprokar_build_cluster_graph <- function(emb, k) {
    n <- nrow(emb)

    if (requireNamespace("RANN", quietly = TRUE)) {
        nn <- RANN::nn2(emb, query = emb, k = k + 1)
        idx <- nn$nn.idx[, -1, drop = FALSE]
        dists <- nn$nn.dists[, -1, drop = FALSE]
    } else {
        dmat <- as.matrix(stats::dist(emb))
        cols <- seq_len(k) + 1L
        idx <- t(apply(dmat, 1, order))[, cols, drop = FALSE]
        dists <- t(apply(dmat, 1, sort, na.last = TRUE))[, cols, drop = FALSE]
    }

    edges <- unique(do.call(
        rbind,
        lapply(seq_len(n), function(i) {
            nbrs <- idx[i, ]
            data.frame(
                from = pmin(i, nbrs),
                to = pmax(i, nbrs),
                weight = 1 / (1 + dists[i, ]),
                stringsAsFactors = FALSE
            )
        })
    ))
    edges <- edges[edges$from != edges$to, , drop = FALSE]

    if (!nrow(edges)) {
        graph <- igraph::make_empty_graph(n = n, directed = FALSE)
        return(list(graph = graph, edges = edges))
    }

    edge_key <- paste(edges$from, edges$to, sep = "_")
    weights <- stats::aggregate(
        edges$weight,
        by = list(edge_key = edge_key),
        FUN = mean
    )
    split_key <- strsplit(weights$edge_key, "_", fixed = TRUE)
    edge_df <- data.frame(
        from = as.integer(vapply(split_key, `[`, character(1), 1)),
        to = as.integer(vapply(split_key, `[`, character(1), 2)),
        weight = weights$x
    )

    graph <- igraph::graph_from_data_frame(
        edge_df,
        directed = FALSE,
        vertices = seq_len(n)
    )
    list(graph = graph, edges = edge_df)
}

#' Assemble the point data frame used by PlotReduction()
#'
#' Combines the first two embedding dimensions with `colData(sce)`, adds the
#' internal `.colour` and `.facet` columns, and optionally shuffles the
#' plotting order.
#'
#' @keywords internal
#' @noRd
.scprokar_reduction_plot_data <- function(
    sce,
    emb,
    colour_by,
    facet_by,
    shuffle
) {
    meta <- as.data.frame(colData(sce))
    plot_df <- data.frame(
        dim1 = emb[, 1],
        dim2 = emb[, 2],
        meta,
        stringsAsFactors = FALSE
    )

    if (!is.null(colour_by)) {
        if (!colour_by %in% colnames(plot_df)) {
            stop(
                "`colour_by` was not found in colData(sce): ",
                colour_by,
                call. = FALSE
            )
        }
        plot_df$.colour <- plot_df[[colour_by]]
    } else {
        plot_df$.colour <- "cells"
    }

    if (!is.null(facet_by)) {
        if (!facet_by %in% colnames(plot_df)) {
            stop(
                "`facet_by` was not found in colData(sce): ",
                facet_by,
                call. = FALSE
            )
        }
        plot_df$.facet <- plot_df[[facet_by]]
    }

    if (isTRUE(shuffle)) {
        plot_df <- plot_df[sample(seq_len(nrow(plot_df))), , drop = FALSE]
    }

    plot_df
}

#' Build the ggplot2 colour scale for PlotReduction()
#'
#' Chooses a fixed grey scale, a continuous gradient, or a discrete manual
#' scale depending on the type of the colouring variable.
#'
#' @keywords internal
#' @noRd
.scprokar_reduction_colour_scale <- function(
    colour_values,
    colour_by,
    palette
) {
    if (is.null(colour_by)) {
        return(ggplot2::scale_color_manual(
            values = c(cells = "grey50"),
            guide = "none"
        ))
    }

    if (is.numeric(colour_values)) {
        return(ggplot2::scale_color_gradientn(
            colours = if (is.null(palette)) c(
                "#2b8cbe", "#fdbb84", "#d7301f"
            ) else palette
        ))
    }

    discrete_values <- unique(as.character(colour_values))
    colors <- if (is.null(palette)) {
        stats::setNames(
            grDevices::hcl.colors(length(discrete_values), "Dark 3"),
            discrete_values
        )
    } else {
        if (is.null(names(palette))) {
            names(palette) <- discrete_values[seq_len(min(
                length(discrete_values), length(palette)
            ))]
        }
        palette[discrete_values]
    }
    ggplot2::scale_color_manual(values = colors)
}

#' Reductions already recorded in the scProkaR integration registry
#'
#' @param sce A `SingleCellExperiment`.
#'
#' @return A character vector of reduced-dimension names, possibly empty.
#' @keywords internal
#' @noRd
.scprokar_registered_reductions <- function(sce) {
    results <- .scprokar_get_metadata(sce)$integration$results
    if (!is.list(results) || length(results) == 0L) {
        return(character(0))
    }
    found <- vapply(
        results,
        function(entry) {
            reduction <- entry[["reduction"]]
            if (is.character(reduction) && length(reduction) == 1L) {
                reduction
            } else {
                NA_character_
            }
        },
        character(1)
    )
    unname(found[!is.na(found)])
}
