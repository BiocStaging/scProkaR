#' Plot the cluster-level TATA graph
#'
#' Plot the directed cluster-level TATA graph. The graph can be drawn in its
#' own abstract layout or positioned on a chosen reduced dimension so that the
#' graph coordinates match the embedding shown to the user.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param dimred Optional reduced dimension used to position cluster nodes. When
#'   supplied together with `layout_mode = "embedding"`, cluster centroids are
#'   computed from this embedding.
#' @param layout_mode Layout mode for the graph. Choices are:
#' - `"graph"`: use the stored abstract graph coordinates, or a force-directed
#'   layout if they are unavailable.
#' - `"embedding"`: place cluster nodes at their centroids in the selected
#'   reduced dimension.
#' @param cluster_col Optional cluster column used when `layout_mode =
#'   "embedding"`. By default this is taken from `tata_result$parameters`.
#' @param show_cells Logical indicating whether to draw cells in the background
#'   when `layout_mode = "embedding"`.
#' @param colour_by Optional `colData` column used to color background cells. If
#'   `NULL` and `show_cells = TRUE`, the cluster column is used when available.
#' @param node_colour_by Node coloring scheme. Choices are:
#' - `"median_time"`: color nodes by cluster median time.
#' - `"cluster"`: color nodes by cluster identity, which is useful for
#'   embedding-aligned graph plots.
#' @param main Plot title.
#' @param vertex_label_size Label size.
#' @param edge_width_scale Relative edge width scale.
#' @param point_size Background point size when cells are drawn.
#' @param point_alpha Background point alpha when cells are drawn.
#' @param arrow_size Arrow size for directed edges.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # Abstract graph layout, nodes coloured by cluster median time.
#' p <- plot_tata_cluster_graph(tata)
#'
#' # The same graph positioned on the UMAP embedding.
#' p_embedded <- plot_tata_cluster_graph(tata,
#'     dimred = "UMAP",
#'     layout_mode = "embedding", show_cells = TRUE,
#'     node_colour_by = "cluster"
#' )
#' p_embedded
plot_tata_cluster_graph <- function(
    tata_result,
    dimred = NULL,
    layout_mode = NULL,
    cluster_col = NULL,
    show_cells = FALSE,
    colour_by = NULL,
    node_colour_by = c("median_time", "cluster"),
    main = "TATA cluster graph",
    vertex_label_size = 3.5,
    edge_width_scale = 5,
    point_size = 0.4,
    point_alpha = 0.45,
    arrow_size = 0.15
) {
    cluster_graph <- tata_result$cluster_graph

    if (!inherits(cluster_graph, "igraph")) {
        stop(
            "`tata_result$cluster_graph` must be an igraph object.",
            call. = FALSE
        )
    }

    if (igraph::vcount(cluster_graph) == 0L) {
        stop("The cluster graph is empty.", call. = FALSE)
    }

    if (is.null(layout_mode)) {
        layout_mode <- if (!is.null(dimred)) "embedding" else "graph"
    }
    layout_mode <- match.arg(layout_mode, c("graph", "embedding"))
    node_colour_by <- match.arg(node_colour_by)

    graph_vertex_df <- data.frame(
        cluster = igraph::V(cluster_graph)$name,
        size = igraph::V(cluster_graph)$size,
        median_time = igraph::V(cluster_graph)$median_time,
        stringsAsFactors = FALSE
    )

    cell_df <- NULL
    if (identical(layout_mode, "embedding")) {
        context <- .tata_cluster_embedding_context(
            tata_result, dimred, cluster_col
        )
        dimred <- context$dimred
        cluster_col <- context$cluster_col
        node_df <- .tata_cluster_embedding_nodes(
            context$meta, graph_vertex_df
        )

        if (isTRUE(show_cells)) {
            if (is.null(colour_by)) {
                colour_by <- cluster_col
            }

            cell_df <- .tata_cluster_cell_table(context$meta, colour_by)
        }
    } else {
        if (isTRUE(show_cells)) {
            warning(
                "`show_cells = TRUE` is only used with `layout_mode = ",
                "\"embedding\"`.",
                call. = FALSE
            )
        }

        node_df <- .tata_cluster_abstract_nodes(
            cluster_graph, graph_vertex_df
        )
    }

    if (any(!is.finite(node_df$x)) || any(!is.finite(node_df$y))) {
        stop(
            "Cluster node positions could not be computed for plotting.",
            call. = FALSE
        )
    }

    edge_df <- .tata_cluster_edge_table(
        cluster_graph, node_df, edge_width_scale
    )

    node_style <- .tata_cluster_node_aesthetics(node_df)
    node_df <- node_style$node_df
    cluster_palette <- node_style$palette

    p <- ggplot2::ggplot()
    p <- .tata_add_cluster_cell_layer(
        p, cell_df, colour_by, cluster_col, cluster_palette,
        levels(node_df$cluster), point_size, point_alpha
    )
    p <- .tata_add_cluster_edge_layers(p, edge_df, arrow_size)
    p <- .tata_add_cluster_node_layers(
        p, node_df, cluster_palette, node_colour_by, vertex_label_size
    )

    p <- p +
        ggplot2::scale_size_identity() +
        ggplot2::scale_linewidth_identity() +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = main,
            x = if (identical(layout_mode, "embedding"))
                paste0(dimred, "_1") else "graph_1",
            y = if (identical(layout_mode, "embedding"))
                paste0(dimred, "_2") else "graph_2"
        )

    p
}


#' Plot TATA trajectories over an embedding
#'
#' Overlay inferred TATA trajectory paths on top of a reduced-dimensional
#' embedding such as UMAP. The plotted paths follow root-to-terminal routes in
#' the cluster graph and can be shown as smoothed curves for a more classical
#' trajectory-style visualization.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param dimred Reduced dimension to plot.
#' @param colour_by Optional `colData` column used to color cells. Set to `NULL`
#'   for a grey background.
#' @param max_paths Optional maximum number of paths to plot after path
#'   selection. `NULL` keeps all available paths.
#' @param point_size Background point size.
#' @param point_alpha Background point alpha.
#' @param smooth_paths Whether to smooth the inferred trajectory paths.
#' @param n_curve_points Number of interpolation points per curve.
#' @param show_centroids Whether to draw cluster centroids.
#' @param show_labels Whether to label cluster centroids.
#' @param curve_colour Colour of the inferred trajectory curves.
#' @param label_size Cluster label size.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # Overlay the three root-to-terminal routes on the UMAP embedding.
#' p <- plot_tata_trajectory_embedding(tata,
#'     dimred = "UMAP",
#'     colour_by = "timepoint", max_paths = 3
#' )
#' p
plot_tata_trajectory_embedding <- function(
    tata_result,
    dimred = "UMAP",
    colour_by = "timepoint",
    max_paths = NULL,
    point_size = 0.6,
    point_alpha = 0.75,
    smooth_paths = TRUE,
    n_curve_points = 100,
    show_centroids = TRUE,
    show_labels = TRUE,
    curve_colour = "black",
    label_size = 3
) {
    sce <- tata_result$sce
    if (!methods::is(sce, "SingleCellExperiment")) {
        stop("`tata_result$sce` must be a SingleCellExperiment.", call. = FALSE)
    }

    max_paths <- .tata_trajectory_max_paths(max_paths)

    cluster_col <- tata_result$parameters$cluster_col
    meta <- .tata_trajectory_meta(sce, dimred, cluster_col)
    cluster_centroids <- .tata_trajectory_centroids(
        meta, tata_result$cluster_pseudotime
    )

    curve_df <- .tata_trajectory_curve_table(
        tata_result, cluster_centroids, max_paths, smooth_paths,
        n_curve_points
    )

    if (!is.null(colour_by) && colour_by %in% colnames(meta)) {
        meta$colour_value <- meta[[colour_by]]
        discrete_colour <- !is.numeric(meta$colour_value)
    } else {
        meta$colour_value <- "cells"
        discrete_colour <- TRUE
    }

    p <- ggplot2::ggplot(meta, ggplot2::aes(x = Dim1, y = Dim2))
    if (!is.null(colour_by) && colour_by %in% colnames(meta)) {
        p <- p + ggplot2::geom_point(
            ggplot2::aes(colour = colour_value),
            size = point_size,
            alpha = point_alpha
        )
    } else {
        p <- p + ggplot2::geom_point(
            colour = "grey80",
            size = point_size,
            alpha = point_alpha
        )
    }

    if (!is.null(curve_df) && nrow(curve_df) > 0L) {
        p <- p + ggplot2::geom_path(
            data = curve_df,
            mapping = ggplot2::aes(x = Dim1, y = Dim2, group = path_id),
            inherit.aes = FALSE,
            colour = curve_colour,
            linewidth = 1.1,
            lineend = "round",
            arrow = grid::arrow(
                length = grid::unit(0.15, "inches"),
                type = "closed"
            )
        )
    }

    if (isTRUE(show_centroids)) {
        p <- p + ggplot2::geom_point(
            data = cluster_centroids,
            mapping = ggplot2::aes(x = Dim1, y = Dim2),
            inherit.aes = FALSE,
            shape = 21,
            size = 2.6,
            fill = "white",
            colour = "black",
            stroke = 0.6
        )
    }

    if (isTRUE(show_labels)) {
        p <- p + ggplot2::geom_text(
            data = cluster_centroids,
            mapping = ggplot2::aes(x = Dim1, y = Dim2, label = cluster),
            inherit.aes = FALSE,
            size = label_size,
            nudge_y = 0.15
        )
    }

    p <- p +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = paste0(dimred, "_1"),
            y = paste0(dimred, "_2"),
            colour = colour_by
        )

    if (!is.null(colour_by) && colour_by %in% colnames(meta) &&
            !discrete_colour) {
        p <- p + ggplot2::scale_colour_viridis_c(option = "plasma")
    }

    p
}


#' Plot TATA edge causality diagnostics
#'
#' Visualize how topology support and temporal evidence combine for each
#' cluster-cluster edge in the abstract TATA graph.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param label_edges Logical indicating whether to label each point by the
#'   cluster-pair direction label.
#' @param kept_only Logical indicating whether to plot only retained edges.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # All candidate edges, showing which ones survived pruning.
#' p <- plot_tata_edge_causality(tata)
#'
#' # Retained edges only, labelled by their inferred direction.
#' p_kept <- plot_tata_edge_causality(tata,
#'     kept_only = TRUE,
#'     label_edges = TRUE
#' )
#' p_kept
plot_tata_edge_causality <- function(
    tata_result,
    label_edges = FALSE,
    kept_only = FALSE
) {
    edge_df <- tata_result$cluster_edge_table
    if (is.null(edge_df) || nrow(edge_df) == 0L) {
        stop(
            "`tata_result` does not contain a cluster edge table.",
            call. = FALSE
        )
    }

    if (isTRUE(kept_only)) {
        edge_df <- edge_df[edge_df$kept, , drop = FALSE]
    }

    edge_df$edge_status <- ifelse(edge_df$kept, "kept", "pruned")
    edge_df$direction_group <- ifelse(
        edge_df$direction == "ambiguous",
        "ambiguous",
        ifelse(edge_df$signed_time_score >= 0, "forward", "backward")
    )

    p <- ggplot2::ggplot(
        edge_df,
        ggplot2::aes(
            x = topology_weight,
            y = signed_time_score,
            colour = edge_status,
            shape = direction_group,
            size = tata_weight
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::scale_colour_manual(
            values = c(kept = "#1b9e77", pruned = "#c94c4c")
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "TATA edge causality diagnostic",
            x = "topology weight",
            y = "signed temporal evidence",
            colour = "edge status",
            shape = "direction",
            size = "TATA weight"
        )

    if (isTRUE(label_edges)) {
        p <- p + ggplot2::geom_text(
            ggplot2::aes(label = direction_label),
            size = 3,
            nudge_y = 0.03,
            show.legend = FALSE
        )
    }

    p
}


#' Plot TATA time calibration
#'
#' Compare inferred TATA pseudotime against the known experimental time labels,
#' optionally split by a branch or condition label.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param time_col `colData` column containing experimental time labels. By
#'   default this is taken from `tata_result$parameters`.
#' @param branch_col Optional `colData` column used to group curves, for
#'   example a true simulated branch or treatment label.
#' @param summary Plot style. Choices are:
#' - `"median_iqr"`: median pseudotime with an IQR ribbon for each group.
#' - `"boxplot"`: per-timepoint boxplots of pseudotime.
#' @param pseudotime_col `colData` column containing the cell-level TATA
#'   pseudotime values to plot.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # Median pseudotime with an IQR ribbon for each treatment condition.
#' p <- plot_tata_time_calibration(tata, branch_col = "condition")
#'
#' # The same comparison shown as per-timepoint boxplots.
#' p_box <- plot_tata_time_calibration(tata,
#'     branch_col = "condition",
#'     summary = "boxplot"
#' )
#' p_box
plot_tata_time_calibration <- function(
    tata_result,
    time_col = NULL,
    branch_col = NULL,
    summary = c("median_iqr", "boxplot"),
    pseudotime_col = "tata_pseudotime_scaled"
) {
    summary <- match.arg(summary)
    sce <- tata_result$sce
    meta <- as.data.frame(SummarizedExperiment::colData(sce))

    if (is.null(time_col)) {
        time_col <- tata_result$parameters$time_col
    }
    if (is.null(time_col) || !time_col %in% colnames(meta)) {
        stop(
            "A valid `time_col` is required for time calibration plotting.",
            call. = FALSE
        )
    }
    if (!pseudotime_col %in% colnames(meta)) {
        stop(
            "`pseudotime_col` is not present in `colData(sce)`.",
            call. = FALSE
        )
    }

    plot_df <- .tata_time_calibration_table(
        meta, time_col, pseudotime_col, branch_col
    )

    if (summary == "boxplot") {
        return(
            ggplot2::ggplot(
                plot_df,
                ggplot2::aes(
                    x = factor(timepoint), y = pseudotime, fill = branch
                )
            ) +
                ggplot2::geom_boxplot(outlier.size = 0.2, alpha = 0.8) +
                ggplot2::theme_classic() +
                ggplot2::labs(
                    title = "TATA time calibration",
                    x = "experimental time",
                    y = "TATA pseudotime",
                    fill = if (!is.null(branch_col)) branch_col else NULL
                )
        )
    }

    summary_df <- .tata_time_calibration_summary(plot_df)

    ggplot2::ggplot(
        summary_df,
        ggplot2::aes(
            x = timepoint, y = median, colour = branch,
            fill = branch, group = branch
        )
    ) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = q25, ymax = q75),
            alpha = 0.18, linewidth = 0, colour = NA
        ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 1.8) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "TATA time calibration",
            x = "experimental time",
            y = "median TATA pseudotime",
            colour = if (!is.null(branch_col)) branch_col else NULL,
            fill = if (!is.null(branch_col)) branch_col else NULL
        )
}


#' Plot truth versus inferred TATA branch assignments
#'
#' For simulated data, compare a ground-truth branch label with the inferred
#' hard TATA branch assignment.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param truth_col `colData` column containing the reference branch labels.
#' @param predicted_col `colData` column containing the inferred TATA branch
#'   assignment.
#' @param normalize Normalization scheme. Choices are:
#' - `"truth"`: normalize each truth branch to sum to 1.
#' - `"predicted"`: normalize each predicted branch to sum to 1.
#' - `"none"`: raw counts.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # `simulated_branch` records the true branch of each simulated cell.
#' p <- plot_tata_branch_confusion(tata,
#'     truth_col = "simulated_branch",
#'     normalize = "truth"
#' )
#' p
plot_tata_branch_confusion <- function(
    tata_result,
    truth_col = "simulated_branch",
    predicted_col = "tata_branch_assignment",
    normalize = c("truth", "predicted", "none")
) {
    normalize <- match.arg(normalize)
    meta <- as.data.frame(SummarizedExperiment::colData(tata_result$sce))
    if (!truth_col %in% colnames(meta) || !predicted_col %in% colnames(meta)) {
        stop(
            "Both `truth_col` and `predicted_col` must be present in ",
            "`colData(sce)`.",
            call. = FALSE
        )
    }

    tab <- as.data.frame(table(
        truth = as.character(meta[[truth_col]]),
        predicted = as.character(meta[[predicted_col]])
    ))
    names(tab)[3] <- "value"

    if (normalize == "truth") {
        tab$value <- tab$value / stats::ave(
            tab$value, tab$truth,
            FUN = function(x) pmax(sum(x), 1)
        )
    } else if (normalize == "predicted") {
        tab$value <- tab$value / stats::ave(
            tab$value, tab$predicted,
            FUN = function(x) pmax(sum(x), 1)
        )
    }

    ggplot2::ggplot(tab, ggplot2::aes(x = predicted, y = truth, fill = value)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::scale_fill_viridis_c(option = "C") +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(
            title = "Truth-aligned TATA branch assignment matrix",
            x = "predicted branch",
            y = "truth branch",
            fill = if (normalize == "none") "count" else "fraction"
        )
}


#' Plot TATA branch probabilities on an embedding
#'
#' Plot one panel per terminal branch, colored by the corresponding TATA branch
#' probability.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param dimred Reduced dimension to plot.
#' @param branches Optional subset of terminal branches to display.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @param ncol Number of facet columns.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # One panel per terminal branch, coloured by branch probability.
#' p <- plot_tata_branch_probabilities(tata, dimred = "UMAP", ncol = 2)
#' p
plot_tata_branch_probabilities <- function(
    tata_result,
    dimred = "UMAP",
    branches = NULL,
    point_size = 0.6,
    point_alpha = 0.8,
    ncol = 2
) {
    sce <- tata_result$sce
    if (!dimred %in% SingleCellExperiment::reducedDimNames(sce)) {
        stop(
            "Reduced dimension `", dimred,
            "` is not present in `tata_result$sce`.",
            call. = FALSE
        )
    }

    branch_prob_df <- tata_result$branch_probabilities
    if (is.null(branch_prob_df)) {
        stop(
            "`tata_result` does not contain branch probabilities.",
            call. = FALSE
        )
    }

    if (is.null(branches)) {
        branches <- tata_result$terminal_clusters
    }
    branch_cols <- .make_tata_probability_colnames(branches)

    embedding <- as.matrix(SingleCellExperiment::reducedDim(sce, dimred))
    plot_df <- data.frame(
        cell_id = colnames(sce),
        Dim1 = embedding[, 1],
        Dim2 = embedding[, 2],
        branch_prob_df[, branch_cols, drop = FALSE],
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    long_df <- stats::reshape(
        plot_df,
        varying = branch_cols,
        v.names = "probability",
        timevar = "branch",
        times = branches,
        direction = "long"
    )

    ggplot2::ggplot(
        long_df,
        ggplot2::aes(x = Dim1, y = Dim2, colour = probability)
    ) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::facet_wrap(~branch, ncol = ncol) +
        ggplot2::scale_colour_viridis_c(option = "magma", limits = c(0, 1)) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "TATA terminal-branch probabilities",
            x = paste0(dimred, "_1"),
            y = paste0(dimred, "_2"),
            colour = "probability"
        )
}


#' Plot selected cells from one TATA branch
#'
#' Highlight cells that belong to one terminal branch using either a hard or a
#' soft branch-selection rule.
#'
#' @inheritParams select_tata_branch_cells
#' @param dimred Reduced dimension to plot.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # Highlight cells with an appreciable probability of one terminal fate.
#' branch <- tata$terminal_clusters[1]
#' p <- plot_tata_branch_selection(tata,
#'     branch = branch, dimred = "UMAP",
#'     selection_mode = "soft", probability_threshold = 0.25
#' )
#' p
plot_tata_branch_selection <- function(
    tata_result,
    branch,
    dimred = "UMAP",
    selection_mode = c("soft", "hard"),
    probability_threshold = 0.25,
    point_size = 0.6,
    point_alpha = 0.8
) {
    selection_df <- select_tata_branch_cells(
        tata_result = tata_result,
        branch = branch,
        selection_mode = selection_mode,
        probability_threshold = probability_threshold
    )
    sce <- tata_result$sce
    embedding <- as.matrix(SingleCellExperiment::reducedDim(sce, dimred))
    plot_df <- data.frame(
        cell_id = colnames(sce),
        Dim1 = embedding[, 1],
        Dim2 = embedding[, 2],
        selection_df[
            match(colnames(sce), selection_df$cell_id),
            c("probability", "selected"),
            drop = FALSE
        ],
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2)) +
        ggplot2::geom_point(
            colour = "grey85", size = point_size,
            alpha = point_alpha * 0.6
        ) +
        ggplot2::geom_point(
            data = plot_df[plot_df$selected, , drop = FALSE],
            mapping = ggplot2::aes(colour = probability),
            size = point_size,
            alpha = point_alpha
        ) +
        ggplot2::scale_colour_viridis_c(option = "plasma", limits = c(0, 1)) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = paste("Cells selected for branch", branch),
            x = paste0(dimred, "_1"),
            y = paste0(dimred, "_2"),
            colour = "probability"
        )
}


#' Plot branch probabilities for representative cells
#'
#' Show a barplot of terminal-branch probabilities for selected cells.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param cell_ids Optional cell ids. If `NULL`, the cells with the highest TATA
#'   branch entropy are shown.
#' @param n_cells Number of representative cells to display when `cell_ids` is
#'   `NULL`.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # The most branch-ambiguous cells have the flattest probability bars.
#' p <- plot_tata_terminal_probabilities(tata, n_cells = 4)
#' p
plot_tata_terminal_probabilities <- function(
    tata_result,
    cell_ids = NULL,
    n_cells = 6
) {
    branch_prob_df <- tata_result$branch_probabilities
    if (is.null(branch_prob_df)) {
        stop(
            "`tata_result` does not contain branch probabilities.",
            call. = FALSE
        )
    }

    if (is.null(cell_ids)) {
        rank_idx <- order(
            branch_prob_df$tata_branch_entropy,
            decreasing = TRUE, na.last = NA
        )
        cell_ids <- branch_prob_df$cell_id[utils::head(rank_idx, n_cells)]
    }

    keep_df <- branch_prob_df[
        match(cell_ids, branch_prob_df$cell_id), ,
        drop = FALSE
    ]
    branch_cols <- grep("^tata_prob_", colnames(keep_df), value = TRUE)
    branch_names <- sub("^tata_prob_", "", branch_cols)

    plot_df <- stats::reshape(
        keep_df[, c("cell_id", branch_cols), drop = FALSE],
        varying = branch_cols,
        v.names = "probability",
        timevar = "branch",
        times = branch_names,
        direction = "long"
    )

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = branch, y = probability, fill = branch)
    ) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::facet_wrap(~cell_id) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(
            title = "Terminal-branch probabilities for representative cells",
            x = "branch",
            y = "probability",
            fill = "branch"
        )
}


#' Plot fitted TATA gene trends
#'
#' Plot fitted branch-specific gene-expression trends over TATA pseudotime.
#'
#' @param trend_result Result list returned by `compute_tata_gene_trends()`.
#' @param genes Optional subset of genes to plot.
#' @param branches Optional subset of branches to plot.
#' @param ncol Number of facet columns.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # Fit trends for a few early- and late-programme marker genes.
#' module <- SummarizedExperiment::rowData(tata$sce)$feature_module
#' genes <- head(rownames(tata$sce)[module %in% c("early", "late")], 4)
#' trends <- compute_tata_gene_trends(tata, genes = genes)
#' p <- plot_tata_gene_trends(trends, ncol = 2)
#' p
plot_tata_gene_trends <- function(
    trend_result,
    genes = NULL,
    branches = NULL,
    ncol = 2
) {
    trend_df <- trend_result$trend_table
    if (!is.null(genes)) {
        trend_df <- trend_df[trend_df$gene %in% genes, , drop = FALSE]
    }
    if (!is.null(branches)) {
        trend_df <- trend_df[trend_df$branch %in% branches, , drop = FALSE]
    }

    ggplot2::ggplot(
        trend_df,
        ggplot2::aes(
            x = pseudotime, y = fitted_expression,
            colour = branch, group = branch
        )
    ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_wrap(~gene, scales = "free_y", ncol = ncol) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "Branch-specific TATA gene trends",
            x = "TATA pseudotime",
            y = "fitted expression",
            colour = "branch"
        )
}


#' Plot a heatmap of TATA gene trends
#'
#' Plot a heatmap of branch-specific fitted gene-expression trends over TATA
#' pseudotime.
#'
#' @param trend_result Result list returned by `compute_tata_gene_trends()`.
#' @param branches Optional subset of branches.
#' @param scale_rows Logical indicating whether to z-score each gene-branch
#'   trend before plotting.
#' @param n_clusters Optional number of k-means clusters used to order genes by
#'   trend shape.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' # Z-scored trends, with genes grouped by k-means on trend shape.
#' trends <- compute_tata_gene_trends(tata,
#'     genes = head(rownames(tata$sce), 12)
#' )
#' p <- plot_tata_gene_trend_heatmap(trends,
#'     scale_rows = TRUE,
#'     n_clusters = 3
#' )
#' p
plot_tata_gene_trend_heatmap <- function(
    trend_result,
    branches = NULL,
    scale_rows = TRUE,
    n_clusters = NULL
) {
    trend_df <- trend_result$trend_table
    if (!is.null(branches)) {
        trend_df <- trend_df[trend_df$branch %in% branches, , drop = FALSE]
    }

    trend_df$gene_branch <- paste(trend_df$gene, trend_df$branch, sep = " | ")

    if (isTRUE(scale_rows)) {
        split_idx <- split(seq_len(nrow(trend_df)), trend_df$gene_branch)
        for (idx in split_idx) {
            value <- trend_df$fitted_expression[idx]
            value_sd <- stats::sd(value, na.rm = TRUE)
            if (is.finite(value_sd) && value_sd > 0) {
                trend_df$fitted_expression[idx] <- as.numeric(scale(value))
            } else {
                trend_df$fitted_expression[idx] <- 0
            }
        }
        trend_df$fitted_expression[is.na(trend_df$fitted_expression)] <- 0
    }

    if (!is.null(n_clusters)) {
        cluster_df <- cluster_tata_gene_trends(
            trend_result = list(trend_table = trend_df),
            n_clusters = n_clusters
        )
        trend_df <- merge(
            trend_df, cluster_df,
            by = "gene", all.x = TRUE, sort = FALSE
        )
        trend_df <- trend_df[
            order(
                trend_df$trend_cluster, trend_df$gene,
                trend_df$branch, trend_df$pseudotime
            ), ,
            drop = FALSE
        ]
    }

    gene_levels <- unique(trend_df$gene_branch)
    trend_df$gene_branch <- factor(
        trend_df$gene_branch,
        levels = rev(gene_levels)
    )

    ggplot2::ggplot(
        trend_df,
        ggplot2::aes(x = pseudotime, y = gene_branch, fill = fitted_expression)
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
            low = "#2166ac", mid = "white", high = "#b2182b"
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "Heatmap of branch-specific TATA gene trends",
            x = "TATA pseudotime",
            y = "gene | branch",
            fill = if (isTRUE(scale_rows)) "z-score" else "expression"
        )
}


#' Plot an overview of TATA results
#'
#' Create a compact set of standard TATA result plots analogous to a trajectory
#' overview dashboard.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param dimred Reduced dimension to use for the embedding plots.
#' @param n_probability_cells Number of representative cells to show in the
#'   terminal-probability barplot.
#'
#' @return A named list of `ggplot2` objects.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#'
#' plots <- plot_tata_results(tata, dimred = "UMAP")
#' names(plots)
#' plots$pseudotime
plot_tata_results <- function(
    tata_result,
    dimred = "UMAP",
    n_probability_cells = 6
) {
    sce <- tata_result$sce
    embedding <- as.matrix(SingleCellExperiment::reducedDim(sce, dimred))
    plot_df <- data.frame(
        Dim1 = embedding[, 1],
        Dim2 = embedding[, 2],
        pseudotime = SummarizedExperiment::colData(sce)$tata_pseudotime_scaled,
        entropy = SummarizedExperiment::colData(sce)$tata_branch_entropy,
        stringsAsFactors = FALSE
    )

    pseudotime_plot <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = Dim1, y = Dim2, colour = pseudotime)
    ) +
        ggplot2::geom_point(size = 0.6, alpha = 0.8) +
        ggplot2::scale_colour_viridis_c(option = "plasma") +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "TATA pseudotime",
            x = paste0(dimred, "_1"),
            y = paste0(dimred, "_2"),
            colour = "pseudotime"
        )

    entropy_plot <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = Dim1, y = Dim2, colour = entropy)
    ) +
        ggplot2::geom_point(size = 0.6, alpha = 0.8) +
        ggplot2::scale_colour_viridis_c(option = "cividis") +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = "TATA branch entropy",
            x = paste0(dimred, "_1"),
            y = paste0(dimred, "_2"),
            colour = "entropy"
        )

    list(
        pseudotime = pseudotime_plot,
        entropy = entropy_plot,
        time_calibration = plot_tata_time_calibration(tata_result),
        edge_causality = plot_tata_edge_causality(tata_result),
        branch_probabilities = plot_tata_branch_probabilities(
            tata_result,
            dimred = dimred
        ),
        representative_cells = plot_tata_terminal_probabilities(
            tata_result,
            n_cells = n_probability_cells
        )
    )
}


#' Resolve the embedding context for the TATA cluster-graph layout
#'
#' Validates the reduced dimension and cluster column, then returns the
#' per-cell metadata used to place cluster nodes on the embedding.
#'
#' @keywords internal
#' @noRd
.tata_cluster_embedding_context <- function(
    tata_result, dimred, cluster_col
) {
    sce <- tata_result$sce
    if (!methods::is(sce, "SingleCellExperiment")) {
        stop(
            "`layout_mode = \"embedding\"` requires `tata_result$sce` to ",
            "be a SingleCellExperiment.",
            call. = FALSE
        )
    }

    if (is.null(dimred)) {
        if ("UMAP" %in% SingleCellExperiment::reducedDimNames(sce)) {
            dimred <- "UMAP"
        } else {
            stop(
                "Please supply `dimred` for `layout_mode = \"embedding\"`.",
                call. = FALSE
            )
        }
    }

    if (!dimred %in% SingleCellExperiment::reducedDimNames(sce)) {
        stop(
            "Reduced dimension `", dimred,
            "` is not present in `tata_result$sce`.",
            call. = FALSE
        )
    }

    if (is.null(cluster_col)) {
        cluster_col <- tata_result$parameters$cluster_col
    }

    if (is.null(cluster_col) || !cluster_col %in%
            colnames(SummarizedExperiment::colData(sce))) {
        stop(
            "A valid `cluster_col` is required for `layout_mode = ",
            "\"embedding\"`.",
            call. = FALSE
        )
    }

    embedding <- as.matrix(SingleCellExperiment::reducedDim(sce, dimred))
    if (ncol(embedding) < 2L) {
        stop(
            "The requested reduced dimension must contain at least two ",
            "columns.",
            call. = FALSE
        )
    }

    meta <- as.data.frame(SummarizedExperiment::colData(sce))
    meta$cluster_plot <- as.character(meta[[cluster_col]])
    meta$Dim1 <- embedding[, 1]
    meta$Dim2 <- embedding[, 2]

    list(meta = meta, dimred = dimred, cluster_col = cluster_col)
}


#' Place cluster nodes at their embedding centroids
#'
#' @keywords internal
#' @noRd
.tata_cluster_embedding_nodes <- function(meta, graph_vertex_df) {
    node_df <- stats::aggregate(
        meta[, c("Dim1", "Dim2")],
        by = list(cluster = meta$cluster_plot),
        FUN = stats::median
    )
    names(node_df)[2:3] <- c("x", "y")
    node_df <- merge(
        graph_vertex_df, node_df,
        by = "cluster", all.x = TRUE, sort = FALSE
    )
    node_df[
        match(graph_vertex_df$cluster, node_df$cluster), ,
        drop = FALSE
    ]
}


#' Build the background cell table for the TATA cluster graph
#'
#' @keywords internal
#' @noRd
.tata_cluster_cell_table <- function(meta, colour_by) {
    cell_df <- meta[, c("Dim1", "Dim2", "cluster_plot"), drop = FALSE]
    if (!is.null(colour_by) && colour_by %in% colnames(meta)) {
        cell_df$colour_value <- meta[[colour_by]]
    } else {
        cell_df$colour_value <- "cells"
    }
    cell_df
}


#' Place cluster nodes using the stored abstract graph coordinates
#'
#' Falls back to a force-directed layout when the stored coordinates are
#' missing or non-finite.
#'
#' @keywords internal
#' @noRd
.tata_cluster_abstract_nodes <- function(cluster_graph, graph_vertex_df) {
    layout_mat <- cbind(
        igraph::V(cluster_graph)$x,
        igraph::V(cluster_graph)$y
    )
    if (any(!is.finite(layout_mat))) {
        layout_mat <- igraph::layout_with_fr(cluster_graph)
    }

    node_df <- graph_vertex_df
    node_df$x <- layout_mat[, 1]
    node_df$y <- layout_mat[, 2]
    node_df
}


#' Assemble the drawable edge table for the TATA cluster graph
#'
#' Attaches endpoint coordinates and the scaled plotting width to each edge.
#'
#' @keywords internal
#' @noRd
.tata_cluster_edge_table <- function(
    cluster_graph, node_df, edge_width_scale
) {
    edge_df <- igraph::as_data_frame(cluster_graph, what = "edges")
    if (nrow(edge_df) > 0L) {
        from_idx <- match(edge_df$from, node_df$cluster)
        to_idx <- match(edge_df$to, node_df$cluster)
        edge_df$x <- node_df$x[from_idx]
        edge_df$y <- node_df$y[from_idx]
        edge_df$xend <- node_df$x[to_idx]
        edge_df$yend <- node_df$y[to_idx]
    } else {
        edge_df <- data.frame(
            from = character(0),
            to = character(0),
            weight = numeric(0),
            ambiguous = logical(0),
            x = numeric(0),
            y = numeric(0),
            xend = numeric(0),
            yend = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    if (!"ambiguous" %in% colnames(edge_df)) {
        edge_df$ambiguous <- FALSE
    }

    if (nrow(edge_df) > 0L && max(edge_df$weight, na.rm = TRUE) > 0) {
        edge_df$plot_width <- 0.5 + edge_width_scale * edge_df$weight /
            max(edge_df$weight, na.rm = TRUE)
    } else {
        edge_df$plot_width <- rep(0.5, nrow(edge_df))
    }

    edge_df
}


#' Derive node sizes and the cluster colour palette
#'
#' @keywords internal
#' @noRd
.tata_cluster_node_aesthetics <- function(node_df) {
    if (length(unique(node_df$size)) > 1L) {
        node_df$plot_size <- 4 + 6 * (node_df$size - min(node_df$size)) /
            diff(range(node_df$size))
    } else {
        node_df$plot_size <- rep(7, nrow(node_df))
    }

    node_df$cluster <- factor(node_df$cluster, levels = node_df$cluster)
    cluster_palette <- grDevices::hcl.colors(nrow(node_df), palette = "Dark 3")
    names(cluster_palette) <- as.character(node_df$cluster)

    list(node_df = node_df, palette = cluster_palette)
}


#' Add the background cell layer to the TATA cluster-graph plot
#'
#' @keywords internal
#' @noRd
.tata_add_cluster_cell_layer <- function(
    p, cell_df, colour_by, cluster_col, cluster_palette, cluster_levels,
    point_size, point_alpha
) {
    if (is.null(cell_df) || nrow(cell_df) == 0L) {
        return(p)
    }

    if (!is.null(colour_by) && colour_by == cluster_col) {
        cell_df$colour_value <- factor(
            as.character(cell_df$colour_value),
            levels = cluster_levels
        )
        p <- p + ggplot2::geom_point(
            data = cell_df,
            mapping = ggplot2::aes(
                x = Dim1, y = Dim2, colour = colour_value
            ),
            size = point_size,
            alpha = point_alpha
        ) +
            ggplot2::scale_colour_manual(
                values = cluster_palette, name = colour_by
            )
    } else if (!is.null(colour_by) && is.numeric(cell_df$colour_value)) {
        p <- p + ggplot2::geom_point(
            data = cell_df,
            mapping = ggplot2::aes(
                x = Dim1, y = Dim2, colour = colour_value
            ),
            size = point_size,
            alpha = point_alpha
        ) +
            ggplot2::scale_colour_viridis_c(
                option = "plasma", name = colour_by
            )
    } else if (!is.null(colour_by) && !is.numeric(cell_df$colour_value)) {
        p <- p + ggplot2::geom_point(
            data = cell_df,
            mapping = ggplot2::aes(
                x = Dim1, y = Dim2, colour = colour_value
            ),
            size = point_size,
            alpha = point_alpha
        )
    } else {
        p <- p + ggplot2::geom_point(
            data = cell_df,
            mapping = ggplot2::aes(x = Dim1, y = Dim2),
            colour = "grey80",
            size = point_size,
            alpha = point_alpha
        )
    }

    p
}


#' Add the ambiguous and directed edge layers to the cluster-graph plot
#'
#' @keywords internal
#' @noRd
.tata_add_cluster_edge_layers <- function(p, edge_df, arrow_size) {
    directed_edges <- edge_df[!edge_df$ambiguous, , drop = FALSE]
    ambiguous_edges <- edge_df[edge_df$ambiguous, , drop = FALSE]

    if (nrow(ambiguous_edges) > 0L) {
        p <- p + ggplot2::geom_segment(
            data = ambiguous_edges,
            mapping = ggplot2::aes(
                x = x, y = y, xend = xend, yend = yend,
                linewidth = plot_width
            ),
            inherit.aes = FALSE,
            colour = "grey55",
            linetype = 2,
            alpha = 0.8
        )
    }

    if (nrow(directed_edges) > 0L) {
        p <- p + ggplot2::geom_segment(
            data = directed_edges,
            mapping = ggplot2::aes(
                x = x, y = y, xend = xend, yend = yend,
                linewidth = plot_width
            ),
            inherit.aes = FALSE,
            colour = "black",
            alpha = 0.85,
            arrow = grid::arrow(
                length = grid::unit(arrow_size, "inches"),
                type = "closed"
            )
        )
    }

    p
}


#' Add the cluster node and label layers to the cluster-graph plot
#'
#' @keywords internal
#' @noRd
.tata_add_cluster_node_layers <- function(
    p, node_df, cluster_palette, node_colour_by, vertex_label_size
) {
    if (identical(node_colour_by, "cluster")) {
        p <- p + ggplot2::geom_point(
            data = node_df,
            mapping = ggplot2::aes(
                x = x, y = y, size = plot_size, fill = cluster
            ),
            inherit.aes = FALSE,
            shape = 21,
            colour = "black",
            stroke = 0.5
        ) +
            ggplot2::scale_fill_manual(
                values = cluster_palette, name = "cluster"
            )
    } else {
        p <- p + ggplot2::geom_point(
            data = node_df,
            mapping = ggplot2::aes(
                x = x, y = y, size = plot_size, fill = median_time
            ),
            inherit.aes = FALSE,
            shape = 21,
            colour = "black",
            stroke = 0.5
        ) +
            ggplot2::scale_fill_gradientn(
                colours = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"),
                name = "median_time"
            )
    }

    p + ggplot2::geom_text(
        data = node_df,
        mapping = ggplot2::aes(x = x, y = y, label = cluster),
        inherit.aes = FALSE,
        size = vertex_label_size,
        nudge_y = 0.12
    )
}


#' Validate the `max_paths` argument of the trajectory embedding plot
#'
#' @keywords internal
#' @noRd
.tata_trajectory_max_paths <- function(max_paths) {
    if (is.null(max_paths)) {
        return(NULL)
    }

    max_paths <- as.integer(max_paths)
    if (length(max_paths) != 1L || is.na(max_paths) || max_paths < 1L) {
        stop(
            "`max_paths` must be NULL or a single positive integer.",
            call. = FALSE
        )
    }

    max_paths
}


#' Build the per-cell metadata used by the trajectory embedding plot
#'
#' @keywords internal
#' @noRd
.tata_trajectory_meta <- function(sce, dimred, cluster_col) {
    if (!dimred %in% SingleCellExperiment::reducedDimNames(sce)) {
        stop(
            "Reduced dimension `", dimred,
            "` is not present in `tata_result$sce`.",
            call. = FALSE
        )
    }

    embedding <- as.matrix(SingleCellExperiment::reducedDim(sce, dimred))
    if (ncol(embedding) < 2L) {
        stop(
            "The requested reduced dimension must contain at least two ",
            "columns.",
            call. = FALSE
        )
    }

    meta <- as.data.frame(SummarizedExperiment::colData(sce))
    meta$Dim1 <- embedding[, 1]
    meta$Dim2 <- embedding[, 2]
    meta$cluster_plot <- meta[[cluster_col]]
    meta
}


#' Compute cluster centroids joined to the cluster pseudotime table
#'
#' @keywords internal
#' @noRd
.tata_trajectory_centroids <- function(meta, cluster_pt) {
    cluster_centroids <- stats::aggregate(
        meta[, c("Dim1", "Dim2")],
        by = list(cluster = as.character(meta$cluster_plot)),
        FUN = stats::median
    )

    if (!is.null(cluster_pt) && nrow(cluster_pt) > 0L) {
        cluster_centroids <- merge(
            cluster_centroids, cluster_pt,
            by = "cluster", all.x = TRUE, sort = FALSE
        )
    } else {
        cluster_centroids$pseudotime <- NA_real_
        cluster_centroids$scaled_pseudotime <- NA_real_
        cluster_centroids$reachable <- TRUE
        cluster_centroids$is_root <- FALSE
    }

    cluster_centroids
}


#' Identify the reachable clusters present in the cluster graph
#'
#' @keywords internal
#' @noRd
.tata_trajectory_reachable_clusters <- function(
    cluster_centroids, cluster_graph
) {
    reachable_clusters <- cluster_centroids$cluster[
        (cluster_centroids$reachable %in% TRUE) |
            is.na(cluster_centroids$reachable)
    ]
    intersect(reachable_clusters, igraph::V(cluster_graph)$name)
}


#' Choose the root cluster used to seed the trajectory paths
#'
#' @keywords internal
#' @noRd
.tata_trajectory_root_cluster <- function(cluster_centroids) {
    if (!any(cluster_centroids$is_root, na.rm = TRUE)) {
        return(cluster_centroids$cluster[
            which.min(cluster_centroids$pseudotime)
        ])
    }

    cluster_centroids$cluster[which(cluster_centroids$is_root)[1]]
}


#' Identify the terminal clusters of the TATA trajectory
#'
#' A reachable cluster is terminal when it is not the root and has no
#' retained outgoing edge to a later-pseudotime cluster.
#'
#' @keywords internal
#' @noRd
.tata_trajectory_terminal_clusters <- function(
    cluster_centroids, directed_edges, reachable_clusters, root_cluster,
    cluster_pt_lookup
) {
    reachable_order <- cluster_centroids$cluster[
        order(cluster_centroids$pseudotime, decreasing = TRUE, na.last = NA)
    ]
    reachable_order <- intersect(reachable_order, reachable_clusters)

    terminal_clusters <- reachable_order[vapply(
        reachable_order,
        FUN.VALUE = logical(1),
        FUN = function(cluster_name) {
            if (identical(cluster_name, root_cluster)) {
                return(FALSE)
            }

            outgoing <- directed_edges$to[directed_edges$from == cluster_name]
            outgoing <- intersect(outgoing, reachable_clusters)
            if (length(outgoing) == 0L) {
                return(TRUE)
            }

            later_outgoing <- outgoing[
                cluster_pt_lookup[outgoing] > cluster_pt_lookup[cluster_name]
            ]
            length(later_outgoing) == 0L
        }
    )]

    if (length(terminal_clusters) == 0L) {
        terminal_clusters <- utils::head(
            setdiff(reachable_order, root_cluster), 3
        )
    }

    terminal_clusters
}


#' Derive the drawable trajectory curves from the TATA cluster graph
#'
#' Restricts the cluster graph to the reachable clusters, finds the
#' root-to-terminal routes and converts the retained routes into curves.
#'
#' @keywords internal
#' @noRd
.tata_trajectory_curve_table <- function(
    tata_result, cluster_centroids, max_paths, smooth_paths, n_curve_points
) {
    cluster_graph <- tata_result$cluster_graph
    reachable_clusters <- .tata_trajectory_reachable_clusters(
        cluster_centroids, cluster_graph
    )
    root_cluster <- .tata_trajectory_root_cluster(cluster_centroids)

    subgraph <- igraph::induced_subgraph(
        cluster_graph,
        vids = reachable_clusters
    )
    undirected_subgraph <- igraph::as_undirected(subgraph, mode = "collapse")

    edge_table <- tata_result$cluster_edge_table
    directed_edges <- edge_table[
        edge_table$kept & edge_table$direction != "ambiguous",
        c("from", "to"),
        drop = FALSE
    ]

    cluster_pt_lookup <- stats::setNames(
        cluster_centroids$pseudotime,
        cluster_centroids$cluster
    )

    terminal_clusters <- .tata_trajectory_terminal_clusters(
        cluster_centroids, directed_edges, reachable_clusters,
        root_cluster, cluster_pt_lookup
    )

    path_info_list <- .tata_trajectory_paths(
        subgraph, undirected_subgraph, terminal_clusters, root_cluster,
        cluster_centroids, cluster_pt_lookup
    )
    path_info_list <- .tata_trajectory_select_paths(
        path_info_list, max_paths
    )

    .tata_trajectory_curves(path_info_list, smooth_paths, n_curve_points)
}


#' Find the vertices of one root-to-terminal route
#'
#' Uses the directed subgraph first and falls back to its undirected
#' collapse when no directed route exists.
#'
#' @keywords internal
#' @noRd
.tata_trajectory_path_vertices <- function(
    subgraph, undirected_subgraph, root_cluster, terminal_cluster
) {
    path_vertices <- .tata_shortest_path_vertices(
        subgraph, root_cluster, terminal_cluster, mode = "out"
    )

    if (length(path_vertices) == 0L) {
        path_vertices <- .tata_shortest_path_vertices(
            undirected_subgraph, root_cluster, terminal_cluster, mode = "all"
        )
    }

    path_vertices
}


#' Extract the root-to-terminal trajectory paths
#'
#' @keywords internal
#' @noRd
.tata_trajectory_paths <- function(
    subgraph, undirected_subgraph, terminal_clusters, root_cluster,
    cluster_centroids, cluster_pt_lookup
) {
    path_info_list <- vector("list", length(terminal_clusters))
    for (i in seq_along(terminal_clusters)) {
        terminal_cluster <- terminal_clusters[i]
        path_vertices <- .tata_trajectory_path_vertices(
            subgraph, undirected_subgraph, root_cluster, terminal_cluster
        )

        path_clusters <- igraph::V(subgraph)$name[as.integer(path_vertices)]
        path_df <- cluster_centroids[
            match(path_clusters, cluster_centroids$cluster),
            c("cluster", "Dim1", "Dim2", "scaled_pseudotime"),
            drop = FALSE
        ]
        path_df <- path_df[
            stats::complete.cases(path_df[, c("Dim1", "Dim2")]), ,
            drop = FALSE
        ]
        if (nrow(path_df) < 2L) {
            next
        }

        terminal_pt <- cluster_pt_lookup[terminal_cluster]
        if (is.na(terminal_pt)) {
            terminal_pt <- nrow(path_df)
        }

        path_info_list[[i]] <- list(
            terminal_cluster = terminal_cluster,
            terminal_pseudotime = terminal_pt,
            path_clusters = path_clusters,
            path_df = path_df
        )
    }

    Filter(Negate(is.null), path_info_list)
}


#' Keep the `max_paths` latest-terminating trajectory paths
#'
#' @keywords internal
#' @noRd
.tata_trajectory_select_paths <- function(path_info_list, max_paths) {
    if (is.null(max_paths) || length(path_info_list) <= max_paths) {
        return(path_info_list)
    }

    path_rank <- order(
        vapply(
            path_info_list, function(x) x$terminal_pseudotime,
            numeric(1)
        ),
        decreasing = TRUE
    )
    path_info_list[path_rank[seq_len(max_paths)]]
}


#' Turn trajectory paths into drawable curves
#'
#' @keywords internal
#' @noRd
.tata_trajectory_curves <- function(
    path_info_list, smooth_paths, n_curve_points
) {
    curve_list <- vector("list", length(path_info_list))
    for (i in seq_along(path_info_list)) {
        path_df <- path_info_list[[i]]$path_df

        if (isTRUE(smooth_paths) && nrow(path_df) >= 4L) {
            spline_index <- seq_len(nrow(path_df))
            x_spline <- stats::smooth.spline(
                spline_index, path_df$Dim1, spar = 0.5
            )
            y_spline <- stats::smooth.spline(
                spline_index, path_df$Dim2, spar = 0.5
            )
            new_index <- seq(
                min(spline_index), max(spline_index),
                length.out = n_curve_points
            )
            curve_df <- data.frame(
                path_id = paste0("path_", i),
                terminal_cluster = path_info_list[[i]]$terminal_cluster,
                Dim1 = stats::predict(x_spline, x = new_index)$y,
                Dim2 = stats::predict(y_spline, x = new_index)$y,
                stringsAsFactors = FALSE
            )
        } else {
            curve_df <- data.frame(
                path_id = paste0("path_", i),
                terminal_cluster = path_info_list[[i]]$terminal_cluster,
                Dim1 = path_df$Dim1,
                Dim2 = path_df$Dim2,
                stringsAsFactors = FALSE
            )
        }

        curve_list[[i]] <- curve_df
    }

    do.call(rbind, Filter(Negate(is.null), curve_list))
}


#' Assemble the time-calibration table of paired time and pseudotime
#'
#' @keywords internal
#' @noRd
.tata_time_calibration_table <- function(
    meta, time_col, pseudotime_col, branch_col
) {
    plot_df <- data.frame(
        timepoint = .coerce_time_to_numeric(meta[[time_col]]),
        pseudotime = as.numeric(meta[[pseudotime_col]]),
        stringsAsFactors = FALSE
    )
    if (!is.null(branch_col) && branch_col %in% colnames(meta)) {
        plot_df$branch <- as.character(meta[[branch_col]])
    } else {
        plot_df$branch <- "all"
    }

    plot_df[
        is.finite(plot_df$timepoint) & is.finite(plot_df$pseudotime), ,
        drop = FALSE
    ]
}


#' Summarise pseudotime by timepoint and branch as median with IQR
#'
#' @keywords internal
#' @noRd
.tata_time_calibration_summary <- function(plot_df) {
    summary_df <- stats::aggregate(
        plot_df$pseudotime,
        by = list(timepoint = plot_df$timepoint, branch = plot_df$branch),
        FUN = function(x) {
            c(
                median = stats::median(x),
                q25 = stats::quantile(x, 0.25),
                q75 = stats::quantile(x, 0.75)
            )
        }
    )
    summary_df <- do.call(data.frame, summary_df)
    colnames(summary_df) <- c("timepoint", "branch", "median", "q25", "q75")
    summary_df
}
