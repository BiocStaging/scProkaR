#' Compute TATA pseudotime
#'
#' Compute a simple cluster-level pseudotime on the directed TATA graph using
#' shortest-path distances from the root cluster with inverse edge weights, then
#' propagate that pseudotime to cells using within-cluster time ranks. If some
#' clusters are not reachable in the directed graph, TATA falls back to the
#' undirected abstract topology for those clusters so late branches are not
#' dropped simply because an intermediate direction is ambiguous.
#'
#' @param cluster_graph Directed cluster-level TATA graph.
#' @param sce A `SingleCellExperiment`.
#' @param cluster_col Column in `colData(sce)` containing inferred clusters.
#' @param time_col Column in `colData(sce)` containing experimental time labels.
#' @param root_cluster Optional root cluster. If `NULL`, the root is chosen as
#'   the cluster with the earliest median timepoint.
#' @param within_cluster_fraction Scale factor used to spread cells within each
#'   cluster after assigning cluster-level pseudotime.
#' @param eps Small constant to avoid division by zero.
#'
#' @return A list with root cluster, cluster pseudotime table, raw cell
#'   pseudotime, scaled cell pseudotime, and the within-cluster step size.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' tata <- run_tata(sce, dimred = "PCA", k = 15)
#'
#' # Recompute pseudotime directly from the abstracted cluster graph.
#' pt <- compute_tata_pseudotime(
#'     cluster_graph = tata$cluster_graph,
#'     sce = tata$sce,
#'     cluster_col = "cluster",
#'     time_col = "timepoint"
#' )
#'
#' pt$root_cluster
#' head(pt$cluster_pseudotime)
#' summary(pt$cell_pseudotime_scaled)
compute_tata_pseudotime <- function(
    cluster_graph,
    sce,
    cluster_col = "cluster",
    time_col = "timepoint",
    root_cluster = NULL,
    within_cluster_fraction = 0.5,
    eps = 1e-8
) {
    if (!inherits(cluster_graph, "igraph")) {
        stop("`cluster_graph` must be an igraph object.", call. = FALSE)
    }

    meta <- as.data.frame(colData(sce))
    .tata_check_pseudotime_columns(meta, cluster_col, time_col)

    inputs <- .tata_pseudotime_inputs(
        meta,
        cluster_graph,
        cluster_col,
        time_col
    )
    clusters <- inputs$clusters
    has_time <- inputs$has_time
    time_numeric <- inputs$time_numeric
    cluster_size <- inputs$cluster_size
    cluster_median_time <- inputs$cluster_median_time

    if (is.null(root_cluster)) {
        root_cluster <- .tata_choose_root_cluster(
            cluster_graph,
            clusters,
            time_numeric,
            cluster_size,
            cluster_median_time,
            has_time,
            eps
        )
    }

    if (!root_cluster %in% igraph::V(cluster_graph)$name) {
        stop(
            "`root_cluster` is not present in the cluster graph.",
            call. = FALSE
        )
    }

    distances <- .tata_cluster_distances(
        cluster_graph,
        root_cluster,
        has_time,
        cluster_median_time,
        time_numeric,
        eps
    )
    directed_distance <- distances$directed
    undirected_distance <- distances$undirected
    cluster_pseudotime <- distances$pseudotime

    use_fallback <- !is.finite(directed_distance) &
        is.finite(undirected_distance)
    cluster_pseudotime <- ifelse(
        is.finite(cluster_pseudotime),
        cluster_pseudotime,
        NA_real_
    )
    local_step <- .tata_within_cluster_step(cluster_pseudotime)

    cell_pseudotime <- .tata_cell_pseudotime(
        clusters,
        cluster_pseudotime,
        time_numeric,
        has_time,
        within_cluster_fraction,
        local_step,
        rownames(meta)
    )
    cell_pseudotime_scaled <- .scale_to_unit(cell_pseudotime)
    names(cell_pseudotime_scaled) <- names(cell_pseudotime)

    cluster_pt_table <- .tata_cluster_pt_table(
        cluster_pseudotime,
        cluster_median_time,
        cluster_size,
        directed_distance,
        use_fallback,
        root_cluster
    )

    list(
        root_cluster = root_cluster,
        cluster_pseudotime = cluster_pt_table,
        cell_pseudotime = cell_pseudotime,
        cell_pseudotime_scaled = cell_pseudotime_scaled,
        within_cluster_step = local_step
    )
}


#' Validate the colData columns required for TATA pseudotime
#'
#' @keywords internal
#' @noRd
.tata_check_pseudotime_columns <- function(meta, cluster_col, time_col) {
    if (!cluster_col %in% colnames(meta)) {
        stop(
            "Cluster column `", cluster_col, "` was not found in colData.",
            call. = FALSE
        )
    }
    if (!is.null(time_col) && !time_col %in% colnames(meta)) {
        stop(
            "Time column `", time_col, "` was not found in colData.",
            call. = FALSE
        )
    }
    invisible(NULL)
}


#' Assemble the cluster and time vectors used by TATA pseudotime
#'
#' @keywords internal
#' @noRd
.tata_pseudotime_inputs <- function(
    meta,
    cluster_graph,
    cluster_col,
    time_col
) {
    clusters <- factor(
        meta[[cluster_col]],
        levels = igraph::V(cluster_graph)$name
    )
    has_time <- !is.null(time_col)
    time_numeric <- if (isTRUE(has_time)) .coerce_time_to_numeric(
        meta[[time_col]]
    ) else rep(NA_real_, nrow(meta))

    cluster_size <- tapply(seq_along(clusters), clusters, length)
    cluster_median_time <- if (isTRUE(has_time)) {
        tapply(time_numeric, clusters, stats::median)
    } else {
        stats::setNames(
            rep(NA_real_, length(levels(clusters))),
            levels(clusters)
        )
    }

    list(
        clusters = clusters,
        has_time = has_time,
        time_numeric = time_numeric,
        cluster_size = cluster_size,
        cluster_median_time = cluster_median_time
    )
}


#' Weighted closeness of every cluster in the collapsed abstract graph
#'
#' @keywords internal
#' @noRd
.tata_cluster_closeness <- function(cluster_graph, eps) {
    undirected_graph <- igraph::as_undirected(
        cluster_graph,
        mode = "collapse"
    )
    edge_weights <- igraph::E(undirected_graph)$weight
    if (length(edge_weights) > 0L) {
        return(
            igraph::closeness(
                undirected_graph,
                mode = "all",
                weights = 1 / pmax(edge_weights, eps),
                normalized = TRUE
            )
        )
    }

    closeness_score <- rep(0, igraph::vcount(undirected_graph))
    names(closeness_score) <- igraph::V(undirected_graph)$name
    closeness_score
}


#' Candidate root clusters holding cells from the earliest timepoint
#'
#' @keywords internal
#' @noRd
.tata_earliest_time_candidates <- function(
    clusters,
    time_numeric,
    cluster_median_time
) {
    earliest_time <- min(time_numeric, na.rm = TRUE)
    earliest_fraction <- tapply(
        time_numeric == earliest_time,
        clusters,
        mean
    )
    earliest_fraction[is.na(earliest_fraction)] <- 0

    candidate_mask <- earliest_fraction > 0
    candidates <- names(earliest_fraction)[candidate_mask]
    if (length(candidates) == 0L) {
        earliest_cluster_time <- min(cluster_median_time, na.rm = TRUE)
        candidates <- names(cluster_median_time)[
            cluster_median_time == earliest_cluster_time
        ]
    }

    list(
        cluster = candidates,
        earliest_fraction = as.numeric(earliest_fraction[candidates])
    )
}


#' Pick the TATA root cluster from time, closeness and size
#'
#' @keywords internal
#' @noRd
.tata_choose_root_cluster <- function(
    cluster_graph,
    clusters,
    time_numeric,
    cluster_size,
    cluster_median_time,
    has_time,
    eps
) {
    closeness_score <- .tata_cluster_closeness(cluster_graph, eps)

    if (isTRUE(has_time)) {
        candidates <- .tata_earliest_time_candidates(
            clusters,
            time_numeric,
            cluster_median_time
        )
        candidate_table <- data.frame(
            cluster = candidates$cluster,
            earliest_fraction = candidates$earliest_fraction,
            closeness = as.numeric(closeness_score[candidates$cluster]),
            size = as.numeric(cluster_size[candidates$cluster]),
            stringsAsFactors = FALSE
        )
        candidate_table <- candidate_table[
            order(
                -candidate_table$earliest_fraction,
                -candidate_table$closeness,
                -candidate_table$size
            ), ,
            drop = FALSE
        ]
    } else {
        candidate_table <- data.frame(
            cluster = names(closeness_score),
            closeness = as.numeric(closeness_score),
            size = as.numeric(cluster_size[names(closeness_score)]),
            stringsAsFactors = FALSE
        )
        candidate_table <- candidate_table[
            order(
                -candidate_table$closeness,
                -candidate_table$size
            ), ,
            drop = FALSE
        ]
    }

    candidate_table$cluster[1]
}


#' Weighted distances from the root on the collapsed undirected graph
#'
#' @keywords internal
#' @noRd
.tata_undirected_distance <- function(cluster_graph, root_cluster, eps) {
    undirected_graph <- igraph::as_undirected(
        cluster_graph,
        mode = "collapse"
    )
    undirected_inverse_weight <- 1 / pmax(
        igraph::E(undirected_graph)$weight,
        eps
    )
    undirected_distance <- as.numeric(
        igraph::distances(
            undirected_graph,
            v = root_cluster,
            to = igraph::V(undirected_graph),
            mode = "all",
            weights = undirected_inverse_weight
        )
    )
    names(undirected_distance) <- igraph::V(undirected_graph)$name
    undirected_distance
}


#' Cluster distances for an abstract graph without any edges
#'
#' @keywords internal
#' @noRd
.tata_empty_graph_distances <- function(cluster_graph, root_cluster) {
    directed_distance <- rep(Inf, igraph::vcount(cluster_graph))
    names(directed_distance) <- igraph::V(cluster_graph)$name
    directed_distance[root_cluster] <- 0

    list(
        directed = directed_distance,
        undirected = directed_distance,
        pseudotime = directed_distance
    )
}


#' Shortest-path cluster distances when no timepoints are available
#'
#' @keywords internal
#' @noRd
.tata_graph_distances <- function(cluster_graph, root_cluster, eps) {
    inverse_weight <- 1 / pmax(igraph::E(cluster_graph)$weight, eps)
    directed_distance <- as.numeric(
        igraph::distances(
            cluster_graph,
            v = root_cluster,
            to = igraph::V(cluster_graph),
            mode = "out",
            weights = inverse_weight
        )
    )
    names(directed_distance) <- igraph::V(cluster_graph)$name

    undirected_distance <- .tata_undirected_distance(
        cluster_graph,
        root_cluster,
        eps
    )

    cluster_pseudotime <- directed_distance
    use_fallback <- !is.finite(cluster_pseudotime) &
        is.finite(undirected_distance)
    cluster_pseudotime[use_fallback] <- undirected_distance[use_fallback]

    list(
        directed = directed_distance,
        undirected = undirected_distance,
        pseudotime = cluster_pseudotime
    )
}


#' Edge table of the abstract graph with default direction annotations
#'
#' @keywords internal
#' @noRd
.tata_annotated_edge_table <- function(cluster_graph) {
    edge_df <- igraph::as_data_frame(cluster_graph, what = "edges")
    if (!"ambiguous" %in% colnames(edge_df)) {
        edge_df$ambiguous <- FALSE
    }
    if (!"time_consistency" %in% colnames(edge_df)) {
        edge_df$time_consistency <- 1
    }
    edge_df
}


#' Order clusters by median time, keeping the root first
#'
#' @keywords internal
#' @noRd
.tata_cluster_time_order <- function(
    cluster_levels,
    median_lookup,
    root_cluster
) {
    order_table <- data.frame(
        cluster = cluster_levels,
        median_time = as.numeric(median_lookup[cluster_levels]),
        stringsAsFactors = FALSE
    )
    order_table$median_time[!is.finite(order_table$median_time)] <- max(
        order_table$median_time[is.finite(order_table$median_time)],
        na.rm = TRUE
    )
    order_table <- order_table[
        order(order_table$median_time, order_table$cluster), ,
        drop = FALSE
    ]

    unique(c(root_cluster, setdiff(order_table$cluster, root_cluster)))
}


#' Pseudotime increment contributed by each incoming abstract edge
#'
#' @keywords internal
#' @noRd
.tata_edge_increment <- function(incoming, typical_step, max_weight, eps) {
    delta_units <- pmax(
        incoming$delta_time / pmax(typical_step, eps),
        0
    )
    base_increment <- pmax(0.35, delta_units)
    weight_penalty <- 1 + 0.35 * (1 - incoming$weight / max_weight)
    ambiguity_penalty <- ifelse(incoming$ambiguous, 1.35, 1)
    consistency_penalty <- 1 + 0.25 * (
        1 - pmin(pmax(incoming$time_consistency, 0), 1)
    )

    base_increment * weight_penalty * ambiguity_penalty *
        consistency_penalty
}


#' Propagate cluster pseudotime forward along time-ordered edges
#'
#' @keywords internal
#' @noRd
.tata_propagate_pseudotime <- function(
    forward_edge_df,
    cluster_order,
    cluster_levels,
    root_cluster,
    typical_step,
    max_weight,
    eps
) {
    cluster_pseudotime <- stats::setNames(
        rep(NA_real_, length(cluster_levels)),
        cluster_levels
    )
    cluster_pseudotime[root_cluster] <- 0

    Reduce(function(cluster_pseudotime, cluster_name) {
        incoming <- forward_edge_df[
            forward_edge_df$to == cluster_name, ,
            drop = FALSE
        ]
        pred_pt <- cluster_pseudotime[incoming$from]
        valid <- is.finite(pred_pt)
        if (!any(valid)) {
            return(cluster_pseudotime)
        }

        pred_pt <- pred_pt[valid]
        edge_increment <- .tata_edge_increment(
            incoming[valid, , drop = FALSE],
            typical_step,
            max_weight,
            eps
        )
        cluster_pseudotime[cluster_name] <- max(
            pred_pt + edge_increment,
            na.rm = TRUE
        )
        cluster_pseudotime
    }, cluster_order[-1], cluster_pseudotime)
}


#' Fill unreachable clusters from the undirected topology and time anchor
#'
#' @keywords internal
#' @noRd
.tata_apply_time_fallback <- function(
    cluster_pseudotime,
    undirected_distance,
    time_anchor
) {
    finite_path <- cluster_pseudotime[is.finite(cluster_pseudotime)]
    max_path <- if (length(finite_path) > 0L) max(finite_path) else 0
    fallback_scale <- .scale_to_unit(undirected_distance)
    fallback_pt <- time_anchor + 0.75 + fallback_scale * (max_path + 0.75)

    use_fallback <- !is.finite(cluster_pseudotime) &
        is.finite(undirected_distance)
    cluster_pseudotime[use_fallback] <- fallback_pt[use_fallback]
    cluster_pseudotime[is.finite(cluster_pseudotime)] <- pmax(
        cluster_pseudotime[is.finite(cluster_pseudotime)],
        0.85 * time_anchor[is.finite(cluster_pseudotime)]
    )
    cluster_pseudotime
}


#' Dispatch cluster distances on graph edges and timepoint availability
#'
#' @keywords internal
#' @noRd
.tata_cluster_distances <- function(
    cluster_graph,
    root_cluster,
    has_time,
    cluster_median_time,
    time_numeric,
    eps
) {
    if (igraph::ecount(cluster_graph) == 0L) {
        return(.tata_empty_graph_distances(cluster_graph, root_cluster))
    }
    if (!isTRUE(has_time)) {
        return(.tata_graph_distances(cluster_graph, root_cluster, eps))
    }

    .tata_time_aware_distances(
        cluster_graph,
        root_cluster,
        cluster_median_time,
        time_numeric,
        eps
    )
}


#' Time-aware cluster pseudotime and its directed/undirected distances
#'
#' @keywords internal
#' @noRd
.tata_time_aware_distances <- function(
    cluster_graph,
    root_cluster,
    cluster_median_time,
    time_numeric,
    eps
) {
    edge_df <- .tata_annotated_edge_table(cluster_graph)
    cluster_levels <- igraph::V(cluster_graph)$name
    typical_step <- .typical_time_step(time_numeric)
    median_lookup <- stats::setNames(
        cluster_median_time[cluster_levels],
        cluster_levels
    )
    max_weight <- max(edge_df$weight, na.rm = TRUE)
    if (!is.finite(max_weight) || max_weight <= 0) {
        max_weight <- 1
    }

    cluster_order <- .tata_cluster_time_order(
        cluster_levels,
        median_lookup,
        root_cluster
    )
    rank_lookup <- stats::setNames(seq_along(cluster_order), cluster_order)

    edge_df$delta_time <- median_lookup[edge_df$to] -
        median_lookup[edge_df$from]
    edge_df$from_rank <- rank_lookup[edge_df$from]
    edge_df$to_rank <- rank_lookup[edge_df$to]
    forward_edge_df <- edge_df[
        edge_df$to_rank > edge_df$from_rank, ,
        drop = FALSE
    ]

    cluster_pseudotime <- .tata_propagate_pseudotime(
        forward_edge_df,
        cluster_order,
        cluster_levels,
        root_cluster,
        typical_step,
        max_weight,
        eps
    )

    time_anchor <- pmax(
        0,
        (median_lookup - median_lookup[root_cluster]) /
            pmax(typical_step, eps)
    )

    directed_distance <- cluster_pseudotime
    directed_distance[!is.finite(directed_distance)] <- Inf
    undirected_distance <- .tata_undirected_distance(
        cluster_graph,
        root_cluster,
        eps
    )

    list(
        directed = directed_distance,
        undirected = undirected_distance,
        pseudotime = .tata_apply_time_fallback(
            cluster_pseudotime,
            undirected_distance,
            time_anchor
        )
    )
}


#' Smallest positive gap between cluster pseudotime values
#'
#' @keywords internal
#' @noRd
.tata_within_cluster_step <- function(cluster_pseudotime) {
    finite_pt <- sort(unique(cluster_pseudotime[!is.na(cluster_pseudotime)]))
    positive_diffs <- diff(finite_pt)
    positive_diffs <- positive_diffs[positive_diffs > 0]

    if (length(positive_diffs) > 0L) min(positive_diffs) else 1
}


#' Rank cells within one cluster by their experimental time
#'
#' @keywords internal
#' @noRd
.tata_within_cluster_rank <- function(time_in_cluster, has_time) {
    n_cells <- length(time_in_cluster)
    finite_times <- time_in_cluster[is.finite(time_in_cluster)]
    if (!isTRUE(has_time) || length(unique(finite_times)) <= 1L) {
        return(rep(0, n_cells))
    }

    (rank(time_in_cluster, ties.method = "average") - 1) / (n_cells - 1)
}


#' Spread cluster-level pseudotime across the cells of each cluster
#'
#' @keywords internal
#' @noRd
.tata_cell_pseudotime <- function(
    clusters,
    cluster_pseudotime,
    time_numeric,
    has_time,
    within_cluster_fraction,
    local_step,
    cell_names
) {
    base_pt <- cluster_pseudotime[as.character(clusters)]
    cell_pseudotime <- base_pt + stats::ave(
        time_numeric,
        replace(clusters, is.na(base_pt), NA),
        FUN = function(x) {
            x[] <- within_cluster_fraction * local_step *
                .tata_within_cluster_rank(x, has_time)
            x
        }
    )
    names(cell_pseudotime) <- cell_names
    cell_pseudotime
}


#' Build the ordered cluster-level pseudotime summary table
#'
#' @keywords internal
#' @noRd
.tata_cluster_pt_table <- function(
    cluster_pseudotime,
    cluster_median_time,
    cluster_size,
    directed_distance,
    use_fallback,
    root_cluster
) {
    cluster_names <- names(cluster_pseudotime)
    cluster_pseudotime_scaled <- .scale_to_unit(cluster_pseudotime)
    names(cluster_pseudotime_scaled) <- cluster_names

    cluster_pt_table <- data.frame(
        cluster = cluster_names,
        pseudotime = as.numeric(cluster_pseudotime),
        scaled_pseudotime = as.numeric(
            cluster_pseudotime_scaled[cluster_names]
        ),
        median_time = as.numeric(cluster_median_time[cluster_names]),
        size = as.integer(cluster_size[cluster_names]),
        reachable = !is.na(cluster_pseudotime),
        reachable_directed = is.finite(directed_distance[cluster_names]),
        used_undirected_fallback = as.logical(use_fallback[cluster_names]),
        is_root = cluster_names == root_cluster,
        stringsAsFactors = FALSE
    )

    cluster_pt_table[
        order(cluster_pt_table$pseudotime, na.last = TRUE), ,
        drop = FALSE
    ]
}
