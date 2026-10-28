#' Compute topology-based cluster connectivity weights
#'
#' For each cluster pair `(a, b)`, count observed inter-cluster edges and
#' compare them to a simple degree-aware expectation.
#'
#' @param adjacency Sparse cell-cell adjacency matrix.
#' @param clusters Cluster labels for cells.
#' @param eps Small constant to avoid division by zero.
#'
#' @return A data frame of observed edges, expected edges, and topology weights.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' knn <- build_knn_graph(sce, dimred = "PCA", k = 12)
#' clusters <- cluster_graph_states(knn$graph, method = "louvain")
#' topology_table <- compute_topology_weights(
#'     adjacency = knn$adjacency,
#'     clusters = clusters
#' )
#' head(topology_table[order(-topology_table$topology_weight), ])
compute_topology_weights <- function(
    adjacency,
    clusters,
    eps = 1e-8
) {
    if (!inherits(adjacency, "sparseMatrix")) {
        adjacency <- as(adjacency, "dgCMatrix")
    }

    clusters <- factor(clusters)
    cluster_levels <- levels(clusters)
    n_clusters <- length(cluster_levels)

    membership_matrix <- Matrix::sparseMatrix(
        i = seq_along(clusters),
        j = as.integer(clusters),
        x = 1,
        dims = c(length(clusters), n_clusters)
    )

    cluster_adjacency <- Matrix::t(membership_matrix) %*%
        adjacency %*% membership_matrix
    dimnames(cluster_adjacency) <- list(cluster_levels, cluster_levels)

    degree_per_cell <- Matrix::rowSums(adjacency)
    degree_by_cluster <- tapply(degree_per_cell, clusters, sum)
    size_by_cluster <- table(clusters)
    total_degree <- sum(degree_per_cell)

    pair_grid <- utils::combn(cluster_levels, 2)
    do.call(rbind, lapply(seq_len(ncol(pair_grid)), function(idx) {
        cluster_a <- pair_grid[1, idx]
        cluster_b <- pair_grid[2, idx]
        observed_edges <- as.numeric(cluster_adjacency[cluster_a, cluster_b])
        expected_edges <- if (total_degree > 0) {
            (degree_by_cluster[[cluster_a]] *
                degree_by_cluster[[cluster_b]]) / total_degree
        } else {
            0
        }

        data.frame(
            cluster_a = cluster_a,
            cluster_b = cluster_b,
            observed_edges = observed_edges,
            degree_a = as.numeric(degree_by_cluster[[cluster_a]]),
            degree_b = as.numeric(degree_by_cluster[[cluster_b]]),
            size_a = as.integer(size_by_cluster[[cluster_a]]),
            size_b = as.integer(size_by_cluster[[cluster_b]]),
            expected_edges = expected_edges,
            topology_weight = observed_edges / (expected_edges + eps),
            stringsAsFactors = FALSE
        )
    }))
}


#' Compute time consistency for cluster-cluster edges
#'
#' For each inter-cluster cell-cell graph edge, compute `sign(t_b - t_a)` after
#' orienting the cluster pair as `(cluster_a, cluster_b)`. The mean sign is the
#' temporal flow score.
#'
#' @param adjacency Sparse cell-cell adjacency matrix.
#' @param clusters Cluster labels for cells.
#' @param timepoint Experimental time labels for cells.
#'
#' @return A data frame of temporal flow scores and time weights.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' knn <- build_knn_graph(sce, dimred = "PCA", k = 12)
#' clusters <- cluster_graph_states(knn$graph, method = "louvain")
#' time_table <- compute_time_weights(
#'     adjacency = knn$adjacency,
#'     clusters = clusters,
#'     timepoint = SummarizedExperiment::colData(sce)$timepoint
#' )
#' head(time_table[, c(
#'     "cluster_a", "cluster_b", "flow_score",
#'     "time_weight"
#' )])
compute_time_weights <- function(
    adjacency,
    clusters,
    timepoint
) {
    if (!inherits(adjacency, "sparseMatrix")) {
        adjacency <- as(adjacency, "dgCMatrix")
    }

    clusters <- factor(clusters)
    cluster_levels <- levels(clusters)
    cluster_id <- as.integer(clusters)
    time_numeric <- .coerce_time_to_numeric(timepoint)
    typical_step <- .typical_time_step(time_numeric)
    cluster_time_list <- split(time_numeric, clusters)
    cluster_median_time <- vapply(cluster_time_list, stats::median, numeric(1))

    pair_grid <- utils::combn(cluster_levels, 2)
    pair_key <- paste(pair_grid[1, ], pair_grid[2, ], sep = "||")
    flow_store <- .tata_pair_time_deltas(
        adjacency = adjacency,
        cluster_id = cluster_id,
        cluster_levels = cluster_levels,
        time_numeric = time_numeric,
        pair_key = pair_key
    )

    do.call(rbind, lapply(seq_len(ncol(pair_grid)), function(idx) {
        cluster_a <- pair_grid[1, idx]
        cluster_b <- pair_grid[2, idx]
        cluster_gap <- as.numeric(
            cluster_median_time[[cluster_b]] -
                cluster_median_time[[cluster_a]]
        )
        time_overlap <- .time_distribution_overlap(
            cluster_time_list[[cluster_a]],
            cluster_time_list[[cluster_b]]
        )

        pair_stats <- .tata_pair_time_stats(
            pair_delta = flow_store[[pair_key[idx]]],
            cluster_gap = cluster_gap,
            time_overlap = time_overlap,
            typical_step = typical_step
        )

        data.frame(
            cluster_a = cluster_a,
            cluster_b = cluster_b,
            n_time_edges = pair_stats$n_time_edges,
            flow_score = pair_stats$flow_score,
            signed_time_score = pair_stats$signed_time_score,
            median_time_gap = cluster_gap,
            edge_median_delta = pair_stats$edge_median_delta,
            time_overlap = time_overlap,
            jump_penalty = pair_stats$jump_penalty,
            time_weight = pair_stats$time_weight,
            time_consistency = pair_stats$time_consistency,
            stringsAsFactors = FALSE
        )
    }))
}


.make_neutral_time_table <- function(clusters) {
    cluster_levels <- levels(factor(clusters))
    pair_grid <- utils::combn(cluster_levels, 2)
    if (ncol(pair_grid) == 0L) {
        return(data.frame(
            cluster_a = character(0),
            cluster_b = character(0),
            n_time_edges = integer(0),
            flow_score = numeric(0),
            signed_time_score = numeric(0),
            median_time_gap = numeric(0),
            edge_median_delta = numeric(0),
            time_overlap = numeric(0),
            jump_penalty = numeric(0),
            time_weight = numeric(0),
            time_consistency = numeric(0),
            stringsAsFactors = FALSE
        ))
    }

    data.frame(
        cluster_a = pair_grid[1, ],
        cluster_b = pair_grid[2, ],
        n_time_edges = 0L,
        flow_score = 0,
        signed_time_score = 0,
        median_time_gap = 0,
        edge_median_delta = 0,
        time_overlap = 1,
        jump_penalty = 1,
        time_weight = 0.5,
        time_consistency = 1,
        stringsAsFactors = FALSE
    )
}


#' Build the TATA cluster graph
#'
#' Combine topology and temporal consistency into final TATA edge weights,
#' assign edge direction, and prune weak edges to obtain the cluster-level TATA
#' graph.
#'
#' @param topology_table Output from `compute_topology_weights()`.
#' @param time_table Output from `compute_time_weights()`.
#' @param clusters Cluster labels for cells.
#' @param embedding Low-dimensional embedding used to compute cluster centroids.
#' @param timepoint Experimental time labels.
#' @param alpha Exponent for topology weights.
#' @param beta Exponent for time weights.
#' @param direction_threshold Threshold for assigning edge direction.
#' @param prune_threshold Minimum TATA weight to retain an edge.
#' @param time_weight_mode How to incorporate time into the final edge weight.
#'   Choices are:
#' - `"directional_confidence"`: reward strong directional evidence regardless
#'   of whether the supported direction is `a -> b` or `b -> a`.
#' - `"ordered"`: use the original ordered-pair formulation
#'   `(flow_score + 1) / 2`, which is stricter but depends on the ordering of
#'   `cluster_a` and `cluster_b`.
#'
#' @return A list containing the directed cluster graph, the cluster edge table,
#'   and cluster vertex metadata.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' knn <- build_knn_graph(sce, dimred = "PCA", k = 12)
#' clusters <- cluster_graph_states(knn$graph, method = "louvain")
#' timepoint <- SummarizedExperiment::colData(sce)$timepoint
#' topology_table <- compute_topology_weights(knn$adjacency, clusters)
#' time_table <- compute_time_weights(knn$adjacency, clusters, timepoint)
#' tata_graph <- build_tata_graph(
#'     topology_table = topology_table,
#'     time_table = time_table,
#'     clusters = clusters,
#'     embedding = knn$embedding,
#'     timepoint = timepoint
#' )
#' igraph::ecount(tata_graph$graph)
#' head(tata_graph$vertex_table)
build_tata_graph <- function(
    topology_table,
    time_table,
    clusters,
    embedding,
    timepoint,
    alpha = 1,
    beta = 1,
    direction_threshold = 0.1,
    prune_threshold = 0.01,
    time_weight_mode = c("directional_confidence", "ordered")
) {
    time_weight_mode <- match.arg(time_weight_mode)

    if (!all(c("cluster_a", "cluster_b") %in% colnames(topology_table))) {
        stop(
            "`topology_table` is missing required cluster columns.",
            call. = FALSE
        )
    }

    if (!all(c("cluster_a", "cluster_b") %in% colnames(time_table))) {
        stop("`time_table` is missing required cluster columns.", call. = FALSE)
    }

    edge_table <- .tata_merge_edge_tables(topology_table, time_table)

    edge_table$time_factor <- if (time_weight_mode == "ordered") {
        pmax(edge_table$time_weight, 0)
    } else {
        pmax(edge_table$time_consistency, 0)
    }

    edge_table$tata_weight <- (pmax(edge_table$topology_weight, 0)^alpha) *
        (edge_table$time_factor^beta)

    edge_table <- .tata_assign_edge_direction(
        edge_table,
        direction_threshold
    )

    vertex_df <- .tata_vertex_table(clusters, embedding, timepoint)

    edge_table$kept <- edge_table$observed_edges > 0 &
        edge_table$tata_weight >= prune_threshold

    kept_edges <- edge_table[edge_table$kept, , drop = FALSE]
    graph_edge_df <- .tata_expand_graph_edges(kept_edges)

    graph_vertex_df <- data.frame(
        name = vertex_df$cluster,
        size = vertex_df$size,
        median_time = vertex_df$median_time,
        x = vertex_df$x,
        y = vertex_df$y,
        stringsAsFactors = FALSE
    )

    cluster_graph <- igraph::graph_from_data_frame(
        d = graph_edge_df,
        directed = TRUE,
        vertices = graph_vertex_df
    )

    list(
        graph = cluster_graph,
        edge_table = edge_table,
        vertex_table = vertex_df
    )
}


#' Collect per-cluster-pair time deltas across inter-cluster graph edges
#'
#' Walks the upper triangle of the adjacency matrix, keeps only edges that
#' connect two different clusters, and records `t_b - t_a` for each edge after
#' orienting the pair as `(cluster_a, cluster_b)`.
#'
#' @keywords internal
#' @noRd
.tata_pair_time_deltas <- function(
    adjacency,
    cluster_id,
    cluster_levels,
    time_numeric,
    pair_key
) {
    flow_store <- stats::setNames(vector("list", length(pair_key)), pair_key)

    upper_adjacency <- Matrix::triu(adjacency, k = 1)
    edge_summary <- Matrix::summary(upper_adjacency)

    if (nrow(edge_summary) > 0L) {
        same_cluster <- cluster_id[edge_summary$i] ==
            cluster_id[edge_summary$j]
        edge_summary <- edge_summary[!same_cluster, , drop = FALSE]
    }

    if (nrow(edge_summary) == 0L) {
        return(flow_store)
    }

    swap <- cluster_id[edge_summary$i] > cluster_id[edge_summary$j]
    if (anyNA(swap)) {
        stop("missing value where TRUE/FALSE needed", call. = FALSE)
    }
    from <- ifelse(swap, edge_summary$j, edge_summary$i)
    to <- ifelse(swap, edge_summary$i, edge_summary$j)
    pair_name <- paste(
        cluster_levels[cluster_id[from]],
        cluster_levels[cluster_id[to]],
        sep = "||"
    )
    deltas <- split(time_numeric[to] - time_numeric[from], pair_name)
    flow_store[names(deltas)] <- deltas
    flow_store
}


#' Summarise the temporal statistics of a single cluster pair
#'
#' Turns the raw per-edge time deltas of one cluster pair into the flow score,
#' directional consistency, jump penalty and derived time weight. Pairs with no
#' inter-cluster edges fall back to a neutral summary.
#'
#' @keywords internal
#' @noRd
.tata_pair_time_stats <- function(
    pair_delta,
    cluster_gap,
    time_overlap,
    typical_step
) {
    if (length(pair_delta) == 0L) {
        return(list(
            flow_score = 0,
            time_weight = 0.5,
            time_consistency = 0,
            signed_time_score = 0,
            edge_median_delta = cluster_gap,
            jump_penalty = 1,
            n_time_edges = 0L
        ))
    }

    flow_score <- mean(sign(pair_delta))
    edge_median_delta <- stats::median(pair_delta)
    lag_strength <- tanh(abs(cluster_gap) / pmax(typical_step, 1e-8))
    overlap_separation <- 1 - time_overlap
    jump_ratio <- abs(cluster_gap) / pmax(typical_step, 1e-8)
    jump_penalty <- 1 / (1 + pmax(jump_ratio - 2, 0) / 2)
    direction_source <- if (abs(flow_score) > 1e-8)
        flow_score else sign(cluster_gap)
    direction_strength <- max(
        abs(flow_score),
        min(1, abs(cluster_gap) / pmax(2 * typical_step, 1e-8))
    )
    time_consistency <- pmin(
        1,
        (0.60 * direction_strength + 0.25 * lag_strength +
            0.15 * overlap_separation) * jump_penalty
    )
    signed_time_score <- sign(direction_source) * time_consistency

    list(
        flow_score = flow_score,
        time_weight = (signed_time_score + 1) / 2,
        time_consistency = time_consistency,
        signed_time_score = signed_time_score,
        edge_median_delta = edge_median_delta,
        jump_penalty = jump_penalty,
        n_time_edges = length(pair_delta)
    )
}


#' Join the topology and time tables and fill missing temporal columns
#'
#' @keywords internal
#' @noRd
.tata_merge_edge_tables <- function(topology_table, time_table) {
    edge_table <- merge(
        topology_table,
        time_table,
        by = c("cluster_a", "cluster_b"),
        all.x = TRUE,
        sort = FALSE
    )

    edge_table$flow_score[is.na(edge_table$flow_score)] <- 0
    if (!"signed_time_score" %in% colnames(edge_table)) {
        edge_table$signed_time_score <- edge_table$flow_score
    }
    edge_table$signed_time_score[is.na(edge_table$signed_time_score)] <- 0
    edge_table$time_weight[is.na(edge_table$time_weight)] <- 0.5
    if (!"time_consistency" %in% colnames(edge_table)) {
        edge_table$time_consistency <- abs(edge_table$flow_score)
    }
    edge_table$time_consistency[is.na(edge_table$time_consistency)] <- 0
    edge_table$n_time_edges[is.na(edge_table$n_time_edges)] <- 0

    edge_table
}


#' Assign edge direction, endpoints and a printable direction label
#'
#' @keywords internal
#' @noRd
.tata_assign_edge_direction <- function(edge_table, direction_threshold) {
    edge_table$direction <- ifelse(
        edge_table$signed_time_score > direction_threshold,
        "a_to_b",
        ifelse(
            edge_table$signed_time_score < -direction_threshold,
            "b_to_a",
            "ambiguous"
        )
    )

    edge_table$from <- ifelse(
        edge_table$direction == "b_to_a",
        edge_table$cluster_b,
        edge_table$cluster_a
    )
    edge_table$to <- ifelse(
        edge_table$direction == "b_to_a",
        edge_table$cluster_a,
        edge_table$cluster_b
    )
    edge_table$direction_label <- ifelse(
        edge_table$direction == "a_to_b",
        paste(edge_table$cluster_a, "->", edge_table$cluster_b),
        ifelse(
            edge_table$direction == "b_to_a",
            paste(edge_table$cluster_b, "->", edge_table$cluster_a),
            paste(edge_table$cluster_a, "--", edge_table$cluster_b)
        )
    )

    edge_table
}


#' Build the cluster vertex table with sizes, median times and centroids
#'
#' @keywords internal
#' @noRd
.tata_vertex_table <- function(clusters, embedding, timepoint) {
    clusters <- factor(clusters)
    cluster_levels <- levels(clusters)
    time_numeric <- .coerce_time_to_numeric(timepoint)
    cluster_size <- table(clusters)
    cluster_median_time <- tapply(time_numeric, clusters, stats::median)
    centroid_df <- .compute_cluster_centroids(embedding, clusters)

    vertex_df <- data.frame(
        cluster = cluster_levels,
        size = as.integer(cluster_size[cluster_levels]),
        median_time = as.numeric(cluster_median_time[cluster_levels]),
        stringsAsFactors = FALSE
    )

    vertex_df <- merge(
        vertex_df,
        centroid_df,
        by = "cluster",
        all.x = TRUE,
        sort = FALSE
    )

    vertex_df[
        match(cluster_levels, vertex_df$cluster), ,
        drop = FALSE
    ]
}


#' Expand retained cluster pairs into directed igraph edge rows
#'
#' Directed pairs contribute one row; ambiguous pairs contribute a reciprocal
#' pair of rows so the undirected relationship survives in a directed graph.
#'
#' @keywords internal
#' @noRd
.tata_expand_graph_edges <- function(kept_edges) {
    if (nrow(kept_edges) == 0L) {
        return(data.frame(
            from = character(0),
            to = character(0),
            weight = numeric(0),
            topology_weight = numeric(0),
            time_weight = numeric(0),
            flow_score = numeric(0),
            direction = character(0),
            ambiguous = logical(0),
            stringsAsFactors = FALSE
        ))
    }

    do.call(rbind, lapply(seq_len(nrow(kept_edges)), function(idx) {
        edge_row <- kept_edges[idx, , drop = FALSE]

        if (edge_row$direction == "ambiguous") {
            data.frame(
                from = c(edge_row$cluster_a, edge_row$cluster_b),
                to = c(edge_row$cluster_b, edge_row$cluster_a),
                weight = rep(edge_row$tata_weight, 2),
                topology_weight = rep(edge_row$topology_weight, 2),
                time_weight = rep(edge_row$time_weight, 2),
                time_consistency = rep(edge_row$time_consistency, 2),
                flow_score = rep(edge_row$flow_score, 2),
                signed_time_score = rep(edge_row$signed_time_score, 2),
                median_time_gap = rep(edge_row$median_time_gap, 2),
                time_overlap = rep(edge_row$time_overlap, 2),
                direction = rep("ambiguous", 2),
                ambiguous = rep(TRUE, 2),
                stringsAsFactors = FALSE
            )
        } else {
            data.frame(
                from = edge_row$from,
                to = edge_row$to,
                weight = edge_row$tata_weight,
                topology_weight = edge_row$topology_weight,
                time_weight = edge_row$time_weight,
                time_consistency = edge_row$time_consistency,
                flow_score = edge_row$flow_score,
                signed_time_score = edge_row$signed_time_score,
                median_time_gap = edge_row$median_time_gap,
                time_overlap = edge_row$time_overlap,
                direction = "directed",
                ambiguous = FALSE,
                stringsAsFactors = FALSE
            )
        }
    }))
}
