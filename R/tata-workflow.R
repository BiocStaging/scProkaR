#' Run the full TATA workflow
#'
#' Run Time Aware Trajectory Analysis (TATA) end-to-end on a
#' `SingleCellExperiment`.
#'
#' @param sce A `SingleCellExperiment`.
#' @param dimred Reduced dimension to use for graph construction.
#' @param dims Optional subset of dimensions.
#' @param time_col Column containing experimental time labels.
#' @param use_time Logical indicating whether known experimental time should be
#'   used during graph abstraction. Set `FALSE` to run a topology-only TATA
#'   variant.
#' @param cluster_col Name of the column where inferred cluster labels will be
#'   stored.
#' @param k Number of neighbours for the cell-cell kNN graph.
#' @param cluster_method Graph clustering method. Choices are:
#' - `"louvain"`: fast modularity-based clustering and the default option.
#' - `"leiden"`: Leiden community detection if available in the installed
#'   `igraph`; otherwise TATA falls back to Louvain.
#' - `"walktrap"`: random-walk clustering that can give slightly coarser states.
#' @param alpha Exponent for topology weights.
#' @param beta Exponent for time weights.
#' @param direction_threshold Threshold for direction assignment.
#' @param prune_threshold Minimum TATA weight to retain an abstract graph edge.
#' @param time_weight_mode How the time signal contributes to the final edge
#'   score. Choices are:
#' - `"directional_confidence"`: use direction strength regardless of whether
#'   the label ordering is `a -> b` or `b -> a`; this is the recommended option
#'   for multi-branch graphs.
#' - `"ordered"`: use the original ordered-pair score `(flow_score + 1) / 2`,
#'   which is stricter but more sensitive to the arbitrary ordering of cluster
#'   labels.
#' @param do_pseudotime Logical indicating whether to compute pseudotime.
#' @param root_cluster Optional root cluster for pseudotime.
#' @param refine_clusters Logical indicating whether to apply a conservative
#'   temporal refinement step after graph clustering.
#' @param refine_min_cells Minimum cluster size considered for temporal
#'   refinement.
#' @param refine_min_split_fraction Minimum fraction of cells required in both
#'   temporal subclusters for a split to be accepted.
#' @param do_branch_probs Logical indicating whether to compute terminal-branch
#'   probabilities and entropy/plasticity scores.
#' @param expected_branches Optional expected number of terminal branches. If
#'   supplied, TATA uses it to guide terminal-state selection for downstream
#'   branch probabilities and branch assignment.
#' @param n_pcs_if_missing Number of PCs to compute if `dimred` is missing.
#'
#' @return A list containing the updated `sce`, cell graph, adjacency matrix,
#'   cluster graph, cluster edge table, cluster pseudotime, cell pseudotime, and
#'   parameters used.
#' @export
#'
#' @examples
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' tata <- run_tata(sce, dimred = "PCA", k = 12)
#' names(tata)
run_tata <- function(
    sce,
    dimred = "PCA",
    dims = NULL,
    time_col = "timepoint",
    use_time = TRUE,
    cluster_col = "cluster",
    k = 15,
    cluster_method = c("louvain", "leiden", "walktrap"),
    alpha = 1,
    beta = 1,
    direction_threshold = 0.1,
    prune_threshold = 0.01,
    time_weight_mode = c("directional_confidence", "ordered"),
    do_pseudotime = TRUE,
    root_cluster = NULL,
    refine_clusters = TRUE,
    refine_min_cells = NULL,
    refine_min_split_fraction = 0.20,
    do_branch_probs = TRUE,
    expected_branches = NULL,
    n_pcs_if_missing = 20
) {
    if (!is(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment.", call. = FALSE)
    }

    cluster_method <- match.arg(cluster_method)
    time_weight_mode <- match.arg(time_weight_mode)

    .validate_tata_workflow_args(sce, time_col, use_time, alpha, beta)
    # Read settings lazily, preserving argument evaluation and RNG stage order.
    parameters <- environment()
    .run_tata_stages(sce, parameters, root_cluster, n_pcs_if_missing)
}


# Prepare the cell states and the weighted abstract graph together, so that
# downstream stages consume the same embedding, labels, and edge tables.
.tata_graph_stage <- function(sce, parameters, n_pcs_if_missing) {
    sce <- .ensure_reduced_dim(
        sce, dimred = parameters$dimred, n_pcs = n_pcs_if_missing
    )
    knn_result <- build_knn_graph(
        x = sce, dimred = parameters$dimred, dims = parameters$dims,
        k = parameters$k
    )
    cluster_stage <- .tata_cluster_stage(
        sce = sce, knn_result = knn_result,
        cluster_method = parameters$cluster_method,
        use_time = parameters$use_time, time_col = parameters$time_col,
        refine_clusters = parameters$refine_clusters,
        refine_min_cells = parameters$refine_min_cells,
        refine_min_split_fraction = parameters$refine_min_split_fraction
    )
    clusters <- cluster_stage$clusters
    colData(sce)[[parameters$cluster_col]] <- clusters
    weight_tables <- .tata_weight_tables(
        sce = sce, knn_result = knn_result, clusters = clusters,
        use_time = parameters$use_time, time_col = parameters$time_col
    )
    graph_result <- build_tata_graph(
        topology_table = weight_tables$topology_table,
        time_table = weight_tables$time_table, clusters = clusters,
        embedding = knn_result$embedding,
        timepoint = .tata_timepoint_vector(
            sce, parameters$use_time, parameters$time_col
        ),
        alpha = parameters$alpha, beta = parameters$beta,
        direction_threshold = parameters$direction_threshold,
        prune_threshold = parameters$prune_threshold,
        time_weight_mode = parameters$time_weight_mode
    )
    list(
        sce = sce, knn = knn_result, graph = graph_result,
        weights = weight_tables,
        refine_min_cells = cluster_stage$refine_min_cells
    )
}


# Run optional inference stages and store their shared outputs in one place.
.run_tata_stages <- function(sce, parameters, root_cluster, n_pcs_if_missing) {
    graph_stage <- .tata_graph_stage(sce, parameters, n_pcs_if_missing)
    pseudotime_stage <- .tata_pseudotime_stage(
        sce = graph_stage$sce, cluster_graph = graph_stage$graph$graph,
        cluster_col = parameters$cluster_col,
        use_time = parameters$use_time, time_col = parameters$time_col,
        root_cluster = root_cluster, do_pseudotime = parameters$do_pseudotime
    )
    pseudotime <- pseudotime_stage$pseudotime
    branch_stage <- .tata_branch_stage(
        sce = pseudotime_stage$sce, cluster_graph = graph_stage$graph$graph,
        pseudotime_result = pseudotime, cluster_col = parameters$cluster_col,
        expected_branches = parameters$expected_branches,
        do_branch_probs = parameters$do_branch_probs
    )
    branch <- branch_stage$branch
    cell_space <- .compute_tata_cell_space(
        embedding = graph_stage$knn$embedding,
        pseudotime = pseudotime$cell_pseudotime_scaled,
        branch_probabilities = branch$branch_probabilities,
        branch_probability_columns = branch$branch_probability_columns
    )
    sce <- branch_stage$sce
    reducedDim(sce, "TATA") <- cell_space
    parameters <- .tata_recorded_parameters(parameters, graph_stage, pseudotime)
    .assemble_tata_result(
        sce = sce, knn_result = graph_stage$knn,
        tata_graph_result = graph_stage$graph, pseudotime_result = pseudotime,
        branch_result = branch, cell_space = cell_space,
        topology_table = graph_stage$weights$topology_table,
        time_table = graph_stage$weights$time_table, parameters = parameters
    )
}


# Convert lazy settings into the public metadata only after all stages finish.
.tata_recorded_parameters <- function(settings, graph_stage, pseudotime) {
    list(
        dimred = settings$dimred,
        dims = settings$dims,
        time_col = settings$time_col,
        use_time = settings$use_time,
        cluster_col = settings$cluster_col,
        k = graph_stage$knn$k,
        knn_engine = graph_stage$knn$knn_engine,
        cluster_method = settings$cluster_method,
        alpha = settings$alpha,
        beta = settings$beta,
        direction_threshold = settings$direction_threshold,
        prune_threshold = settings$prune_threshold,
        time_weight_mode = settings$time_weight_mode,
        do_pseudotime = settings$do_pseudotime,
        refine_clusters = settings$refine_clusters,
        refine_min_cells = graph_stage$refine_min_cells,
        refine_min_split_fraction = settings$refine_min_split_fraction,
        do_branch_probs = settings$do_branch_probs,
        expected_branches = settings$expected_branches,
        root_cluster = pseudotime$root_cluster
    )
}


#' Validate the time and edge-weight arguments of the TATA workflow
#'
#' @keywords internal
#' @noRd
.validate_tata_workflow_args <- function(
    sce,
    time_col,
    use_time,
    alpha,
    beta
) {
    if (isTRUE(use_time) && !time_col %in% colnames(colData(sce))) {
        stop(
            "Time column `", time_col,
            "` is not present in colData(sce).",
            call. = FALSE
        )
    }

    if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0) {
        stop(
            "`alpha` must be a single non-negative numeric value.",
            call. = FALSE
        )
    }

    if (!is.numeric(beta) || length(beta) != 1L || beta < 0) {
        stop(
            "`beta` must be a single non-negative numeric value.",
            call. = FALSE
        )
    }

    invisible(TRUE)
}


#' Experimental time vector used for graph abstraction
#'
#' Returns the stored time labels when known time is used, and an all-missing
#' vector of the right length otherwise.
#'
#' @keywords internal
#' @noRd
.tata_timepoint_vector <- function(sce, use_time, time_col) {
    if (isTRUE(use_time)) {
        colData(sce)[[time_col]]
    } else {
        rep(NA_real_, ncol(sce))
    }
}


#' Infer cell states and optionally refine them with experimental time
#'
#' @return A list with the cluster labels and the minimum cluster size that
#'   was actually used for temporal refinement.
#'
#' @keywords internal
#' @noRd
.tata_cluster_stage <- function(
    sce,
    knn_result,
    cluster_method,
    use_time,
    time_col,
    refine_clusters,
    refine_min_cells,
    refine_min_split_fraction
) {
    inferred_cluster <- cluster_graph_states(
        cell_graph = knn_result$graph,
        method = cluster_method
    )

    if (isTRUE(refine_clusters) && isTRUE(use_time)) {
        if (is.null(refine_min_cells)) {
            refine_min_cells <- max(40L, 2L * knn_result$k)
        }

        inferred_cluster <- .refine_clusters_temporally(
            clusters = inferred_cluster,
            timepoint = colData(sce)[[time_col]],
            embedding = knn_result$embedding,
            min_cluster_size = refine_min_cells,
            min_split_fraction = refine_min_split_fraction
        )
    }

    list(
        clusters = inferred_cluster,
        refine_min_cells = refine_min_cells
    )
}


#' Compute the topology and time edge-weight tables
#'
#' @return A list with the topology table and the time table.
#'
#' @keywords internal
#' @noRd
.tata_weight_tables <- function(
    sce,
    knn_result,
    clusters,
    use_time,
    time_col
) {
    topology_table <- compute_topology_weights(
        adjacency = knn_result$adjacency,
        clusters = clusters
    )

    time_table <- if (isTRUE(use_time)) {
        compute_time_weights(
            adjacency = knn_result$adjacency,
            clusters = clusters,
            timepoint = colData(sce)[[time_col]]
        )
    } else {
        .make_neutral_time_table(clusters)
    }

    list(
        topology_table = topology_table,
        time_table = time_table
    )
}


#' Compute TATA pseudotime and store it on the object
#'
#' @return A list with the updated `sce` and the pseudotime result.
#'
#' @keywords internal
#' @noRd
.tata_pseudotime_stage <- function(
    sce,
    cluster_graph,
    cluster_col,
    use_time,
    time_col,
    root_cluster,
    do_pseudotime
) {
    if (isTRUE(do_pseudotime)) {
        pseudotime_result <- compute_tata_pseudotime(
            cluster_graph = cluster_graph,
            sce = sce,
            cluster_col = cluster_col,
            time_col = if (isTRUE(use_time)) time_col else NULL,
            root_cluster = root_cluster
        )
colData(sce)$tata_pseudotime <- pseudotime_result$cell_pseudotime
        colData(sce)$tata_pseudotime_scaled <-
            pseudotime_result$cell_pseudotime_scaled
    } else {
        pseudotime_result <- list(
            root_cluster = NULL,
            cluster_pseudotime = NULL,
            cell_pseudotime = NULL,
            cell_pseudotime_scaled = NULL,
            within_cluster_step = NULL
        )
    }

    list(sce = sce, pseudotime = pseudotime_result)
}


#' Compute terminal-branch probabilities and store them on the object
#'
#' @return A list with the updated `sce` and the branch result.
#'
#' @keywords internal
#' @noRd
.tata_branch_stage <- function(
    sce,
    cluster_graph,
    pseudotime_result,
    cluster_col,
    expected_branches,
    do_branch_probs
) {
    if (!isTRUE(do_branch_probs)) {
        branch_result <- list(
            terminal_clusters = NULL,
            transition_matrix = NULL,
            cluster_branch_probabilities = NULL,
            branch_probabilities = NULL,
            branch_probability_columns = character(0)
        )
        return(list(sce = sce, branch = branch_result))
    }

    branch_result <- compute_tata_branch_probabilities(
        tata_result = list(
            sce = sce,
            cluster_graph = cluster_graph,
            cluster_pseudotime = pseudotime_result$cluster_pseudotime,
            parameters = list(
                cluster_col = cluster_col,
                root_cluster = pseudotime_result$root_cluster
            )
        ),
        expected_branches = expected_branches
    )

    sce <- .store_tata_branch_coldata(sce = sce, branch_result = branch_result)

    list(sce = sce, branch = branch_result)
}


#' Copy branch probabilities and branch summaries into `colData`
#'
#' @keywords internal
#' @noRd
.store_tata_branch_coldata <- function(sce, branch_result) {
    prob_cols <- branch_result$branch_probability_columns
    prob_values <-
        branch_result$branch_probabilities[, prob_cols, drop = FALSE]
    rownames(prob_values) <- branch_result$branch_probabilities$cell_id
    prob_values <- prob_values[colnames(sce), , drop = FALSE]
    colData(sce) <- Reduce(function(cell_data, col_name) {
        cell_data[[col_name]] <- prob_values[[col_name]]
        cell_data
    }, prob_cols, colData(sce))
    colData(sce)$tata_branch_entropy <-
        branch_result$branch_probabilities$tata_branch_entropy
    colData(sce)$tata_branch_plasticity <-
        branch_result$branch_probabilities$tata_branch_plasticity
    colData(sce)$tata_branch_commitment <-
        branch_result$branch_probabilities$tata_branch_commitment
    colData(sce)$tata_max_branch_probability <-
        branch_result$branch_probabilities$tata_max_branch_probability
    colData(sce)$tata_branch_assignment <-
        branch_result$branch_probabilities$tata_branch_assignment

    sce
}


#' Store the TATA metadata and assemble the workflow return value
#'
#' @keywords internal
#' @noRd
.assemble_tata_result <- function(
    sce,
    knn_result,
    tata_graph_result,
    pseudotime_result,
    branch_result,
    cell_space,
    topology_table,
    time_table,
    parameters
) {
    metadata(sce)$TATA <- list(
        parameters = parameters,
        cluster_edge_table = tata_graph_result$edge_table,
        cluster_vertex_table = tata_graph_result$vertex_table,
        cluster_pseudotime = pseudotime_result$cluster_pseudotime,
        terminal_clusters = branch_result$terminal_clusters,
        cluster_branch_probabilities =
            branch_result$cluster_branch_probabilities
    )

    list(
        sce = sce,
        cell_graph = knn_result$graph,
        adjacency_matrix = knn_result$adjacency,
        cluster_graph = tata_graph_result$graph,
        cluster_edge_table = tata_graph_result$edge_table,
        cluster_vertex_table = tata_graph_result$vertex_table,
        cluster_pseudotime = pseudotime_result$cluster_pseudotime,
        cell_pseudotime = pseudotime_result$cell_pseudotime,
        cell_pseudotime_scaled = pseudotime_result$cell_pseudotime_scaled,
        cell_space = cell_space,
        terminal_clusters = branch_result$terminal_clusters,
        cluster_branch_probabilities =
            branch_result$cluster_branch_probabilities,
        branch_probabilities = branch_result$branch_probabilities,
        transition_matrix = branch_result$transition_matrix,
        parameters = parameters,
        topology_table = topology_table,
        time_table = time_table
    )
}
