# Internal helpers for TATA branch probabilities, uncertainty, and gene trends.

.make_tata_probability_colnames <- function(branches) {
    paste0("tata_prob_", make.names(branches))
}


.compute_normalized_entropy <- function(prob_matrix) {
    if (is.null(dim(prob_matrix))) {
        prob_matrix <- matrix(prob_matrix, nrow = 1L)
    }

    n_branches <- ncol(prob_matrix)
    if (n_branches <= 1L) {
        return(rep(0, nrow(prob_matrix)))
    }

    prob_matrix <- prob_matrix / pmax(rowSums(prob_matrix), 1e-12)
    log_prob <- prob_matrix
    keep <- prob_matrix > 0
    log_prob[keep] <- log(prob_matrix[keep])
    log_prob[!keep] <- 0
    -rowSums(prob_matrix * log_prob) / log(n_branches)
}


.identify_tata_terminal_clusters <- function(
    cluster_graph,
    cluster_pseudotime = NULL,
    root_cluster = NULL,
    expected_branches = NULL
) {
    cluster_names <- igraph::V(cluster_graph)$name
    if (length(cluster_names) == 0L) {
        return(character(0))
    }

    pt_root <- .tata_resolve_pseudotime_root(
        cluster_graph = cluster_graph,
        cluster_names = cluster_names,
        cluster_pseudotime = cluster_pseudotime,
        root_cluster = root_cluster
    )
    pt_lookup <- pt_root$pt_lookup
    root_cluster <- pt_root$root_cluster

    terminal_clusters <- .tata_initial_terminal_clusters(
        cluster_graph = cluster_graph,
        cluster_names = cluster_names,
        pt_lookup = pt_lookup,
        root_cluster = root_cluster
    )

    if (!is.null(expected_branches)) {
        return(.tata_expected_branch_terminals(
            cluster_graph = cluster_graph,
            cluster_names = cluster_names,
            pt_lookup = pt_lookup,
            root_cluster = root_cluster,
            terminal_clusters = terminal_clusters,
            expected_branches = expected_branches
        ))
    }

    if (length(terminal_clusters) > 1L &&
            !is.null(root_cluster) &&
            root_cluster %in% cluster_names) {
        terminal_clusters <- .tata_collapse_terminal_clusters(
            cluster_graph = cluster_graph,
            terminal_clusters = terminal_clusters,
            root_cluster = root_cluster,
            pt_lookup = pt_lookup
        )
    }

    unique(as.character(terminal_clusters))
}


.build_tata_transition_matrix <- function(
    cluster_graph,
    terminal_clusters,
    cluster_pseudotime = NULL,
    backward_penalty = 0.15
) {
    cluster_names <- as.vector(igraph::V(cluster_graph)$name)
    n_clusters <- length(cluster_names)
    transition <- matrix(
        0,
        nrow = n_clusters,
        ncol = n_clusters,
        dimnames = list(cluster_names, cluster_names)
    )

    edge_df <- igraph::as_data_frame(cluster_graph, what = "edges")
    if (nrow(edge_df) == 0L) {
        diag(transition) <- 1
        return(transition)
    }

    pt_lookup <- .tata_cluster_pt_lookup(cluster_pseudotime)

    transition <- Reduce(function(transition, cluster_name) {
        outgoing <- edge_df[edge_df$from == cluster_name, , drop = FALSE]
        weights <- outgoing$weight
        if (!is.null(pt_lookup) &&
                all(c(cluster_name, outgoing$to) %in% names(pt_lookup))) {
            delta <- pt_lookup[outgoing$to] - pt_lookup[cluster_name]
            weights <- weights *
                ifelse(is.na(delta) | delta >= -1e-8, 1, backward_penalty)
        }
        if (cluster_name %in% terminal_clusters ||
                sum(weights, na.rm = TRUE) <= 0) {
            transition[cluster_name, cluster_name] <- 1
            return(transition)
        }
        prob <- tapply(weights / sum(weights), outgoing$to, sum)
        transition[cluster_name, names(prob)] <- as.numeric(prob)
        residual <- 1 - sum(transition[cluster_name, ])
        if (residual > 1e-8) {
            transition[cluster_name, cluster_name] <- residual
        }
        transition
    }, cluster_names, transition)

    transition / pmax(rowSums(transition), 1e-12)
}


.tata_cluster_pt_lookup <- function(cluster_pseudotime, fallback = NULL) {
    if (!is.null(cluster_pseudotime) && nrow(cluster_pseudotime) > 0L) {
        return(stats::setNames(
            cluster_pseudotime$pseudotime,
            cluster_pseudotime$cluster
        ))
    }
    fallback
}


.solve_absorbing_probabilities <- function(transition, terminal_clusters) {
    cluster_names <- rownames(transition)
    transient_clusters <- setdiff(cluster_names, terminal_clusters)
    n_terminals <- length(terminal_clusters)

    branch_prob <- matrix(
        0,
        nrow = length(cluster_names),
        ncol = n_terminals,
        dimnames = list(cluster_names, terminal_clusters)
    )

    if (n_terminals == 0L) {
        return(branch_prob)
    }

    branch_prob <- Reduce(function(bp, i) {
        bp[terminal_clusters[i], terminal_clusters[i]] <- 1
        bp
    }, seq_along(terminal_clusters), branch_prob)

    if (length(transient_clusters) == 0L) {
        return(branch_prob)
    }

    q_mat <- transition[transient_clusters, transient_clusters, drop = FALSE]
    r_mat <- transition[transient_clusters, terminal_clusters, drop = FALSE]
    identity_mat <- diag(nrow(q_mat))

    transient_prob <- tryCatch(
        {
            solve(identity_mat - q_mat, r_mat)
        },
        error = function(e) {
            power_mat <- Reduce(`%*%`, rep(list(transition), 201L))
            power_mat[transient_clusters, terminal_clusters, drop = FALSE]
        }
    )

    branch_prob[transient_clusters, terminal_clusters] <- transient_prob
    branch_prob / pmax(rowSums(branch_prob), 1e-12)
}


#' Compute TATA branch probabilities and entropy
#'
#' Estimate terminal-branch absorption probabilities for clusters and cells on
#' the TATA abstract graph, then derive an entropy-style plasticity score for
#' each cell.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param terminal_clusters Optional character vector of terminal clusters. If
#'   `NULL`, TATA identifies terminal clusters from the directed abstract graph
#'   and cluster pseudotime.
#' @param assignment_threshold Minimum maximum branch probability required for a
#'   hard branch assignment. Cells below this threshold are labelled
#'   `"ambiguous"`.
#' @param expected_branches Optional expected number of terminal branches. If
#'   supplied, TATA keeps up to that many high-pseudotime terminal candidates
#'   rather than relying only on automatically collapsed terminal states.
#' @param backward_penalty Multiplicative penalty applied to transitions that
#'   move backward in cluster pseudotime when constructing the transition
#'   matrix.
#'
#' @return A list with terminal clusters, cluster-level branch probabilities,
#'   cell-level branch probabilities, and entropy/plasticity summaries.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#' branch_res <- compute_tata_branch_probabilities(
#'     tata,
#'     expected_branches = 3
#' )
#' branch_res$terminal_clusters
#' head(branch_res$cluster_branch_probabilities)
#' table(branch_res$branch_probabilities$tata_branch_assignment)
compute_tata_branch_probabilities <- function(
    tata_result,
    terminal_clusters = NULL,
    assignment_threshold = 0.5,
    expected_branches = NULL,
    backward_penalty = 0.15
) {
    if (!is.list(tata_result) ||
            is.null(tata_result$cluster_graph) ||
            is.null(tata_result$sce)) {
        stop(
            "`tata_result` must be a list returned by `run_tata()`.",
            call. = FALSE
        )
    }

    cluster_graph <- tata_result$cluster_graph
    sce <- tata_result$sce
    cluster_pt <- tata_result$cluster_pseudotime
    cluster_col <- tata_result$parameters$cluster_col

    if (is.null(terminal_clusters)) {
        terminal_clusters <- .identify_tata_terminal_clusters(
            cluster_graph = cluster_graph,
            cluster_pseudotime = cluster_pt,
            root_cluster = tata_result$parameters$root_cluster,
            expected_branches = expected_branches
        )
    }

    transition <- .build_tata_transition_matrix(
        cluster_graph = cluster_graph,
        terminal_clusters = terminal_clusters,
        cluster_pseudotime = cluster_pt,
        backward_penalty = backward_penalty
    )
    cluster_prob_mat <- .solve_absorbing_probabilities(
        transition,
        terminal_clusters
    )

    prob_colnames <- .make_tata_probability_colnames(terminal_clusters)
    cluster_prob_df <- .tata_cluster_branch_table(
        cluster_prob_mat = cluster_prob_mat,
        terminal_clusters = terminal_clusters,
        prob_colnames = prob_colnames,
        assignment_threshold = assignment_threshold
    )
    cell_prob_df <- .tata_cell_branch_table(
        sce = sce,
        cluster_col = cluster_col,
        cluster_prob_df = cluster_prob_df,
        prob_colnames = prob_colnames
    )

    list(
        terminal_clusters = terminal_clusters,
        transition_matrix = transition,
        cluster_branch_probabilities = cluster_prob_df,
        branch_probabilities = cell_prob_df,
        branch_probability_columns = prob_colnames
    )
}


#' Select cells from a TATA branch
#'
#' Select cells that belong to a chosen terminal branch using either hard branch
#' assignment or a probability threshold.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param branch Terminal branch or terminal cluster to select.
#' @param selection_mode How to define selected cells. Choices are:
#' - `"soft"`: select cells whose probability for the branch is at least
#'   `probability_threshold`.
#' - `"hard"`: select cells whose hard TATA branch assignment equals `branch`.
#' @param probability_threshold Minimum branch probability required for
#'   `selection_mode = "soft"`.
#'
#' @return A data frame with cell ids, branch probabilities, and a logical
#'   `selected` column.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#' branch <- tata$terminal_clusters[1]
#' soft <- select_tata_branch_cells(tata,
#'     branch = branch,
#'     selection_mode = "soft", probability_threshold = 0.3
#' )
#' head(soft)
#' table(soft$selected)
#' hard <- select_tata_branch_cells(tata,
#'     branch = branch,
#'     selection_mode = "hard"
#' )
#' table(hard$selected)
select_tata_branch_cells <- function(
    tata_result,
    branch,
    selection_mode = c("soft", "hard"),
    probability_threshold = 0.25
) {
    selection_mode <- match.arg(selection_mode)
    branch_prob_df <- tata_result$branch_probabilities
    if (is.null(branch_prob_df)) {
        stop(
            "`tata_result` does not contain branch probabilities.",
            call. = FALSE
        )
    }

    branch_col <- .make_tata_probability_colnames(branch)
    if (!branch_col %in% colnames(branch_prob_df)) {
        stop(
            "Branch `", branch,
            "` is not present in the branch probability table.",
            call. = FALSE
        )
    }

    selected <- if (selection_mode == "hard") {
        branch_prob_df$tata_branch_assignment == branch
    } else {
        branch_prob_df[[branch_col]] >= probability_threshold
    }

    data.frame(
        cell_id = branch_prob_df$cell_id,
        cluster = branch_prob_df$cluster,
        branch = branch,
        probability = branch_prob_df[[branch_col]],
        tata_branch_assignment = branch_prob_df$tata_branch_assignment,
        selected = selected,
        stringsAsFactors = FALSE
    )
}


#' Compute branch-specific TATA gene trends
#'
#' Fit smooth branch-weighted gene-expression trends along TATA pseudotime for
#' one or more genes.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param genes Character vector of genes to model.
#' @param branches Optional subset of terminal branches.
#' @param assay_name Assay to use. If `NULL`, TATA prefers `logcounts` and
#'   falls back to `counts`.
#' @param probability_threshold Minimum branch probability used to retain cells
#'   when fitting a branch-specific trend.
#' @param grid_length Number of pseudotime grid points for fitted curves.
#' @param span Smoothing span passed to `stats::loess()`.
#'
#' @return A list containing the fitted trend table, the selected branches, and
#'   the genes used.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#' trends <- compute_tata_gene_trends(
#'     tata,
#'     genes = head(rownames(tata$sce), 6),
#'     grid_length = 30
#' )
#' trends$branches
#' head(trends$trend_table)
#' with(trends$trend_table, tapply(fitted_expression, branch, range))
compute_tata_gene_trends <- function(
    tata_result,
    genes,
    branches = NULL,
    assay_name = NULL,
    probability_threshold = 0.05,
    grid_length = 100,
    span = 0.5
) {
    if (missing(genes) || length(genes) == 0L) {
        stop("Please provide at least one gene in `genes`.", call. = FALSE)
    }

    sce <- tata_result$sce
    assay_name <- .scprokar_pick_assay(sce, assay_name, stop(
        "No suitable assay was found for gene-trend fitting.", call. = FALSE
    ))
    genes <- .tata_validate_gene_trend_inputs(
        sce = sce,
        genes = genes,
        branch_prob_df = tata_result$branch_probabilities
    )
    branch_prob_df <- tata_result$branch_probabilities

    if (is.null(branches)) {
        branches <- tata_result$terminal_clusters
    }
    branch_cols <- .make_tata_probability_colnames(branches)

    expr <- assay(sce, assay_name)
    pseudotime <- colData(sce)$tata_pseudotime_scaled
    trend_grid <- seq(0, 1, length.out = grid_length)

    trend_rows <- .tata_fit_gene_trend_rows(
        expr = expr,
        pseudotime = pseudotime,
        branch_prob_df = branch_prob_df,
        genes = genes,
        branches = branches,
        branch_cols = branch_cols,
        trend_grid = trend_grid,
        probability_threshold = probability_threshold,
        span = span
    )

    trend_rows <- Filter(Negate(is.null), trend_rows)
    if (length(trend_rows) == 0L) {
        stop(
            "No branch-specific gene trends could be fit with the current ",
            "settings.",
            call. = FALSE
        )
    }

    list(
        trend_table = do.call(rbind, trend_rows),
        genes = genes,
        branches = branches,
        assay_name = assay_name,
        probability_threshold = probability_threshold,
        grid = trend_grid
    )
}


#' Cluster TATA gene trends by shape
#'
#' Cluster branch-specific fitted gene trends by shape using k-means.
#'
#' @param trend_result Result list returned by `compute_tata_gene_trends()`.
#' @param n_clusters Number of trend clusters.
#'
#' @return A data frame mapping genes to trend clusters.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' tata <- run_tata(sce,
#'     dimred = "PCA", time_col = "timepoint", k = 15,
#'     expected_branches = 3
#' )
#' trends <- compute_tata_gene_trends(
#'     tata,
#'     genes = head(rownames(tata$sce), 12),
#'     grid_length = 30
#' )
#' set.seed(2)
#' trend_clusters <- cluster_tata_gene_trends(trends, n_clusters = 3)
#' head(trend_clusters)
#' table(trend_clusters$trend_cluster)
cluster_tata_gene_trends <- function(
    trend_result,
    n_clusters = 4
) {
    trend_table <- trend_result$trend_table
    split_key <- paste(trend_table$gene, trend_table$branch, sep = "||")
    trend_mat <- stats::reshape(
        trend_table[, c("gene", "branch", "pseudotime", "fitted_expression")],
        idvar = c("gene", "branch"),
        timevar = "pseudotime",
        direction = "wide"
    )

    rownames(trend_mat) <- paste(trend_mat$gene, trend_mat$branch, sep = "||")
    numeric_mat <- as.matrix(trend_mat[
        , grepl("^fitted_expression\\.", colnames(trend_mat)),
        drop = FALSE
    ])
    numeric_mat <- t(scale(t(numeric_mat)))
    numeric_mat[is.na(numeric_mat)] <- 0

    km <- stats::kmeans(
        numeric_mat,
        centers = min(n_clusters, nrow(numeric_mat))
    )
    gene_cluster <- tapply(km$cluster, trend_mat$gene, function(x) x[1])

    data.frame(
        gene = names(gene_cluster),
        trend_cluster = as.integer(gene_cluster),
        stringsAsFactors = FALSE
    )
}


#' Resolve the cluster pseudotime lookup and root cluster for TATA
#'
#' @keywords internal
#' @noRd
.tata_resolve_pseudotime_root <- function(
    cluster_graph,
    cluster_names,
    cluster_pseudotime,
    root_cluster
) {
    pt_lookup <- .tata_cluster_pt_lookup(
        cluster_pseudotime,
        stats::setNames(igraph::V(cluster_graph)$median_time, cluster_names)
    )
    if (is.null(root_cluster) && NROW(cluster_pseudotime) > 0L &&
            any(cluster_pseudotime$is_root %in% TRUE)) {
        root_cluster <- cluster_pseudotime$cluster[
            which(cluster_pseudotime$is_root)[1]
        ]
    }

    if (is.null(root_cluster) && length(pt_lookup) > 0L) {
        root_cluster <- names(which.min(pt_lookup))[1]
    }

    list(pt_lookup = pt_lookup, root_cluster = root_cluster)
}


#' Identify candidate terminal clusters from the directed abstract graph
#'
#' @keywords internal
#' @noRd
.tata_initial_terminal_clusters <- function(
    cluster_graph,
    cluster_names,
    pt_lookup,
    root_cluster
) {
    edge_df <- igraph::as_data_frame(cluster_graph, what = "edges")
    if (!"ambiguous" %in% colnames(edge_df)) {
        edge_df$ambiguous <- FALSE
    }
    directed_edges <- edge_df[!edge_df$ambiguous, , drop = FALSE]

    terminal_clusters <- cluster_names[vapply(
        cluster_names,
        FUN.VALUE = logical(1),
        FUN = function(cluster_name) {
            if (!is.null(root_cluster) &&
                    identical(cluster_name, root_cluster)) {
                return(FALSE)
            }

            outgoing <- directed_edges$to[directed_edges$from == cluster_name]
            if (length(outgoing) == 0L) {
                return(TRUE)
            }

            if (all(is.na(pt_lookup[c(cluster_name, outgoing)]))) {
                return(FALSE)
            }

            later_outgoing <- outgoing[
                pt_lookup[outgoing] > (pt_lookup[cluster_name] + 1e-8)
            ]
            length(later_outgoing) == 0L
        }
    )]

    if (length(terminal_clusters) == 0L) {
        component_id <- igraph::components(
            igraph::as_undirected(cluster_graph, mode = "collapse")
        )$membership
        split_clusters <- split(cluster_names, component_id)
        terminal_clusters <- vapply(
            split_clusters,
            FUN.VALUE = character(1),
            FUN = function(cluster_set) {
                cluster_set[which.max(pt_lookup[cluster_set])]
            }
        )
    }

    unique(as.character(terminal_clusters))
}


#' Trace weighted shortest paths from the root cluster to target clusters
#'
#' @keywords internal
#' @noRd
.tata_paths_from_root <- function(cluster_graph, root_cluster, targets) {
    undirected_graph <- igraph::as_undirected(cluster_graph, mode = "collapse")
    lapply(targets, function(target_cluster) {
        path <- .tata_trajectory_path_vertices(
            cluster_graph, undirected_graph, root_cluster, target_cluster
        )
        igraph::V(cluster_graph)$name[as.integer(path)]
    })
}

.tata_shortest_path_vertices <- function(graph, from, to, mode) {
    ## Check topology before requesting a weighted route. An unreachable target
    ## is an expected empty path, while warnings about a reachable route remain
    ## visible to the caller.
    reachable <- igraph::distances(
        graph, v = from, to = to, mode = mode, weights = NA
    )[1]
    if (!is.finite(reachable)) {
        return(igraph::V(graph)[integer()])
    }
    igraph::shortest_paths(
        graph,
        from = from,
        to = to,
        mode = mode,
        weights = 1 / pmax(igraph::E(graph)$weight, 1e-8)
    )$vpath[[1]]
}



#' Length of the cluster prefix shared by every root-to-target path
#'
#' @keywords internal
#' @noRd
.tata_common_prefix_length <- function(path_list) {
    min_path_len <- min(vapply(path_list, length, integer(1)))
    shared <- vapply(seq_len(min_path_len), function(pos) {
        length(unique(vapply(path_list, `[`, character(1), pos))) == 1L
    }, logical(1))
    match(FALSE, c(shared, FALSE)) - 1L
}


#' Label each root-to-target path by the cluster where it leaves the trunk
#'
#' @keywords internal
#' @noRd
.tata_path_branch_keys <- function(path_list, common_prefix_len) {
    vapply(
        path_list,
        FUN.VALUE = character(1),
        FUN = function(path_clusters) {
            if (length(path_clusters) > common_prefix_len) {
                path_clusters[common_prefix_len + 1L]
            } else {
                utils::tail(path_clusters, 1)
            }
        }
    )
}


#' Rank one representative terminal cluster per branch of the abstract graph
#'
#' @keywords internal
#' @noRd
.tata_rank_branch_representatives <- function(
    path_candidates,
    candidate_clusters,
    pt_lookup,
    expected_branches
) {
    common_prefix_len <- .tata_common_prefix_length(path_candidates)
    branch_key <- .tata_path_branch_keys(path_candidates, common_prefix_len)

    grouped_candidates <- split(names(path_candidates), branch_key)
    branch_representatives <- vapply(
        grouped_candidates,
        FUN.VALUE = character(1),
        FUN = function(cluster_set) {
            cluster_set[which.max(pt_lookup[cluster_set])]
        }
    )
    branch_representatives <- unique(as.character(branch_representatives))
    ordered_representatives <- names(sort(
        pt_lookup[branch_representatives],
        decreasing = TRUE
    ))
    ordered_representatives <- ordered_representatives[
        !is.na(ordered_representatives)
    ]

    if (length(ordered_representatives) >= expected_branches) {
        return(utils::head(ordered_representatives, expected_branches))
    }

    remaining_candidates <- setdiff(
        candidate_clusters,
        ordered_representatives
    )
    ordered_remaining <- names(sort(
        pt_lookup[remaining_candidates],
        decreasing = TRUE
    ))
    ordered_remaining <- ordered_remaining[!is.na(ordered_remaining)]

    utils::head(
        unique(c(ordered_representatives, ordered_remaining)),
        expected_branches
    )
}


#' Select terminal clusters when an expected branch count is supplied
#'
#' @keywords internal
#' @noRd
.tata_expected_branch_terminals <- function(
    cluster_graph,
    cluster_names,
    pt_lookup,
    root_cluster,
    terminal_clusters,
    expected_branches
) {
    expected_branches <- max(1L, as.integer(expected_branches))

    candidate_clusters <- setdiff(cluster_names, root_cluster)
    candidate_clusters <- candidate_clusters[
        !is.na(pt_lookup[candidate_clusters])
    ]

    if (!is.null(root_cluster) &&
            root_cluster %in% cluster_names &&
            length(candidate_clusters) > 0L) {
        path_candidates <- .tata_paths_from_root(
            cluster_graph,
            root_cluster,
            candidate_clusters
        )
        names(path_candidates) <- candidate_clusters
        path_candidates <- Filter(
            function(path_clusters) length(path_clusters) > 1L,
            path_candidates
        )

        if (length(path_candidates) > 0L) {
            return(.tata_rank_branch_representatives(
                path_candidates = path_candidates,
                candidate_clusters = candidate_clusters,
                pt_lookup = pt_lookup,
                expected_branches = expected_branches
            ))
        }
    }

    ordered_candidates <- names(sort(
        pt_lookup[terminal_clusters],
        decreasing = TRUE
    ))
    ordered_candidates <- ordered_candidates[!is.na(ordered_candidates)]
    utils::head(
        unique(as.character(ordered_candidates)),
        expected_branches
    )
}


#' Collapse terminal clusters that share the same branch of the graph
#'
#' @keywords internal
#' @noRd
.tata_collapse_terminal_clusters <- function(
    cluster_graph,
    terminal_clusters,
    root_cluster,
    pt_lookup
) {
    path_list <- .tata_paths_from_root(
        cluster_graph,
        root_cluster,
        terminal_clusters
    )
    path_list <- Filter(function(x) length(x) > 0L, path_list)
    if (length(path_list) <= 1L) {
        return(terminal_clusters)
    }

    common_prefix_len <- .tata_common_prefix_length(path_list)
    branch_key <- .tata_path_branch_keys(path_list, common_prefix_len)

    grouped_paths <- split(seq_along(path_list), branch_key)
    keep_idx <- vapply(
        grouped_paths,
        FUN.VALUE = integer(1),
        FUN = function(idx) {
            candidate_terminal <- terminal_clusters[idx]
            idx[which.max(pt_lookup[candidate_terminal])]
        }
    )
    terminal_clusters[unname(keep_idx)]
}


#' Assemble the cluster-level TATA branch probability table
#'
#' @keywords internal
#' @noRd
.tata_cluster_branch_table <- function(
    cluster_prob_mat,
    terminal_clusters,
    prob_colnames,
    assignment_threshold
) {
    cluster_prob_df <- data.frame(
        cluster = rownames(cluster_prob_mat),
        cluster_prob_mat,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    colnames(cluster_prob_df)[-1] <- prob_colnames

    cluster_entropy <- .compute_normalized_entropy(cluster_prob_mat)
    cluster_prob_df$tata_branch_entropy <- cluster_entropy
    cluster_prob_df$tata_branch_plasticity <- cluster_entropy
    cluster_prob_df$tata_branch_commitment <- 1 - cluster_entropy
    cluster_prob_df$tata_max_branch_probability <-
        apply(cluster_prob_mat, 1, max)
    branch_assignment <- terminal_clusters[
        max.col(cluster_prob_mat, ties.method = "first")
    ]
    branch_assignment[
        cluster_prob_df$tata_max_branch_probability < assignment_threshold
    ] <- "ambiguous"
    cluster_prob_df$tata_branch_assignment <- branch_assignment
    cluster_prob_df
}


#' Expand the cluster-level branch table to one row per cell
#'
#' @keywords internal
#' @noRd
.tata_cell_branch_table <- function(
    sce,
    cluster_col,
    cluster_prob_df,
    prob_colnames
) {
    cell_clusters <- as.character(colData(sce)[[cluster_col]])
    cluster_lookup <- match(cell_clusters, cluster_prob_df$cluster)

    data.frame(
        cell_id = colnames(sce),
        cluster = cell_clusters,
        cluster_prob_df[cluster_lookup, prob_colnames, drop = FALSE],
        tata_branch_entropy =
            cluster_prob_df$tata_branch_entropy[cluster_lookup],
        tata_branch_plasticity =
            cluster_prob_df$tata_branch_plasticity[cluster_lookup],
        tata_branch_commitment =
            cluster_prob_df$tata_branch_commitment[cluster_lookup],
        tata_max_branch_probability =
            cluster_prob_df$tata_max_branch_probability[cluster_lookup],
        tata_branch_assignment =
            cluster_prob_df$tata_branch_assignment[cluster_lookup],
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}



#' Validate gene-trend inputs and return the usable gene set
#'
#' @keywords internal
#' @noRd
.tata_validate_gene_trend_inputs <- function(sce, genes, branch_prob_df) {
    genes <- intersect(unique(genes), rownames(sce))
    if (length(genes) == 0L) {
        stop("None of the requested genes are present in `sce`.", call. = FALSE)
    }
if (!"tata_pseudotime_scaled" %in% colnames(colData(sce))) {
        stop(
            "`tata_pseudotime_scaled` was not found in `colData(sce)`.",
            call. = FALSE
        )
    }

    if (is.null(branch_prob_df)) {
        stop(
            "`tata_result` does not contain branch probabilities.",
            call. = FALSE
        )
    }

    genes
}


#' Build the weighted fitting frame for one gene and one branch
#'
#' @keywords internal
#' @noRd
.tata_trend_fit_frame <- function(
    pseudotime,
    gene_expr,
    branch_prob,
    probability_threshold
) {
    keep <- !is.na(pseudotime) & !is.na(branch_prob) &
        branch_prob >= probability_threshold

    if (sum(keep) < 10L ||
            length(unique(round(pseudotime[keep], 3))) < 4L) {
        return(NULL)
    }

    fit_df <- data.frame(
        pseudotime = pseudotime[keep],
        expression = gene_expr[keep],
        weight = branch_prob[keep]
    )
    fit_df <- fit_df[
        is.finite(fit_df$pseudotime) &
            is.finite(fit_df$expression) &
            is.finite(fit_df$weight) &
            fit_df$weight > 0, ,
        drop = FALSE
    ]
    fit_df <- fit_df[order(fit_df$pseudotime), , drop = FALSE]

    if (nrow(fit_df) < 10L ||
            length(unique(round(fit_df$pseudotime, 3))) < 4L) {
        return(NULL)
    }

    list(fit_df = fit_df, n_cells = sum(keep))
}


#' Fit one weighted expression curve on the TATA pseudotime grid
#'
#' Uses a weighted loess fit, falling back to a smoothing spline and then to
#' the weighted mean when the loess fit cannot be evaluated.
#'
#' @keywords internal
#' @noRd
.tata_fit_trend_curve <- function(fit_df, trend_grid, span) {
    fallback_value <- rep(
        stats::weighted.mean(fit_df$expression, fit_df$weight),
        length(trend_grid)
    )

    fit_obj <- .scprokar_attempt(
        stats::loess(
            expression ~ pseudotime,
            data = fit_df,
            weights = weight,
            span = span
        ),
        fallback = NULL
    )

    fitted_value <- if (is.null(fit_obj)) {
        rep(NA_real_, length(trend_grid))
    } else {
        .scprokar_attempt(
            stats::predict(
                fit_obj,
                newdata = data.frame(pseudotime = trend_grid)
            ),
            fallback = rep(NA_real_, length(trend_grid))
        )
    }

    if (all(is.na(fitted_value))) {
        spline_fit <- .scprokar_attempt(
            stats::smooth.spline(
                x = fit_df$pseudotime,
                y = fit_df$expression,
                w = fit_df$weight,
                spar = 0.6
            ),
            fallback = NULL
        )

        if (is.null(spline_fit)) {
            fitted_value <- fallback_value
        } else {
            fitted_value <- .scprokar_attempt(
                stats::predict(spline_fit, x = trend_grid)$y,
                fallback = fallback_value
            )
        }
    } else if (anyNA(fitted_value)) {
        non_missing <- which(!is.na(fitted_value))
        fitted_value <- stats::approx(
            x = trend_grid[non_missing],
            y = fitted_value[non_missing],
            xout = trend_grid,
            rule = 2
        )$y
    }

    fitted_value
}


#' Fit branch-weighted expression trends for every branch and gene pair
#'
#' @keywords internal
#' @noRd
.tata_fit_gene_trend_rows <- function(
    expr,
    pseudotime,
    branch_prob_df,
    genes,
    branches,
    branch_cols,
    trend_grid,
    probability_threshold,
    span
) {
    Reduce(c, lapply(seq_along(branches), function(branch_idx) {
        branch_prob <- branch_prob_df[[branch_cols[branch_idx]]]

        lapply(as.vector(genes), function(gene) {
            fit_data <- .tata_trend_fit_frame(
                pseudotime = pseudotime,
                gene_expr = as.numeric(expr[gene, ]),
                branch_prob = branch_prob,
                probability_threshold = probability_threshold
            )
            if (is.null(fit_data)) {
                return(NULL)
            }

            fitted_value <- .tata_fit_trend_curve(
                fit_df = fit_data$fit_df,
                trend_grid = trend_grid,
                span = span
            )

            data.frame(
                gene = gene,
                branch = branches[branch_idx],
                pseudotime = trend_grid,
                fitted_expression = as.numeric(fitted_value),
                n_cells = fit_data$n_cells,
                stringsAsFactors = FALSE
            )
        })
    }), list())
}
