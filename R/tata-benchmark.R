.tata_correlation <- function(x, y, method) {
    ## Correlation is undefined for a constant vector or fewer than two pairs.
    if (length(x) < 2L || length(unique(x)) < 2L ||
            length(unique(y)) < 2L) {
        return(NA_real_)
    }
    stats::cor(x, y, method = method)
}

#' Orient pseudotime to match increasing real time
#'
#' Flip an inferred pseudotime vector when needed so that it best agrees with
#' increasing experimental time. This is useful when benchmarking methods whose
#' pseudotime direction is arbitrary up to reversal.
#'
#' @param pseudotime Numeric vector of inferred pseudotime values.
#' @param timepoint Experimental time labels for the same cells.
#' @param method Correlation used to choose the orientation. Choices are:
#' - `"spearman"`: rank-based and the default choice for monotone trajectories.
#' - `"kendall"`: more conservative rank correlation.
#' - `"pearson"`: linear correlation on the scaled values.
#'
#' @return A numeric vector scaled to `[0, 1]` and oriented so that larger
#'   values agree as well as possible with later experimental time.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce, dimred = "PCA", time_col = "timepoint", k = 15)
#'
#' oriented <- orient_pseudotime_to_time(
#'     pseudotime = tata$cell_pseudotime,
#'     timepoint = tata$sce$timepoint
#' )
#' summary(oriented)
#'
#' # Larger oriented values should now fall on later experimental times.
#' tapply(oriented, tata$sce$timepoint, median)
orient_pseudotime_to_time <- function(
    pseudotime,
    timepoint,
    method = c("spearman", "kendall", "pearson")
) {
    method <- match.arg(method)

    pseudotime <- as.numeric(pseudotime)
    time_numeric <- .coerce_time_to_numeric(timepoint)
    ok <- is.finite(pseudotime) & is.finite(time_numeric)

    oriented <- rep(NA_real_, length(pseudotime))
    if (sum(ok) < 3L) {
        return(oriented)
    }

    forward <- .scale_to_unit(pseudotime[ok])
    flipped_values <- max(pseudotime[ok], na.rm = TRUE) +
        min(pseudotime[ok], na.rm = TRUE) - pseudotime[ok]
    flipped <- .scale_to_unit(flipped_values)
    truth <- .scale_to_unit(time_numeric[ok])

    cor_forward <- .tata_correlation(forward, truth, method)
    cor_flipped <- .tata_correlation(flipped, truth, method)

    if (is.na(cor_forward)) {
        cor_forward <- -Inf
    }
    if (is.na(cor_flipped)) {
        cor_flipped <- -Inf
    }

    oriented[ok] <- if (cor_flipped > cor_forward) flipped else forward
    oriented
}


#' Compute local temporal order metrics on a cell graph
#'
#' Summarize how often pseudotime agrees with the known direction of time on
#' local cell-cell graph edges.
#'
#' @param pseudotime Numeric vector of inferred pseudotime values.
#' @param adjacency Sparse adjacency matrix for the cell-cell graph.
#' @param timepoint Experimental time labels for the same cells.
#'
#' @return A one-row `data.frame` with the number of compared edges, local time
#'   concordance, local temporal inversion rate, a higher-is-better temporal
#'   direction score, and the fraction of edges tied in pseudotime.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce, dimred = "PCA", time_col = "timepoint", k = 15)
#'
#' oriented <- orient_pseudotime_to_time(
#'     pseudotime = tata$cell_pseudotime,
#'     timepoint = tata$sce$timepoint
#' )
#' local_temporal_order_metrics(
#'     pseudotime = oriented,
#'     adjacency = tata$adjacency_matrix,
#'     timepoint = tata$sce$timepoint
#' )
local_temporal_order_metrics <- function(
    pseudotime,
    adjacency,
    timepoint
) {
    if (!inherits(adjacency, "sparseMatrix")) {
        adjacency <- as(adjacency, "dgCMatrix")
    }

    pseudotime <- as.numeric(pseudotime)
    time_numeric <- .coerce_time_to_numeric(timepoint)
    edge_df <- Matrix::summary(Matrix::triu(adjacency, k = 1))

    if (nrow(edge_df) == 0L) {
        return(data.frame(
            n_compared_edges = 0L,
            local_time_concordance = NA_real_,
            local_temporal_inversion_rate = NA_real_,
            pseudotime_tie_rate = NA_real_,
            stringsAsFactors = FALSE
        ))
    }

    valid <- is.finite(pseudotime[edge_df$i]) &
        is.finite(pseudotime[edge_df$j]) &
        is.finite(time_numeric[edge_df$i]) &
        is.finite(time_numeric[edge_df$j])
    edge_df <- edge_df[valid, , drop = FALSE]

    if (nrow(edge_df) == 0L) {
        return(data.frame(
            n_compared_edges = 0L,
            local_time_concordance = NA_real_,
            local_temporal_inversion_rate = NA_real_,
            pseudotime_tie_rate = NA_real_,
            stringsAsFactors = FALSE
        ))
    }

    true_sign <- sign(time_numeric[edge_df$j] - time_numeric[edge_df$i])
    inferred_sign <- sign(pseudotime[edge_df$j] - pseudotime[edge_df$i])
    keep <- true_sign != 0

    if (!any(keep)) {
        return(data.frame(
            n_compared_edges = 0L,
            local_time_concordance = NA_real_,
            local_temporal_inversion_rate = NA_real_,
            pseudotime_tie_rate = NA_real_,
            stringsAsFactors = FALSE
        ))
    }

    true_sign <- true_sign[keep]
    inferred_sign <- inferred_sign[keep]

    data.frame(
        n_compared_edges = length(true_sign),
        local_time_concordance = mean(inferred_sign == true_sign),
        local_temporal_inversion_rate = mean(inferred_sign == -true_sign),
        local_temporal_direction_score = 1 - mean(inferred_sign == -true_sign),
        pseudotime_tie_rate = mean(inferred_sign == 0),
        stringsAsFactors = FALSE
    )
}


#' Compute adjacent-timepoint AUC for a pseudotime ordering
#'
#' For each adjacent pair of experimental timepoints, compute the probability
#' that a cell from the later timepoint has larger pseudotime than a cell from
#' the earlier timepoint. This is equivalent to the Mann-Whitney AUC.
#'
#' @param pseudotime Numeric vector of inferred pseudotime values.
#' @param timepoint Experimental time labels for the same cells.
#' @param group Optional grouping variable, for example treatment or condition.
#'   If provided, adjacent-timepoint AUC is computed within each group before
#'   averaging.
#' @param average How to summarize the per-pair AUC values. Choices are:
#' - `"weighted"`: weight each comparison by the number of earlier-later cell
#'   pairs.
#' - `"macro"`: give every valid comparison equal weight.
#' @param min_cells Minimum number of cells required in each timepoint group for
#'   a comparison to be included.
#'
#' @return A list with:
#' - `summary`: a one-row `data.frame` containing the overall adjacent-timepoint
#'   AUC and the number of valid comparisons.
#' - `pair_table`: a `data.frame` with AUC values for each adjacent timepoint
#'   comparison.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce, dimred = "PCA", time_col = "timepoint", k = 15)
#'
#' oriented <- orient_pseudotime_to_time(
#'     pseudotime = tata$cell_pseudotime,
#'     timepoint = tata$sce$timepoint
#' )
#'
#' # Pooled across all cells.
#' auc_all <- adjacent_timepoint_auc(
#'     pseudotime = oriented,
#'     timepoint = tata$sce$timepoint
#' )
#' auc_all$summary
#' head(auc_all$pair_table)
#'
#' # Computed within each drug condition before averaging.
#' auc_by_drug <- adjacent_timepoint_auc(
#'     pseudotime = oriented,
#'     timepoint = tata$sce$timepoint,
#'     group = tata$sce$condition,
#'     min_cells = 5
#' )
#' auc_by_drug$summary
adjacent_timepoint_auc <- function(
    pseudotime,
    timepoint,
    group = NULL,
    average = c("weighted", "macro"),
    min_cells = 10L
) {
    average <- match.arg(average)
    pseudotime <- as.numeric(pseudotime)
    time_numeric <- .coerce_time_to_numeric(timepoint)

    if (is.null(group)) {
        group <- rep("all", length(pseudotime))
    }
    group <- as.character(group)

    time_levels <- sort(unique(time_numeric[is.finite(time_numeric)]))
    if (length(time_levels) < 2L) {
        return(.empty_auc_result(average))
    }

    pair_rows <- .auc_pair_rows(
        pseudotime = pseudotime,
        time_numeric = time_numeric,
        group = group,
        time_levels = time_levels,
        min_cells = min_cells
    )

    if (length(pair_rows) == 0L) {
        return(.empty_auc_result(average))
    }

    pair_table <- do.call(rbind, pair_rows)
    overall_auc <- if (average == "weighted") {
        stats::weighted.mean(pair_table$auc, w = pair_table$weight)
    } else {
        mean(pair_table$auc)
    }

    list(
        summary = data.frame(
            adjacent_timepoint_auc = overall_auc,
            n_valid_comparisons = nrow(pair_table),
            averaging = average,
            stringsAsFactors = FALSE
        ),
        pair_table = pair_table
    )
}


#' Measure branch-emergence timing from pseudotime
#'
#' Compare the known order or timing of branch emergence with the inferred
#' pseudotime of sufficiently committed cells in each branch.
#'
#' @param pseudotime Numeric vector of inferred pseudotime values.
#' @param branch Branch or condition label for each cell.
#' @param branch_activation Numeric activation score indicating how committed a
#'   cell is to its branch.
#' @param true_onset Named numeric vector giving the known onset time of each
#'   branch. The names must match a subset of the values in `branch`.
#' @param activation_threshold Minimum activation value required for a cell
#'   to be treated as branch-committed when estimating the branch onset.
#'
#' @return A list with:
#' - `summary`: a one-row `data.frame` containing Kendall tau, Spearman rho,
#'   normalized branch-onset mean absolute error, a higher-is-better branch
#'   onset order score, a higher-is-better branch onset accuracy score, and the
#'   number of branches used.
#' - `onset_table`: a `data.frame` with the true and inferred onset values for
#'   each branch.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce, dimred = "PCA", time_col = "timepoint", k = 15)
#'
#' oriented <- orient_pseudotime_to_time(
#'     pseudotime = tata$cell_pseudotime,
#'     timepoint = tata$sce$timepoint
#' )
#'
#' # The simulation activates drug_A first, then drug_C, then drug_B.
#' onset <- branch_onset_metrics(
#'     pseudotime = oriented,
#'     branch = tata$sce$simulated_branch,
#'     branch_activation = tata$sce$branch_activation,
#'     true_onset = c(drug_A = 24, drug_B = 96, drug_C = 48)
#' )
#' onset$summary
#' onset$onset_table
branch_onset_metrics <- function(
    pseudotime,
    branch,
    branch_activation,
    true_onset,
    activation_threshold = 0.25
) {
    pseudotime <- as.numeric(pseudotime)
    branch <- as.character(branch)
    branch_activation <- as.numeric(branch_activation)

    if (is.null(names(true_onset)) || anyNA(names(true_onset))) {
        stop("`true_onset` must be a named numeric vector.", call. = FALSE)
    }

    branch_names <- names(true_onset)
    keep <- branch %in% branch_names &
        is.finite(pseudotime) &
        is.finite(branch_activation) &
        branch_activation >= activation_threshold

    onset_table <- data.frame(
        branch = branch_names,
        true_onset = as.numeric(true_onset[branch_names]),
        inferred_onset = NA_real_,
        n_committed_cells = 0L,
        stringsAsFactors = FALSE
    )

    committed <- lapply(branch_names, function(b) {
        pseudotime[which(keep & branch == b)]
    })
    onset_table$inferred_onset <- vapply(
        committed, stats::median, numeric(1), na.rm = TRUE
    )
    onset_table$n_committed_cells <- lengths(committed)

    valid <- is.finite(onset_table$true_onset) &
        is.finite(onset_table$inferred_onset)
    if (sum(valid) < 2L) {
        summary_df <- data.frame(
            n_branches_used = sum(valid),
            kendall_tau = NA_real_,
            spearman_rho = NA_real_,
            normalized_onset_mae = NA_real_,
            branch_onset_order_score = NA_real_,
            branch_onset_accuracy = NA_real_,
            activation_threshold = activation_threshold,
            stringsAsFactors = FALSE
        )
        return(list(summary = summary_df, onset_table = onset_table))
    }

    inferred_scaled <- .scale_to_unit(onset_table$inferred_onset[valid])
    truth_scaled <- .scale_to_unit(onset_table$true_onset[valid])
    kendall_tau <- .tata_correlation(
        onset_table$true_onset[valid],
        onset_table$inferred_onset[valid],
        method = "kendall"
    )

    summary_df <- data.frame(
        n_branches_used = sum(valid),
        kendall_tau = kendall_tau,
        spearman_rho = .tata_correlation(
            onset_table$true_onset[valid],
            onset_table$inferred_onset[valid],
            method = "spearman"
        ),
        normalized_onset_mae = mean(abs(inferred_scaled - truth_scaled)),
        branch_onset_order_score = (kendall_tau + 1) / 2,
        branch_onset_accuracy = 1 - mean(abs(inferred_scaled - truth_scaled)),
        activation_threshold = activation_threshold,
        stringsAsFactors = FALSE
    )

    list(summary = summary_df, onset_table = onset_table)
}


#' Quantify how well a directed cluster graph follows real time
#'
#' For ordered cluster pairs with increasing median experimental time, measure
#' whether the directed cluster graph supports reachability in the forward
#' direction.
#'
#' @param cluster_graph Directed cluster-level graph as an `igraph`.
#' @param clusters Cluster label for each cell.
#' @param timepoint Experimental time labels for the same cells.
#' @param min_time_gap Minimum difference in median cluster time required before
#'   a pair is considered a meaningful forward-time comparison.
#' @param weight_by How to weight cluster-pair comparisons. Choices are:
#' - `"cells"`: weight each ordered pair by the product of the two cluster
#'   sizes.
#' - `"pairs"`: weight every ordered pair equally.
#'
#' @return A list with:
#' - `summary`: a one-row `data.frame` containing the weighted forward
#'   reachability concordance and the number of ordered cluster pairs.
#' - `pair_table`: a `data.frame` describing each ordered cluster pair.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce, dimred = "PCA", time_col = "timepoint", k = 15)
#'
#' reachability <- directed_reachability_concordance(
#'     cluster_graph = tata$cluster_graph,
#'     clusters = tata$sce$cluster,
#'     timepoint = tata$sce$timepoint
#' )
#' reachability$summary
#'
#' # Cluster pairs separated by at least one full sampling interval.
#' head(reachability$pair_table[reachability$pair_table$time_gap >= 24, ])
directed_reachability_concordance <- function(
    cluster_graph,
    clusters,
    timepoint,
    min_time_gap = 0,
    weight_by = c("cells", "pairs")
) {
    if (!inherits(cluster_graph, "igraph")) {
        stop("`cluster_graph` must be an igraph object.", call. = FALSE)
    }

    weight_by <- match.arg(weight_by)
    clusters <- factor(clusters, levels = igraph::V(cluster_graph)$name)
    time_numeric <- .coerce_time_to_numeric(timepoint)

    cluster_median_time <- tapply(time_numeric, clusters, stats::median)
    cluster_size <- tapply(seq_along(clusters), clusters, length)
    cluster_names <- names(cluster_median_time)

    pair_rows <- .reachability_pair_rows(
        cluster_graph = cluster_graph,
        cluster_names = cluster_names,
        cluster_median_time = cluster_median_time,
        cluster_size = cluster_size,
        min_time_gap = min_time_gap,
        weight_by = weight_by
    )

    if (length(pair_rows) == 0L) {
        return(.empty_reachability_result(weight_by))
    }

    pair_table <- do.call(rbind, pair_rows)

    list(
        summary = data.frame(
            n_ordered_cluster_pairs = nrow(pair_table),
            directed_reachability_concordance = stats::weighted.mean(
                pair_table$reachability_score,
                w = pair_table$weight
            ),
            weighting = weight_by,
            stringsAsFactors = FALSE
        ),
        pair_table = pair_table
    )
}


#' Measure the weight of anti-causal edges in a directed graph
#'
#' Compare each directed graph edge with the median experimental time of its
#' source and target clusters, and quantify how much total edge weight points
#' backward in time.
#'
#' @param cluster_graph Directed cluster-level graph as an `igraph`.
#' @param clusters Cluster label for each cell.
#' @param timepoint Experimental time labels for the same cells.
#' @param tolerance A non-negative tolerance for small median-time differences.
#'   Directed edges with time differences in `[-tolerance, tolerance]` are
#'   treated as time-neutral rather than anti-causal.
#' @param weight_attr Edge attribute to use as the edge weight.
#'
#' @return A list with:
#' - `summary`: a one-row `data.frame` containing the total edge weight, the
#'   anti-causal edge mass, its proportion, and a higher-is-better causal edge
#'   score.
#' - `edge_table`: a `data.frame` describing each directed edge.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60,
#'     n_pcs = 10
#' )
#' tata <- run_tata(sce, dimred = "PCA", time_col = "timepoint", k = 15)
#'
#' causality <- anticausal_edge_mass(
#'     cluster_graph = tata$cluster_graph,
#'     clusters = tata$sce$cluster,
#'     timepoint = tata$sce$timepoint
#' )
#' causality$summary
#' table(causality$edge_table$temporal_class)
#'
#' # Treat median-time gaps below one sampling interval as time-neutral.
#' anticausal_edge_mass(
#'     cluster_graph = tata$cluster_graph,
#'     clusters = tata$sce$cluster,
#'     timepoint = tata$sce$timepoint,
#'     tolerance = 24
#' )$summary
anticausal_edge_mass <- function(
    cluster_graph,
    clusters,
    timepoint,
    tolerance = 0,
    weight_attr = "weight"
) {
    if (!inherits(cluster_graph, "igraph")) {
        stop("`cluster_graph` must be an igraph object.", call. = FALSE)
    }

    clusters <- factor(clusters, levels = igraph::V(cluster_graph)$name)
    time_numeric <- .coerce_time_to_numeric(timepoint)
    cluster_median_time <- tapply(time_numeric, clusters, stats::median)

    edge_df <- igraph::as_data_frame(cluster_graph, what = "edges")
    if (!nrow(edge_df)) {
        empty_edges <- data.frame(
            from = character(0),
            to = character(0),
            weight = numeric(0),
            median_time_from = numeric(0),
            median_time_to = numeric(0),
            time_delta = numeric(0),
            temporal_class = character(0),
            stringsAsFactors = FALSE
        )
        return(list(
            summary = data.frame(
                total_edge_weight = 0,
                anticausal_edge_mass = 0,
                anticausal_edge_mass_fraction = NA_real_,
                causal_edge_score = NA_real_,
                stringsAsFactors = FALSE
            ),
            edge_table = empty_edges
        ))
    }

    if (!weight_attr %in% colnames(edge_df)) {
        edge_df[[weight_attr]] <- 1
    }

    edge_df$median_time_from <- cluster_median_time[edge_df$from]
    edge_df$median_time_to <- cluster_median_time[edge_df$to]
    edge_df$time_delta <- edge_df$median_time_to - edge_df$median_time_from
    edge_df$temporal_class <- ifelse(
        edge_df$time_delta < -tolerance,
        "anticausal",
        ifelse(abs(edge_df$time_delta) <= tolerance, "neutral", "forward")
    )

    total_weight <- sum(edge_df[[weight_attr]], na.rm = TRUE)
    anticausal_mass <- sum(
        edge_df[[weight_attr]][edge_df$temporal_class == "anticausal"],
        na.rm = TRUE
    )

    out_edge_table <- data.frame(
        from = edge_df$from,
        to = edge_df$to,
        weight = edge_df[[weight_attr]],
        median_time_from = edge_df$median_time_from,
        median_time_to = edge_df$median_time_to,
        time_delta = edge_df$time_delta,
        temporal_class = edge_df$temporal_class,
        stringsAsFactors = FALSE
    )

    list(
        summary = data.frame(
            total_edge_weight = total_weight,
            anticausal_edge_mass = anticausal_mass,
            anticausal_edge_mass_fraction =
                if (total_weight > 0) anticausal_mass / total_weight
                else NA_real_,
            causal_edge_score =
                if (total_weight > 0) 1 - (anticausal_mass / total_weight)
                else NA_real_,
            stringsAsFactors = FALSE
        ),
        edge_table = out_edge_table
    )
}


#' Empty adjacent-timepoint AUC pair table
#'
#' @keywords internal
#' @noRd
.empty_auc_pair_table <- function() {
    data.frame(
        group = character(0),
        time_earlier = numeric(0),
        time_later = numeric(0),
        n_earlier = integer(0),
        n_later = integer(0),
        auc = numeric(0),
        weight = numeric(0),
        stringsAsFactors = FALSE
    )
}


#' Empty adjacent-timepoint AUC result for a given averaging scheme
#'
#' @keywords internal
#' @noRd
.empty_auc_result <- function(average) {
    list(
        summary = data.frame(
            adjacent_timepoint_auc = NA_real_,
            n_valid_comparisons = 0L,
            averaging = average,
            stringsAsFactors = FALSE
        ),
        pair_table = .empty_auc_pair_table()
    )
}


#' Mann-Whitney AUC of later versus earlier pseudotime values
#'
#' @keywords internal
#' @noRd
.mann_whitney_auc <- function(values_earlier, values_later) {
    combined <- c(values_earlier, values_later)
    ranks <- rank(combined, ties.method = "average")
    n_earlier <- length(values_earlier)
    n_later <- length(values_later)
    rank_later <- ranks[(n_earlier + 1L):(n_earlier + n_later)]
    u_stat <- sum(rank_later) - n_later * (n_later + 1) / 2
    u_stat / (n_earlier * n_later)
}


#' Build one adjacent-timepoint AUC row per group and timepoint pair
#'
#' @keywords internal
#' @noRd
.auc_pair_rows <- function(
    pseudotime,
    time_numeric,
    group,
    time_levels,
    min_cells
) {
    pair_rows <- lapply(unique(group), function(group_name) {
        idx_group <- which(
            group == group_name &
                is.finite(pseudotime) &
                is.finite(time_numeric)
        )
        if (length(idx_group) == 0L) {
            return(NULL)
        }

        lapply(seq_len(length(time_levels) - 1L), function(k) {
            earlier <- time_levels[k]
            later <- time_levels[k + 1L]
            idx_earlier <- idx_group[time_numeric[idx_group] == earlier]
            idx_later <- idx_group[time_numeric[idx_group] == later]

            if (length(idx_earlier) < min_cells ||
                    length(idx_later) < min_cells) {
                return(NULL)
            }

            n_earlier <- length(idx_earlier)
            n_later <- length(idx_later)
            data.frame(
                group = group_name,
                time_earlier = earlier,
                time_later = later,
                n_earlier = n_earlier,
                n_later = n_later,
                auc = .mann_whitney_auc(
                    pseudotime[idx_earlier],
                    pseudotime[idx_later]
                ),
                weight = n_earlier * n_later,
                stringsAsFactors = FALSE
            )
        })
    })

    Filter(Negate(is.null), Reduce(c, pair_rows, list()))
}


#' Empty directed reachability pair table
#'
#' @keywords internal
#' @noRd
.empty_reachability_pair_table <- function() {
    data.frame(
        cluster_earlier = character(0),
        cluster_later = character(0),
        median_time_earlier = numeric(0),
        median_time_later = numeric(0),
        time_gap = numeric(0),
        size_earlier = integer(0),
        size_later = integer(0),
        forward_reachable = logical(0),
        backward_reachable = logical(0),
        reachability_score = numeric(0),
        weight = numeric(0),
        stringsAsFactors = FALSE
    )
}


#' Empty directed reachability result for a given weighting scheme
#'
#' @keywords internal
#' @noRd
.empty_reachability_result <- function(weight_by) {
    list(
        summary = data.frame(
            n_ordered_cluster_pairs = 0L,
            directed_reachability_concordance = NA_real_,
            weighting = weight_by,
            stringsAsFactors = FALSE
        ),
        pair_table = .empty_reachability_pair_table()
    )
}


#' Test whether one cluster reaches another along directed graph edges
#'
#' @keywords internal
#' @noRd
.cluster_reachable <- function(cluster_graph, from, to) {
    is.finite(igraph::distances(
        cluster_graph,
        v = from,
        to = to,
        mode = "out"
    )[1, 1])
}


#' Score a single ordered cluster pair by directed reachability
#'
#' @keywords internal
#' @noRd
.reachability_score <- function(forward_reachable, backward_reachable) {
    if (forward_reachable && !backward_reachable) {
        1
    } else if (forward_reachable && backward_reachable) {
        0.5
    } else {
        0
    }
}


#' Build one row per forward-in-time ordered cluster pair
#'
#' @keywords internal
#' @noRd
.reachability_pair_rows <- function(
    cluster_graph,
    cluster_names,
    cluster_median_time,
    cluster_size,
    min_time_gap,
    weight_by
) {
    pair_rows <- lapply(cluster_names, function(cluster_a) {
        lapply(setdiff(cluster_names, cluster_a), function(cluster_b) {
            time_gap <- cluster_median_time[[cluster_b]] -
                cluster_median_time[[cluster_a]]
            if (!is.finite(time_gap) || time_gap <= min_time_gap) {
                return(NULL)
            }

            forward_reachable <- .cluster_reachable(
                cluster_graph, cluster_a, cluster_b
            )
            backward_reachable <- .cluster_reachable(
                cluster_graph, cluster_b, cluster_a
            )

            data.frame(
                cluster_earlier = cluster_a,
                cluster_later = cluster_b,
                median_time_earlier = cluster_median_time[[cluster_a]],
                median_time_later = cluster_median_time[[cluster_b]],
                time_gap = time_gap,
                size_earlier = cluster_size[[cluster_a]],
                size_later = cluster_size[[cluster_b]],
                forward_reachable = forward_reachable,
                backward_reachable = backward_reachable,
                reachability_score = .reachability_score(
                    forward_reachable, backward_reachable
                ),
                weight = if (weight_by == "cells") {
                    cluster_size[[cluster_a]] * cluster_size[[cluster_b]]
                } else {
                    1
                },
                stringsAsFactors = FALSE
            )
        })
    })

    Filter(Negate(is.null), Reduce(c, pair_rows, list()))
}
