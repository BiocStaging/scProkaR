test_that("run_tata returns graph, pseudotime, and plotting outputs", {
    sce <- make_toy_multidrug_sce(n_cells = 900, n_features = 70, seed = 31)

    set.seed(31)
    tata <- run_tata(
        sce = sce,
        dimred = "PCA",
        time_col = "timepoint",
        cluster_col = "tata_cluster",
        k = 20,
        cluster_method = "louvain",
        prune_threshold = 0.001,
        time_weight_mode = "directional_confidence",
        expected_branches = 3
    )

    expect_true(all(c(
        "sce",
        "cell_graph",
        "adjacency_matrix",
        "cluster_graph",
        "cluster_edge_table",
        "cluster_pseudotime",
        "cell_pseudotime"
    ) %in% names(tata)))
    expect_s4_class(tata$sce, "SingleCellExperiment")
    expect_true(inherits(tata$cell_graph, "igraph"))
    expect_true(inherits(tata$cluster_graph, "igraph"))
    expect_true(nrow(tata$cluster_edge_table) > 0)
    expect_true(
        "tata_pseudotime" %in%
            colnames(SummarizedExperiment::colData(tata$sce))
    )
    expect_true(length(tata$terminal_clusters) >= 1)
    expect_true(
        "tata_branch_entropy" %in%
            colnames(SummarizedExperiment::colData(tata$sce))
    )
    expect_true(
        "tata_branch_assignment" %in%
            colnames(SummarizedExperiment::colData(tata$sce))
    )

    graph_plot <- plot_tata_cluster_graph(
        tata,
        dimred = "UMAP",
        layout_mode = "embedding",
        show_cells = TRUE,
        colour_by = "timepoint",
        node_colour_by = "cluster"
    )
    abstract_plot <- plot_tata_cluster_graph(
        tata,
        layout_mode = "graph",
        node_colour_by = "median_time"
    )
    trajectory_plot <- plot_tata_trajectory_embedding(
        tata,
        dimred = "UMAP",
        colour_by = "timepoint",
        max_paths = 3
    )
    branch_prob_plot <- plot_tata_branch_probabilities(tata, dimred = "UMAP")
    branch_sel_plot <- plot_tata_branch_selection(
        tata,
        branch = tata$terminal_clusters[1],
        dimred = "UMAP"
    )
    branch_bar_plot <- plot_tata_terminal_probabilities(tata, n_cells = 4)

    feature_modules <- SummarizedExperiment::rowData(tata$sce)$feature_module
    trend_genes <- unique(c(
        rownames(tata$sce)[feature_modules == "branch_A"][1:2],
        rownames(tata$sce)[feature_modules == "branch_B"][1:2],
        rownames(tata$sce)[feature_modules == "condition"][1:2]
    ))
    trend_genes <- unique(stats::na.omit(trend_genes))
    if (length(trend_genes) == 0L) {
        trend_genes <- rownames(tata$sce)[seq_len(min(4, nrow(tata$sce)))]
    }
    ## The tied pseudotimes in this simulation can yield singular local fits.
    ## These diagnostics are expected here and remain visible in package code.
    trend_result <- withCallingHandlers(
        compute_tata_gene_trends(tata, genes = trend_genes),
        warning = function(w) {
            expect_match(conditionMessage(w), paste(
                "pseudoinverse|neighborhood radius|reciprocal condition",
                "near singularities|span too small", sep = "|"
            ))
            invokeRestart("muffleWarning")
        }
    )
    trend_plot <- plot_tata_gene_trends(trend_result, ncol = 2)
    heatmap_plot <- plot_tata_gene_trend_heatmap(trend_result, n_clusters = 2)
    overview_plot <- plot_tata_results(tata, dimred = "UMAP")

    expect_s3_class(graph_plot, "ggplot")
    expect_s3_class(abstract_plot, "ggplot")
    expect_s3_class(trajectory_plot, "ggplot")
    expect_s3_class(branch_prob_plot, "ggplot")
    expect_s3_class(branch_sel_plot, "ggplot")
    expect_s3_class(branch_bar_plot, "ggplot")
    expect_s3_class(trend_plot, "ggplot")
    expect_s3_class(heatmap_plot, "ggplot")
    expect_type(overview_plot, "list")
    expect_s3_class(overview_plot$pseudotime, "ggplot")
})

test_that("TATA benchmark helpers return finite summary statistics", {
    sce <- make_toy_multidrug_sce(n_cells = 900, n_features = 70, seed = 32)
    set.seed(32)
    tata <- run_tata(
        sce = sce,
        dimred = "PCA",
        time_col = "timepoint",
        cluster_col = "tata_cluster",
        k = 20,
        cluster_method = "louvain",
        prune_threshold = 0.001,
        time_weight_mode = "directional_confidence",
        expected_branches = 3
    )

    local_metrics <- local_temporal_order_metrics(
        pseudotime = tata$cell_pseudotime,
        adjacency = tata$adjacency_matrix,
        timepoint = tata$sce$timepoint
    )
    auc_metrics <- adjacent_timepoint_auc(
        pseudotime = tata$cell_pseudotime,
        timepoint = tata$sce$timepoint,
        group = tata$sce$condition,
        min_cells = 5
    )
    onset_metrics <- branch_onset_metrics(
        pseudotime = tata$cell_pseudotime,
        branch = as.character(tata$sce$condition),
        branch_activation = tata$sce$branch_activation,
        true_onset = c(drug_A = 24, drug_C = 48, drug_B = 96)
    )
    reachability <- directed_reachability_concordance(
        tata$cluster_graph,
        clusters = tata$sce$tata_cluster,
        timepoint = tata$sce$timepoint
    )
    causal_mass <- anticausal_edge_mass(
        tata$cluster_graph,
        clusters = tata$sce$tata_cluster,
        timepoint = tata$sce$timepoint
    )

    expect_true(is.data.frame(local_metrics))
    expect_true(is.list(auc_metrics))
    expect_true(is.data.frame(onset_metrics$summary))
    expect_true(is.data.frame(reachability$summary))
    expect_true(is.data.frame(causal_mass$summary))
    expect_true(is.finite(local_metrics$local_temporal_direction_score))
    expect_true(is.finite(auc_metrics$summary$adjacent_timepoint_auc))
})

test_that("expected_branches guides TATA terminal-state selection", {
    sce <- make_toy_multidrug_sce(n_cells = 1200, n_features = 80, seed = 33)

    set.seed(33)
    tata <- run_tata(
        sce = sce,
        dimred = "PCA",
        time_col = "timepoint",
        cluster_col = "tata_cluster",
        k = 20,
        cluster_method = "louvain",
        prune_threshold = 0.001,
        time_weight_mode = "directional_confidence",
        expected_branches = 3
    )

    expect_equal(length(tata$terminal_clusters), 3)
})

test_that("leiden clustering groups cells rather than isolating them", {
    skip_if_not_installed("igraph")
    skip_if_not(
        "cluster_leiden" %in% getNamespaceExports("igraph"),
        "igraph does not provide cluster_leiden()"
    )

    sce <- make_toy_multidrug_sce(n_cells = 900, n_features = 70, seed = 41)
    knn <- build_knn_graph(sce, dimred = "PCA", k = 20)

    set.seed(41)
    leiden_clusters <- cluster_graph_states(knn$graph, method = "leiden")
    set.seed(41)
    louvain_clusters <- cluster_graph_states(knn$graph, method = "louvain")

    n_leiden <- nlevels(factor(leiden_clusters))
    n_louvain <- nlevels(factor(louvain_clusters))

    ## Regression guard: cluster_leiden() defaults to the CPM objective, which
    ## puts every cell in its own community on an unweighted kNN graph.
    expect_length(leiden_clusters, ncol(sce))
    expect_lt(n_leiden, ncol(sce) / 10)
    expect_gt(n_leiden, 1)

    ## Leiden and Louvain optimise the same objective, so they should land in
    ## the same ballpark rather than differing by orders of magnitude.
    expect_lt(n_leiden, 5 * n_louvain)
})
