test_that("disabled TATA stages preserve existing cell annotations", {
    sce <- make_toy_multidrug_sce(n_cells = 240, n_features = 60, seed = 71)
    sce$timepoint <- NULL
    sce$tata_pseudotime <- seq_len(ncol(sce))
    sce$tata_branch_assignment <- rep("previous", ncol(sce))
    S4Vectors::metadata(sce)$project <- "preserve me"
    original_embedding <- SingleCellExperiment::reducedDim(sce, "PCA")

    result <- run_tata(
        sce, use_time = FALSE, do_pseudotime = FALSE,
        do_branch_probs = FALSE, cluster_method = "walktrap", k = 12,
        root_cluster = stop("Skipped pseudotime must not evaluate the root."),
        n_pcs_if_missing = stop("An existing embedding must be reused.")
    )

    expect_null(result$cell_pseudotime)
    expect_null(result$cluster_pseudotime)
    expect_null(result$terminal_clusters)
    expect_null(result$branch_probabilities)
    expect_null(result$transition_matrix)
    expect_null(result$parameters$root_cluster)
    expect_null(result$parameters$refine_min_cells)
    expect_identical(result$sce$tata_pseudotime, sce$tata_pseudotime)
    expect_identical(result$sce$tata_branch_assignment, sce$tata_branch_assignment)
    expect_identical(S4Vectors::metadata(result$sce)$project, "preserve me")
    expect_identical(
        SingleCellExperiment::reducedDim(result$sce, "PCA"), original_embedding
    )
    expect_identical(S4Vectors::metadata(result$sce)$TATA$parameters, result$parameters)
    expect_identical(
        SingleCellExperiment::reducedDim(result$sce, "TATA"), result$cell_space
    )
    expect_false("TATA" %in% SingleCellExperiment::reducedDimNames(sce))
})


test_that("TATA records realized settings and aligned inference outputs", {
    sce <- make_toy_multidrug_sce(n_cells = 240, n_features = 60, seed = 72)
    set.seed(72)
    result <- run_tata(sce, k = 12.9, prune_threshold = 0.001)

    expect_equal(result$parameters$k, 12L)
    expect_equal(result$parameters$refine_min_cells, 40L)
    expect_identical(
        result$sce$tata_pseudotime, result$cell_pseudotime
    )
    expect_identical(
        result$sce$tata_pseudotime_scaled, result$cell_pseudotime_scaled
    )
    expect_identical(names(result$cell_pseudotime), colnames(sce))
    expect_identical(rownames(result$cell_space), colnames(sce))
    expect_identical(result$branch_probabilities$cell_id, colnames(sce))
    stored <- S4Vectors::metadata(result$sce)$TATA
    expect_identical(stored$parameters, result$parameters)
    expect_identical(stored$cluster_edge_table, result$cluster_edge_table)
    expect_identical(stored$cluster_pseudotime, result$cluster_pseudotime)
    expect_identical(
        stored$cluster_branch_probabilities, result$cluster_branch_probabilities
    )
    expect_identical(stored$terminal_clusters, result$terminal_clusters)
})


test_that("topology-only TATA can infer pseudotime and branches without time labels", {
    sce <- make_toy_multidrug_sce(n_cells = 240, n_features = 60, seed = 73)
    sce$timepoint <- NULL
    set.seed(73)
    result <- run_tata(
        sce, use_time = FALSE, do_pseudotime = TRUE,
        do_branch_probs = TRUE, k = 12, expected_branches = 2
    )

    expect_length(result$cell_pseudotime, ncol(sce))
    expect_true("tata_pseudotime" %in% names(SummarizedExperiment::colData(result$sce)))
    expect_length(result$terminal_clusters, 2L)
    expect_identical(result$branch_probabilities$cell_id, colnames(sce))
    expect_true(all(is.finite(result$cell_space)))
    expect_true(all(result$time_table$flow_score == 0))
    expect_true(all(result$time_table$time_consistency == 1))
    expect_true("tata_branch_assignment" %in% names(SummarizedExperiment::colData(result$sce)))
})


test_that("TATA keeps validation errors before graph construction", {
    sce <- make_toy_multidrug_sce(n_cells = 120, n_features = 60, seed = 74)
    expect_error(run_tata(matrix(1)), "`sce` must be a SingleCellExperiment.", fixed = TRUE)
    expect_error(
        run_tata(sce, time_col = "missing", alpha = -1),
        "Time column `missing` is not present in colData(sce).", fixed = TRUE
    )
    expect_error(
        run_tata(sce, alpha = -1),
        "`alpha` must be a single non-negative numeric value.", fixed = TRUE
    )
    expect_error(
        run_tata(sce, beta = -1),
        "`beta` must be a single non-negative numeric value.", fixed = TRUE
    )
})


test_that("sampled TATA settings are evaluated after cell clustering", {
    sce <- make_toy_multidrug_sce(n_cells = 240, n_features = 60, seed = 71)
    set.seed(75)
    fixed <- run_tata(
        sce, k = 12, refine_clusters = FALSE, refine_min_cells = 40L,
        do_pseudotime = FALSE, do_branch_probs = FALSE
    )
    # Later stages are deterministic here, so this draw represents evaluation
    # of the unused refinement threshold after cell clustering has finished.
    sample(c(40L, 40L), 1L)
    expected_rng <- .Random.seed

    set.seed(75)
    sampled <- run_tata(
        sce, k = 12, refine_clusters = FALSE,
        refine_min_cells = sample(c(40L, 40L), 1L),
        do_pseudotime = FALSE, do_branch_probs = FALSE
    )

    expect_identical(sampled$sce$cluster, fixed$sce$cluster)
    expect_identical(sampled$cluster_edge_table, fixed$cluster_edge_table)
    expect_identical(sampled$parameters, fixed$parameters)
    expect_identical(sampled$cell_space, fixed$cell_space)
    expect_identical(.Random.seed, expected_rng)
})
