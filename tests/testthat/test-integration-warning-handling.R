mock_harmony_runner <- function(exports, functions) {
    runner <- .scprokar_run_harmony
    shim <- new.env(parent = environment(runner))
    shim$getNamespaceExports <- function(ns) {
        expect_identical(ns, "harmony")
        exports
    }
    shim$getExportedValue <- function(ns, name) {
        expect_identical(ns, "harmony")
        functions[[name]]
    }
    environment(runner) <- shim
    runner
}

with_warning_errors <- function(code) {
    previous <- options(warn = 2)
    on.exit(options(previous))
    force(code)
}

test_that("Harmony exposes warnings from its current backend", {
    calls <- 0L
    backend <- function(data_mat, meta_data, vars_use, nclust, verbose) {
        calls <<- calls + 1L
        expect_identical(vars_use, "batch")
        expect_identical(nclust, 3L)
        expect_true(verbose)
        warning("Harmony did not converge", call. = FALSE)
        data_mat + 1
    }
    runner <- mock_harmony_runner("RunHarmony", list(RunHarmony = backend))
    x <- matrix(seq_len(12), nrow = 6)
    meta <- data.frame(batch = rep(1:2, each = 3))

    expect_warning(actual <- runner(x, meta, "batch", nclust = 3L,
        verbose = TRUE), "Harmony did not converge")
    expect_equal(actual, x + 1)
    expect_identical(calls, 1L)
})

test_that("Harmony retains signature retries and legacy backend warnings", {
    calls <- 0L
    backend <- function(
        data_mat, meta_data, vars_use, do_pca, nclust, verbose
    ) {
        calls <<- calls + 1L
        if (missing(do_pca)) {
            stop("Legacy backend requires do_pca")
        }
        expect_false(do_pca)
        warning("Legacy Harmony iteration limit reached", call. = FALSE)
        data_mat + 2
    }
    runner <- mock_harmony_runner("RunHarmony", list(RunHarmony = backend))
    x <- matrix(seq_len(12), nrow = 6)
    meta <- data.frame(batch = rep(1:2, each = 3))

    expect_warning(actual <- runner(x, meta, "batch"),
        "Legacy Harmony iteration limit reached")
    expect_equal(actual, x + 2)
    expect_identical(calls, 2L)
})

test_that("Harmony exposes matrix fallback warnings after backend errors", {
    run_calls <- matrix_calls <- 0L
    current <- function(...) {
        run_calls <<- run_calls + 1L
        stop("Current backend is unavailable")
    }
    legacy <- function(
        data_mat, meta_data, vars_use, do_pca, nclust, verbose
    ) {
        matrix_calls <<- matrix_calls + 1L
        expect_false(do_pca)
        warning("Matrix Harmony did not converge", call. = FALSE)
        data_mat + 3
    }
    runner <- mock_harmony_runner(c("RunHarmony", "HarmonyMatrix"),
        list(RunHarmony = current, HarmonyMatrix = legacy))
    x <- matrix(seq_len(12), nrow = 6)
    meta <- data.frame(batch = rep(1:2, each = 3))

    expect_warning(actual <- runner(x, meta, "batch"),
        "Matrix Harmony did not converge")
    expect_equal(actual, x + 3)
    expect_identical(run_calls, 2L)
    expect_identical(matrix_calls, 1L)
})

test_that("Harmony warning errors propagate without compatibility retries", {
    calls <- 0L
    backend <- function(...) {
        calls <<- calls + 1L
        warning("Harmony did not converge", call. = FALSE)
    }
    runner <- mock_harmony_runner(c("RunHarmony", "HarmonyMatrix"),
        list(RunHarmony = backend, HarmonyMatrix = backend))
    x <- matrix(seq_len(12), nrow = 6)
    meta <- data.frame(batch = rep(1:2, each = 3))

    with_warning_errors(expect_error(runner(x, meta, "batch"),
        "converted from warning.*Harmony did not converge"))
    expect_identical(calls, 1L)
})

test_that("Harmony retries genuine errors after a caller handles warnings", {
    calls <- 0L
    backend <- function(
        data_mat, meta_data, vars_use, do_pca, nclust, verbose
    ) {
        calls <<- calls + 1L
        if (missing(do_pca)) {
            warning("Backend diagnostic", call. = FALSE)
            stop("Legacy backend requires do_pca")
        }
        data_mat + 4
    }
    runner <- mock_harmony_runner("RunHarmony", list(RunHarmony = backend))
    x <- matrix(seq_len(12), nrow = 6)
    meta <- data.frame(batch = rep(1:2, each = 3))

    actual <- with_warning_errors(withCallingHandlers(
        runner(x, meta, "batch"),
        warning = function(w) invokeRestart("muffleWarning")
    ))
    expect_equal(actual, x + 4)
    expect_identical(calls, 2L)
})

make_clustering_warning_sce <- function() {
    sce <- make_toy_sce()
    embedding <- matrix(c(
        0, 0.1, 0.2, 3, 3.1, 3.2,
        0, 0.2, 0.1, 3, 3.2, 3.1,
        1, 1.1, 1.2, 2, 2.1, 2.2
    ), nrow = 6, dimnames = list(colnames(sce), paste0("latent", 1:3)))
    SingleCellExperiment::reducedDim(sce, "integrated_warning_test") <-
        embedding
    sce
}

test_that("Clustering preserves edgeless results and selected dimensions", {
    sce <- make_clustering_warning_sce()
    original <- SingleCellExperiment::reducedDim(
        sce, "integrated_warning_test"
    )
    runner <- RunIntegratedClustering
    shim <- new.env(parent = environment(runner))
    shim$.scprokar_build_cluster_graph <- function(emb, k) {
        expect_equal(emb, original[, c(3, 1), drop = FALSE])
        expect_identical(k, 5L)
        list(graph = igraph::make_empty_graph(n = nrow(emb), directed = FALSE))
    }
    environment(runner) <- shim

    expect_warning(actual <- runner(sce, dims = c(3, 1), k = 99,
        algorithm = "walktrap"), NA)
    expect_identical(actual$integrated_clusters, factor(rep(1L, 6)))
    expect_identical(SingleCellExperiment::reducedDim(
        actual, "integrated_warning_test"), original)
    expect_identical(S4Vectors::metadata(actual)$scProkaR$integration$clustering,
        list(reduction = "integrated_warning_test",
            cluster_col = "integrated_clusters", dims = c(3, 1), k = 5L,
            algorithm = "walktrap", resolution = 0.8))
})

test_that("Walktrap preserves its ignored-resolution warning and results", {
    sce <- make_clustering_warning_sce()
    expect_warning(with_resolution <- RunIntegratedClustering(sce,
        algorithm = "walktrap", k = 2, resolution = 0.8),
        "`resolution` is ignored")
    expect_warning(without_resolution <- RunIntegratedClustering(sce,
        algorithm = "walktrap", k = 2, resolution = NULL), NA)
    expect_identical(with_resolution$integrated_clusters,
        without_resolution$integrated_clusters)
    expect_null(S4Vectors::metadata(without_resolution)$scProkaR$
        integration$clustering$resolution)
    expect_error(RunIntegratedClustering(sce, dims = 1),
        "must have at least two dimensions")
})
