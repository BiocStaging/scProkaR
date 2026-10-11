test_that("IntegrateBacData validates batch column", {
    sce <- make_toy_sce()
    expect_error(
        IntegrateBacData(sce, batch_col = "missing_batch", method = "mnn"),
        "Missing batch column"
    )
})

test_that("IntegrateBacData runs optional backends when installed", {
    sce <- make_toy_sce()

    if (requireNamespace("batchelor", quietly = TRUE)) {
        out <- IntegrateBacData(sce, batch_col = "batch", method = "mnn", dims = 1:2)
        expect_true(
            "integrated_mnn" %in% SingleCellExperiment::reducedDimNames(out)
        )
    } else {
        skip("batchelor not installed")
    }

    if (requireNamespace("harmony", quietly = TRUE)) {
        ## Current Harmony adjusts its block size for this tiny fixture.
        ## Accept that documented backend diagnostic, while rejecting others.
        out <- withCallingHandlers(
            IntegrateBacData(
                sce, batch_col = "batch", method = "harmony", dims = 1:2
            ),
            warning = function(w) {
                expect_match(conditionMessage(w),
                    "Too few cells.*Setting block_size")
                invokeRestart("muffleWarning")
            }
        )
        expect_true(
            "integrated_harmony" %in%
                SingleCellExperiment::reducedDimNames(out)
        )
    } else {
        skip("harmony not installed")
    }
})

test_that("RegisterIntegrationEmbedding records custom embeddings", {
    sce <- make_toy_sce()
    embedding <- matrix(
        seq_len(ncol(sce) * 2),
        ncol = 2,
        dimnames = list(colnames(sce), c("dim1", "dim2"))
    )

    out <- RegisterIntegrationEmbedding(sce, embedding = embedding, method_name = "custom")

    expect_true(
        "integrated_custom" %in% SingleCellExperiment::reducedDimNames(out)
    )
    expect_equal(
        S4Vectors::metadata(out)$scProkaR$integration$results$custom$source,
        "external"
    )
})

test_that("PlotReduction returns a ggplot object", {
    sce <- make_toy_sce()
    embedding <- matrix(
        seq_len(ncol(sce) * 2),
        ncol = 2,
        dimnames = list(colnames(sce), c("dim1", "dim2"))
    )
    sce <- RegisterIntegrationEmbedding(sce, embedding = embedding, method_name = "custom")

    expect_s3_class(
        PlotReduction(
            sce,
            reduction = "integrated_custom",
            colour_by = "batch"
        ),
        "ggplot"
    )
})

test_that("PlotIntegrationOverview returns Seurat-style integration panels", {
    sce <- make_toy_sce()
    SummarizedExperiment::colData(sce)$seurat_clusters <- SummarizedExperiment::colData(sce)$cluster

    unintegrated_umap <- matrix(
        seq_len(ncol(sce) * 2),
        ncol = 2,
        dimnames = list(colnames(sce), c("UMAP_1", "UMAP_2"))
    )
    integrated_umap <- matrix(
        seq_len(ncol(sce) * 2) / 10,
        ncol = 2,
        dimnames = list(colnames(sce), c("UMAP_1", "UMAP_2"))
    )

    SingleCellExperiment::reducedDim(sce, "umap") <- unintegrated_umap
    sce <- RegisterIntegrationEmbedding(sce, embedding = integrated_umap, method_name = "custom")
    SingleCellExperiment::reducedDim(sce, "umap_integrated_custom") <- integrated_umap

    plots <- PlotIntegrationOverview(
        sce,
        unintegrated_reduction = "umap",
        integrated_reduction = "umap_integrated_custom",
        batch_col = "batch",
        label_col = "seurat_clusters",
        split_by = "batch"
    )

    expect_true(all(c(
        "unintegrated_by_batch",
        "integrated_by_batch",
        "unintegrated_by_label",
        "integrated_by_label",
        "integrated_split"
    ) %in% names(plots)))
    expect_s3_class(plots$unintegrated_by_batch, "ggplot")
    expect_s3_class(plots$integrated_by_batch, "ggplot")
    expect_s3_class(plots$integrated_split, "ggplot")
})

test_that("RunIntegratedUMAP stores a new UMAP reduction when uwot is available", {
    if (!requireNamespace("uwot", quietly = TRUE)) {
        skip("uwot not installed")
    }

    sce <- make_toy_sce()
    embedding <- matrix(
        rnorm(ncol(sce) * 4),
        ncol = 4,
        dimnames = list(colnames(sce), paste0("latent_", 1:4))
    )
    sce <- RegisterIntegrationEmbedding(sce, embedding = embedding, method_name = "custom")
    sce <- RunIntegratedUMAP(sce, reduction = "integrated_custom")

    expect_true(
        "umap_integrated_custom" %in%
            SingleCellExperiment::reducedDimNames(sce)
    )
})

test_that("RunIntegratedClustering stores cluster labels from an integrated embedding", {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        skip("igraph not installed")
    }

    sce <- make_toy_sce()
    embedding <- matrix(
        c(
            0, 0,
            0.1, 0.1,
            0.2, 0.1,
            5, 5,
            5.1, 5.1,
            5.2, 5.0
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(colnames(sce), c("latent_1", "latent_2"))
    )
    sce <- RegisterIntegrationEmbedding(sce, embedding = embedding, method_name = "custom")
    sce <- RunIntegratedClustering(
        sce,
        reduction = "integrated_custom",
        k = 2,
        resolution = 0.8
    )

    expect_true(
        "integrated_clusters" %in% colnames(SummarizedExperiment::colData(sce))
    )
    expect_equal(
        length(SummarizedExperiment::colData(sce)$integrated_clusters),
        ncol(sce)
    )
    expect_true(length(unique(SummarizedExperiment::colData(sce)$integrated_clusters)) >= 2)
    expect_equal(
        S4Vectors::metadata(sce)$scProkaR$integration$clustering$resolution,
        0.8
    )
})

test_that("RunIntegratedClustering leiden respects resolution", {
    skip_if_not_installed("igraph")
    skip_if_not(
        "cluster_leiden" %in% getNamespaceExports("igraph"),
        "igraph does not provide cluster_leiden()"
    )

    sce <- make_toy_multidrug_sce(n_cells = 900, n_features = 70, seed = 42)
    embedding <- SingleCellExperiment::reducedDim(sce, "PCA")
    sce <- RegisterIntegrationEmbedding(
        sce,
        embedding = embedding,
        method_name = "custom"
    )

    cluster_count <- function(resolution) {
        set.seed(42)
        out <- RunIntegratedClustering(
            sce,
            reduction = "integrated_custom",
            cluster_col = "leiden_clusters",
            algorithm = "leiden",
            k = 20,
            resolution = resolution
        )
        nlevels(factor(
            SummarizedExperiment::colData(out)$leiden_clusters
        ))
    }

    coarse <- cluster_count(0.4)
    fine <- cluster_count(1.6)

    ## Regression guard: with the CPM objective this returned hundreds of
    ## singleton clusters and ignored the resolution entirely.
    expect_lt(coarse, ncol(sce) / 10)
    expect_lt(fine, ncol(sce) / 10)
    expect_gt(coarse, 1)
    expect_gt(fine, coarse)
})
