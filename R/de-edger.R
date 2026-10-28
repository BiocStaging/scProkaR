# Internal helpers for the edgeR workflows.

.scprokar_numeric_labels <- function(values) {
    ## Use R's parser to check the whole vector without coercion warnings.
    ## Logical and complex labels must retain their categorical ordering.
    converted <- utils::type.convert(as.character(values), as.is = TRUE)
    if (is.numeric(converted) && !is.complex(converted)) {
        return(as.numeric(converted))
    }
    rep(NA_real_, length(values))
}

.most_frequent_value <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) {
        return(NA)
    }
    tab <- table(x)
    names(tab)[which.max(tab)]
}


.sanitize_levels <- function(x) {
    stats::setNames(make.names(x), x)
}


.make_default_reference <- function(values) {
    if (is.numeric(values) || is.integer(values)) {
        unique_vals <- sort(unique(as.numeric(values)))
        return(as.character(unique_vals[1]))
    }

    if (is.factor(values)) {
        return(as.character(levels(values)[1]))
    }

    unique_vals <- unique(as.character(values))
    numeric_vals <- .scprokar_numeric_labels(unique_vals)
    if (!anyNA(numeric_vals)) {
        return(as.character(unique_vals[order(numeric_vals)][1]))
    }
    sort(unique_vals)[1]
}


.build_group_factor <- function(meta, group_col, combine_group_cols = NULL) {
    if (is.null(combine_group_cols)) {
        values <- meta[[group_col]]
        if (is.numeric(values) || is.integer(values)) {
            ordered_levels <- as.character(sort(unique(values)))
            return(factor(as.character(values), levels = ordered_levels))
        }
        if (is.factor(values)) {
            return(factor(as.character(values), levels = levels(values)))
        }
        values_chr <- as.character(values)
        values_num <- .scprokar_numeric_labels(values_chr)
        if (!anyNA(values_num)) {
            ordered_levels <- unique(values_chr[order(values_num)])
            return(factor(values_chr, levels = ordered_levels))
        }
        return(factor(values_chr, levels = sort(unique(values_chr))))
    }

    combined <- interaction(
        meta[, combine_group_cols, drop = FALSE],
        drop = TRUE,
        lex.order = TRUE,
        sep = "__"
    )
    factor(as.character(combined), levels = levels(combined))
}


.make_pairwise_pairs <- function(
    group_levels,
    contrast_type,
    reference_level = NULL,
    manual_pairs = NULL
) {
    contrast_type <- match.arg(contrast_type, c("reference", "all", "manual"))

    if (contrast_type == "manual") {
        if (is.null(manual_pairs)) {
            stop(
                "`manual_pairs` must be supplied when `contrast_type = ",
                "\"manual\"`.",
                call. = FALSE
            )
        }
        if (!is.matrix(manual_pairs) && !is.data.frame(manual_pairs)) {
            stop(
                "`manual_pairs` must be a two-column matrix or data.frame.",
                call. = FALSE
            )
        }
        manual_pairs <- as.data.frame(manual_pairs, stringsAsFactors = FALSE)
        if (ncol(manual_pairs) != 2L) {
            stop("`manual_pairs` must have exactly two columns.", call. = FALSE)
        }
        colnames(manual_pairs) <- c("group1", "group2")
        return(manual_pairs)
    }

    if (contrast_type == "reference") {
        if (is.null(reference_level)) {
            reference_level <- group_levels[1]
        }
        if (!reference_level %in% group_levels) {
            stop(
                "`reference_level` is not among the group levels.",
                call. = FALSE
            )
        }
        out <- data.frame(
            group1 = setdiff(group_levels, reference_level),
            group2 = reference_level,
            stringsAsFactors = FALSE
        )
        return(out)
    }

    combn_out <- utils::combn(group_levels, 2)
    data.frame(
        group1 = combn_out[2, ],
        group2 = combn_out[1, ],
        stringsAsFactors = FALSE
    )
}


.prepare_edger_dge <- function(counts, meta = NULL, normalization = "TMM") {
    dge <- edgeR::DGEList(counts = counts, samples = meta)
    if (!is.null(normalization) && normalization != "none") {
        dge <- edgeR::calcNormFactors(dge, method = normalization)
    } else {
        dge$samples$norm.factors <- rep(1, ncol(counts))
    }
    dge
}


.resolve_pairwise_de_table <- function(de_result, contrast) {
    if (is.null(de_result$tables) || !length(de_result$tables)) {
        stop("`de_result` does not contain pairwise DE tables.", call. = FALSE)
    }

    if (is.null(contrast)) {
        contrast <- names(de_result$tables)[1]
    }

    if (!contrast %in% names(de_result$tables)) {
        stop(
            "Contrast `",
            contrast,
            "` is not present in `de_result$tables`.",
            call. = FALSE
        )
    }

    tab <- de_result$tables[[contrast]]
    if (!all(c("gene", "logFC", "FDR") %in% colnames(tab))) {
        stop(
            "The selected DE table does not contain the required columns.",
            call. = FALSE
        )
    }

    tab
}


.rbind_fill_data_frames <- function(df_list) {
    df_list <- Filter(function(x) !is.null(x), df_list)
    if (!length(df_list)) {
        return(data.frame())
    }

    all_cols <- unique(unlist(lapply(df_list, colnames), use.names = FALSE))
    df_list <- lapply(df_list, function(df) {
        df[setdiff(all_cols, colnames(df))] <- list(rep(NA, nrow(df)))
        df[, all_cols, drop = FALSE]
    })

    do.call(rbind, df_list)
}


.predict_curves_from_spline_result <- function(one_result, genes = NULL) {
    if (is.null(one_result$fit) ||
        is.null(one_result$curve_design) ||
        !length(one_result$time_grid)) {
        return(data.frame())
    }

    available_genes <- rownames(stats::coef(one_result$fit))
    if (is.null(genes)) {
        genes <- available_genes
    }
    .build_spline_curve_table(
        fit = one_result$fit,
        curve_genes = intersect(as.character(genes), available_genes),
        curve_design = one_result$curve_design,
        time_grid = one_result$time_grid,
        label = unique(one_result$table$condition)[1],
        tab = one_result$table
    )
}


.predict_curves_from_pseudobulk <- function(
    pb,
    genes,
    time_col,
    condition_col,
    curve_grid_length = 200,
    assay_name = NULL
) {
    if (!is(pb, "SingleCellExperiment")) {
        return(data.frame())
    }

    if (is.null(condition_col) || !nzchar(condition_col)) {
        return(data.frame())
    }

    assay_name <- .scprokar_pick_assay(pb, assay_name, "counts")

    genes <- intersect(as.character(genes), rownames(pb))
    if (!length(genes)) {
        return(data.frame())
    }

    meta <- as.data.frame(colData(pb))
    if (!all(c(time_col, condition_col) %in% colnames(meta))) {
        return(data.frame())
    }

    time_numeric <- .coerce_time_to_numeric(meta[[time_col]])
    condition_values <- as.character(meta[[condition_col]])
    expr_mat <- as.matrix(assay(pb, assay_name)[genes, , drop = FALSE])

    out <- lapply(sort(unique(condition_values)), function(cond) {
        idx_cond <- which(condition_values == cond & is.finite(time_numeric))
        if (!length(idx_cond)) {
            return(NULL)
        }

        cond_time <- time_numeric[idx_cond]
        time_grid <- seq(
            min(cond_time),
            max(cond_time),
            length.out = curve_grid_length
        )

        lapply(genes, function(gene) {
            pred <- .interpolate_gene_time_curve(
                cond_time = cond_time,
                expr_values = expr_mat[gene, idx_cond],
                time_grid = time_grid
            )
            if (is.null(pred)) {
                return(NULL)
            }

            data.frame(
                gene = gene,
                time = time_grid,
                logCPM = pred,
                condition = cond,
                status = "NotSig",
                stringsAsFactors = FALSE
            )
        })
    })

    .rbind_fill_data_frames(unlist(out, recursive = FALSE))
}


#' Run pairwise pseudobulk differential expression with edgeR
#'
#' Fit a quasi-likelihood negative-binomial model with edgeR and test pairwise
#' contrasts between groups such as timepoints, treatments, or composite
#' treatment-timepoint groups.
#'
#' @param pb A pseudobulk `SingleCellExperiment`.
#' @param assay_name Assay containing pseudobulk counts.
#' @param group_col Sample metadata column defining the primary grouping, for
#'   example `"timepoint"`.
#' @param covariates Optional metadata columns to include as additive covariates
#'   in the design matrix, for example treatment, batch, or donor.
#' @param combine_group_cols Optional metadata columns to combine into a single
#'   composite group via `interaction()`, for example
#'   `c("treatment", "timepoint")`.
#' @param contrast_type Contrast construction strategy. Choices are:
#' - `"reference"`: compare every non-reference group against one reference
#'   level. This is the default.
#' - `"all"`: test all pairwise group contrasts.
#' - `"manual"`: use `manual_pairs` exactly as supplied.
#' @param reference_level Optional reference level used when
#'   `contrast_type = "reference"`. By default the earliest or first group is
#'   used.
#' @param manual_pairs Optional two-column matrix or data frame of manual
#'   contrasts with `group1` and `group2` values. Each test is
#'   `group1 - group2`.
#' @param normalization Library-size normalization method. Choices are:
#' - `"TMM"`: edgeR TMM normalization.
#' - `"upperquartile"`: edgeR upper-quartile normalization.
#' - `"none"`: no normalization factors beyond raw library sizes.
#' @param filter_by_expr Logical indicating whether to run
#'   `edgeR::filterByExpr()`.
#' @param robust Logical passed to `edgeR::glmQLFit()`.
#' @param abundance_trend Logical passed to `edgeR::glmQLFit()`.
#' @param fdr_cutoff FDR cutoff used to mark significant genes.
#'
#' @return A list containing the fitted edgeR objects, contrast matrix, result
#'   tables, and analysis settings.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' sce$replicate <- rep(paste0("rep", 1:3), length.out = ncol(sce))
#'
#' pb <- aggregate_pseudobulk(
#'     sce,
#'     sample_cols = c("condition", "timepoint", "replicate"),
#'     aggregation = "sum",
#'     min_cells = 5
#' )
#' pb <- normalize_pseudobulk(pb, method = "TMM")
#'
#' # Compare every timepoint against the baseline, adjusting for drug.
#' de <- run_edger_pairwise_de(
#'     pb,
#'     group_col = "timepoint",
#'     covariates = "condition",
#'     contrast_type = "reference"
#' )
#' names(de$tables)
#' head(de$tables[["24_vs_0"]][, c("gene", "logFC", "FDR", "significant")])
#' table(de$combined_table$contrast, de$combined_table$significant)
run_edger_pairwise_de <- function(
    pb,
    assay_name = "counts",
    group_col = "timepoint",
    covariates = NULL,
    combine_group_cols = NULL,
    contrast_type = c("reference", "all", "manual"),
    reference_level = NULL,
    manual_pairs = NULL,
    normalization = c("TMM", "upperquartile", "none"),
    filter_by_expr = TRUE,
    robust = TRUE,
    abundance_trend = FALSE,
    fdr_cutoff = 0.05
) {
    contrast_type <- match.arg(contrast_type)
    normalization <- match.arg(normalization)

    complete <- .filter_complete_de_samples(
        pb = pb,
        assay_name = assay_name,
        needed_cols = c(group_col, combine_group_cols, covariates),
        "All grouping and covariate columns must exist in `colData(pb)`."
    )
    meta <- complete$meta
    counts <- complete$counts

    groups <- .resolve_pairwise_groups(
        meta = meta,
        group_col = group_col,
        combine_group_cols = combine_group_cols,
        contrast_type = contrast_type,
        reference_level = reference_level
    )
    group_factor <- groups$group_factor
    group_levels <- groups$group_levels
    reference_level <- groups$reference_level

    safe_group_levels <- .sanitize_levels(group_levels)
    design <- .build_pairwise_design(
        meta = meta,
        group_factor = group_factor,
        covariates = covariates,
        safe_group_levels = safe_group_levels
    )

    edger_fit <- .fit_edger_qlf(
        counts = counts,
        meta = meta,
        design = design,
        normalization = normalization,
        filter_by_expr = filter_by_expr,
        robust = robust,
        abundance_trend = abundance_trend,
        group = group_factor
    )

    pair_df <- .make_pairwise_pairs(
        group_levels = group_levels,
        contrast_type = contrast_type,
        reference_level = reference_level,
        manual_pairs = manual_pairs
    )

    contrasts <- .run_pairwise_contrasts(
        fit = edger_fit$fit,
        design = design,
        pair_df = pair_df,
        group_levels = group_levels,
        safe_group_levels = safe_group_levels,
        fdr_cutoff = fdr_cutoff
    )

    list(
        dge = edger_fit$dge,
        design = design,
        fit = edger_fit$fit,
        group_factor = group_factor,
        contrast_matrix = contrasts$contrast_matrix,
        tests = contrasts$tests,
        tables = contrasts$tables,
        combined_table = do.call(rbind, contrasts$tables),
        settings = list(
            assay_name = assay_name,
            group_col = group_col,
            covariates = covariates,
            combine_group_cols = combine_group_cols,
            contrast_type = contrast_type,
            reference_level = reference_level,
            normalization = normalization,
            filter_by_expr = filter_by_expr,
            robust = robust,
            abundance_trend = abundance_trend,
            fdr_cutoff = fdr_cutoff
        )
    )
}


#' Plot a volcano plot for a pairwise edgeR contrast
#'
#' Visualize one pairwise contrast from [run_edger_pairwise_de()] as a volcano
#' plot using log-fold change and FDR.
#'
#' @param de_result Result list returned by `run_edger_pairwise_de()`.
#' @param contrast Name of the contrast to plot. If `NULL`, the first available
#'   contrast is used.
#' @param fdr_cutoff FDR threshold used to define significant genes.
#' @param lfc_cutoff Absolute log-fold-change threshold used to highlight larger
#'   effects.
#' @param top_n_labels Number of top genes to label by FDR. Set to `0` to omit
#'   labels.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' sce$replicate <- rep(paste0("rep", 1:3), length.out = ncol(sce))
#'
#' pb <- aggregate_pseudobulk(
#'     sce,
#'     sample_cols = c("condition", "timepoint", "replicate"),
#'     aggregation = "sum",
#'     min_cells = 5
#' )
#' pb <- normalize_pseudobulk(pb, method = "TMM")
#'
#' de <- run_edger_pairwise_de(
#'     pb,
#'     group_col = "timepoint",
#'     covariates = "condition",
#'     contrast_type = "reference"
#' )
#'
#' p <- plot_pairwise_de_volcano(
#'     de,
#'     contrast = "120_vs_0",
#'     fdr_cutoff = 0.05,
#'     lfc_cutoff = 0.25,
#'     top_n_labels = 5
#' )
#' p
plot_pairwise_de_volcano <- function(
    de_result,
    contrast = NULL,
    fdr_cutoff = 0.05,
    lfc_cutoff = 0,
    top_n_labels = 10L,
    point_size = 1.5,
    point_alpha = 0.75
) {
    tab <- .resolve_pairwise_de_table(de_result, contrast = contrast)
    contrast_name <- unique(tab$contrast)[1]

    tab$neg_log10_fdr <- -log10(pmax(tab$FDR, .Machine$double.xmin))
    tab$de_class <- ifelse(
        tab$FDR <= fdr_cutoff & tab$logFC >= lfc_cutoff,
        "Up",
        ifelse(
            tab$FDR <= fdr_cutoff & tab$logFC <= -lfc_cutoff,
            "Down",
            "NotSig"
        )
    )
    tab$de_class <- factor(tab$de_class, levels = c("Down", "NotSig", "Up"))

    p <- ggplot2::ggplot(
        tab,
        ggplot2::aes(x = logFC, y = neg_log10_fdr, colour = de_class)
    ) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::geom_vline(
            xintercept = c(-lfc_cutoff, lfc_cutoff),
            linetype = 2,
            colour = "grey60"
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(fdr_cutoff),
            linetype = 2,
            colour = "grey60"
        ) +
        ggplot2::scale_colour_manual(
            values = c(Down = "#2c7bb6", NotSig = "grey70", Up = "#d7191c")
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = paste("Volcano plot:", contrast_name),
            x = "logFC",
            y = expression(-log[10](FDR)),
            colour = "status"
        )

    top_n_labels <- as.integer(top_n_labels)
    if (!is.na(top_n_labels) && top_n_labels > 0L) {
        label_df <- tab[
            order(tab$FDR, -abs(tab$logFC)),
            c("gene", "logFC", "neg_log10_fdr"),
            drop = FALSE
        ]
        label_df <- utils::head(label_df, top_n_labels)
        p <- p + ggplot2::geom_text(
            data = label_df,
            mapping = ggplot2::aes(x = logFC, y = neg_log10_fdr, label = gene),
            inherit.aes = FALSE,
            size = 3,
            vjust = -0.4,
            check_overlap = TRUE
        )
    }

    p
}


#' Plot an MA plot for a pairwise edgeR contrast
#'
#' Visualize one pairwise contrast from [run_edger_pairwise_de()] as an MA plot
#' using average abundance and log-fold change.
#'
#' @param de_result Result list returned by `run_edger_pairwise_de()`.
#' @param contrast Name of the contrast to plot. If `NULL`, the first available
#'   contrast is used.
#' @param fdr_cutoff FDR threshold used to define significant genes.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' sce$replicate <- rep(paste0("rep", 1:3), length.out = ncol(sce))
#'
#' pb <- aggregate_pseudobulk(
#'     sce,
#'     sample_cols = c("condition", "timepoint", "replicate"),
#'     aggregation = "sum",
#'     min_cells = 5
#' )
#' pb <- normalize_pseudobulk(pb, method = "TMM")
#'
#' de <- run_edger_pairwise_de(
#'     pb,
#'     group_col = "timepoint",
#'     covariates = "condition",
#'     contrast_type = "reference"
#' )
#'
#' p <- plot_pairwise_de_ma(de, contrast = "120_vs_0", fdr_cutoff = 0.05)
#' p
plot_pairwise_de_ma <- function(
    de_result,
    contrast = NULL,
    fdr_cutoff = 0.05,
    point_size = 1.5,
    point_alpha = 0.75
) {
    tab <- .resolve_pairwise_de_table(de_result, contrast = contrast)
    contrast_name <- unique(tab$contrast)[1]

    if (!"logCPM" %in% colnames(tab)) {
        stop(
            "The selected DE table does not contain a `logCPM` column for MA ",
            "plotting.",
            call. = FALSE
        )
    }

    tab$de_class <- ifelse(
        tab$FDR <= fdr_cutoff & tab$logFC > 0,
        "Up",
        ifelse(tab$FDR <= fdr_cutoff & tab$logFC < 0, "Down", "NotSig")
    )
    tab$de_class <- factor(tab$de_class, levels = c("Down", "NotSig", "Up"))

    ggplot2::ggplot(
        tab,
        ggplot2::aes(x = logCPM, y = logFC, colour = de_class)
    ) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
        ggplot2::scale_colour_manual(
            values = c(Down = "#2c7bb6", NotSig = "grey70", Up = "#d7191c")
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            title = paste("MA plot:", contrast_name),
            x = "average logCPM",
            y = "logFC",
            colour = "status"
        )
}


.fit_one_spline_de <- function(
    counts,
    meta,
    time_col,
    covariates,
    normalization,
    df,
    filter_by_expr,
    robust,
    abundance_trend,
    fdr_cutoff,
    return_curve_for,
    curve_grid_length,
    label = "all"
) {
    spline_model <- .build_spline_de_design(
        meta = meta,
        time_col = time_col,
        covariates = covariates,
        df = df
    )
    design <- spline_model$design

    edger_fit <- .fit_edger_qlf(
        counts = counts,
        meta = meta,
        design = design,
        normalization = normalization,
        filter_by_expr = filter_by_expr,
        robust = robust,
        abundance_trend = abundance_trend
    )
    dge <- edger_fit$dge
    fit <- edger_fit$fit

    spline_coef <- seq_len(ncol(spline_model$spline_df)) + 1L
    qlf <- edgeR::glmQLFTest(fit, coef = spline_coef)

    tab <- edgeR::topTags(qlf, n = Inf)$table
    tab$gene <- rownames(tab)
    tab$condition <- label
    tab$significant <- tab$FDR <= fdr_cutoff

    if (return_curve_for == "significant") {
        curve_genes <- tab$gene[tab$significant]
    } else {
        curve_genes <- tab$gene
    }

    curve_model <- .build_spline_curve_design(
        meta = meta,
        covariates = covariates,
        spline_model = spline_model,
        curve_grid_length = curve_grid_length
    )

    curve_df <- .build_spline_curve_table(
        fit = fit,
        curve_genes = curve_genes,
        curve_design = curve_model$curve_design,
        time_grid = curve_model$time_grid,
        label = label,
        tab = tab
    )

    list(
        dge = dge,
        design = design,
        fit = fit,
        test = qlf,
        table = tab,
        curve_data = curve_df,
        time_grid = curve_model$time_grid,
        curve_design = curve_model$curve_design,
        effective_df = spline_model$effective_df
    )
}


#' Run pseudobulk spline or time-series differential expression with edgeR
#'
#' Fit spline-based edgeR models to identify genes that change over time within
#' each condition, then optionally return fitted expression curves for
#' visualization.
#'
#' @param pb A pseudobulk `SingleCellExperiment`.
#' @param time_col Sample metadata column containing numeric or ordered time.
#' @param condition_col Optional sample metadata column defining conditions
#'   to be analysed separately, for example treatment.
#' @param covariates Optional additive covariates to include in each
#'   condition-specific model.
#' @param assay_name Assay containing pseudobulk counts.
#' @param normalization Library-size normalization method. Choices are:
#' - `"TMM"`: edgeR TMM normalization.
#' - `"upperquartile"`: edgeR upper-quartile normalization.
#' - `"none"`: no normalization factors beyond raw library sizes.
#' @param df Degrees of freedom for the natural spline basis.
#' @param filter_by_expr Logical indicating whether to run
#'   `edgeR::filterByExpr()`.
#' @param robust Logical passed to `edgeR::glmQLFit()`.
#' @param abundance_trend Logical passed to `edgeR::glmQLFit()`.
#' @param fdr_cutoff FDR cutoff used to define significant genes.
#' @param return_curve_for Which genes should receive fitted spline curves.
#'   Choices are:
#' - `"significant"`: return curves only for significant genes.
#' - `"all"`: return curves for all tested genes.
#' @param curve_grid_length Number of points used when drawing fitted spline
#'   curves.
#'
#' @return A list containing condition-specific edgeR fits, tables of temporal
#'   DE genes, and fitted curve data for plotting.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' sce$replicate <- rep(paste0("rep", 1:3), length.out = ncol(sce))
#'
#' pb <- aggregate_pseudobulk(
#'     sce,
#'     sample_cols = c("condition", "timepoint", "replicate"),
#'     aggregation = "sum",
#'     min_cells = 5
#' )
#' pb <- normalize_pseudobulk(pb, method = "TMM")
#'
#' # Six timepoints per drug, so a natural spline in time is estimable.
#' spline_de <- run_edger_spline_de(
#'     pb,
#'     time_col = "timepoint",
#'     condition_col = "condition",
#'     df = 3,
#'     curve_grid_length = 50
#' )
#' names(spline_de$by_condition)
#' head(spline_de$combined_table[, c("gene", "FDR", "condition")])
#' table(
#'     spline_de$combined_table$condition,
#'     spline_de$combined_table$significant
#' )
run_edger_spline_de <- function(
    pb,
    time_col = "timepoint",
    condition_col = NULL,
    covariates = NULL,
    assay_name = "counts",
    normalization = c("TMM", "upperquartile", "none"),
    df = 3,
    filter_by_expr = TRUE,
    robust = TRUE,
    abundance_trend = FALSE,
    fdr_cutoff = 0.05,
    return_curve_for = c("significant", "all"),
    curve_grid_length = 200
) {
    normalization <- match.arg(normalization)
    return_curve_for <- match.arg(return_curve_for)

    complete <- .filter_complete_de_samples(
        pb = pb,
        assay_name = assay_name,
        needed_cols = c(time_col, condition_col, covariates),
        "All requested time, condition, and covariate columns must exist ",
        "in `colData(pb)`."
    )
    keep <- complete$keep
    meta <- complete$meta
    counts <- complete$counts

    settings <- list(
        time_col = time_col,
        condition_col = condition_col,
        covariates = covariates,
        assay_name = assay_name,
        normalization = normalization,
        df = df,
        filter_by_expr = filter_by_expr,
        robust = robust,
        abundance_trend = abundance_trend,
        fdr_cutoff = fdr_cutoff,
        return_curve_for = return_curve_for,
        curve_grid_length = curve_grid_length
    )

    if (is.null(condition_col)) {
        res <- .fit_spline_de_with_settings(
            counts = counts,
            meta = meta,
            settings = settings,
            label = "all"
        )
        return(list(
            by_condition = list(all = res),
            combined_table = res$table,
            curve_data = res$curve_data,
            pseudobulk = pb[, keep, drop = FALSE],
            settings = settings
        ))
    }

    res_list <- .fit_spline_de_by_condition(
        counts = counts,
        meta = meta,
        condition_col = condition_col,
        settings = settings
    )

    if (!length(res_list)) {
        stop(
            "No condition had enough pseudobulk samples for spline DE.",
            call. = FALSE
        )
    }

    res_list <- .add_spline_condition_column(res_list, condition_col)

    combined_table <- .rbind_fill_data_frames(lapply(res_list, `[[`, "table"))
    curve_data <- .rbind_fill_data_frames(lapply(res_list, `[[`, "curve_data"))

    list(
        by_condition = res_list,
        combined_table = combined_table,
        curve_data = curve_data,
        pseudobulk = pb[, keep, drop = FALSE],
        settings = settings
    )
}


#' Plot fitted spline expression curves from time-series DE
#'
#' Plot fitted spline curves for selected genes, optionally overlaid across
#' conditions. When a requested gene was significant in one condition but not
#' another, the non-significant condition is still drawn using a dotted line.
#'
#' @param spline_result Result list returned by `run_edger_spline_de()`.
#' @param genes Optional character vector of genes to plot. If `NULL`, the top
#'   genes ranked by FDR are used.
#' @param top_n Number of genes to plot when `genes = NULL`.
#' @param scales Facet scale setting passed to `ggplot2::facet_wrap()`.
#' @param ncol Optional number of columns used by `ggplot2::facet_wrap()`.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' sce$replicate <- rep(paste0("rep", 1:3), length.out = ncol(sce))
#'
#' pb <- aggregate_pseudobulk(
#'     sce,
#'     sample_cols = c("condition", "timepoint", "replicate"),
#'     aggregation = "sum",
#'     min_cells = 5
#' )
#' pb <- normalize_pseudobulk(pb, method = "TMM")
#'
#' spline_de <- run_edger_spline_de(
#'     pb,
#'     time_col = "timepoint",
#'     condition_col = "condition",
#'     df = 3,
#'     curve_grid_length = 50
#' )
#'
#' # Top temporally variable genes, overlaid across the three drugs.
#' p <- plot_time_series_deg_curves(spline_de, top_n = 4, ncol = 2)
#' p
plot_time_series_deg_curves <- function(
    spline_result,
    genes = NULL,
    top_n = 12L,
    scales = "free_y",
    ncol = NULL
) {
    table_df <- spline_result$combined_table

    if (is.null(genes)) {
        rank_df <- stats::aggregate(FDR ~ gene, data = table_df, FUN = min)
        rank_df <- rank_df[order(rank_df$FDR, rank_df$gene), , drop = FALSE]
        genes <- utils::head(rank_df$gene, as.integer(top_n))
    }

    curve_df <- .rbind_fill_data_frames(lapply(
        spline_result$by_condition,
        .predict_curves_from_spline_result,
        genes = genes
    ))

    fallback_df <- .predict_curves_from_pseudobulk(
        pb = spline_result$pseudobulk,
        genes = genes,
        time_col = spline_result$settings$time_col,
        condition_col = spline_result$settings$condition_col,
        curve_grid_length = spline_result$settings$curve_grid_length
    )

    if (nrow(fallback_df) > 0L) {
        existing_keys <- if (nrow(curve_df) > 0L) {
            paste(curve_df$gene, curve_df$condition, sep = "||")
        } else {
            character(0)
        }
        fallback_keys <- paste(
            fallback_df$gene,
            fallback_df$condition,
            sep = "||"
        )
        fallback_df <- fallback_df[
            !fallback_keys %in% existing_keys, ,
            drop = FALSE
        ]
        curve_df <- .rbind_fill_data_frames(list(curve_df, fallback_df))
    }

    if (!nrow(curve_df)) {
        stop(
            "No fitted spline curves are available for the requested genes.",
            call. = FALSE
        )
    }

    if (!is.null(spline_result$settings$condition_col) &&
        spline_result$settings$condition_col %in% colnames(curve_df)) {
        curve_df$condition <- curve_df[[spline_result$settings$condition_col]]
    }

    curve_df$gene <- factor(curve_df$gene, levels = genes)
    curve_df$status <- factor(curve_df$status, levels = c("Sig", "NotSig"))

    colour_label <- if (is.null(spline_result$settings$condition_col)) {
        "condition"
    } else {
        spline_result$settings$condition_col
    }

    ggplot2::ggplot(
        curve_df,
        ggplot2::aes(
            x = time,
            y = logCPM,
            colour = condition,
            linetype = status,
            group = interaction(gene, condition)
        )
    ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::scale_linetype_manual(
            values = c(Sig = "solid", NotSig = "dotted")
        ) +
        ggplot2::facet_wrap(~gene, scales = scales, ncol = ncol) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = spline_result$settings$time_col,
            y = "fitted logCPM",
            colour = colour_label,
            linetype = "status"
        )
}


#' Run a full scProkaR differential expression workflow
#'
#' Aggregate single-cell counts into pseudobulk samples and then dispatch to a
#' pairwise or time-series edgeR workflow.
#'
#' @param sce A `SingleCellExperiment`.
#' @param mode Workflow mode. Choices are:
#' - `"pairwise"`: pseudobulk followed by pairwise edgeR contrasts.
#' - `"time_series"`: pseudobulk followed by spline-based time-series edgeR.
#' @param sample_cols Metadata columns used to define pseudobulk samples.
#' @param assay_name Assay containing input counts.
#' @param aggregation Pseudobulk aggregation method. Choices are:
#' - `"sum"`: sum counts within each pseudobulk sample.
#' - `"mean"`: average counts within each pseudobulk sample.
#' @param min_cells Minimum number of cells required per pseudobulk sample.
#' @param normalization Pseudobulk normalization method. Choices are:
#' - `"TMM"`: edgeR TMM normalization.
#' - `"upperquartile"`: edgeR upper-quartile normalization.
#' - `"none"`: no normalization factors beyond raw library sizes.
#' @param ... Additional arguments passed to `run_edger_pairwise_de()` or
#'   `run_edger_spline_de()`, depending on `mode`.
#'
#' @return A list with the pseudobulk object, the DE result, and workflow
#'   settings.
#' @export
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 600, n_features = 60, n_pcs = 10
#' )
#' sce$replicate <- rep(paste0("rep", 1:3), length.out = ncol(sce))
#'
#' # Pseudobulk aggregation plus pairwise edgeR contrasts in one call.
#' pairwise_res <- run_scproka_de(
#'     sce,
#'     mode = "pairwise",
#'     sample_cols = c("condition", "timepoint", "replicate"),
#'     min_cells = 5,
#'     group_col = "timepoint",
#'     covariates = "condition"
#' )
#' dim(pairwise_res$pseudobulk)
#' names(pairwise_res$de_result$tables)
#'
#' # The same pseudobulk design analysed as a spline time series.
#' ts_res <- run_scproka_de(
#'     sce,
#'     mode = "time_series",
#'     sample_cols = c("condition", "timepoint", "replicate"),
#'     min_cells = 5,
#'     time_col = "timepoint",
#'     condition_col = "condition",
#'     curve_grid_length = 50
#' )
#' head(ts_res$de_result$combined_table[, c("gene", "FDR", "condition")])
run_scproka_de <- function(
    sce,
    mode = c("pairwise", "time_series"),
    sample_cols,
    assay_name = "counts",
    aggregation = c("sum", "mean"),
    min_cells = 10L,
    normalization = c("TMM", "upperquartile", "none"),
    ...
) {
    mode <- match.arg(mode)
    aggregation <- match.arg(aggregation)
    normalization <- match.arg(normalization)

    pb <- aggregate_pseudobulk(
        sce = sce,
        sample_cols = sample_cols,
        assay_name = assay_name,
        aggregation = aggregation,
        min_cells = min_cells,
        add_logcounts = TRUE
    )

    pb <- normalize_pseudobulk(
        pb,
        assay_name = "counts",
        method = normalization,
        output_assay = "logcounts"
    )

    de_result <- if (mode == "pairwise") {
        run_edger_pairwise_de(
            pb,
            assay_name = "counts",
            normalization = normalization,
            ...
        )
    } else {
        run_edger_spline_de(
            pb,
            assay_name = "counts",
            normalization = normalization,
            ...
        )
    }

    list(
        pseudobulk = pb,
        de_result = de_result,
        settings = list(
            mode = mode,
            sample_cols = sample_cols,
            assay_name = assay_name,
            aggregation = aggregation,
            min_cells = min_cells,
            normalization = normalization
        )
    )
}


#' Validate that a pseudobulk object carries the requested assay
#'
#' @keywords internal
#' @noRd
.check_pb_assay <- function(pb, assay_name) {
    if (!is(pb, "SingleCellExperiment")) {
        stop("`pb` must be a SingleCellExperiment.", call. = FALSE)
    }
    if (!assay_name %in% assayNames(pb)) {
        stop("Assay `", assay_name, "` is not present in `pb`.", call. = FALSE)
    }
    invisible(TRUE)
}


#' Keep pseudobulk samples with complete values for the modelled columns
#'
#' @keywords internal
#' @noRd
.filter_complete_de_samples <- function(pb, assay_name, needed_cols, ...) {
    .check_pb_assay(pb, assay_name)
    meta <- as.data.frame(colData(pb))
    needed_cols <- unique(needed_cols)
    if (!all(needed_cols %in% colnames(meta))) {
        stop(..., call. = FALSE)
    }

    keep <- stats::complete.cases(meta[, needed_cols, drop = FALSE])
    if (!any(keep)) {
        stop(
            "No pseudobulk samples have complete values for the requested ",
            "model.",
            call. = FALSE
        )
    }

    list(
        keep = keep,
        meta = meta[keep, , drop = FALSE],
        counts = assay(pb, assay_name)[, keep, drop = FALSE]
    )
}


#' Build the pairwise group factor and resolve its reference level
#'
#' @keywords internal
#' @noRd
.resolve_pairwise_groups <- function(
    meta,
    group_col,
    combine_group_cols,
    contrast_type,
    reference_level
) {
    group_factor <- .build_group_factor(
        meta,
        group_col = group_col,
        combine_group_cols = combine_group_cols
    )
    group_levels <- levels(group_factor)
    if (length(group_levels) < 2L) {
        stop("At least two groups are required for pairwise DE.", call. = FALSE)
    }

    if (is.null(reference_level) && contrast_type == "reference") {
        reference_level <- .resolve_default_reference_level(
            meta = meta,
            group_col = group_col,
            combine_group_cols = combine_group_cols,
            group_levels = group_levels
        )
    }

    list(
        group_factor = group_factor,
        group_levels = group_levels,
        reference_level = reference_level
    )
}


#' Choose the default reference level for pairwise contrasts
#'
#' @keywords internal
#' @noRd
.resolve_default_reference_level <- function(
    meta,
    group_col,
    combine_group_cols,
    group_levels
) {
    source_col <- if (is.null(combine_group_cols)) {
        group_col
    } else {
        combine_group_cols[1]
    }
    reference_level <- .make_default_reference(meta[[source_col]])
    if (!reference_level %in% group_levels) {
        reference_level <- group_levels[1]
    }
    reference_level
}


#' Build the pairwise design matrix with sanitized group column names
#'
#' @keywords internal
#' @noRd
.build_pairwise_design <- function(
    meta,
    group_factor,
    covariates,
    safe_group_levels
) {
    design_df <- data.frame(
        group_factor = group_factor,
        stringsAsFactors = FALSE
    )
    if (!is.null(covariates) && length(covariates) > 0L) {
        design_df <- cbind(design_df, meta[, covariates, drop = FALSE])
    }

    design <- stats::model.matrix(~ 0 + ., data = design_df)
    group_cols <- seq_along(safe_group_levels)
    colnames(design)[group_cols] <- paste0("group__", unname(safe_group_levels))
    colnames(design) <- make.names(colnames(design))
    design
}


#' Normalize, filter, and fit a quasi-likelihood edgeR model
#'
#' @keywords internal
#' @noRd
.fit_edger_qlf <- function(
    counts,
    meta,
    design,
    normalization,
    filter_by_expr,
    robust,
    abundance_trend,
    group = NULL
) {
    dge <- .prepare_edger_dge(
        counts = as.matrix(counts),
        meta = meta,
        normalization = normalization
    )
    if (isTRUE(filter_by_expr)) {
        keep_genes <- edgeR::filterByExpr(
            dge,
            design = design,
            group = group
        )
        dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
    }

    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmQLFit(
        dge,
        design,
        robust = robust,
        abundance.trend = abundance_trend
    )
    list(dge = dge, fit = fit)
}


#' Build the contrast vector testing `g1 - g2` in a pairwise design
#'
#' @keywords internal
#' @noRd
.pairwise_contrast_vector <- function(
    design,
    safe_group_levels,
    group_levels,
    g1,
    g2
) {
    if (!g1 %in% group_levels || !g2 %in% group_levels) {
        stop(
            "All manual contrast groups must appear in the design.",
            call. = FALSE
        )
    }

    contrast_vec <- rep(0, ncol(design))
    contrast_vec[
        match(paste0("group__", safe_group_levels[[g1]]), colnames(design))
    ] <- 1
    contrast_vec[
        match(paste0("group__", safe_group_levels[[g2]]), colnames(design))
    ] <- -1
    contrast_vec
}


#' Test every requested pairwise contrast and collect the result tables
#'
#' @keywords internal
#' @noRd
.run_pairwise_contrasts <- function(
    fit,
    design,
    pair_df,
    group_levels,
    safe_group_levels,
    fdr_cutoff
) {
    contrast_names <- paste0(pair_df$group1, "_vs_", pair_df$group2)
    contrast_matrix <- matrix(
        0,
        nrow = ncol(design),
        ncol = nrow(pair_df),
        dimnames = list(colnames(design), contrast_names)
    )
    res <- lapply(seq_len(nrow(pair_df)), function(i) {
        g1 <- pair_df$group1[i]
        g2 <- pair_df$group2[i]
        contrast_vec <- .pairwise_contrast_vector(
            design = design,
            safe_group_levels = safe_group_levels,
            group_levels = group_levels,
            g1 = g1,
            g2 = g2
        )

        qlf <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
        tab <- edgeR::topTags(qlf, n = Inf)$table
        tab$gene <- rownames(tab)
        tab$contrast <- contrast_names[i]
        tab$group1 <- g1
        tab$group2 <- g2
        tab$significant <- tab$FDR <= fdr_cutoff
        list(contrast_vec = contrast_vec, test = qlf, table = tab)
    })
    names(res) <- colnames(contrast_matrix)
    contrast_matrix[] <- vapply(
        res, `[[`, numeric(nrow(contrast_matrix)), "contrast_vec"
    )

    list(
        contrast_matrix = contrast_matrix,
        tests = lapply(res, `[[`, "test"),
        tables = lapply(res, `[[`, "table")
    )
}


#' Build the natural-spline design used by the time-series edgeR fit
#'
#' @keywords internal
#' @noRd
.build_spline_de_design <- function(meta, time_col, covariates, df) {
    time_numeric <- .coerce_time_to_numeric(meta[[time_col]])
    if (length(unique(time_numeric)) < 3L) {
        stop(
            "Spline DE requires at least three distinct time values.",
            call. = FALSE
        )
    }

    effective_df <- min(as.integer(df), length(unique(time_numeric)) - 1L)
    if (effective_df < 1L) {
        stop(
            "Spline DE requires at least one estimable spline degree of ",
            "freedom.",
            call. = FALSE
        )
    }

    spline_basis <- splines::ns(time_numeric, df = effective_df)
    spline_df <- as.data.frame(spline_basis)
    colnames(spline_df) <- paste0("spline_", seq_len(ncol(spline_df)))

    design_df <- spline_df
    if (!is.null(covariates) && length(covariates) > 0L) {
        design_df <- cbind(design_df, meta[, covariates, drop = FALSE])
    }
    design_formula <- stats::as.formula(
        paste("~", paste(colnames(design_df), collapse = " + "))
    )
    design <- as.matrix(stats::model.matrix(design_formula, data = design_df))

    list(
        time_numeric = time_numeric,
        spline_basis = spline_basis,
        spline_df = spline_df,
        design_formula = design_formula,
        design = design,
        effective_df = effective_df
    )
}


#' Build a covariate template holding each covariate at a typical value
#'
#' @keywords internal
#' @noRd
.representative_covariate_frame <- function(meta, covariates, n_rows) {
    template <- meta[rep(1, n_rows), covariates, drop = FALSE]
    Reduce(function(template, cov_name) {
        x <- meta[[cov_name]]
        template[[cov_name]] <- if (is.numeric(x) || is.integer(x)) {
            mean(x, na.rm = TRUE)
        } else if (is.factor(x)) {
            factor(rep(.most_frequent_value(x), n_rows), levels = levels(x))
        } else {
            rep(.most_frequent_value(x), n_rows)
        }
        template
    }, as.vector(covariates), template)
}


#' Build the prediction grid and design used to draw fitted spline curves
#'
#' @keywords internal
#' @noRd
.build_spline_curve_design <- function(
    meta,
    covariates,
    spline_model,
    curve_grid_length
) {
    time_grid <- seq(
        min(spline_model$time_numeric),
        max(spline_model$time_numeric),
        length.out = curve_grid_length
    )
    new_basis <- as.data.frame(
        stats::predict(spline_model$spline_basis, newx = time_grid)
    )
    colnames(new_basis) <- colnames(spline_model$spline_df)

    if (!is.null(covariates) && length(covariates) > 0L) {
        template <- .representative_covariate_frame(
            meta = meta,
            covariates = covariates,
            n_rows = curve_grid_length
        )
        new_design_df <- cbind(new_basis, template)
    } else {
        new_design_df <- new_basis
    }

    list(
        time_grid = time_grid,
        curve_design = stats::model.matrix(
            spline_model$design_formula,
            data = new_design_df
        )
    )
}


#' Turn fitted spline coefficients into a long curve table
#'
#' @keywords internal
#' @noRd
.build_spline_curve_table <- function(
    fit,
    curve_genes,
    curve_design,
    time_grid,
    label,
    tab
) {
    if (!length(curve_genes)) {
        return(data.frame())
    }

    coef_mat <- stats::coef(fit)[curve_genes, , drop = FALSE]
    pred_eta <- coef_mat %*% t(curve_design)
    pred_lcpm <- (pred_eta + log(1e6)) / log(2)

    curve_df <- data.frame(
        gene = rep(rownames(pred_lcpm), each = ncol(pred_lcpm)),
        time = rep(time_grid, times = nrow(pred_lcpm)),
        logCPM = as.vector(t(pred_lcpm)),
        condition = label,
        stringsAsFactors = FALSE
    )
    sig_flag <- stats::setNames(tab$significant, tab$gene)[curve_df$gene]
    sig_flag[is.na(sig_flag)] <- FALSE
    curve_df$status <- ifelse(sig_flag, "Sig", "NotSig")
    curve_df[
        order(curve_df$gene, curve_df$condition, curve_df$time), ,
        drop = FALSE
    ]
}


#' Fit one spline DE model using a `run_edger_spline_de()` settings list
#'
#' @keywords internal
#' @noRd
.fit_spline_de_with_settings <- function(counts, meta, settings, label) {
    .fit_one_spline_de(
        counts = counts,
        meta = meta,
        time_col = settings$time_col,
        covariates = settings$covariates,
        normalization = settings$normalization,
        df = settings$df,
        filter_by_expr = settings$filter_by_expr,
        robust = settings$robust,
        abundance_trend = settings$abundance_trend,
        fdr_cutoff = settings$fdr_cutoff,
        return_curve_for = settings$return_curve_for,
        curve_grid_length = settings$curve_grid_length,
        label = label
    )
}


#' Fit one spline DE model per condition, skipping under-powered conditions
#'
#' @keywords internal
#' @noRd
.fit_spline_de_by_condition <- function(counts, meta, condition_col, settings) {
    condition_values <- unique(as.character(meta[[condition_col]]))
    condition_values <- condition_values[order(condition_values)]
    res_list <- lapply(condition_values, function(cond) {
        idx <- which(as.character(meta[[condition_col]]) == cond)
        if (length(idx) < max(settings$df + 1L, 4L)) {
            return(NULL)
        }
        .fit_spline_de_with_settings(
            counts = counts[, idx, drop = FALSE],
            meta = meta[idx, , drop = FALSE],
            settings = settings,
            label = cond
        )
    })
    names(res_list) <- condition_values

    Filter(Negate(is.null), res_list)
}


#' Mirror the internal `condition` column under the user's condition column
#'
#' @keywords internal
#' @noRd
.add_spline_condition_column <- function(res_list, condition_col) {
    if (is.null(condition_col) || !nzchar(condition_col) ||
        condition_col == "condition") {
        return(res_list)
    }

    res_list[] <- lapply(res_list, function(res) {
        res$table[[condition_col]] <- res$table$condition
        if (nrow(res$curve_data) > 0L) {
            res$curve_data[[condition_col]] <- res$curve_data$condition
        }
        res
    })
    res_list
}



#' Interpolate one gene's mean expression onto a time grid
#'
#' @keywords internal
#' @noRd
.interpolate_gene_time_curve <- function(cond_time, expr_values, time_grid) {
    summary_df <- stats::aggregate(
        x = expr_values,
        by = list(time = cond_time),
        FUN = mean
    )
    summary_df <- summary_df[order(summary_df$time), , drop = FALSE]

    if (nrow(summary_df) == 0L) {
        return(NULL)
    }
    if (nrow(summary_df) == 1L) {
        return(rep(summary_df$x[1], length(time_grid)))
    }
    if (nrow(summary_df) == 2L) {
        return(stats::approx(
            x = summary_df$time,
            y = summary_df$x,
            xout = time_grid,
            rule = 2
        )$y)
    }

    interp_fun <- stats::splinefun(
        x = summary_df$time,
        y = summary_df$x,
        method = "natural"
    )
    interp_fun(time_grid)
}
