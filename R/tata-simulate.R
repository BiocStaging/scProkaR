#' Simulate a divergent TATA dataset
#'
#' Simulate a `SingleCellExperiment` containing single-bacterium observations on
#' a divergent two-way trajectory. Time `0` is centred at the origin of the
#' latent trajectory, and later timepoints move away from that central root into
#' two branches. This makes the simulated dataset useful for demonstrating TATA
#' on a divergent trajectory where the earliest state appears in the middle of a
#' UMAP.
#'
#' The returned object contains:
#' - `counts` and `logcounts`
#' - `colData` with `cell_id`, `timepoint`, `condition`, `simulated_state`, and
#'   `simulated_branch`
#' - `reducedDims(sce)$PCA`
#' - `reducedDims(sce)$UMAP`
#' - `reducedDims(sce)$truth`
#'
#' @param n_cells Number of cells to simulate.
#' @param n_features Number of features to simulate.
#' @param n_pcs Number of principal components to store.
#'
#' @return A `SingleCellExperiment`.
#' @export
#'
#' @examples
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' dim(sce)
simulate_tata_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20
) {
    simulate_tata_divergent_sce(
        n_cells = n_cells,
        n_features = n_features,
        n_pcs = n_pcs
    )
}


#' Simulate an extreme spiral-plus-divergent TATA dataset
#'
#' Simulate a stress-test `SingleCellExperiment` whose `UMAP` embedding follows
#' a spiral-like manifold that transitions into two divergent tails. This is
#' useful for visually stress-testing TATA in an intentionally extreme geometry.
#'
#' The embedding is constructed directly in UMAP-like space and paired with
#' smooth expression programs so that neighboring cells along the path remain
#' transcriptionally similar.
#'
#' @inheritParams simulate_tata_sce
#'
#' @return A `SingleCellExperiment`.
#' @export
#'
#' @examples
#' sce <- simulate_tata_spiral_sce(n_cells = 1000, n_features = 80)
#' SingleCellExperiment::reducedDimNames(sce)
simulate_tata_spiral_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20
) {
    if (n_features < 60) {
        stop(
            "Please simulate at least 60 features for the spiral stress test.",
            call. = FALSE
        )
    }

    walk <- .tata_spiral_progress(n_cells)
    geom <- .tata_spiral_coordinates(walk$branch_label, walk$progress)
    labelled <- .tata_spiral_metadata(geom, n_cells)
    meta <- labelled$meta
    geom$radius <- geom$radius / max(geom$radius)
    sig <- .tata_spiral_signals(geom, labelled$condition)
    modules <- .tata_spiral_modules(n_features)
    assay_list <- .tata_simulate_counts(
        sig, modules$loads, n_cells, n_features, meta$cell_id,
        noise_sd = 0.07, eta_max = 4.2, size_scale = 6.5, nb_size = 10
    )
    branch_codes <- c(spiral = 0, tail_up = 1, tail_down = -1)

    .tata_assemble_sce(
        assay_list,
        meta,
        x = cbind(
            progress = geom$progress,
            radius = geom$radius,
            do.call(cbind, sig[c(
                "angle_sin", "angle_cos", "tail_up", "tail_down", "condition"
            )]),
            umap1 = scale(geom$umap1)[, 1],
            umap2 = scale(geom$umap2)[, 1]
        ),
        truth = cbind(
            latent_progress = geom$progress,
            branch_code = branch_codes[geom$branch_label]
        ),
        umap = cbind(geom$umap1, geom$umap2),
        n_pcs = n_pcs,
        simulation = list(
            trajectory = "spiral_divergent",
            feature_module = modules$feature_module
        )
    )
}


#' Simulate a multi-drug branching TATA dataset
#'
#' Simulate a `SingleCellExperiment` with one shared origin trajectory that can
#' diverge into three drug-associated fates. The three drug-associated branches
#' emerge at different real treatment times, making the dataset useful for
#' checking whether TATA uses time information rather than feature similarity
#' alone when orienting cluster-level transitions.
#'
#' The simulated design uses treatment times in hours:
#' - `drug_A` branch emerges early, around 24 hours
#' - `drug_C` branch emerges around 48 hours
#' - `drug_B` branch emerges late, around 96 hours
#'
#' The geometry is intentionally center-rooted rather than tip-rooted:
#' `drug_A` and `drug_B` extend in opposite directions from a shared middle
#' state, and `drug_C` diverges along a third arm. This is meant to stress-test
#' whether time labels can recover direction when the root is not visually
#' obvious from topology alone.
#'
#' @inheritParams simulate_tata_sce
#'
#' @return A `SingleCellExperiment`.
#' @export
#'
#' @examples
#' sce <- simulate_tata_multidrug_sce(n_cells = 1000, n_features = 80)
#' table(sce$condition, sce$timepoint)
simulate_tata_multidrug_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20
) {
    if (n_features < 60) {
        stop(
"Please simulate at least 60 features for the multi-drug branching example.",
            call. = FALSE
        )
    }

    condition_levels <- c("drug_A", "drug_B", "drug_C")
    time_levels <- c(0, 24, 48, 72, 96, 120)
    branch_onset <- c(drug_A = 24, drug_B = 96, drug_C = 48)

    meta <- .tata_multidrug_design(n_cells, condition_levels, time_levels)
    activation <- .tata_multidrug_activation(meta, branch_onset, time_levels)
    geom <- .tata_multidrug_coordinates(meta, activation)
    meta <- .tata_multidrug_metadata(meta, geom, condition_levels)
    sig <- .tata_multidrug_signals(meta, geom)
    modules <- .tata_multidrug_modules(n_features)
    assay_list <- .tata_simulate_counts(
        sig, modules$loads, nrow(meta), n_features, meta$cell_id,
        noise_sd = 0.08, eta_max = 4.2, size_scale = 6.8, nb_size = 10
    )
    branch_codes <- c(origin = 0, drug_A = 1, drug_B = 2, drug_C = 3)

    .tata_assemble_sce(
        assay_list,
        meta,
        x = cbind(
            progress = geom$progress,
            do.call(cbind, sig[c(
                "branch_A", "branch_B", "branch_C", "transition",
                "condition", "umap1", "umap2"
            )])
        ),
        truth = cbind(
            progress = geom$progress,
            branch_activation = geom$branch_activation,
            branch_code = branch_codes[as.character(meta$simulated_branch)]
        ),
        umap = cbind(geom$umap1, geom$umap2),
        n_pcs = n_pcs,
        simulation = list(
            trajectory = "multidrug_branching",
            branch_onset_hours = branch_onset,
            feature_module = modules$feature_module
        )
    )
}


#' @rdname simulate_tata_sce
#' @export
simulate_tata_divergent_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20
) {
    if (n_cells < 1000) {
        warning(
            "The divergent simulation is designed for larger datasets; ",
            "consider at least 1000 cells for clearer branch structure.",
            call. = FALSE
        )
    }

    if (n_features < 50) {
        stop(
            "Please simulate at least 50 features for this example.",
            call. = FALSE
        )
    }

    allocation <- .tata_divergent_allocation(n_cells)
    meta <- .tata_divergent_metadata(allocation)
    tau <- .tata_divergent_tau(meta)
    coords <- .tata_divergent_coordinates(tau$branch_sign, tau$tau_continuous)
    sig <- .tata_divergent_signals(tau, coords)
    modules <- .tata_divergent_modules(n_features)
    assay_list <- .tata_simulate_counts(
        sig, modules$loads, nrow(meta), n_features, meta$cell_id,
        noise_sd = 0.08, eta_max = 4.5, size_scale = 7, nb_size = 8
    )

    .tata_assemble_sce(
        assay_list,
        meta,
        x = t(assay_list$logcounts),
        truth = cbind(latent_x = coords$latent_x, latent_y = coords$latent_y),
        umap = cbind(
            coords$latent_x + 0.15 * sin(coords$latent_y),
            0.85 * coords$latent_y + 0.10 * sin(coords$latent_x / 2)
        ),
        n_pcs = n_pcs,
        simulation = list(
            allocation = allocation,
            feature_module = modules$feature_module,
            trajectory = "divergent"
        )
    )
}


# ---------------------------------------------------------------------------
# Internal helpers shared by the TATA simulators
# ---------------------------------------------------------------------------

#' Split a feature budget into named expression modules
#'
#' @keywords internal
#' @noRd
.tata_feature_modules <- function(module_names, proportions, n_features) {
    module_sizes <- round(proportions * n_features)
    module_sizes[1] <- module_sizes[1] + (n_features - sum(module_sizes))
    rep(module_names, times = module_sizes)
}


#' Draw the baseline expression level and the per-module loading vectors
#'
#' `ranges` is an ordered named list of `c(min, max)` uniform bounds. The
#' loadings are drawn in list order, so the random number stream follows the
#' order in which the modules are declared.
#'
#' @keywords internal
#' @noRd
.tata_module_loads <- function(
    feature_module,
    baseline_mean,
    baseline_sd,
    ranges
) {
    n_features <- length(feature_module)
    baseline <- stats::rnorm(
        n_features,
        mean = baseline_mean,
        sd = baseline_sd
    )
    loads <- lapply(names(ranges), function(module) {
        load <- rep(0, n_features)
        in_module <- feature_module == module
        load[in_module] <- stats::runif(
            sum(in_module),
            ranges[[module]][1],
            ranges[[module]][2]
        )
        load
    })
    names(loads) <- names(ranges)
    c(list(baseline = baseline), loads)
}


.tata_simulate_counts <- function(
    sig,
    loads,
    n_obs,
    n_features,
    cell_ids,
    noise_sd,
    eta_max,
    size_scale,
    nb_size
) {
    eta <- Reduce(
        function(acc, module) acc + outer(sig[[module]], loads[[module]]),
        names(sig),
        matrix(rep(loads$baseline, each = n_obs), nrow = n_obs)
    ) + matrix(
        stats::rnorm(n_obs * n_features, mean = 0, sd = noise_sd),
        nrow = n_obs
    )
    mu <- t(exp(pmax(pmin(eta, eta_max), -1.5)))
    size_factor <- stats::rgamma(n_obs, shape = 10, rate = 10)
    mu <- sweep(mu, 2, size_factor * size_scale, `*`)

    counts <- matrix(
        stats::rnbinom(n = length(mu), mu = as.vector(mu), size = nb_size),
        nrow = n_features,
        ncol = n_obs,
        dimnames = list(paste0("gene_", seq_len(n_features)), cell_ids)
    )
    library_size <- colSums(counts)
    norm_factor <- library_size / stats::median(library_size)

    list(
        counts = counts,
        logcounts = log1p(sweep(counts, 2, norm_factor, `/`))
    )
}


#' Assemble the simulated pieces into a SingleCellExperiment
#'
#' @keywords internal
#' @noRd
.tata_assemble_sce <- function(
    assay_list,
    meta,
    x,
    truth,
    umap,
    n_pcs,
    simulation
) {
    pca <- stats::prcomp(x, center = TRUE, scale. = TRUE)
    n_keep <- min(n_pcs, ncol(pca$x))
    pcs <- pca$x[, seq_len(n_keep), drop = FALSE]
    colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
    colnames(umap) <- c("UMAP1", "UMAP2")
    rownames(pcs) <- rownames(umap) <- rownames(truth) <- meta$cell_id

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = assay_list,
        colData = DataFrame(meta)
    )

    reducedDim(sce, "PCA") <- pcs
    reducedDim(sce, "UMAP") <- umap
    reducedDim(sce, "truth") <- truth
    rowData(sce)$feature_module <- simulation$feature_module

    metadata(sce)$simulation <- simulation

    sce
}


# ---------------------------------------------------------------------------
# Spiral simulation helpers
# ---------------------------------------------------------------------------

#' Assign each spiral cell to a branch and a position along the walk
#'
#' @keywords internal
#' @noRd
.tata_spiral_progress <- function(n_cells) {
    branch_label <- sample(
        x = c("spiral", "tail_up", "tail_down"),
        size = n_cells,
        replace = TRUE,
        prob = c(0.72, 0.14, 0.14)
    )

    progress <- numeric(n_cells)
    progress[branch_label == "spiral"] <- stats::rbeta(
        sum(branch_label == "spiral"),
        shape1 = 1.3,
        shape2 = 1.5
    ) * 0.82
    progress[branch_label == "tail_up"] <- 0.82 + stats::rbeta(
        sum(branch_label == "tail_up"),
        shape1 = 1.2,
        shape2 = 1.1
    ) * 0.18
    progress[branch_label == "tail_down"] <- 0.82 + stats::rbeta(
        sum(branch_label == "tail_down"),
        shape1 = 1.2,
        shape2 = 1.1
    ) * 0.18

    list(branch_label = branch_label, progress = progress)
}


#' Place the spiral walk into a UMAP-like plane with two divergent tails
#'
#' @keywords internal
#' @noRd
.tata_spiral_coordinates <- function(branch_label, progress) {
    angle <- 1.35 * pi + 3.7 * pi * pmin(progress, 0.82) / 0.82
    radius <- 1.2 + 9.2 * pmin(progress, 0.82) / 0.82

    umap1 <- radius * cos(angle)
    umap2 <- radius * sin(angle)

    end_angle <- 1.35 * pi + 3.7 * pi
    end_radius <- 1.2 + 9.2
    end_x <- end_radius * cos(end_angle)
    end_y <- end_radius * sin(end_angle)

    tail_progress <- pmax((progress - 0.82) / 0.18, 0)
    up_idx <- branch_label == "tail_up"
    down_idx <- branch_label == "tail_down"

    if (any(up_idx)) {
        umap1[up_idx] <- end_x + 5.6 * tail_progress[up_idx] +
            stats::rnorm(sum(up_idx), sd = 0.40)
        umap2[up_idx] <- end_y + 0.8 + 7.5 * tail_progress[up_idx] +
            stats::rnorm(sum(up_idx), sd = 0.42)
    }

    if (any(down_idx)) {
        umap1[down_idx] <- end_x + 5.2 * tail_progress[down_idx] +
            stats::rnorm(sum(down_idx), sd = 0.42)
        umap2[down_idx] <- end_y - 0.8 - 7.2 * tail_progress[down_idx] +
            stats::rnorm(sum(down_idx), sd = 0.42)
    }

    spiral_noise <- 0.25 + 0.15 * (1 - pmin(progress, 0.82) / 0.82)
    spiral_idx <- branch_label == "spiral"
    umap1[spiral_idx] <- umap1[spiral_idx] +
        stats::rnorm(sum(spiral_idx), sd = spiral_noise[spiral_idx])
    umap2[spiral_idx] <- umap2[spiral_idx] +
        stats::rnorm(sum(spiral_idx), sd = spiral_noise[spiral_idx])

    list(
        branch_label = branch_label,
        progress = progress,
        angle = angle,
        radius = radius,
        tail_progress = tail_progress,
        umap1 = umap1,
        umap2 = umap2
    )
}


#' Bin the spiral walk into discrete timepoints and named states
#'
#' @keywords internal
#' @noRd
.tata_spiral_state_labels <- function(progress, branch_label) {
    timepoint <- cut(
        progress,
        breaks = c(-Inf, 0.16, 0.34, 0.54, 0.74, Inf),
        labels = c(0, 1, 2, 3, 4),
        right = TRUE
    )
    timepoint <- as.integer(as.character(timepoint))

    state_band <- cut(
        progress,
        breaks = c(-Inf, 0.12, 0.28, 0.46, 0.66, 0.82, 0.90, Inf),
        labels = c(
            "core", "inner", "mid", "outer", "pre_tail", "tail_early",
            "tail_late"
        ),
        right = TRUE
    )

    simulated_state <- ifelse(
        branch_label == "spiral",
        paste0("spiral_", state_band),
        ifelse(
            branch_label == "tail_up",
            paste0("tail_up_", state_band),
            paste0("tail_down_", state_band)
        )
    )

    list(timepoint = timepoint, simulated_state = simulated_state)
}


#' Derive the spiral cell annotation table from the latent geometry
#'
#' @keywords internal
#' @noRd
.tata_spiral_metadata <- function(geom, n_cells) {
    branch_label <- geom$branch_label
    labels <- .tata_spiral_state_labels(geom$progress, branch_label)

    condition <- ifelse(
        branch_label == "tail_down",
        ifelse(stats::runif(n_cells) < 0.85, "drug", "control"),
        ifelse(
            branch_label == "tail_up",
            ifelse(stats::runif(n_cells) < 0.80, "control", "drug"),
            ifelse(stats::runif(n_cells) < 0.55, "control", "drug")
        )
    )

    meta <- data.frame(
        cell_id = paste0("bacterium_", seq_len(n_cells)),
        timepoint = labels$timepoint,
        condition = factor(condition, levels = c("control", "drug")),
        simulated_state = factor(labels$simulated_state),
        simulated_branch = factor(
            branch_label,
            levels = c("spiral", "tail_up", "tail_down")
        ),
        progress = geom$progress,
        stringsAsFactors = FALSE
    )
    rownames(meta) <- meta$cell_id

    list(meta = meta, condition = condition)
}


#' Smooth expression programs driving the spiral simulation
#'
#' @keywords internal
#' @noRd
.tata_spiral_signals <- function(geom, condition) {
    progress <- geom$progress
    tail_down <- as.numeric(geom$branch_label == "tail_down") *
        geom$tail_progress

    list(
        root = exp(-((progress - 0.06)^2) / (2 * 0.06^2)),
        inner = exp(-((progress - 0.20)^2) / (2 * 0.08^2)),
        mid = exp(-((progress - 0.42)^2) / (2 * 0.09^2)),
        outer = exp(-((progress - 0.66)^2) / (2 * 0.08^2)) +
            pmax((progress - 0.78) / 0.22, 0),
        progress = progress + geom$radius,
        angle_sin = sin(geom$angle),
        angle_cos = cos(geom$angle),
        tail_up = as.numeric(geom$branch_label == "tail_up") *
            geom$tail_progress,
        tail_down = tail_down,
        condition = 0.20 * as.numeric(condition == "drug") +
            0.22 * tail_down
    )
}


#' Declare the spiral feature modules and draw their loadings
#'
#' @keywords internal
#' @noRd
.tata_spiral_modules <- function(n_features) {
    module_names <- c(
        "root", "inner", "mid", "outer", "progress", "angle_sin",
        "angle_cos", "tail_up", "tail_down", "condition", "noise"
    )
    feature_module <- .tata_feature_modules(
        module_names,
        c(0.10, 0.10, 0.12, 0.10, 0.16, 0.08, 0.08, 0.08, 0.08, 0.06, 0.04),
        n_features
    )

    loads <- .tata_module_loads(
        feature_module,
        baseline_mean = 0.78,
        baseline_sd = 0.12,
        ranges = list(
            root = c(0.70, 0.95),
            inner = c(0.60, 0.85),
            mid = c(0.60, 0.85),
            outer = c(0.60, 0.85),
            progress = c(0.60, 0.90),
            angle_sin = c(0.35, 0.55),
            angle_cos = c(0.35, 0.55),
            tail_up = c(0.55, 0.85),
            tail_down = c(0.55, 0.85),
            condition = c(0.18, 0.32)
        )
    )

    list(feature_module = feature_module, loads = loads)
}


# ---------------------------------------------------------------------------
# Multi-drug simulation helpers
# ---------------------------------------------------------------------------

#' Allocate cells across the multi-drug condition-by-time design
#'
#' @keywords internal
#' @noRd
.tata_multidrug_design <- function(n_cells, condition_levels, time_levels) {
    design <- expand.grid(
        condition = condition_levels,
        timepoint = time_levels,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    design$weight <- c(
        420, 520, 480, 430, 380, 320,
        420, 520, 500, 470, 430, 360,
        420, 520, 500, 430, 360, 320
    )
    design$n <- round(n_cells * design$weight / sum(design$weight))
    delta <- n_cells - sum(design$n)
    if (delta != 0L) {
        adjust_idx <- seq_len(abs(delta))
        design$n[adjust_idx] <- design$n[adjust_idx] + sign(delta)
    }

    meta <- design[
        rep(seq_len(nrow(design)), design$n),
        c("condition", "timepoint"),
        drop = FALSE
    ]
    rownames(meta) <- NULL
    meta <- meta[sample(seq_len(nrow(meta))), , drop = FALSE]

    meta$cell_id <- paste0("bacterium_", seq_len(nrow(meta)))
    meta
}


#' Draw the continuous progression and per-condition branch activation
#'
#' @keywords internal
#' @noRd
.tata_multidrug_activation <- function(meta, branch_onset, time_levels) {
    n_obs <- nrow(meta)
    time_scaled <- meta$timepoint / max(time_levels)
    onset_scaled <- branch_onset[meta$condition] / max(time_levels)

    # Continuous progression within each labeled timepoint creates many
    # transitional cells instead of isolated condition-time islands.
    progress <- pmin(
        1,
        pmax(
            0,
            time_scaled + stats::rnorm(n_obs, mean = 0, sd = 0.07)
        )
    )

    branch_activation <- pmax(
        0,
        progress - onset_scaled + stats::rnorm(n_obs, mean = 0, sd = 0.020)
    )
    branch_activation <- pmin(0.48, branch_activation) / 0.48
    branch_activation <- pmin(1, pmax(0, branch_activation))

    list(
        progress = progress,
        branch_activation = branch_activation,
        onset_scaled = onset_scaled
    )
}


#' Place the shared trunk and the three drug arms in a UMAP-like plane
#'
#' @keywords internal
#' @noRd
.tata_multidrug_coordinates <- function(meta, activation) {
    n_obs <- nrow(meta)
    progress <- activation$progress
    branch_activation <- activation$branch_activation
    onset_scaled <- activation$onset_scaled

    branch_direction <- rbind(
        drug_A = c(-1.18, 0.12),
        drug_B = c(1.18, -0.10),
        drug_C = c(0.05, 1.22)
    )

    dir_x <- branch_direction[meta$condition, 1]
    dir_y <- branch_direction[meta$condition, 2]

    central_radius <- 0.25 + 1.05 * pmin(progress, 0.35) / 0.35
    central_angle <- 1.6 * pi * progress +
        stats::rnorm(n_obs, mean = 0, sd = 0.08)
    trunk_x <- central_radius * cos(central_angle) * (1 - branch_activation)
    trunk_y <- -0.8 + 1.8 * progress +
        0.7 * sin(1.2 * pi * progress) * (1 - 0.4 * branch_activation)

    branch_length <- 7.0 * branch_activation^1.10
    branch_bend <- 1.2 * branch_activation * (progress - onset_scaled)
    bend_sign <- c(drug_A = 1, drug_B = -1, drug_C = 1)[meta$condition]

    umap1 <- trunk_x +
        branch_length * dir_x +
        0.9 * branch_bend * ifelse(meta$condition == "drug_C", 0.35, 1) +
        stats::rnorm(n_obs, sd = 0.18 + 0.16 * (1 - branch_activation))

    umap2 <- trunk_y +
        branch_length * dir_y +
        0.8 * bend_sign * branch_bend *
            ifelse(meta$condition == "drug_C", 0.15, 0.55) +
        stats::rnorm(n_obs, sd = 0.18 + 0.16 * (1 - branch_activation))

    list(
        progress = progress,
        branch_activation = branch_activation,
        umap1 = umap1,
        umap2 = umap2
    )
}


#' Attach fate labels and factor levels to the multi-drug annotation table
#'
#' @keywords internal
#' @noRd
.tata_multidrug_metadata <- function(meta, geom, condition_levels) {
    branch_activation <- geom$branch_activation

    fate_label <- ifelse(
        branch_activation < 0.10,
        "origin",
        ifelse(branch_activation < 0.55, "transition", "fate")
    )
    simulated_state <- paste(
        meta$condition,
        fate_label,
        paste0("h", meta$timepoint),
        sep = "_"
    )

    meta$condition <- factor(meta$condition, levels = condition_levels)
    meta$timepoint <- as.numeric(meta$timepoint)
    meta$simulated_state <- factor(simulated_state)
    meta$simulated_branch <- factor(
        ifelse(
            branch_activation < 0.10, "origin", as.character(meta$condition)
        ),
        levels = c("origin", condition_levels)
    )
    meta$branch_activation <- branch_activation
    meta$progress <- geom$progress
    meta <- meta[, c(
        "cell_id", "timepoint", "condition", "simulated_state",
        "simulated_branch", "progress", "branch_activation"
    )]
    rownames(meta) <- meta$cell_id
    meta
}


#' Smooth expression programs driving the multi-drug simulation
#'
#' @keywords internal
#' @noRd
.tata_multidrug_signals <- function(meta, geom) {
    progress <- geom$progress
    branch_activation <- geom$branch_activation

    # Condition-specific offsets should only appear once a branch is activated;
    # otherwise the early shared origin would split by treatment too soon.
    condition_signal <- branch_activation * (
        0.12 * as.numeric(meta$condition == "drug_B") +
            0.08 * as.numeric(meta$condition == "drug_C")
    )

    list(
        root = exp(-((progress - 0.08)^2) / (2 * 0.06^2)),
        early = exp(-((progress - 0.24)^2) / (2 * 0.08^2)),
        mid = exp(-((progress - 0.50)^2) / (2 * 0.11^2)),
        late = exp(-((progress - 0.78)^2) / (2 * 0.12^2)),
        center = exp(-((progress - 0.14)^2) / (2 * 0.10^2)) *
            (1 - branch_activation),
        transition = pmax(0, 1 - abs(branch_activation - 0.28) / 0.28),
        branch_A = branch_activation *
            as.numeric(meta$condition == "drug_A"),
        branch_B = branch_activation *
            as.numeric(meta$condition == "drug_B"),
        branch_C = branch_activation *
            as.numeric(meta$condition == "drug_C"),
        condition = condition_signal,
        umap1 = scale(geom$umap1)[, 1],
        umap2 = scale(geom$umap2)[, 1]
    )
}


#' Declare the multi-drug feature modules and draw their loadings
#'
#' @keywords internal
#' @noRd
.tata_multidrug_modules <- function(n_features) {
    module_names <- c(
        "root", "early", "mid", "late", "center", "transition",
        "branch_A", "branch_B", "branch_C", "condition", "umap1",
        "umap2", "noise"
    )
    feature_module <- .tata_feature_modules(
        module_names,
        c(
            0.09, 0.09, 0.10, 0.10, 0.08, 0.10, 0.09, 0.09, 0.09, 0.07,
            0.06, 0.07, 0.03
        ),
        n_features
    )

    loads <- .tata_module_loads(
        feature_module,
        baseline_mean = 0.82,
        baseline_sd = 0.13,
        ranges = list(
            root = c(0.75, 1.00),
            early = c(0.60, 0.88),
            mid = c(0.60, 0.88),
            late = c(0.60, 0.88),
            center = c(0.40, 0.65),
            transition = c(0.45, 0.70),
            branch_A = c(0.45, 0.72),
            branch_B = c(0.45, 0.72),
            branch_C = c(0.45, 0.72),
            condition = c(0.16, 0.30),
            umap1 = c(0.18, 0.35),
            umap2 = c(0.18, 0.35)
        )
    )

    list(feature_module = feature_module, loads = loads)
}


# ---------------------------------------------------------------------------
# Divergent simulation helpers
# ---------------------------------------------------------------------------

#' Allocate cells across the divergent timepoint-by-branch design
#'
#' @keywords internal
#' @noRd
.tata_divergent_allocation <- function(n_cells) {
    allocation <- data.frame(
        timepoint = c(0, 1, 1, 2, 2, 3, 3, 4, 4),
        simulated_branch = c(
            "root", "adaptive", "collapse", "adaptive", "collapse",
            "adaptive", "collapse", "adaptive", "collapse"
        ),
        weight = c(900, 700, 700, 600, 600, 500, 500, 250, 250),
        stringsAsFactors = FALSE
    )

    allocation$n <- round(n_cells * allocation$weight / sum(allocation$weight))
    delta <- n_cells - sum(allocation$n)
    if (delta != 0L) {
        adjust_idx <- seq_len(abs(delta))
        allocation$n[adjust_idx] <- allocation$n[adjust_idx] + sign(delta)
    }
    allocation
}


#' Expand the divergent allocation into a labelled cell annotation table
#'
#' @keywords internal
#' @noRd
.tata_divergent_metadata <- function(allocation) {
    meta <- allocation[
        rep(seq_len(nrow(allocation)), allocation$n),
        c("timepoint", "simulated_branch"),
        drop = FALSE
    ]
    meta <- meta[sample(seq_len(nrow(meta))), , drop = FALSE]
    rownames(meta) <- NULL

    control_probability <- ifelse(
        meta$simulated_branch == "adaptive",
        0.75,
        ifelse(meta$simulated_branch == "collapse", 0.25, 0.50)
    )

    meta$condition <- ifelse(
        stats::runif(nrow(meta)) < control_probability,
        "control",
        "drug"
    )
    meta$simulated_state <- ifelse(
        meta$simulated_branch == "root",
        "root",
        paste(meta$simulated_branch, paste0("t", meta$timepoint), sep = "_")
    )
    meta$cell_id <- paste0("bacterium_", seq_len(nrow(meta)))

    state_levels <- c(
        "root",
        "adaptive_t1", "adaptive_t2", "adaptive_t3", "adaptive_t4",
        "collapse_t1", "collapse_t2", "collapse_t3", "collapse_t4"
    )

    meta$condition <- factor(meta$condition, levels = c("control", "drug"))
    meta$simulated_branch <- factor(
        meta$simulated_branch,
        levels = c("root", "adaptive", "collapse")
    )
    meta$simulated_state <- factor(meta$simulated_state, levels = state_levels)
    meta <- meta[, c(
        "cell_id", "timepoint", "condition", "simulated_state",
        "simulated_branch"
    )]
    rownames(meta) <- meta$cell_id
    meta
}


#' Spread the labelled timepoints onto a continuous latent progression axis
#'
#' @keywords internal
#' @noRd
.tata_divergent_tau <- function(meta) {
    n_obs <- nrow(meta)

    branch_sign <- ifelse(
        meta$simulated_branch == "adaptive",
        -1,
        ifelse(meta$simulated_branch == "collapse", 1, 0)
    )
    is_drug <- as.numeric(meta$condition == "drug")

    # Each experimental timepoint occupies a broad, overlapping interval on a
    # continuous latent progression axis. This creates many intermediate cells
    # between labelled timepoints, which helps UMAP recover a trajectory rather
    # than disconnected blobs.
    tau_min <- c(`0` = 0.00, `1` = 0.08, `2` = 0.24, `3` = 0.44, `4` = 0.64)
    tau_max <- c(`0` = 0.18, `1` = 0.34, `2` = 0.58, `3` = 0.80, `4` = 1.00)
    tau_continuous <- stats::runif(
        n_obs,
        min = tau_min[as.character(meta$timepoint)],
        max = tau_max[as.character(meta$timepoint)]
    )
    tau_continuous <- pmin(
        1,
        pmax(
            0,
            tau_continuous + stats::rnorm(n_obs, mean = 0, sd = 0.025)
        )
    )

    list(
        branch_sign = branch_sign,
        is_drug = is_drug,
        tau_continuous = tau_continuous
    )
}


#' Turn branch commitment into the divergent latent coordinates
#'
#' @keywords internal
#' @noRd
.tata_divergent_coordinates <- function(branch_sign, tau_continuous) {
    n_obs <- length(tau_continuous)

    # Branch commitment turns on gradually after the root. Near the split, cells
    # from the two future branches still lie close together, generating visible
    # transition states in PCA/UMAP space.
    branch_strength <- ifelse(
        branch_sign == 0,
        0,
        stats::plogis((tau_continuous - 0.34) / 0.09) +
            stats::rnorm(n_obs, mean = 0, sd = 0.035)
    )
    branch_strength <- pmin(1, pmax(0, branch_strength))
    branch_effect <- branch_sign * branch_strength

    latent_x <- 4.20 * branch_effect * (0.55 + 0.55 * tau_continuous) +
        stats::rnorm(
            n_obs,
            mean = 0,
            sd = 0.11 + 0.05 * (1 - branch_strength)
        )

    latent_y <- 5.00 * tau_continuous - 1.15 +
        0.45 * tau_continuous^2 -
        0.15 * branch_strength^2 +
        stats::rnorm(n_obs, mean = 0, sd = 0.10)

    list(
        branch_effect = branch_effect,
        latent_x = latent_x,
        latent_y = latent_y
    )
}


#' Smooth expression programs driving the divergent simulation
#'
#' @keywords internal
#' @noRd
.tata_divergent_signals <- function(tau, coords) {
    tau_continuous <- tau$tau_continuous
    is_drug <- tau$is_drug
    branch_signal <- coords$branch_effect
    collapse_signal <- pmax(branch_signal, 0)

    list(
        root = exp(-((tau_continuous - 0.05)^2) / (2 * 0.08^2)),
        early = exp(-((tau_continuous - 0.22)^2) / (2 * 0.10^2)),
        mid = exp(-((tau_continuous - 0.48)^2) / (2 * 0.14^2)),
        late = exp(-((tau_continuous - 0.78)^2) / (2 * 0.14^2)),
        branch = branch_signal,
        adaptive = pmax(-branch_signal, 0),
        collapse = collapse_signal,
        stress = 0.25 * is_drug + 0.30 * collapse_signal +
            0.10 * tau_continuous +
            (0.18 * is_drug + 0.08 * tau_continuous * is_drug),
        latent_y = coords$latent_y
    )
}


#' Declare the divergent feature modules and draw their loadings
#'
#' The `branch` module uses a fixed alternating +/- loading rather than a
#' random one, so it is filled in after the uniform draws and consumes no
#' random numbers.
#'
#' @keywords internal
#' @noRd
.tata_divergent_modules <- function(n_features) {
    module_names <- c(
        "root", "early", "mid", "late", "branch", "adaptive",
        "collapse", "stress", "latent_y", "noise"
    )
    feature_module <- .tata_feature_modules(
        module_names,
        c(0.14, 0.14, 0.16, 0.12, 0.12, 0.08, 0.08, 0.07, 0.05, 0.04),
        n_features
    )

    loads <- .tata_module_loads(
        feature_module,
        baseline_mean = 0.82,
        baseline_sd = 0.14,
        ranges = list(
            root = c(0.75, 1.05),
            early = c(0.65, 0.95),
            mid = c(0.65, 0.95),
            late = c(0.65, 0.95),
            adaptive = c(0.30, 0.55),
            collapse = c(0.30, 0.55),
            stress = c(0.18, 0.35),
            latent_y = c(0.22, 0.40)
        )
    )

    loads$branch <- rep(0, n_features)
    loads$branch[feature_module == "branch"] <- rep(
        c(-0.50, 0.50), length.out = sum(feature_module == "branch")
    )

    list(feature_module = feature_module, loads = loads)
}
