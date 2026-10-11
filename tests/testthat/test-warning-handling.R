test_that("categorical and numeric text labels keep their original ordering", {
    labels <- c("late", "early", "late")
    expect_warning(reference <- .make_default_reference(labels), NA)
    expect_identical(reference, "early")
    expect_warning(groups <- .build_group_factor(
        data.frame(time = labels), "time"
    ), NA)
    expect_identical(levels(groups), c("early", "late"))
    expect_warning(times <- .coerce_time_to_numeric(labels), NA)
    expect_identical(times, c(2L, 1L, 2L))

    ordered_labels <- factor(labels, levels = c("late", "early"))
    expect_warning(times <- .coerce_time_to_numeric(ordered_labels), NA)
    expect_identical(times, c(1, 2, 1))

    numeric_labels <- c("10", "2", "0x3", "1e0")
    expect_warning(times <- .coerce_time_to_numeric(numeric_labels), NA)
    expect_identical(times, c(10, 2, 3, 1))
    expect_identical(.make_default_reference(numeric_labels), "1e0")
    groups <- .build_group_factor(data.frame(time = numeric_labels), "time")
    expect_identical(levels(groups), c("1e0", "2", "0x3", "10"))

    ## Logical and complex-looking text are categorical labels, rather than
    ## numeric time values; missing annotations remain missing.
    expect_identical(.coerce_time_to_numeric(c("TRUE", "FALSE")), c(2L, 1L))
    expect_identical(.coerce_time_to_numeric(c("2i", "1i")), c(2L, 1L))
    expect_warning(times <- .coerce_time_to_numeric(c("late", NA, "early")),
        NA)
    expect_identical(times, c(2L, NA_integer_, 1L))
    expect_identical(.coerce_time_to_numeric(
        factor(c("10", NA, "2"), levels = c("2", "10"))
    ), c(2, NA_real_, 1))
    expect_identical(.coerce_time_to_numeric(character()), numeric())
    expect_identical(.make_default_reference(character()), NA_character_)
    expect_identical(levels(.build_group_factor(
        data.frame(time = character()), "time"
    )), character())
})

test_that("undefined time correlations preserve NA metrics without warnings", {
    expect_warning(constant <- orient_pseudotime_to_time(
        rep(4, 3), c(0, 1, 2)
    ), NA)
    expect_identical(constant, c(0, 0, 0))
    expect_warning(uninformative <- orient_pseudotime_to_time(
        c(2, 3, 4), rep("same", 3)
    ), NA)
    expect_identical(uninformative, c(0, 0.5, 1))
    expect_identical(orient_pseudotime_to_time(c(1, 2), c(0, 1)),
        rep(NA_real_, 2))

    for (constant_truth in c(TRUE, FALSE)) {
        inferred <- if (constant_truth) c(0.2, 0.8) else c(0.5, 0.5)
        truth <- if (constant_truth) c(A = 1, B = 1) else c(A = 1, B = 2)
        expect_warning(onsets <- branch_onset_metrics(
            pseudotime = inferred,
            branch = c("A", "B"),
            branch_activation = c(1, 1),
            true_onset = truth
        ), NA)
        expect_true(all(is.na(onsets$summary[
            c("kendall_tau", "spearman_rho", "branch_onset_order_score")
        ])))
        expect_equal(onsets$summary$n_branches_used, 2)
    }
})

test_that("weighted routes prefer direction and handle disconnected targets", {
    skip_if_not_installed("igraph")
    graph <- igraph::graph_from_data_frame(
        data.frame(
            from = c("root", "middle", "back", "reverse"),
            to = c("middle", "back", "root", "root"),
            weight = c(1, 1, 10, 1)
        ),
        directed = TRUE,
        vertices = c("root", "middle", "back", "reverse", "orphan")
    )
    expect_warning(paths <- .tata_paths_from_root(
        graph, "root", c("back", "reverse", "orphan")
    ), NA)
    expect_identical(paths, list(
        c("root", "middle", "back"), c("root", "reverse"), character()
    ))

    undirected <- igraph::as_undirected(graph, mode = "collapse")
    expect_warning(path <- .tata_trajectory_path_vertices(
        graph, undirected, "root", "reverse"
    ), NA)
    expect_identical(igraph::V(graph)$name[as.integer(path)],
        c("root", "reverse"))
    expect_warning(path <- .tata_trajectory_path_vertices(
        graph, undirected, "root", "orphan"
    ), NA)
    expect_length(path, 0)
    expect_error(.tata_paths_from_root(graph, "root", "unknown"))
})

test_that("ill-conditioned expression fits report their numerical warnings", {
    fit_df <- data.frame(
        pseudotime = rep(c(0, 0.25, 0.5, 0.75, 1), each = 4),
        expression = seq_len(20),
        weight = 1
    )
    messages <- character()
    curve <- withCallingHandlers(
        .tata_fit_trend_curve(fit_df, seq(0, 1, length.out = 10), span = 0.1),
        warning = function(w) {
            messages <<- c(messages, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_true(any(grepl("span too small|singular|pseudoinverse", messages)))
    expect_length(curve, 10)
    expect_true(all(is.finite(curve)))

    previous <- options(warn = 2)
    on.exit(options(previous))
    expect_error(
        .tata_fit_trend_curve(fit_df, seq(0, 1, length.out = 10), span = 0.1),
        "converted from warning.*span too small"
    )
})
