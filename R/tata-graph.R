#' Build a cell-cell kNN graph
#'
#' Build an undirected cell-cell k-nearest-neighbour graph from a
#' low-dimensional representation.
#'
#' @param x A `SingleCellExperiment` or a numeric matrix with cells in rows.
#' @param dimred Name of the reduced dimension to use when `x` is a
#'   `SingleCellExperiment`.
#' @param dims Optional subset of dimensions to use.
#' @param k Number of neighbours.
#'
#' @return A list containing an `igraph` object, a sparse adjacency matrix, the
#'   embedding used, neighbour indices, neighbour distances, and kNN metadata.
#' @export
#'
#' @examples
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' knn <- build_knn_graph(sce, dimred = "PCA", k = 10)
#' igraph::vcount(knn$graph)
build_knn_graph <- function(
    x,
    dimred = "PCA",
    dims = NULL,
    k = 15
) {
    embedding <- .get_embedding_matrix(x, dimred = dimred, dims = dims)
    k <- max(1, min(as.integer(k), nrow(embedding) - 1L))
    knn <- .compute_knn(embedding, k = k)

    edge_i <- rep(seq_len(nrow(embedding)), each = k)
    edge_j <- as.vector(t(knn$nn.index))

    adjacency <- Matrix::sparseMatrix(
        i = edge_i,
        j = edge_j,
        x = 1,
        dims = c(nrow(embedding), nrow(embedding)),
        dimnames = list(rownames(embedding), rownames(embedding))
    )

    ## Symmetrise, then binarise. Summing the matrix with its transpose gives a
    ## weight of 2 to mutual neighbours, so collapse every non-zero entry to 1
    ## rather than reaching into the Matrix object's slots.
    adjacency <- adjacency + Matrix::t(adjacency)
    adjacency <- as(adjacency > 0, "dMatrix")
    diag(adjacency) <- 0
    adjacency <- Matrix::drop0(adjacency)

    cell_graph <- igraph::graph_from_adjacency_matrix(
        adjacency,
        mode = "undirected",
        diag = FALSE
    )

    list(
        graph = cell_graph,
        adjacency = adjacency,
        embedding = embedding,
        knn_index = knn$nn.index,
        knn_distance = knn$nn.dist,
        k = k,
        knn_engine = knn$engine
    )
}


#' Cluster graph-connected cell states
#'
#' Cluster the cell graph into coarse states using a graph-based community
#' detection method.
#'
#' @param cell_graph An undirected cell-level `igraph`.
#' @param method Clustering method. Choices are:
#' - `"louvain"`: fast modularity-based community detection.
#' - `"leiden"`: Leiden clustering if supported by the installed `igraph`.
#' - `"walktrap"`: random-walk clustering that can produce broader partitions.
#'
#' @return A factor of inferred cluster labels.
#' @export
#'
#' @examples
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' knn <- build_knn_graph(sce, dimred = "PCA", k = 10)
#' clusters <- cluster_graph_states(knn$graph, method = "louvain")
#' table(clusters)
cluster_graph_states <- function(
    cell_graph,
    method = c("louvain", "leiden", "walktrap")
) {
    method <- match.arg(method)

    if (!inherits(cell_graph, "igraph")) {
        stop("`cell_graph` must be an igraph object.", call. = FALSE)
    }

    if (igraph::ecount(cell_graph) == 0L) {
        membership <- rep(1L, igraph::vcount(cell_graph))
    } else if (method == "leiden") {
        if ("cluster_leiden" %in% getNamespaceExports("igraph")) {
            ## `cluster_leiden()` defaults to the CPM objective, which on an
            ## unweighted kNN graph places every cell in its own community.
            ## Optimise modularity instead, so that Leiden is comparable with
            ## the Louvain and walktrap branches below.
            membership <- igraph::membership(igraph::cluster_leiden(
                cell_graph,
                objective_function = "modularity"
            ))
        } else {
            warning(
                "`cluster_leiden()` is not available; falling back to Louvain.",
                call. = FALSE
            )
            membership <- igraph::membership(
                igraph::cluster_louvain(cell_graph)
            )
        }
    } else if (method == "walktrap") {
        membership <- igraph::membership(igraph::cluster_walktrap(cell_graph))
    } else {
        membership <- igraph::membership(igraph::cluster_louvain(cell_graph))
    }

    membership_levels <- sort(unique(membership))
    cluster_labels <- paste0("C", match(membership, membership_levels))
    factor(cluster_labels, levels = paste0("C", seq_along(membership_levels)))
}
