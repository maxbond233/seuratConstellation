#' @title Assign Cluster Groups
#' @description Assign metadata-based groups to clusters for coloring or hull grouping
#' @name assign_cluster_groups
NULL

#' Assign Cluster Groups
#'
#' @description Maps clusters to metadata groups based on majority voting or other
#'   strategies. This function separates the grouping logic from the graph building
#'   and color resolution steps, making dependencies explicit.
#'
#' @param knn_graph A constellation_knn object from build_knn_graph()
#' @param seu A Seurat object (must be the same one used in build_knn_graph)
#' @param group_by Character. Column name in meta.data for grouping
#' @param strategy Character. Grouping strategy: "majority" (default), "strict", or "manual"
#' @param group_type Character. Type of grouping: "color" or "hull"
#' @param manual_mapping Optional. Named vector for manual cluster-to-group mapping
#'
#' @return A constellation_knn object with group assignments added to cluster_info
#'
#' @details
#' Strategies:
#' \itemize{
#'   \item \strong{majority}: Assign each cluster to the most common group among its cells
#'   \item \strong{strict}: Only assign if >90\% of cells belong to one group, else NA
#'   \item \strong{manual}: Use the provided manual_mapping (named vector)
#' }
#'
#' Group types:
#' \itemize{
#'   \item \strong{color}: Adds "color_group" column for color resolution
#'   \item \strong{hull}: Adds "hull_group" column for hull computation
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' knn_graph <- build_knn_graph(seu, cluster_col = "celltype")
#' knn_graph <- assign_cluster_groups(knn_graph, seu,
#'                                    group_by = "cell_class",
#'                                    group_type = "color")
#' knn_graph <- assign_cluster_groups(knn_graph, seu,
#'                                    group_by = "tissue",
#'                                    group_type = "hull")
#' }
assign_cluster_groups <- function(knn_graph,
                                  seu,
                                  group_by,
                                  strategy = "majority",
                                  group_type = "color",
                                  manual_mapping = NULL) {

  if (!inherits(knn_graph, "constellation_knn")) {
    stop("Input must be a constellation_knn object from build_knn_graph()")
  }

  if (!inherits(seu, "Seurat")) {
    stop("seu must be a Seurat object")
  }

  if (!group_by %in% colnames(seu@meta.data)) {
    stop(paste("Column", group_by, "not found in meta.data"))
  }

  if (!group_type %in% c("color", "hull")) {
    stop("group_type must be 'color' or 'hull'")
  }

  if (!strategy %in% c("majority", "strict", "manual")) {
    stop("strategy must be 'majority', 'strict', or 'manual'")
  }

  # Extract data
  cluster_col <- knn_graph$params$cluster_col
  cluster_vec <- seu@meta.data[[cluster_col]]
  group_vec <- seu@meta.data[[group_by]]
  cluster_levels <- knn_graph$cluster_info$cluster_label

  # Compute group assignment based on strategy
  if (strategy == "manual") {
    if (is.null(manual_mapping)) {
      stop("manual_mapping must be provided when strategy='manual'")
    }
    cluster_groups <- manual_mapping[as.character(cluster_levels)]
  } else if (strategy == "majority") {
    cluster_groups <- .compute_majority_category(cluster_vec, group_vec, cluster_levels)
  } else if (strategy == "strict") {
    cluster_groups <- .compute_strict_category(cluster_vec, group_vec, cluster_levels, threshold = 0.9)
  }

  # Add to cluster_info based on group_type
  col_name <- if (group_type == "color") "color_group" else "hull_group"
  knn_graph$cluster_info[[col_name]] <- cluster_groups

  # If centers already computed, keep cl.center.df in sync
  if (!is.null(knn_graph$cl.center.df)) {
    center_labels <- knn_graph$cl.center.df$cluster_label
    map_idx <- match(center_labels, knn_graph$cluster_info$cluster_label)
    knn_graph$cl.center.df[[col_name]] <- knn_graph$cluster_info[[col_name]][map_idx]
  }

  # Store parameter for reference
  if (is.null(knn_graph$params$group_by)) {
    knn_graph$params$group_by <- list()
  }
  knn_graph$params$group_by[[group_type]] <- group_by

  return(knn_graph)
}

#' Compute strict category assignment
#' @keywords internal
.compute_strict_category <- function(cluster_vec, category_vec, cluster_levels, threshold = 0.9) {
  strict <- sapply(cluster_levels, function(cl) {
    idx <- cluster_vec == cl
    if (sum(idx) == 0) return(NA)
    tab <- table(category_vec[idx])
    max_frac <- max(tab) / sum(tab)
    if (max_frac >= threshold) {
      names(tab)[which.max(tab)]
    } else {
      NA
    }
  })
  names(strict) <- cluster_levels
  strict
}
