#' @title ConstellationPlot - One-stop Convenience Function
#' @description All-in-one function for creating constellation plots from Seurat objects
#' @name ConstellationPlot
NULL

#' ConstellationPlot
#'
#' @description A convenience wrapper that combines all steps of constellation
#'   plot creation into a single function call. This is the recommended entry
#'   point for most users. The function now uses a clearer pipeline with
#'   explicit color resolution.
#'
#' @param seu A Seurat object (V5 compatible)
#' @param cluster_col Character. Column name in meta.data for cluster identity
#' @param reduction Character. Name of dimensional reduction. Default "umap"
#' @param k Integer. Number of nearest neighbors. Default 15
#' @param frac_th Numeric. Minimum fraction threshold for edges. Default 0.05
#' @param colors Optional. Named vector of colors by cluster_label, or NULL for auto-generation
#' @param color_by Character. Column name in meta.data for node coloring. Default NULL
#' @param color_mode Character. Color mode: "group" or "gradient". Default "group"
#' @param hull_by Character. Column name in meta.data for hull grouping. Default NULL
#' @param hull_type Character. Hull type: "none", "convex", or "concave". Default "none"
#' @param hull_alpha Numeric. Hull transparency (0-1). Default 0.2
#' @param hull_expand Numeric. Expansion factor for hulls. Default 0.1
#' @param node.label Character. Column for node labels. Default "cluster_label"
#' @param exxageration Numeric. Edge width exaggeration. Default 2
#' @param curved Logical. Curved edges. Default TRUE
#' @param node.dodge Logical. Dodge overlapping nodes. Default FALSE
#' @param label.size Numeric. Label size. Default 5
#' @param max_size Numeric. Maximum node size. Default 10
#' @param label_repel Logical. Use ggrepel. Default FALSE
#' @param node_trans Character. Size transformation. Default "sqrt"
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' ConstellationPlot(seu, cluster_col = "celltype")
#'
#' # With group coloring and hulls
#' ConstellationPlot(
#'   seu,
#'   cluster_col = "celltype",
#'   color_by = "cell_class",
#'   hull_by = "cell_class",
#'   hull_type = "concave",
#'   label_repel = TRUE
#' )
#'
#' # Custom cluster colors
#' ConstellationPlot(
#'   seu,
#'   cluster_col = "celltype",
#'   colors = c("T cell" = "red", "B cell" = "blue")
#' )
#' }
ConstellationPlot <- function(seu,
                              cluster_col,
                              reduction = "umap",
                              k = 15,
                              frac_th = 0.05,
                              colors = NULL,
                              color_by = NULL,
                              color_mode = "group",
                              hull_by = NULL,
                              hull_type = "none",
                              hull_alpha = 0.2,
                              hull_expand = 0.1,
                              node.label = "cluster_label",
                              exxageration = 2,
                              curved = TRUE,
                              node.dodge = FALSE,
                              label.size = 5,
                              max_size = 10,
                              label_repel = FALSE,
                              node_trans = "sqrt") {

  # Step 1: Build KNN graph (pure data, no styling)
  # Note: color_by/hull_by still accepted here for backward compatibility
  # but will trigger deprecation warnings from build_knn_graph()
  knn_graph <- build_knn_graph(
    seu = seu,
    cluster_col = cluster_col,
    reduction = reduction,
    k = k,
    color_by = color_by,
    hull_by = hull_by
  )

  # Step 2: Filter edges
  knn_graph <- filter_knn_edges(knn_graph, frac_th = frac_th)

  # Step 3: Compute cluster centers (pure geometry)
  knn_graph <- compute_cluster_centers(knn_graph)

  # Step 4: Assign groups (new explicit approach)
  # If color_by/hull_by were passed, they were already handled by build_knn_graph
  # for backward compatibility. For new code, users should call assign_cluster_groups directly.

  # If color_by was NOT passed to build_knn_graph but user wants group coloring,
  # they should use the new workflow. Here we only handle the new explicit case
  # when color_by is NULL but user might have pre-assigned groups.

  # Step 5: Resolve colors with explicit precedence
  knn_graph <- resolve_cluster_colors(
    knn_graph,
    cluster_colors = colors,
    color_mode = color_mode
  )

  # Step 6: Plot
  plot_constellation(
    knn_graph,
    node.label = node.label,
    exxageration = exxageration,
    curved = curved,
    node.dodge = node.dodge,
    label.size = label.size,
    max_size = max_size,
    label_repel = label_repel,
    node_trans = node_trans,
    hull_type = hull_type,
    hull_alpha = hull_alpha,
    hull_expand = hull_expand
  )
}
