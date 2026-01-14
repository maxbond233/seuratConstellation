#' @title ConstellationPlot - One-stop Convenience Function
#' @description All-in-one function for creating constellation plots from Seurat objects
#' @name ConstellationPlot
NULL

#' ConstellationPlot
#'
#' @description A convenience wrapper that combines all steps of constellation
#'   plot creation into a single function call. This is the recommended entry
#'   point for most users.
#'
#' @param seu A Seurat object (V5 compatible)
#' @param cluster_col Character. Column name in meta.data for cluster identity
#' @param reduction Character. Name of dimensional reduction. Default "umap"
#' @param k Integer. Number of nearest neighbors. Default 15
#' @param frac_th Numeric. Minimum fraction threshold for edges. Default 0.05
#' @param colors Optional. Named vector of colors, or NULL for auto-generation
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
#' # With custom options
#' ConstellationPlot(
#'   seu,
#'   cluster_col = "celltype",
#'   reduction = "umap",
#'   frac_th = 0.05,
#'   label_repel = TRUE,
#'   node.dodge = TRUE
#' )
#' }
ConstellationPlot <- function(seu,
                              cluster_col,
                              reduction = "umap",
                              k = 15,
                              frac_th = 0.05,
                              colors = NULL,
                              node.label = "cluster_label",
                              exxageration = 2,
                              curved = TRUE,
                              node.dodge = FALSE,
                              label.size = 5,
                              max_size = 10,
                              label_repel = FALSE,
                              node_trans = "sqrt") {

  knn_graph <- build_knn_graph(
    seu = seu,
    cluster_col = cluster_col,
    reduction = reduction,
    k = k
  )

  knn_graph <- filter_knn_edges(knn_graph, frac_th = frac_th)

  knn_graph <- compute_cluster_centers(knn_graph, colors = colors)

  plot_constellation(
    knn_graph,
    node.label = node.label,
    exxageration = exxageration,
    curved = curved,
    node.dodge = node.dodge,
    label.size = label.size,
    max_size = max_size,
    label_repel = label_repel,
    node_trans = node_trans
  )
}
