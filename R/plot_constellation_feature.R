#' @title Plot Constellation Feature
#' @description Plot a constellation graph with nodes colored by gene expression
#' @name plot_constellation_feature
NULL

#' Plot Constellation Feature
#'
#' @description Creates a constellation plot where each node is colored by
#'   cluster-level expression of a feature (e.g., average expression or percent
#'   expressing). Intended to be used after build_knn_graph(),
#'   filter_knn_edges(), and compute_cluster_centers().
#'
#' @param knn_graph A constellation_knn object with cluster centers computed
#' @param seu A Seurat object used to compute feature expression
#' @param feature Character. Gene/feature name to visualize
#' @param assay Character. Assay to use (defaults to DefaultAssay)
#' @param layer Character. Layer to use (default "data"). For Seurat V5 compatibility.
#' @param feature_stat Character. "avg", "median", or "pct" (percent expressing)
#' @param pct_threshold Numeric. Threshold for percent expressing (default 0)
#' @param node.label Character. Column name for node labels. Default "cluster_label"
#' @param exxageration Numeric. Edge width exaggeration factor. Default 2
#' @param curved Logical. Whether edges should be curved. Default TRUE
#' @param node.dodge Logical. Whether to dodge overlapping nodes. Default FALSE
#' @param label.size Numeric. Size of node labels. Default 5
#' @param max_size Numeric. Maximum node size. Default 10
#' @param label_repel Logical. Use ggrepel for labels. Default FALSE
#' @param node_trans Character. Size transformation. Default "sqrt"
#' @param hull_type Character. Hull type: "none", "convex", or "concave". Default "none"
#' @param hull_alpha Numeric. Hull transparency (0-1). Default 0.2
#' @param hull_expand Numeric. Expansion factor for hulls. Default 0.1
#'
#' @return A ggplot object
#'
#' @export
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom ggplot2 scale_color_viridis_c
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#' seu %>%
#'   build_knn_graph(cluster_col = "celltype") %>%
#'   filter_knn_edges(frac_th = 0.05) %>%
#'   compute_cluster_centers() %>%
#'   plot_constellation_feature(seu, feature = "MS4A1")
#' }
plot_constellation_feature <- function(knn_graph,
                                       seu,
                                       feature,
                                       assay = NULL,
                                       layer = "data",
                                       feature_stat = c("avg", "median", "pct"),
                                       pct_threshold = 0,
                                       node.label = "cluster_label",
                                       exxageration = 2,
                                       curved = TRUE,
                                       node.dodge = FALSE,
                                       label.size = 5,
                                       max_size = 10,
                                       label_repel = FALSE,
                                       node_trans = "sqrt",
                                       hull_type = "none",
                                       hull_alpha = 0.2,
                                       hull_expand = 0.1) {

  if (!inherits(knn_graph, "constellation_knn")) {
    stop("Input must be a constellation_knn object")
  }

  if (is.null(knn_graph$cl.center.df)) {
    stop("Cluster centers not computed. Run compute_cluster_centers() first.")
  }

  if (!inherits(seu, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (length(feature) != 1) {
    stop("feature must be a single gene/feature name")
  }

  feature_stat <- match.arg(feature_stat)

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seu)
  }

  data_mat <- Seurat::GetAssayData(seu, assay = assay, layer = layer)
  if (!feature %in% rownames(data_mat)) {
    stop(paste("Feature", feature, "not found in assay", assay))
  }

  expr_vec <- as.numeric(data_mat[feature, , drop = TRUE])
  names(expr_vec) <- colnames(data_mat)

  cells <- intersect(names(knn_graph$cl), names(expr_vec))
  if (length(cells) == 0) {
    stop("No overlapping cell names between knn_graph and seu")
  }

  if (length(cells) < length(knn_graph$cl)) {
    warning("Some cells in knn_graph are missing from seu; using overlap only")
  }

  cl_vec <- knn_graph$cl[cells]
  expr_vec <- expr_vec[cells]

  stat_fun <- switch(
    feature_stat,
    avg = function(x) mean(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    pct = function(x) 100 * mean(x > pct_threshold, na.rm = TRUE)
  )

  feature_by_cl <- tapply(expr_vec, cl_vec, stat_fun)

  cl.center.df <- knn_graph$cl.center.df
  cl.center.df$feature_value <- feature_by_cl[as.character(cl.center.df$cl)]

  missing <- is.na(cl.center.df$feature_value)
  if (any(missing)) {
    warning("Missing feature values for some clusters; setting to 0")
    cl.center.df$feature_value[missing] <- 0
  }

  knn.cl.df <- knn_graph$knn.cl.df

  nodes <- .prepare_nodes(cl.center.df, knn.cl.df, node.label, max_size, node_trans)

  if (node.dodge) {
    nodes <- .dodge_nodes(nodes)
  }

  nodes <- nodes[order(nodes$cl), ]

  edge_data <- .prepare_edges(knn.cl.df, nodes, exxageration)
  knn.cl.lines <- edge_data$knn.cl.lines
  nodes <- edge_data$nodes

  if (!is.null(knn.cl.lines) && nrow(knn.cl.lines) > 0) {
    line.segments <- .create_line_segments(knn.cl.lines, nodes, exxageration)
    poly.Edges <- .create_poly_edges(line.segments, curved)
  } else {
    poly.Edges <- data.frame(x = numeric(), y = numeric(), Group = character())
  }

  hull_data <- NULL
  if (hull_type != "none" && "hull_group" %in% colnames(nodes)) {
    hull_data <- .compute_hulls(nodes, hull_type, hull_expand)
  }

  stat_label <- switch(
    feature_stat,
    avg = "Average expression",
    median = "Median expression",
    pct = "Percent expressing"
  )
  legend_title <- paste0(feature, " (", stat_label, ")")

  plot.all <- .build_plot_feature(
    nodes,
    poly.Edges,
    node.label,
    label.size,
    max_size,
    node_trans,
    label_repel,
    hull_data,
    hull_alpha,
    legend_title,
    feature_stat
  )

  return(plot.all)
}

#' Build the final plot for feature coloring
#' @keywords internal
.build_plot_feature <- function(nodes, poly.Edges, node.label, label.size,
                                max_size, node_trans, label_repel,
                                hull_data = NULL, hull_alpha = 0.2,
                                legend_title = NULL,
                                feature_stat = "avg") {

  plot.all <- ggplot2::ggplot()

  if (!is.null(hull_data) && nrow(hull_data) > 0) {
    plot.all <- plot.all +
      ggplot2::geom_polygon(
        data = hull_data,
        ggplot2::aes(x = x, y = y, group = hull_group, fill = hull_color),
        alpha = hull_alpha
      ) +
      ggplot2::scale_fill_identity()
  }

  if (nrow(poly.Edges) > 0) {
    plot.all <- plot.all +
      ggplot2::geom_polygon(
        data = poly.Edges,
        alpha = 0.2,
        ggplot2::aes(x = x, y = y, group = Group)
      )
  }

  plot.all <- plot.all +
    ggplot2::geom_point(
      data = nodes,
      alpha = 0.9,
      shape = 19,
      ggplot2::aes(x = x, y = y, size = cluster_size, color = feature_value)
    ) +
    ggplot2::scale_size_area(
      trans = node_trans,
      max_size = max_size,
      breaks = c(100, 1000, 10000, 100000),
      guide = "none"
    )

  if (label_repel) {
    plot.all <- plot.all +
      ggrepel::geom_text_repel(
        data = nodes,
        ggplot2::aes(x = x, y = y, label = .data[[node.label]]),
        size = label.size,
        min.segment.length = Inf
      )
  } else {
    plot.all <- plot.all +
      ggplot2::geom_text(
        data = nodes,
        ggplot2::aes(x = x, y = y, label = .data[[node.label]]),
        size = label.size
      )
  }

  color_limits <- NULL
  if (feature_stat == "pct") {
    color_limits <- c(0, 100)
  }

  plot.all <- plot.all +
    ggplot2::scale_color_viridis_c(
      name = legend_title,
      limits = color_limits,
      na.value = "grey80"
    ) +
    ggplot2::theme_void()

  return(plot.all)
}
