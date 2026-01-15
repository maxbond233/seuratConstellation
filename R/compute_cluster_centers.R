#' @title Compute Cluster Centers
#' @description Calculate cluster centroid coordinates for constellation plot
#' @name compute_cluster_centers
NULL

#' Compute Cluster Centers
#'
#' @description Calculates the centroid coordinates for each cluster based on
#'   dimensional reduction embeddings. Supports custom colors or auto-generation
#'   with group-based coloring modes.
#'
#' @param knn_graph A constellation_knn object from build_knn_graph()
#' @param colors Optional. Named vector of colors for clusters, or NULL for auto
#' @param color_mode Character. Color mode: "group" for uniform group colors,
#'   "gradient" for similar shades within groups. Default "group"
#'
#' @return A constellation_knn object with cluster center coordinates added
#'
#' @export
#' @importFrom dplyr group_by summarise
#' @importFrom scales hue_pal
#' @importFrom grDevices colorRampPalette col2rgb rgb
#'
#' @examples
#' \dontrun{
#' knn_graph %>% compute_cluster_centers()
#' knn_graph %>% compute_cluster_centers(color_mode = "gradient")
#' knn_graph %>% compute_cluster_centers(colors = c("A" = "red", "B" = "blue"))
#' }
compute_cluster_centers <- function(knn_graph, colors = NULL, color_mode = "group") {

  if (!inherits(knn_graph, "constellation_knn")) {
    stop("Input must be a constellation_knn object")
  }

  rd.dat <- knn_graph$rd.dat
  cl <- knn_graph$cl
  cl.df <- knn_graph$cluster_info

  cl_coor <- data.frame(
    x = rd.dat[, 1],
    y = rd.dat[, 2],
    cl = cl
  )

  cl_coor <- cl_coor %>%
    dplyr::group_by(cl) %>%
    dplyr::summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
    as.data.frame()

  rownames(cl_coor) <- cl_coor$cl

  cl.center.df <- merge(cl_coor, cl.df, by.x = "cl", by.y = "cluster_id")

  # Determine colors based on color_group or custom colors
  if (!is.null(colors)) {
    # Use provided colors
    cl.center.df$cluster_color <- colors[cl.center.df$cluster_label]
    if (any(is.na(cl.center.df$cluster_color))) {
      missing <- cl.center.df$cluster_label[is.na(cl.center.df$cluster_color)]
      n_missing <- length(missing)
      cl.center.df$cluster_color[is.na(cl.center.df$cluster_color)] <-
        scales::hue_pal()(n_missing)
    }
  } else if ("color_group" %in% colnames(cl.center.df)) {
    # Use group-based coloring
    colors <- .assign_group_colors(cl.center.df, color_mode)
    cl.center.df$cluster_color <- colors
  } else {
    # Default: one color per cluster
    n_clusters <- nrow(cl.center.df)
    colors <- scales::hue_pal()(n_clusters)
    names(colors) <- cl.center.df$cluster_label
    cl.center.df$cluster_color <- colors[cl.center.df$cluster_label]
  }

  knn_graph$cl.center.df <- cl.center.df
  knn_graph$colors <- colors


  return(knn_graph)
}
