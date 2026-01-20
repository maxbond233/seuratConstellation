#' @title Compute Cluster Centers
#' @description Calculate cluster centroid coordinates for constellation plot
#' @name compute_cluster_centers
NULL

#' Compute Cluster Centers
#'
#' @description Calculates the centroid coordinates for each cluster based on
#'   dimensional reduction embeddings. This is a pure geometry computation step
#'   that does not handle colors. Use resolve_cluster_colors() to assign colors.
#'
#' @param knn_graph A constellation_knn object from build_knn_graph()
#' @param colors DEPRECATED. Use resolve_cluster_colors() instead
#' @param color_mode DEPRECATED. Use resolve_cluster_colors() instead
#'
#' @return A constellation_knn object with cluster center coordinates added
#'
#' @export
#' @importFrom dplyr group_by summarise
#'
#' @examples
#' \dontrun{
#' # New workflow (recommended)
#' knn_graph %>%
#'   compute_cluster_centers() %>%
#'   resolve_cluster_colors()
#'
#' # With group-based coloring
#' knn_graph %>%
#'   compute_cluster_centers() %>%
#'   assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
#'   resolve_cluster_colors(color_mode = "gradient")
#' }
compute_cluster_centers <- function(knn_graph, colors = NULL, color_mode = "group") {

  if (!inherits(knn_graph, "constellation_knn")) {
    stop("Input must be a constellation_knn object")
  }

  # Deprecation warnings
  if (!is.null(colors)) {
    warning(
      "The 'colors' parameter is deprecated and will be ignored.\n",
      "Please use resolve_cluster_colors(knn_graph, cluster_colors = ...) instead.",
      call. = FALSE
    )
  }

  if (!missing(color_mode) && color_mode != "group") {
    warning(
      "The 'color_mode' parameter is deprecated and will be ignored.\n",
      "Please use resolve_cluster_colors(knn_graph, color_mode = ...) instead.",
      call. = FALSE
    )
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

  knn_graph$cl.center.df <- cl.center.df

  return(knn_graph)
}
