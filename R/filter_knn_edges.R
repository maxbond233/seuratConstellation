#' @title Filter KNN Edges
#' @description Filter KNN graph edges by fraction threshold
#' @name filter_knn_edges
NULL

#' Filter KNN Edges
#'
#' @description Filters edges in a KNN graph based on the fraction of
#'   connections between clusters. This is useful for removing weak
#'   connections before visualization.
#'
#' @param knn_graph A constellation_knn object from build_knn_graph()
#' @param frac_th Numeric. Minimum fraction threshold for edges. Default 0.05
#'
#' @return A constellation_knn object with filtered edges
#'
#' @export
#' @importFrom dplyr filter mutate
#'
#' @examples
#' \dontrun{
#' knn_graph %>% filter_knn_edges(frac_th = 0.05)
#' }
filter_knn_edges <- function(knn_graph, frac_th = 0.05) {

  if (!inherits(knn_graph, "constellation_knn")) {
    stop("Input must be a constellation_knn object from build_knn_graph()")
  }

  knn_graph$knn.cl.df <- knn_graph$knn.cl.df %>%
    dplyr::filter(frac >= frac_th) %>%
    dplyr::mutate(
      cl.from = as.numeric(as.character(cl.from)),
      cl.to = as.numeric(as.character(cl.to))
    )

  knn_graph$params$frac_th <- frac_th

  return(knn_graph)
}
