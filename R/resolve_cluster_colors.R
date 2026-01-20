#' @title Resolve Cluster Colors
#' @description Resolve final cluster colors with explicit precedence rules
#' @name resolve_cluster_colors
NULL

#' Resolve Cluster Colors
#'
#' @description Determines the final color for each cluster based on explicit precedence
#'   rules. This function separates color resolution from geometry computation,
#'   making color assignment logic transparent and predictable.
#'
#' @param knn_graph A constellation_knn object with cluster centers computed
#' @param cluster_colors Optional. Named vector of colors by cluster_label
#' @param group_colors Optional. Named vector of colors by group (requires color_group)
#' @param color_mode Character. Color mode: "group" (uniform) or "gradient" (shades). Default "group"
#'
#' @return A constellation_knn object with cluster_color added to cl.center.df
#'
#' @details
#' Color resolution precedence (highest to lowest):
#' \enumerate{
#'   \item \strong{cluster_colors}: If provided and covers all clusters, use directly
#'   \item \strong{group_colors + color_group}: If color_group exists and group_colors provided, map via groups
#'   \item \strong{color_group + auto palette}: If color_group exists, auto-generate group palette
#'   \item \strong{auto cluster palette}: One unique color per cluster
#' }
#'
#' Color modes (when using color_group):
#' \itemize{
#'   \item \strong{group}: All clusters in same group get identical color
#'   \item \strong{gradient}: Clusters in same group get similar shades
#' }
#'
#' @export
#' @importFrom scales hue_pal
#'
#' @examples
#' \dontrun{
#' # Auto colors per cluster
#' knn_graph %>% resolve_cluster_colors()
#'
#' # Custom cluster colors
#' knn_graph %>% resolve_cluster_colors(
#'   cluster_colors = c("A" = "red", "B" = "blue", "C" = "green")
#' )
#'
#' # Group-based coloring (requires assign_cluster_groups first)
#' knn_graph %>%
#'   assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
#'   resolve_cluster_colors(color_mode = "gradient")
#'
#' # Custom group colors
#' knn_graph %>%
#'   assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
#'   resolve_cluster_colors(
#'     group_colors = c("Neuron" = "#E41A1C", "Glia" = "#377EB8"),
#'     color_mode = "group"
#'   )
#' }
resolve_cluster_colors <- function(knn_graph,
                                   cluster_colors = NULL,
                                   group_colors = NULL,
                                   color_mode = "group") {

  if (!inherits(knn_graph, "constellation_knn")) {
    stop("Input must be a constellation_knn object")
  }

  if (is.null(knn_graph$cl.center.df)) {
    stop("Cluster centers not computed. Run compute_cluster_centers() first.")
  }

  if (!color_mode %in% c("group", "gradient")) {
    stop("color_mode must be 'group' or 'gradient'")
  }

  cl.center.df <- knn_graph$cl.center.df
  cluster_labels <- cl.center.df$cluster_label
  n_clusters <- length(cluster_labels)

  # Precedence 1: User-provided cluster_colors (highest priority)
  if (!is.null(cluster_colors)) {
    if (all(cluster_labels %in% names(cluster_colors))) {
      # All clusters covered
      final_colors <- cluster_colors[as.character(cluster_labels)]
      names(final_colors) <- cluster_labels
      message("Using provided cluster_colors")
    } else {
      # Partial coverage - warn and fill missing
      missing <- setdiff(cluster_labels, names(cluster_colors))
      warning(
        "cluster_colors missing for: ",
        paste(missing, collapse = ", "),
        ". Filling with auto-generated colors."
      )
      final_colors <- cluster_colors[as.character(cluster_labels)]
      n_missing <- sum(is.na(final_colors))
      final_colors[is.na(final_colors)] <- scales::hue_pal()(n_missing)
      names(final_colors) <- cluster_labels
    }

  # Precedence 2 & 3: Group-based coloring (if color_group exists)
  } else if ("color_group" %in% colnames(cl.center.df)) {
    groups <- unique(cl.center.df$color_group)
    n_groups <- length(groups)

    # Generate or use provided group_colors
    if (is.null(group_colors)) {
      group_colors <- scales::hue_pal()(n_groups)
      names(group_colors) <- groups
      message("Auto-generating group palette for color_group")
    } else {
      if (!all(groups %in% names(group_colors))) {
        missing_groups <- setdiff(groups, names(group_colors))
        warning(
          "group_colors missing for groups: ",
          paste(missing_groups, collapse = ", "),
          ". Filling with auto-generated colors."
        )
        n_missing_groups <- length(missing_groups)
        missing_cols <- scales::hue_pal()(n_missing_groups)
        names(missing_cols) <- missing_groups
        group_colors <- c(group_colors, missing_cols)
      }
    }

    # Apply color_mode
    if (color_mode == "group") {
      # Uniform color per group
      final_colors <- group_colors[as.character(cl.center.df$color_group)]
      names(final_colors) <- cluster_labels
      message("Using group colors (mode: uniform)")
    } else if (color_mode == "gradient") {
      # Generate gradient shades within each group
      final_colors <- character(n_clusters)
      for (grp in groups) {
        idx <- which(cl.center.df$color_group == grp)
        n_in_group <- length(idx)
        base_col <- group_colors[grp]
        if (n_in_group == 1) {
          final_colors[idx] <- base_col
        } else {
          shades <- .generate_shades(base_col, n_in_group)
          final_colors[idx] <- shades
        }
      }
      names(final_colors) <- cluster_labels
      message("Using group colors (mode: gradient)")
    }

  # Precedence 4: Auto cluster palette (lowest priority)
  } else {
    final_colors <- scales::hue_pal()(n_clusters)
    names(final_colors) <- cluster_labels
    message("Auto-generating colors per cluster")
  }

  # Add to cl.center.df
  knn_graph$cl.center.df$cluster_color <- final_colors[as.character(cluster_labels)]
  knn_graph$colors <- final_colors  # Keep for backward compatibility

  return(knn_graph)
}
