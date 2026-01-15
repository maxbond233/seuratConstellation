#' @title Utility Functions for Constellation Plot
#' @description Internal utility functions for edge drawing and geometric calculations
#' @name utils
#' @keywords internal
NULL

#' Calculate angle between two points
#' @param x numeric vector of length 2
#' @param y numeric vector of length 2
#' @return angle in radians
#' @keywords internal
.angle <- function(x, y) {
  atan2(y[2] - y[1], x[2] - x[1])
}

#' Calculate average angle for three points
#' @param x numeric vector of length 3
#' @param y numeric vector of length 3
#' @return average angle in radians
#' @keywords internal
.avgangle <- function(x, y) {
  a1 <- .angle(x[1:2], y[1:2])
  a2 <- .angle(x[2:3], y[2:3])
  atan2(sin(a1) + sin(a2), cos(a1) + cos(a2))
}

#' Calculate perpendicular points
#' @param x numeric vector
#' @param y numeric vector
#' @param len half-width
#' @param a angle
#' @param mid midpoint index
#' @return matrix with upper and lower points
#' @keywords internal
.perp <- function(x, y, len, a, mid) {
  dx <- len * cos(a + pi/2)
  dy <- len * sin(a + pi/2)
  upper <- c(x[mid] + dx, y[mid] + dy)
  lower <- c(x[mid] - dx, y[mid] - dy)
  rbind(upper, lower)
}

#' Perpendicular at start point
#' @keywords internal
.perpStart <- function(x, y, len) {
  .perp(x, y, len, .angle(x, y), 1)
}

#' Perpendicular at midpoint
#' @keywords internal
.perpMid <- function(x, y, len) {
  .perp(x, y, len, .avgangle(x, y), 2)
}

#' Perpendicular at end point
#' @keywords internal
.perpEnd <- function(x, y, len) {
  .perp(x, y, len, .angle(x, y), 2)
}

#' Predict cluster membership using KNN
#' @param knn KNN matrix
#' @param ref.cells reference cell names
#' @param cl cluster assignments
#' @return list with pred.df and pred.prob
#' @keywords internal
.predict_knn <- function(knn, ref.cells, cl) {
  knn.cl <- matrix(cl[knn], nrow = nrow(knn))
  rownames(knn.cl) <- rownames(knn)

  cl.levels <- sort(unique(cl))
  pred.prob <- t(apply(knn.cl, 1, function(x) {
    tab <- table(factor(x, levels = cl.levels))
    tab / sum(tab)
  }))
  colnames(pred.prob) <- cl.levels

  pred.cl <- apply(pred.prob, 1, which.max)
  pred.cl <- cl.levels[pred.cl]
  names(pred.cl) <- rownames(pred.prob)

  list(pred.cl = pred.cl, pred.prob = pred.prob)
}

#' Compute majority category for each cluster
#' @param cluster_vec vector of cluster assignments
#' @param category_vec vector of category values
#' @param cluster_levels unique cluster levels
#' @return named vector of majority categories per cluster
#' @keywords internal
.compute_majority_category <- function(cluster_vec, category_vec, cluster_levels) {
  majority <- sapply(cluster_levels, function(cl) {
    idx <- cluster_vec == cl
    if (sum(idx) == 0) return(NA)
    tab <- table(category_vec[idx])
    names(tab)[which.max(tab)]
  })
  names(majority) <- cluster_levels
  majority
}

#' Assign colors based on group membership
#' @param cl.center.df data frame with cluster info including color_group
#' @param color_mode "group" for uniform colors, "gradient" for shades
#' @return vector of colors for each cluster
#' @keywords internal
#' @importFrom scales hue_pal
#' @importFrom grDevices colorRampPalette col2rgb rgb
.assign_group_colors <- function(cl.center.df, color_mode = "group") {
  groups <- unique(cl.center.df$color_group)
  n_groups <- length(groups)
  group_colors <- scales::hue_pal()(n_groups)
  names(group_colors) <- groups

  if (color_mode == "group") {
    # Uniform color per group
    colors <- group_colors[cl.center.df$color_group]
  } else if (color_mode == "gradient") {
    # Generate gradient shades within each group
    colors <- character(nrow(cl.center.df))
    for (grp in groups) {
      idx <- which(cl.center.df$color_group == grp)
      n_in_group <- length(idx)
      if (n_in_group == 1) {
        colors[idx] <- group_colors[grp]
      } else {
        base_col <- group_colors[grp]
        shades <- .generate_shades(base_col, n_in_group)
        colors[idx] <- shades
      }
    }
  } else {
    stop("color_mode must be 'group' or 'gradient'")
  }

  names(colors) <- cl.center.df$cluster_label
  colors
}

#' Generate color shades from a base color
#' @param base_color base color
#' @param n number of shades to generate
#' @return vector of color shades
#' @keywords internal
#' @importFrom grDevices col2rgb rgb
.generate_shades <- function(base_color, n) {
  rgb_base <- grDevices::col2rgb(base_color) / 255
  shades <- sapply(seq(0.6, 1, length.out = n), function(factor) {
    adjusted <- rgb_base * factor + (1 - factor) * 0.3
    adjusted <- pmin(pmax(adjusted, 0), 1)
    grDevices::rgb(adjusted[1], adjusted[2], adjusted[3])
  })
  shades
}

#' Compute hulls for each group
#' @param nodes data frame with node coordinates and hull_group
#' @param hull_type "convex" or "concave"
#' @param expand expansion factor
#' @return data frame with hull polygon coordinates
#' @keywords internal
#' @importFrom grDevices chull col2rgb rgb
.compute_hulls <- function(nodes, hull_type = "convex", expand = 0.1) {
  groups <- unique(nodes$hull_group)
  hull_list <- list()

  for (grp in groups) {
    grp_nodes <- nodes[nodes$hull_group == grp, ]
    if (nrow(grp_nodes) < 3) next

    pts <- as.matrix(grp_nodes[, c("x", "y")])
    grp_color <- grp_nodes$cluster_color[1]
    hull_color <- .lighten_color(grp_color, 0.5)

    if (hull_type == "convex") {
      hull_idx <- grDevices::chull(pts)
      hull_pts <- pts[hull_idx, ]
    } else if (hull_type == "concave") {
      if (requireNamespace("concaveman", quietly = TRUE)) {
        hull_pts <- concaveman::concaveman(pts)
      } else {
        warning("concaveman not installed, using convex hull")
        hull_idx <- grDevices::chull(pts)
        hull_pts <- pts[hull_idx, ]
      }
    }

    # Expand hull
    if (expand > 0) {
      hull_pts <- .expand_hull(hull_pts, expand)
    }

    hull_df <- data.frame(
      x = hull_pts[, 1],
      y = hull_pts[, 2],
      hull_group = grp,
      hull_color = hull_color
    )
    hull_list[[grp]] <- hull_df
  }

  do.call(rbind, hull_list)
}

#' Lighten a color
#' @param color input color
#' @param amount amount to lighten (0-1)
#' @return lightened color
#' @keywords internal
#' @importFrom grDevices col2rgb rgb
.lighten_color <- function(color, amount = 0.5) {
  rgb_val <- grDevices::col2rgb(color) / 255
  lightened <- rgb_val + (1 - rgb_val) * amount
  grDevices::rgb(lightened[1], lightened[2], lightened[3])
}

#' Expand hull outward from centroid
#' @param hull_pts matrix of hull points
#' @param expand expansion factor
#' @return expanded hull points
#' @keywords internal
.expand_hull <- function(hull_pts, expand) {
  centroid <- colMeans(hull_pts)
  expanded <- t(apply(hull_pts, 1, function(pt) {
    direction <- pt - centroid
    pt + direction * expand
  }))
  expanded
}
