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
