#' @title Edge Maker Functions
#' @description Functions for creating curved edges between nodes
#' @name edge_maker
#' @keywords internal
NULL

#' Create a curved edge between two points
#'
#' @param whichRow row index in line.segments
#' @param len number of points along the edge
#' @param line.segments data frame with edge information
#' @param curved logical, whether to curve the edge
#' @return data frame with edge coordinates
#' @keywords internal
.edgeMaker <- function(whichRow, len = 100, line.segments, curved = TRUE) {
  fromC <- unlist(line.segments[whichRow, c(3, 4)])
  toC <- unlist(line.segments[whichRow, c(5, 6)])


  graphCenter <- colMeans(line.segments[, c(3, 4)])
  bezierMid <- c(fromC[1], toC[2])
  distance1 <- sum((graphCenter - bezierMid)^2)

  if (distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)) {
    bezierMid <- c(toC[1], fromC[2])
  }

  bezierMid <- (fromC + toC + bezierMid) / 3

  if (!curved) {
    bezierMid <- (fromC + toC) / 2
  }

  edge <- data.frame(
    Hmisc::bezier(
      c(fromC[1], bezierMid[1], toC[1]),
      c(fromC[2], bezierMid[2], toC[2]),
      evaluation = len
    )
  )

  edge$fraction <- seq(
    line.segments$ex.line.from[whichRow],
    line.segments$ex.line.to[whichRow],
    length.out = len
  )

  edge$Group <- paste(line.segments[whichRow, 1:2], collapse = ">")
  return(edge)
}
