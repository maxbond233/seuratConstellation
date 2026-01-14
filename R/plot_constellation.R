#' @title Plot Constellation
#' @description Main plotting function for constellation visualization
#' @name plot_constellation
NULL

#' Prepare nodes for plotting
#' @keywords internal
.prepare_nodes <- function(cl.center.df, knn.cl.df, node.label, max_size, node_trans) {

  knn.cl.same <- knn.cl.df[knn.cl.df$cl.from == knn.cl.df$cl.to, ]
  cl.center.df$edge.frac.within <- knn.cl.same$frac[
    match(cl.center.df$cl, knn.cl.same$cl.from)
  ]

  labels <- cl.center.df[[node.label]]

  p.nodes <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = cl.center.df,
      shape = 19,
      ggplot2::aes(x = x, y = y, size = cluster_size,
                   color = alpha(cluster_color, 0.8))
    ) +
    ggplot2::scale_size_area(
      trans = node_trans,
      max_size = max_size,
      breaks = c(100, 1000, 10000, 100000)
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_text(
      data = cl.center.df,
      ggplot2::aes(x = x, y = y, label = labels),
      size = 3
    )

  g <- ggplot2::ggplot_build(p.nodes)
  dots <- g[["data"]][[1]]

  nodes <- dplyr::left_join(cl.center.df, dots, by = c("x", "y")) %>%
    dplyr::ungroup()

  return(nodes)
}

#' Dodge overlapping nodes
#' @keywords internal
.dodge_nodes <- function(nodes) {
  nodes$r <- (nodes$size / 10) / 2

  x.list <- c(mean(nodes$x), nodes$x)
  y.list <- c(mean(nodes$y), nodes$y)
  dist.test <- as.matrix(dist(cbind(x.list, y.list)))
  nodes$distance <- dist.test[2:nrow(dist.test), 1]
  nodes <- nodes[order(nodes$distance), ]

  for (pass in 1:2) {
    for (d1 in 1:(nrow(nodes) - 1)) {
      j <- d1 + 1
      for (d2 in j:nrow(nodes)) {
        distSq <- sqrt(
          ((nodes$x[d1] - nodes$x[d2])^2) +
          ((nodes$y[d1] - nodes$y[d2])^2)
        )
        radSumSq <- (nodes$r[d1] * 1.25) + (nodes$r[d2] * 1.25)

        if (distSq < radSumSq) {
          angsk <- seq(0, 2 * pi, length.out = 2)
          nodes$x[d2] <- nodes$x[d2] + cos(angsk[1]) * (nodes$r[d1] + nodes$r[d2] + 0.5)
          nodes$y[d2] <- nodes$y[d2] + sin(angsk[1]) * (nodes$r[d1] + nodes$r[d2] + 0.5)
        }
      }
    }
  }

  return(nodes)
}

#' Find bidirectional edges
#' @keywords internal
.find_bidirectional <- function(knn.cl.d) {
  knn.cl.bid <- NULL
  for (i in seq_len(nrow(knn.cl.d))) {
    line <- knn.cl.d[i, ]
    r <- knn.cl.d[i:nrow(knn.cl.d), ]
    r <- r[(line$cl.from == r$cl.to & line$cl.to == r$cl.from), ]
    if (nrow(r) > 0) {
      line$Freq.to <- r$Freq[1]
      line$node.pt.to <- r$node.pt.from[1]
      line$frac.to <- r$frac[1]
      knn.cl.bid <- rbind(knn.cl.bid, line)
    }
  }
  return(knn.cl.bid)
}

#' Find unidirectional edges
#' @keywords internal
.find_unidirectional <- function(knn.cl.d) {
  knn.cl.uni <- NULL
  for (i in seq_len(nrow(knn.cl.d))) {
    line <- knn.cl.d[i, ]
    r <- knn.cl.d[(line$cl.from == knn.cl.d$cl.to & line$cl.to == knn.cl.d$cl.from), ]
    if (nrow(r) == 0) {
      knn.cl.uni <- rbind(knn.cl.uni, line)
    }
  }
  return(knn.cl.uni)
}

#' Prepare edge data
#' @keywords internal
.prepare_edges <- function(knn.cl.df, nodes, exxageration) {
  knn.cl <- knn.cl.df
  knn.cl.d <- knn.cl[!(knn.cl$cl.from == knn.cl$cl.to), ]

  nodes$cl <- as.numeric(as.character(nodes$cl))
  knn.cl.d$cl.from <- as.numeric(as.character(knn.cl.d$cl.from))
  knn.cl.d$cl.to <- as.numeric(as.character(knn.cl.d$cl.to))

  nodes$node.width <- nodes$size
  knn.cl.d <- dplyr::left_join(
    knn.cl.d,
    dplyr::select(nodes, cl, node.width),
    by = c("cl.from" = "cl")
  )
  colnames(knn.cl.d)[colnames(knn.cl.d) == "node.width"] <- "node.pt.from"
  knn.cl.d$node.pt.to <- ""
  knn.cl.d$Freq.to <- ""
  knn.cl.d$frac.to <- ""

  knn.cl.bid <- .find_bidirectional(knn.cl.d)
  knn.cl.uni <- .find_unidirectional(knn.cl.d)

  if (!is.null(knn.cl.uni) && nrow(knn.cl.uni) > 0) {
    knn.cl.uni$node.pt.to <- nodes$node.width[match(knn.cl.uni$cl.to, nodes$cl)]
    knn.cl.uni$Freq.to <- 1
    knn.cl.uni$frac.to <- 0.01
  }

  knn.cl.lines <- rbind(knn.cl.bid, knn.cl.uni)
  return(list(knn.cl.lines = knn.cl.lines, nodes = nodes))
}

#' Create line segments
#' @keywords internal
.create_line_segments <- function(knn.cl.lines, nodes, exxageration) {
  line.segments <- knn.cl.lines[, c("cl.from", "cl.to")]
  nodes$cl <- as.numeric(as.character(nodes$cl))

  line.segments <- dplyr::left_join(
    line.segments,
    dplyr::select(nodes, x, y, cl),
    by = c("cl.from" = "cl")
  )
  line.segments <- dplyr::left_join(
    line.segments,
    dplyr::select(nodes, x, y, cl),
    by = c("cl.to" = "cl")
  )
  colnames(line.segments) <- c("cl.from", "cl.to", "x.from", "y.from", "x.to", "y.to")

  line.segments <- data.frame(
    line.segments,
    freq.from = knn.cl.lines$Freq,
    freq.to = as.numeric(knn.cl.lines$Freq.to),
    frac.from = knn.cl.lines$frac,
    frac.to = as.numeric(knn.cl.lines$frac.to),
    node.pt.from = as.numeric(knn.cl.lines$node.pt.from),
    node.pt.to = as.numeric(knn.cl.lines$node.pt.to)
  )

  line.segments$node.size.from <- line.segments$node.pt.from / 10
  line.segments$node.size.to <- line.segments$node.pt.to / 10

  max_frac <- max(line.segments$frac.from, line.segments$frac.to, na.rm = TRUE)
  line.segments$line.width.from <- (line.segments$frac.from / max_frac) * line.segments$node.size.from
  line.segments$line.width.to <- (line.segments$frac.to / max_frac) * line.segments$node.size.to

  line.segments$ex.line.from <- pmin(
    line.segments$line.width.from * exxageration,
    line.segments$node.size.from
  )
  line.segments$ex.line.to <- pmin(
    line.segments$line.width.to * exxageration,
    line.segments$node.size.to
  )

  line.segments <- stats::na.omit(line.segments)
  return(line.segments)
}

#' Create polygon edges
#' @keywords internal
.create_poly_edges <- function(line.segments, curved = TRUE) {
  if (nrow(line.segments) == 0) {
    return(data.frame(x = numeric(), y = numeric(), Group = character()))
  }

  allEdges <- lapply(
    seq_len(nrow(line.segments)),
    .edgeMaker,
    len = 500,
    curved = curved,
    line.segments = line.segments
  )
  allEdges <- do.call(rbind, allEdges)

  groups <- unique(allEdges$Group)
  poly.Edges <- data.frame(x = numeric(), y = numeric(), Group = character())

  for (i in seq_along(groups)) {
    select.group <- groups[i]
    select.edge <- allEdges[allEdges$Group == select.group, ]

    x <- select.edge$x
    y <- select.edge$y
    w <- select.edge$fraction
    N <- length(x)

    leftx <- lefty <- rightx <- righty <- numeric(N)

    perps <- .perpStart(x[1:2], y[1:2], w[1] / 2)
    leftx[1] <- perps[1, 1]; lefty[1] <- perps[1, 2]
    rightx[1] <- perps[2, 1]; righty[1] <- perps[2, 2]

    for (ii in 2:(N - 1)) {
      seq_idx <- (ii - 1):(ii + 1)
      perps <- .perpMid(as.numeric(x[seq_idx]), as.numeric(y[seq_idx]), w[ii] / 2)
      leftx[ii] <- perps[1, 1]; lefty[ii] <- perps[1, 2]
      rightx[ii] <- perps[2, 1]; righty[ii] <- perps[2, 2]
    }

    perps <- .perpEnd(x[(N - 1):N], y[(N - 1):N], w[N] / 2)
    leftx[N] <- perps[1, 1]; lefty[N] <- perps[1, 2]
    rightx[N] <- perps[2, 1]; righty[N] <- perps[2, 2]

    lineleft <- data.frame(x = leftx, y = lefty)
    lineright <- data.frame(x = rightx, y = righty)
    lineright <- lineright[nrow(lineright):1, ]
    lines.lr <- rbind(lineleft, lineright)
    lines.lr$Group <- select.group

    poly.Edges <- rbind(poly.Edges, lines.lr)
  }

  return(poly.Edges)
}

#' Plot Constellation
#'
#' @description Creates a constellation plot visualization from a KNN graph.
#'   This function is designed to work in a pipeline with build_knn_graph(),
#'   filter_knn_edges(), and compute_cluster_centers().
#'
#' @param knn_graph A constellation_knn object with cluster centers computed
#' @param node.label Character. Column name for node labels. Default "cluster_label"
#' @param exxageration Numeric. Edge width exaggeration factor. Default 2
#' @param curved Logical. Whether edges should be curved. Default TRUE
#' @param node.dodge Logical. Whether to dodge overlapping nodes. Default FALSE
#' @param label.size Numeric. Size of node labels. Default 5
#' @param max_size Numeric. Maximum node size. Default 10
#' @param label_repel Logical. Use ggrepel for labels. Default FALSE
#' @param node_trans Character. Size transformation. Default "sqrt"
#'
#' @return A ggplot object
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_polygon geom_text
#' @importFrom ggplot2 scale_size_area scale_color_identity theme_void theme
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#' seu %>%
#'   build_knn_graph(cluster_col = "celltype") %>%
#'   filter_knn_edges(frac_th = 0.05) %>%
#'   compute_cluster_centers() %>%
#'   plot_constellation(label_repel = TRUE)
#' }
plot_constellation <- function(knn_graph,
                               node.label = "cluster_label",
                               exxageration = 2,
                               curved = TRUE,
                               node.dodge = FALSE,
                               label.size = 5,
                               max_size = 10,
                               label_repel = FALSE,
                               node_trans = "sqrt") {

  if (!inherits(knn_graph, "constellation_knn")) {
    stop("Input must be a constellation_knn object")
  }

  if (is.null(knn_graph$cl.center.df)) {
    stop("Cluster centers not computed. Run compute_cluster_centers() first.")
  }

  knn.cl.df <- knn_graph$knn.cl.df
  cl.center.df <- knn_graph$cl.center.df

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

  plot.all <- .build_plot(nodes, poly.Edges, node.label, label.size,
                          max_size, node_trans, label_repel)

  return(plot.all)
}

#' Build the final plot
#' @keywords internal
.build_plot <- function(nodes, poly.Edges, node.label, label.size,
                        max_size, node_trans, label_repel) {

  plot.all <- ggplot2::ggplot()

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
      alpha = 0.8,
      shape = 19,
      ggplot2::aes(x = x, y = y, size = cluster_size, color = cluster_color)
    ) +
    ggplot2::scale_size_area(
      trans = node_trans,
      max_size = max_size,
      breaks = c(100, 1000, 10000, 100000)
    ) +
    ggplot2::scale_color_identity()

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

  plot.all <- plot.all +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  return(plot.all)
}
