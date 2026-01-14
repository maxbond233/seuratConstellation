#' @title Build KNN Graph from Seurat Object
#' @description Constructs a KNN graph for constellation plot visualization
#' @name build_knn_graph
NULL

#' Build KNN Graph
#'
#' @description Builds a K-nearest neighbor graph from a Seurat object for
#'   constellation plot visualization. This function extracts cell embeddings
#'   and cluster information to compute inter-cluster connectivity.
#'
#' @param seu A Seurat object (V5 compatible)
#' @param cluster_col Character. Column name in meta.data for cluster identity
#' @param reduction Character. Name of dimensional reduction to use. Default "umap"
#' @param k Integer. Number of nearest neighbors. Default 15
#' @param knn.outlier.th Numeric. Threshold for outlier detection. Default 2
#' @param outlier.frac.th Numeric. Fraction threshold for outlier cells. Default 0.5
#'
#' @return A list with class "constellation_knn" containing:
#'   \item{knn.result}{Raw KNN results from RANN::nn2}
#'   \item{knn.cl.df}{Data frame of cluster-to-cluster edge information}
#'   \item{cluster_info}{Data frame with cluster metadata}
#'   \item{params}{List of parameters used}
#'
#' @export
#' @importFrom Seurat Embeddings
#' @importFrom RANN nn2
#' @importFrom dplyr %>% mutate arrange
#'
#' @examples
#' \dontrun{
#' knn_graph <- build_knn_graph(seu, cluster_col = "celltype", reduction = "umap")
#' }
build_knn_graph <- function(seu,
                            cluster_col,
                            reduction = "umap",
                            k = 15,
                            knn.outlier.th = 2,
                            outlier.frac.th = 0.5) {


  if (!inherits(seu, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(paste("Column", cluster_col, "not found in meta.data"))
  }

  if (!reduction %in% names(seu@reductions)) {
    stop(paste("Reduction", reduction, "not found. Available:",
               paste(names(seu@reductions), collapse = ", ")))
  }

  rd.dat <- Seurat::Embeddings(seu, reduction = reduction)

  cl_raw <- seu@meta.data[[cluster_col]]
  cl_levels <- unique(cl_raw)
  cl <- as.numeric(factor(cl_raw, levels = cl_levels))
  names(cl) <- colnames(seu)

  cl.df <- data.frame(
    cluster_id = seq_along(cl_levels),
    cluster_label = cl_levels,
    cluster_size = as.numeric(table(factor(cl_raw, levels = cl_levels))),
    row.names = seq_along(cl_levels)
  )

  knn.result <- RANN::nn2(rd.dat, k = k)
  rownames(knn.result[[1]]) <- rownames(knn.result[[2]]) <- rownames(rd.dat)
  knn <- knn.result[[1]]
  knn.dist <- knn.result[[2]]

  cl.knn.dist.mean <- tapply(names(cl), cl, function(x) mean(knn.dist[x, -1]))
  cl.knn.dist.sd <- tapply(names(cl), cl, function(x) sd(knn.dist[x, -1]))
  cl.knn.dist.th <- cl.knn.dist.mean + knn.outlier.th * cl.knn.dist.sd

  knn.dist.th <- cl.knn.dist.th[as.character(cl[rownames(knn)])]
  outlier <- apply(knn.dist, 2, function(x) x > knn.dist.th)
  rownames(outlier) <- rownames(knn.dist)
  knn[outlier] <- NA
  select.cells <- rownames(outlier)[rowMeans(outlier) < outlier.frac.th]

  pred.result <- .predict_knn(knn[select.cells, ], rownames(rd.dat), cl)
  pred.prob <- pred.result$pred.prob
  knn.cell.cl.counts <- round(pred.prob * ncol(knn))

  knn.cl.cl.counts <- do.call("rbind", tapply(
    rownames(pred.prob),
    cl[rownames(pred.prob)],
    function(x) colSums(knn.cell.cl.counts[x, , drop = FALSE])
  ))

  knn.cl.df <- as.data.frame(as.table(knn.cl.cl.counts))
  colnames(knn.cl.df)[1:2] <- c("cl.from", "cl.to")

  from.size <- rowSums(knn.cl.cl.counts)
  to.size <- colSums(knn.cl.cl.counts)
  total <- sum(knn.cl.cl.counts)

  knn.cl.df$cl.from.total <- from.size[as.character(knn.cl.df$cl.from)]
  knn.cl.df$cl.to.total <- to.size[as.character(knn.cl.df$cl.to)]
  knn.cl.df <- knn.cl.df[knn.cl.df$Freq > 0, ]
  knn.cl.df$pval.log <- knn.cl.df$odds <- 0

  for (i in seq_len(nrow(knn.cl.df))) {
    q <- knn.cl.df$Freq[i] - 1
    k_val <- knn.cl.df$cl.from.total[i]
    m <- knn.cl.df$cl.to.total[i]
    n <- total - m
    knn.cl.df$pval.log[i] <- phyper(q, m = m, n = n, k = k_val,
                                     lower.tail = FALSE, log.p = TRUE)
    knn.cl.df$odds[i] <- (q + 1) / (k_val * m / total)
  }

  knn.cl.df$frac <- knn.cl.df$Freq / knn.cl.df$cl.from.total
  knn.cl.df$cl.from.label <- cl.df[as.character(knn.cl.df$cl.from), "cluster_label"]
  knn.cl.df$cl.to.label <- cl.df[as.character(knn.cl.df$cl.to), "cluster_label"]

  result <- list(
    knn.result = knn.result,
    knn.cl.df = knn.cl.df,
    cluster_info = cl.df,
    rd.dat = rd.dat,
    cl = cl,
    params = list(
      cluster_col = cluster_col,
      reduction = reduction,
      k = k
    )
  )

  class(result) <- c("constellation_knn", "list")
  return(result)
}
