#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "seuratConstellation v", utils::packageVersion(pkgname), "\n",
    "Constellation plots for Seurat V5 objects"
  )
}

.onLoad <- function(libname, pkgname) {
  # Set default options if needed
  op <- options()
  op.constellation <- list(
    seuratConstellation.verbose = TRUE
  )
  toset <- !(names(op.constellation) %in% names(op))

  if (any(toset)) options(op.constellation[toset])

  invisible()
}

# Global variables to avoid R CMD check notes
utils::globalVariables(c(

  "x", "y", "cl", "cluster_size", "cluster_color",
  "cluster_label", "cluster_id", "node.width",
  "frac", "cl.from", "cl.to", "Group", ".data"
))
