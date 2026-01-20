# seuratConstellation

<p align="center">
  <strong>Constellation Plot Visualization for Seurat V5 Objects</strong>
</p>

A pipeline-style R package for creating constellation plots from Seurat V5 single-cell RNA-seq data. Constellation plots visualize cell type relationships based on K-nearest neighbor (KNN) graph analysis, showing cluster connectivity and relative sizes.

## Acknowledgments

This package is inspired by and adapted from the constellation plot functionality in [**scrattch.hicat**](https://github.com/AllenInstitute/scrattch.hicat) developed by the **Allen Institute for Brain Science**. The original implementation was designed for their internal data structures. This package refactors and extends that functionality to work seamlessly with **Seurat V5** objects and provides a modern, pipeline-friendly API.

**Original source**: https://github.com/AllenInstitute/scrattch.hicat

AI contributors (tooling support): Codex CLI, Gemini CLI.

## Features (v0.2.1)

- **Seurat V5 Compatible**: Direct input from Seurat objects
- **Pipeline-style API**: Supports `%>%` pipe operations for flexible workflows
- **Self-contained**: Implements KNN prediction internally (no scrattch.hicat dependency)
- **Customizable**: Extensive options for colors, labels, node sizes, and edge styling
- **ggplot2 Output**: Returns ggplot objects for further customization
- **Explicit Grouping**: Separate group assignment from graph building
- **Explicit Color Resolution**: Predictable color precedence rules
- **Hull Visualization**: Draw convex/concave hulls to group related clusters
- **Feature Expression**: Visualize gene expression on constellation plots

## Installation

### From GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install seuratConstellation
devtools::install_github("maxbond233/seuratConstellation")
```

### From Local Source

```r
devtools::install("/path/to/seuratConstellation")

# Or load for development
devtools::load_all("/path/to/seuratConstellation")
```

## Dependencies

- R (>= 4.0.0)
- Seurat (>= 5.0.0)
- ggplot2
- ggrepel
- dplyr
- RANN
- Hmisc
- scales
- grDevices

**Optional** (for concave hulls):
- concaveman

## Quick Start

```r
library(seuratConstellation)
library(Seurat)

# Load your Seurat object
seu <- readRDS("your_seurat_object.rds")

# One-line constellation plot
p <- ConstellationPlot(seu, cluster_col = "celltype")
print(p)

# Save the plot
ggsave("constellation_plot.pdf", p, width = 12, height = 12)
```

## Usage

### Method 1: Pipeline Style (Recommended)

The pipeline approach gives you full control over each step:

```r
library(seuratConstellation)
library(dplyr)

p <- seu %>%
  # Step 1: Build KNN graph
  build_knn_graph(
    cluster_col = "celltype",
    reduction = "umap",
    k = 15
  ) %>%
  # Step 2: Filter weak edges
  filter_knn_edges(frac_th = 0.05) %>%
  # Step 3: Compute cluster centers
  compute_cluster_centers() %>%
  # Step 4: Resolve colors
  resolve_cluster_colors() %>%
  # Step 5: Create plot
  plot_constellation(
    node.label = "cluster_label",
    label_repel = TRUE,
    node.dodge = TRUE
  )
```

### Method 2: All-in-One Function

For quick visualization with sensible defaults:

```r
p <- ConstellationPlot(
  seu,
  cluster_col = "celltype",
  reduction = "umap",
  k = 15,
  frac_th = 0.05,
  label_repel = TRUE,
  node.dodge = TRUE
)
```

### Custom Cluster Colors

```r
# Define custom colors for each cluster (names must match cluster_label)
my_colors <- c(
  "Excitatory" = "#E41A1C",
  "Inhibitory" = "#377EB8",
  "Astrocyte" = "#4DAF4A",
  "Microglia" = "#984EA3",
  "Oligodendrocyte" = "#FF7F00"
)

p <- seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  filter_knn_edges(frac_th = 0.05) %>%
  compute_cluster_centers() %>%
  resolve_cluster_colors(cluster_colors = my_colors) %>%
  plot_constellation(label_repel = TRUE)
```

### Group-Based Coloring

Color nodes based on another metadata column:

```r
p <- seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  filter_knn_edges(frac_th = 0.05) %>%
  compute_cluster_centers() %>%
  assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
  resolve_cluster_colors(color_mode = "group") %>%
  plot_constellation(label_repel = TRUE)
```

Use `color_mode = "gradient"` to generate shades within each group.

### Hull Visualization

Draw background hulls to group related clusters:

```r
p <- seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  filter_knn_edges(frac_th = 0.05) %>%
  compute_cluster_centers() %>%
  assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
  assign_cluster_groups(seu, group_by = "cell_class", group_type = "hull") %>%
  resolve_cluster_colors(color_mode = "group") %>%
  plot_constellation(hull_type = "concave", hull_alpha = 0.2)
```

### Gene Expression Visualization

```r
# Color nodes by gene expression
p <- seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  filter_knn_edges(frac_th = 0.05) %>%
  compute_cluster_centers() %>%
  resolve_cluster_colors() %>%
  plot_constellation_feature(seu, feature = "MS4A1")

# Different statistics
plot_constellation_feature(seu, feature = "CD3E", feature_stat = "avg")
plot_constellation_feature(seu, feature = "CD3E", feature_stat = "pct")
```

## Function Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `ConstellationPlot()` | All-in-one convenience function |
| `build_knn_graph()` | Build KNN graph from Seurat object |
| `filter_knn_edges()` | Filter edges by fraction threshold |
| `compute_cluster_centers()` | Calculate cluster centroid coordinates |
| `assign_cluster_groups()` | Map clusters to metadata groups |
| `resolve_cluster_colors()` | Resolve cluster colors with explicit precedence |
| `plot_constellation()` | Generate the constellation plot |
| `plot_constellation_feature()` | Visualize gene expression on constellation plot |

### Key Parameters

#### `build_knn_graph()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seu` | Seurat | required | Seurat V5 object |
| `cluster_col` | character | required | Column name in meta.data for clusters |
| `reduction` | character | "umap" | Dimensional reduction to use |
| `k` | integer | 15 | Number of nearest neighbors |
| `knn.outlier.th` | numeric | 2 | Outlier detection threshold |
| `outlier.frac.th` | numeric | 0.5 | Fraction threshold for outlier cells |
| `color_by` | character | NULL | Deprecated; use `assign_cluster_groups()` |
| `hull_by` | character | NULL | Deprecated; use `assign_cluster_groups()` |

#### `assign_cluster_groups()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seu` | Seurat | required | Seurat V5 object (same one used in build) |
| `group_by` | character | required | Metadata column to group by |
| `group_type` | character | "color" | Grouping target: "color" or "hull" |
| `strategy` | character | "majority" | Grouping strategy: "majority", "strict", or "manual" |

#### `resolve_cluster_colors()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_colors` | named vector | NULL | Colors by cluster_label |
| `group_colors` | named vector | NULL | Colors by group label |
| `color_mode` | character | "group" | "group" or "gradient" |

## Migration Notes (v0.2.1)

- `compute_cluster_centers(colors = ..., color_mode = ...)` is deprecated. Use `resolve_cluster_colors()` instead.
- `build_knn_graph(color_by = ..., hull_by = ...)` is deprecated. Use `assign_cluster_groups()` instead.
- `plot_constellation()` now requires colors to be resolved beforehand with `resolve_cluster_colors()`.

## Changelog (short)

- v0.2.1: Split grouping and color resolution into explicit steps (`assign_cluster_groups()`, `resolve_cluster_colors()`), added clearer pipeline examples, and documented deprecations.
