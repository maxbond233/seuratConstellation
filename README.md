# seuratConstellation

<p align="center">
  <strong>Constellation Plot Visualization for Seurat V5 Objects</strong>
</p>

A pipeline-style R package for creating constellation plots from Seurat V5 single-cell RNA-seq data. Constellation plots visualize cell type relationships based on K-nearest neighbor (KNN) graph analysis, showing cluster connectivity and relative sizes.

## Acknowledgments

This package is inspired by and adapted from the constellation plot functionality in [**scrattch.hicat**](https://github.com/AllenInstitute/scrattch.hicat) developed by the **Allen Institute for Brain Science**. The original implementation was designed for their internal data structures. This package refactors and extends that functionality to work seamlessly with **Seurat V5** objects and provides a modern, pipeline-friendly API.

**Original source**: https://github.com/AllenInstitute/scrattch.hicat

## Features

- **Seurat V5 Compatible**: Direct input from Seurat objects
- **Pipeline-style API**: Supports `%>%` pipe operations for flexible workflows
- **Self-contained**: Implements KNN prediction internally (no scrattch.hicat dependency)
- **Customizable**: Extensive options for colors, labels, node sizes, and edge styling
- **ggplot2 Output**: Returns ggplot objects for further customization
- **Metadata-based Coloring**: Color nodes by any metadata column (v0.2.0+)
- **Hull Visualization**: Draw convex/concave hulls to group related clusters (v0.2.0+)
- **Feature Expression**: Visualize gene expression on constellation plots (v0.2.0+)

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
  compute_cluster_centers(colors = NULL) %>%
  # Step 4: Create plot
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

### Custom Colors

```r
# Define custom colors for each cluster
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
  compute_cluster_centers(colors = my_colors) %>%
  plot_constellation(label_repel = TRUE)
```

### Metadata-based Coloring (v0.2.0+)

Color nodes based on another metadata column:

```r
# Color by cell_class, group clusters by the same class
p <- ConstellationPlot(
  seu,
  cluster_col = "celltype",
  color_by = "cell_class",
  color_mode = "group"  # or "gradient" for color shades
)
```

### Hull Visualization (v0.2.0+)

Draw background hulls to group related clusters:

```r
p <- ConstellationPlot(
  seu,
  cluster_col = "celltype",
  color_by = "cell_class",
  hull_by = "cell_class",
  hull_type = "concave",  # or "convex"
  hull_alpha = 0.2
)
```

### Gene Expression Visualization (v0.2.0+)

```r
# Color nodes by gene expression
p <- seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  filter_knn_edges(frac_th = 0.05) %>%
  compute_cluster_centers() %>%
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
| `plot_constellation()` | Generate the constellation plot |
| `plot_constellation_feature()` | Visualize gene expression on constellation plot (v0.2.0+) |

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
| `color_by` | character | NULL | Metadata column for node coloring (v0.2.0+) |
| `hull_by` | character | NULL | Metadata column for hull grouping (v0.2.0+) |

#### `filter_knn_edges()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `knn_graph` | constellation_knn | required | Output from `build_knn_graph()` |
| `frac_th` | numeric | 0.05 | Minimum edge fraction to keep |

#### `compute_cluster_centers()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `knn_graph` | constellation_knn | required | Output from previous step |
| `colors` | named vector | NULL | Custom colors; NULL for auto-generation |
| `color_mode` | character | "group" | Color mode: "group" (uniform) or "gradient" (shades) (v0.2.0+) |

#### `plot_constellation()`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `knn_graph` | constellation_knn | required | Output from previous step |
| `node.label` | character | "cluster_label" | Column for node labels |
| `exxageration` | numeric | 2 | Edge width exaggeration factor |
| `curved` | logical | TRUE | Whether edges should be curved |
| `node.dodge` | logical | FALSE | Dodge overlapping nodes |
| `label.size` | numeric | 5 | Size of node labels |
| `max_size` | numeric | 10 | Maximum node size |
| `label_repel` | logical | FALSE | Use ggrepel for non-overlapping labels |
| `node_trans` | character | "sqrt" | Size transformation ("sqrt", "identity", "log10") |
| `hull_type` | character | "convex" | Hull type: "convex" or "concave" (v0.2.0+) |
| `hull_alpha` | numeric | 0.2 | Hull fill transparency (v0.2.0+) |
| `hull_expand` | numeric | 0.1 | Hull expansion factor (v0.2.0+) |

#### `plot_constellation_feature()` (v0.2.0+)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `knn_graph` | constellation_knn | required | Output from previous step |
| `seu` | Seurat | required | Seurat V5 object |
| `feature` | character | required | Gene name to visualize |
| `feature_stat` | character | "mean" | Statistic: "mean", "median", "avg", or "pct" |
| `low_color` | character | "lightgrey" | Color for low expression |
| `high_color` | character | "darkred" | Color for high expression |

## Understanding the Plot

A constellation plot visualizes:

- **Nodes**: Each node represents a cell cluster
  - **Size**: Proportional to the number of cells in the cluster
  - **Color**: User-defined or auto-generated

- **Edges**: Connections between clusters based on KNN relationships
  - **Width**: Proportional to the fraction of KNN connections between clusters
  - **Direction**: Edges can be bidirectional or unidirectional

## Advanced Usage

### Combining with ggplot2

Since the output is a ggplot object, you can further customize it:

```r
p <- ConstellationPlot(seu, cluster_col = "celltype")

# Add title and customize theme
p +
  ggtitle("Cell Type Constellation") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))
```

### Using Different Reductions

```r
# Use tSNE instead of UMAP
p_tsne <- ConstellationPlot(
  seu,
  cluster_col = "celltype",
  reduction = "tsne"
)

# Use PCA
p_pca <- ConstellationPlot(
  seu,
  cluster_col = "celltype",
  reduction = "pca"
)
```

## Citation

If you use this package, please cite:

1. **This package**:
   ```
   seuratConstellation: Constellation Plot Visualization for Seurat V5 Objects
   https://github.com/username/seuratConstellation
   ```

2. **Original scrattch.hicat package**:
   ```
   Tasic B, et al. (2018). Shared and distinct transcriptomic cell types across
   neocortical areas. Nature, 563(7729), 72-78.
   
   GitHub: https://github.com/AllenInstitute/scrattch.hicat
   ```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

If you encounter any problems or have suggestions, please open an issue on GitHub.
