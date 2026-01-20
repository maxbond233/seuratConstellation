# Refactor Plan: Reduce Coupling in Constellation Plot Pipeline

## Purpose
The current pipeline has high coupling and hidden parameter dependencies across
`build_knn_graph()`, `compute_cluster_centers()`, and `plot_constellation()`.
This document proposes a refactor that separates data preparation, grouping,
color resolution, and plotting into explicit, testable steps while preserving
backward compatibility via `ConstellationPlot()`.

## Current Pain Points
- Hidden precedence: `compute_cluster_centers(colors=...)` silently overrides
  `color_by` derived grouping.
- `build_knn_graph()` injects `color_group` and `hull_group` into cluster info,
  creating implicit dependencies for later steps.
- `plot_constellation()` assumes fields (`cluster_color`, `hull_group`) exist,
  forcing a strict call order.
- Parameter interactions are implicit and difficult to debug for users.

## Goals
- Make dependencies explicit and local to each function.
- Separate geometry from styling and from metadata grouping.
- Allow independent use of each function with clear inputs/outputs.
- Preserve the current all-in-one workflow (with deprecation warnings).

## Proposed Architecture

### 1) Data Layer: Build graph only
**Function:** `build_knn_graph()`
- Input: Seurat object, `cluster_col`, `reduction`, `k`, outlier thresholds.
- Output: `knn.cl.df`, `cluster_info`, `rd.dat`, `cl`.
- Remove `color_by` and `hull_by` logic from this function (or keep and warn).

### 2) Grouping Layer: Map clusters to groups
**New function:** `derive_cluster_groups()`
- Input:
  - `cluster_info` (cluster labels and sizes)
  - `cluster_col` vector (cluster assignment per cell)
  - `group_by` vector (metadata per cell)
  - `strategy`: `"majority"` (default) or `"strict"` or `"manual"`
- Output: `cluster_group` mapping merged into `cluster_info`.

### 3) Color Layer: Resolve final cluster colors
**New function:** `resolve_cluster_colors()`
- Input:
  - `cluster_info`
  - `cluster_colors` (named by cluster_label, optional)
  - `group_colors` (named by group label, optional)
  - `color_mode`: `"group"` or `"gradient"`
- Output: `cluster_color` vector aligned to `cluster_label`.
- Explicit precedence:
  1. `cluster_colors` if provided and complete.
  2. Else derive from `group_colors` + `color_mode`.
  3. Else fallback to auto palette per cluster.

### 4) Geometry Layer: Cluster centers only
**Function:** `compute_cluster_centers()`
- Input: `knn_graph`.
- Output: `cl.center.df` with coordinates and cluster sizes.
- Remove all color logic from this function.

### 5) Plot Layer: Consume explicit fields
**Function:** `plot_constellation()`
- Require `cluster_color` and optional `hull_group` already present in `nodes`.
- `plot_constellation()` no longer infers or computes colors.

### 6) Hull Layer: Independent hull computation
**Function:** `compute_hulls()` (existing internal can stay)
- Input: `nodes` + `hull_group` + `hull_type`, `hull_expand`.
- Output: `hull_data` used by `plot_constellation()`.

## Backward Compatibility Plan
**Keep** `ConstellationPlot()` as one-stop entry:
1. `build_knn_graph()`
2. `compute_cluster_centers()`
3. If `color_by` or `hull_by` provided:
   - call `derive_cluster_groups()`
4. Call `resolve_cluster_colors()`
5. Call `plot_constellation()`

Deprecate (with warnings):
- `build_knn_graph(color_by=..., hull_by=...)`
- `compute_cluster_centers(colors=..., color_mode=...)`

## Suggested Function Signatures (Draft)
```r
derive_cluster_groups <- function(cluster_vec, group_vec, cluster_levels,
                                  strategy = "majority") { ... }

resolve_cluster_colors <- function(cluster_info,
                                   cluster_colors = NULL,
                                   group_colors = NULL,
                                   color_mode = "group") { ... }
```

## Benefits
- Clear separation of concerns.
- Fewer implicit side effects.
- Easier testing and debugging.
- Users can customize colors/hulls without breaking call order.

## Decisions Needed
1. Preserve 100% API compatibility, or allow mild breaking changes?
2. Should `color_by/hull_by` remain in `ConstellationPlot()` only, or remain
   available in lower-level calls with deprecation warnings?

