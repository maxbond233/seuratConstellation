# Migration Guide: v0.1.x to v0.2.x

## Overview

Version 0.2.0 introduces a major refactoring to reduce coupling and make parameter dependencies explicit. This document guides you through migrating your code to the new architecture.

## What Changed?

### Key Changes

1. **Separated color resolution from geometry computation**
   - `compute_cluster_centers()` now only computes coordinates
   - New function `resolve_cluster_colors()` handles all color logic

2. **Explicit grouping logic**
   - New function `assign_cluster_groups()` for metadata-based grouping
   - Clearer precedence rules for color assignment

3. **Deprecated parameters**
   - `build_knn_graph(color_by=..., hull_by=...)` → Use `assign_cluster_groups()`
   - `compute_cluster_centers(colors=..., color_mode=...)` → Use `resolve_cluster_colors()`

## Migration Examples

### Example 1: Basic Usage (No Changes Needed)

**Old code:**
```r
ConstellationPlot(seu, cluster_col = "celltype")
```

**New code:**
```r
# Still works exactly the same!
ConstellationPlot(seu, cluster_col = "celltype")
```

**Status:** ✅ No changes required

---

### Example 2: Custom Cluster Colors

**Old code:**
```r
ConstellationPlot(
  seu,
  cluster_col = "celltype",
  colors = c("T cell" = "red", "B cell" = "blue")
)
```

**New code:**
```r
# Option 1: Still works with ConstellationPlot (recommended for simple cases)
ConstellationPlot(
  seu,
  cluster_col = "celltype",
  colors = c("T cell" = "red", "B cell" = "blue")
)

# Option 2: Explicit pipeline (recommended for complex cases)
seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  filter_knn_edges(frac_th = 0.05) %>%
  compute_cluster_centers() %>%
  resolve_cluster_colors(
    cluster_colors = c("T cell" = "red", "B cell" = "blue")
  ) %>%
  plot_constellation()
```

**Status:** ✅ No changes required (but new explicit option available)

---

### Example 3: Group-Based Coloring

**Old code:**
```r
ConstellationPlot(
  seu,
  cluster_col = "celltype",
  color_by = "cell_class",
  color_mode = "gradient"
)
```

**New code (Option 1: Quick fix):**
```r
# Still works but shows deprecation warning
ConstellationPlot(
  seu,
  cluster_col = "celltype",
  color_by = "cell_class",
  color_mode = "gradient"
)
```

**New code (Option 2: Recommended):**
```r
# Explicit pipeline (clearer and more flexible)
ConstellationPlot(
  seu,
  cluster_col = "celltype",
  color_by = "cell_class",
  color_mode = "gradient"
)

# Or use the step-by-step approach:
knn_graph <- build_knn_graph(seu, cluster_col = "celltype")
knn_graph <- filter_knn_edges(knn_graph, frac_th = 0.05)
knn_graph <- compute_cluster_centers(knn_graph)
knn_graph <- assign_cluster_groups(knn_graph, seu,
                                   group_by = "cell_class",
                                   group_type = "color")
knn_graph <- resolve_cluster_colors(knn_graph, color_mode = "gradient")
plot_constellation(knn_graph)
```

**Status:** ⚠️ Deprecation warning (still works)

---

### Example 4: Hulls with Group Coloring

**Old code:**
```r
ConstellationPlot(
  seu,
  cluster_col = "celltype",
  color_by = "cell_class",
  hull_by = "tissue",
  hull_type = "concave"
)
```

**New code (Option 1: Quick fix):**
```r
# Still works but shows deprecation warnings
ConstellationPlot(
  seu,
  cluster_col = "celltype",
  color_by = "cell_class",
  hull_by = "tissue",
  hull_type = "concave"
)
```

**New code (Option 2: Recommended):**
```r
# Explicit pipeline showing all steps
seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  filter_knn_edges(frac_th = 0.05) %>%
  compute_cluster_centers() %>%
  assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
  assign_cluster_groups(seu, group_by = "tissue", group_type = "hull") %>%
  resolve_cluster_colors(color_mode = "group") %>%
  plot_constellation(hull_type = "concave")
```

**Status:** ⚠️ Deprecation warnings (still works)

---

### Example 5: Custom Group Colors

**Old code:**
```r
# This was NOT directly supported before
# Users had to manually create cluster_colors from groups
```

**New code:**
```r
# Now directly supported!
seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  compute_cluster_centers() %>%
  assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
  resolve_cluster_colors(
    group_colors = c("Neuron" = "#E41A1C", "Glia" = "#377EB8"),
    color_mode = "group"
  ) %>%
  plot_constellation()
```

**Status:** ✨ New feature!

---

## New Features in v0.2.0

### 1. Explicit Color Precedence

Colors are now resolved with clear precedence:

```r
resolve_cluster_colors(
  knn_graph,
  cluster_colors = c("A" = "red", "B" = "blue"),  # Priority 1: Direct cluster colors
  group_colors = c("Group1" = "green"),           # Priority 2: Group-level colors
  color_mode = "gradient"                         # Priority 3: Auto-generate from groups
)
```

### 2. Multiple Grouping Strategies

```r
# Majority voting (default)
assign_cluster_groups(knn_graph, seu, group_by = "cell_class", strategy = "majority")

# Strict mode (>90% purity required)
assign_cluster_groups(knn_graph, seu, group_by = "cell_class", strategy = "strict")

# Manual mapping
assign_cluster_groups(
  knn_graph, seu,
  group_by = "cell_class",
  strategy = "manual",
  manual_mapping = c("Cluster1" = "GroupA", "Cluster2" = "GroupA")
)
```

### 3. Independent Color and Hull Groups

```r
# Color by cell class, hull by tissue
knn_graph %>%
  assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
  assign_cluster_groups(seu, group_by = "tissue", group_type = "hull") %>%
  resolve_cluster_colors() %>%
  plot_constellation(hull_type = "concave")
```

## Deprecation Timeline

| Version | Status | Action |
|---------|--------|--------|
| v0.2.x | Deprecated with warnings | Old parameters still work but show warnings |
| v0.3.x | Deprecated with warnings | Warnings become more prominent |
| v0.4.0 | Removed | Old parameters will be removed |

## Benefits of the New Architecture

1. **Explicit Dependencies**: No more hidden parameter interactions
2. **Better Testability**: Each function has a single, clear responsibility
3. **More Flexible**: Mix and match grouping strategies and color modes
4. **Easier Debugging**: Clear error messages guide you to the right solution
5. **Backward Compatible**: Existing code continues to work (with warnings)

## Getting Help

If you encounter issues during migration:

1. Check this migration guide
2. Read the function documentation: `?resolve_cluster_colors`
3. Report issues at: https://github.com/your-repo/seuratConstellation/issues

## Quick Reference

### Old Workflow
```r
ConstellationPlot(seu, cluster_col = "celltype", color_by = "cell_class")
```

### New Workflow (Explicit)
```r
seu %>%
  build_knn_graph(cluster_col = "celltype") %>%
  compute_cluster_centers() %>%
  assign_cluster_groups(seu, group_by = "cell_class", group_type = "color") %>%
  resolve_cluster_colors() %>%
  plot_constellation()
```

### Convenience Workflow (Still Works)
```r
ConstellationPlot(seu, cluster_col = "celltype", color_by = "cell_class")
# Shows deprecation warning but works fine
```
