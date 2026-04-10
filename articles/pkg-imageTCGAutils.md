# Package: imageTCGAutils

## Overview

**GitHub:** <https://github.com/waldronlab/imageTCGAutils>

`imageTCGAutils` provides utility functions for integrating and
analyzing multi-modal whole-slide image (WSI) data from TCGA. It is
designed to work alongside `imageFeatureTCGA`, which handles import of
pre-computed features from HoVerNet and Prov-GigaPath.

``` r
library(imageFeatureTCGA)
library(imageTCGAutils)
library(ggplot2)
library(dplyr)
library(sfdep)
library(spdep)
library(SpatialExperiment)
library(data.table)
```

## Import Prov-GigaPath Tile-Level Embeddings

``` r
# Select an ovarian cancer slide as example
tile_prov_url <- paste0(
    "https://store.cancerdatasci.org/provgigapath/tile_level/",
    "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D.csv.gz"
)

example_slide <- ProvGiga(tile_prov_url) |>
    import()
```

## PCA on Tile Embeddings

``` r
# Extract embedding columns (named with numbers: "0", "1", ...)
embedding_cols <- grep("^[0-9]+$", names(example_slide), value = TRUE)

# Run PCA
pca_res <- prcomp(example_slide[, embedding_cols], scale. = TRUE)

pca_example_slide <- bind_cols(
    example_slide,
    as_tibble(pca_res$x)[, 1:2]
)

# PCA scatter plot
ggplot(pca_example_slide, aes(PC1, PC2)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_minimal() +
    labs(title = "Tile-level PCA Ovarian Cancer Embedding: Single Slide")

# Spatial layout coloured by PC1
ggplot(pca_example_slide, aes(tile_x, tile_y, color = PC1)) +
    geom_point(size = 1) +
    scale_color_viridis_c() +
    coord_equal() +
    theme_minimal() +
    labs(title = "Tissue layout colored by PC1")
```

## Spatial Autocorrelation

To investigate spatial patterns, we build a k-nearest neighbour graph
from tile coordinates and compute global and local spatial
autocorrelation on the PCA-reduced embeddings.

``` r
coords <- pca_example_slide[, c("tile_x", "tile_y")]
nb <- knn2nb(knearneigh(coords, k = 6))
lw <- nb2listw(nb, style = "W")
```

### Global Moran’s I and Geary’s C

``` r
mi <- moran.test(pca_example_slide$PC1, lw)
gc <- geary.test(pca_example_slide$PC1, lw)

mi
gc
```

Moran’s I ∈ \[−1, 1\]: positive values indicate spatial clustering.
Geary’s C ∈ \[0, 2\]: values \< 1 indicate positive spatial
autocorrelation.

### Local Moran’s I (LISA)

``` r
lisa <- localmoran(pca_example_slide$PC1, lw)
pca_example_slide$localI      <- lisa[, "Ii"]
pca_example_slide$localI_pval <- lisa[, "Pr(z != E(Ii))"]

# Moran scatterplot
moran.plot(pca_example_slide$PC1, lw, labels = FALSE,
           main = "Moran scatterplot of PC1")

# LISA spatial map
ggplot(pca_example_slide, aes(x = tile_x, y = tile_y, color = localI)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c() +
    coord_equal() +
    theme_minimal() +
    ggtitle("Local Moran's I (LISA) for PC1")
```

## Integration with HoVerNet Nuclei

`matchHoverNetToTiles()` computes the scaling factor between HoVerNet
nucleus coordinates and Prov-GigaPath tile coordinates, then assigns
each tile a dominant cell type and nucleus count based on the
overlapping nuclei.

``` r
hov_file <- paste0(
    "https://store.cancerdatasci.org/hovernet/h5ad/",
    "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D.h5ad.gz"
)

hn_spe <- HoverNet(hov_file, outClass = "SpatialExperiment") |>
    import()
```

``` r
# Visualise nucleus vs tile coordinates — they use different scales
cell_coords <- spatialCoords(hn_spe)

plot(cell_coords[, 1], cell_coords[, 2],
     pch = 16, col = "#0000FF20",
     main = "HoVerNet nuclei (blue) vs GigaPath tiles (red)")
points(pca_example_slide$tile_x, pca_example_slide$tile_y,
       pch = 16, col = "#FF000020")
```

``` r
pg_spe <- ProvGiga(tile_prov_url) |> import()

match_hv_pg <- matchHoverNetToTiles(hn_spe, pg_spe)

# Dominant cell type per tile
ggplot(match_hv_pg$tiles_with_nuclei,
       aes(tile_x, tile_y, color = dominant_cell_type)) +
    geom_point(size = 2) +
    coord_equal() +
    theme_minimal() +
    labs(title = "Per-tile dominant HoVerNet cell type")

# All cell types per tile (sized by count)
ggplot(match_hv_pg$tiles_with_nuclei,
       aes(tile_x, tile_y, color = cell_type_label, size = N)) +
    geom_point(alpha = 0.7) +
    coord_equal() +
    theme_minimal() +
    labs(title = "All HoVerNet cell types per tile")
```

## Session Information

``` r
sessioninfo::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.5.3 (2026-03-11)
#>  os       Ubuntu 24.04.4 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language en
#>  collate  C.UTF-8
#>  ctype    C.UTF-8
#>  tz       UTC
#>  date     2026-04-10
#>  pandoc   3.1.11 @ /opt/hostedtoolcache/pandoc/3.1.11/x64/ (via rmarkdown)
#>  quarto   NA
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package     * version date (UTC) lib source
#>  BiocManager   1.30.27 2025-11-14 [1] RSPM (R 4.5.0)
#>  BiocStyle   * 2.38.0  2025-10-29 [1] Bioconductor 3.22 (R 4.5.3)
#>  bookdown      0.46    2025-12-05 [1] RSPM (R 4.5.0)
#>  bslib         0.10.0  2026-01-26 [1] RSPM (R 4.5.0)
#>  cachem        1.1.0   2024-05-16 [1] RSPM (R 4.5.0)
#>  cli           3.6.6   2026-04-09 [1] RSPM (R 4.5.0)
#>  desc          1.4.3   2023-12-10 [1] RSPM (R 4.5.0)
#>  digest        0.6.39  2025-11-19 [1] RSPM (R 4.5.0)
#>  evaluate      1.0.5   2025-08-27 [1] RSPM (R 4.5.0)
#>  fastmap       1.2.0   2024-05-15 [1] RSPM (R 4.5.0)
#>  fs            2.0.1   2026-03-24 [1] RSPM (R 4.5.0)
#>  htmltools     0.5.9   2025-12-04 [1] RSPM (R 4.5.0)
#>  jquerylib     0.1.4   2021-04-26 [1] RSPM (R 4.5.0)
#>  jsonlite      2.0.0   2025-03-27 [1] RSPM (R 4.5.0)
#>  knitr         1.51    2025-12-20 [1] RSPM (R 4.5.0)
#>  lifecycle     1.0.5   2026-01-08 [1] RSPM (R 4.5.0)
#>  pkgdown       2.2.0   2025-11-06 [1] RSPM (R 4.5.0)
#>  R6            2.6.1   2025-02-15 [1] RSPM (R 4.5.0)
#>  ragg          1.5.2   2026-03-23 [1] RSPM (R 4.5.0)
#>  rlang         1.2.0   2026-04-06 [1] RSPM (R 4.5.0)
#>  rmarkdown     2.31    2026-03-26 [1] RSPM (R 4.5.0)
#>  sass          0.4.10  2025-04-11 [1] RSPM (R 4.5.0)
#>  sessioninfo   1.2.3   2025-02-05 [1] RSPM (R 4.5.0)
#>  systemfonts   1.3.2   2026-03-05 [1] RSPM (R 4.5.0)
#>  textshaping   1.0.5   2026-03-06 [1] RSPM (R 4.5.0)
#>  xfun          0.57    2026-03-20 [1] RSPM (R 4.5.0)
#>  yaml          2.3.12  2025-12-10 [1] RSPM (R 4.5.0)
#> 
#>  [1] /home/runner/work/_temp/Library
#>  [2] /opt/R/4.5.3/lib/R/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
