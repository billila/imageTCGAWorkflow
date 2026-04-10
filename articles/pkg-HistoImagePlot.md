# Package: HistoImagePlot

## Overview

**GitHub:** <https://github.com/waldronlab/HistoImagePlot>

`HistoImagePlot` provides visualization tools for histopathology
research, working alongside `imageFeatureTCGA` to import and display
HoVerNet nuclei segmentation results overlaid on tissue thumbnail
images.

Supported HoVerNet input formats:

| Format | Extension            | Content                                                                          |
|--------|----------------------|----------------------------------------------------------------------------------|
| H5AD   | `.h5ad` / `.h5ad.gz` | Feature matrix + spatial statistics (mean intensity, nearest-neighbour distance) |
| JSON   | `.json` / `.json.gz` | Nucleus centroids, contours, type, probabilities                                 |

``` r
library(imageFeatureTCGA)
library(HistoImagePlot)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
library(cowplot)
```

## Import HoVerNet H5AD Data

Import a HoVerNet H5AD file into a `SpatialExperiment` object. The file
is automatically cached via `BiocFileCache`:

``` r
hov_file <- paste0(
    "https://store.cancerdatasci.org/hovernet/h5ad/",
    "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D.h5ad.gz"
)

hn_spe <- HoverNet(hov_file, outClass = "SpatialExperiment") |>
    import()

hn_spe
```

The corresponding tissue thumbnail is also available from the same
store:

``` r
thumb_path <- paste0(
    "https://store.cancerdatasci.org/hovernet/thumb/",
    "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D.png"
)
```

## Overlaying Segmentation on Tissue Thumbnails

### Basic Overlay

`plotHoverNetH5ADOverlay()` renders nucleus centroids coloured by cell
type on top of the tissue thumbnail:

``` r
plotHoverNetH5ADOverlay(hn_spe, thumb_path)
```

### Customised Overlay

Point size and legend appearance can be adjusted:

``` r
plotHoverNetH5ADOverlay(
    hn_spe,
    thumb_path,
    title             = "Ovarian Cancer Tissue - Cell Segmentation",
    point_size        = 0.02,
    legend_point_size = 3
)
```

### Custom Colour Palette

Pass a named character vector to override the default cell-type colours:

``` r
custom_colors <- c(
    "no label"          = "#808080",
    "neoplastic"        = "#E31A1C",
    "inflammatory"      = "#1F78B4",
    "stromal"           = "#33A02C",
    "necrotic"          = "#FF7F00",
    "benign epithelial" = "#6A3D9A"
)

plotHoverNetH5ADOverlay(
    hn_spe,
    thumb_path,
    color_palette = custom_colors,
    title         = "Custom Color Scheme"
)
```

## Visualising H5AD Features

H5AD files include additional per-nucleus statistics stored in
`colData()`: `mean_intensity` and `nearest_neighbor_distance`. Spatial
coordinates are available via `spatialCoords()` as `x_centroid` and
`y_centroid`.

``` r
h5ad_coords <- data.frame(spatialCoords(hn_spe), colData(hn_spe))

p1 <- ggplot(h5ad_coords, aes(x = x_centroid, y = y_centroid,
                               color = mean_intensity)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Mean Intensity", color = "Intensity")

p2 <- ggplot(h5ad_coords, aes(x = x_centroid, y = y_centroid,
                               color = nearest_neighbor_distance)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Nearest Neighbor Distance", color = "Distance")

cowplot::plot_grid(p1, p2, ncol = 2)
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
