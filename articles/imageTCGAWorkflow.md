# Histopathology Image Analysis with TCGA Data: an end-to-end workflow

## Introduction

This workflow demonstrates a complete analysis pipeline for
histopathology image data from The Cancer Genome Atlas (TCGA), using a
suite of Bioconductor packages developed around the `imageTCGA`
ecosystem:

| Package                                                            | Role                                                                 |
|--------------------------------------------------------------------|----------------------------------------------------------------------|
| [imageTCGA](https://github.com/billila/imageTCGA)                  | Interactive Shiny app to explore and select TCGA diagnostic images   |
| [imageFeatureTCGA](https://github.com/waldronlab/imageFeatureTCGA) | Import HoVerNet and Prov-GigaPath features into Bioconductor objects |
| [imageTCGAutils](https://github.com/waldronlab/imageTCGAutils)     | Spatial statistics and dimensionality reduction on tile/cell data    |
| [HistoImagePlot](https://github.com/waldronlab/HistoImagePlot)     | Overlay cell segmentation on tissue thumbnails                       |

The TCGA image database contains ~11,765 diagnostic whole-slide images
(WSI) from ~9,640 patients across many cancer types. Pre-computed
features are available for 54,253 files: 33,177 from HoVerNet (nuclei
segmentation/ classification) and 21,076 from Prov-GigaPath (slide- and
tile-level deep learning embeddings).

### Installation

These packages are available from GitHub during active development and
will be submitted to Bioconductor. Install them as follows:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Bioconductor dependencies
BiocManager::install(c(
    "SpatialExperiment", "SingleCellExperiment", "MultiAssayExperiment",
    "BiocFileCache", "BiocIO", "BumpyMatrix", "GenomicDataCommons"
))

# imageTCGA ecosystem (GitHub)
BiocManager::install(c(
    "billila/imageTCGA",
    "waldronlab/imageFeatureTCGA",
    "waldronlab/imageTCGAutils",
    "waldronlab/HistoImagePlot"
))
```

### Load packages

``` r
library(imageTCGA)
library(imageFeatureTCGA)
library(imageTCGAutils)
library(HistoImagePlot)
library(MultiAssayExperiment)
library(SpatialExperiment)
library(dplyr)
```

------------------------------------------------------------------------

## Step 1: Data Discovery with imageTCGA

The `imageTCGA` Shiny application provides an interactive interface to
the TCGA diagnostic image database. It allows filtering by cancer type,
clinical parameters (histological type, grade, stage), and patient
geography, and generates ready-to-run `GenomicDataCommons` download
code.

``` r
# Launch the Shiny app
imageTCGA::imageTCGA()
```

The app displays:

- **Dot plots** of image counts stratified by clinical variables for
  gynaecological tumours (BRCA, OV, UCS, UCEC).
- **Interactive maps** showing the origin institution of each image at
  country, US-state, and centre level.
- A **code generation panel** that produces the `GenomicDataCommons`
  query corresponding to the current filter selection.

Once you have identified the case(s) of interest, copy the generated
code for use in the next step.

------------------------------------------------------------------------

## Step 2: Retrieving Pre-computed Features with imageFeatureTCGA

`imageFeatureTCGA` provides access to a catalogue of 54,253 pre-computed
feature files hosted on the Waldron Lab server. Two computational
pipelines are covered:

- **HoVerNet** – nuclei detection, segmentation, and classification at
  single-nucleus resolution. Output formats: JSON (coordinates +
  contours), GeoJSON, H5AD (features), and PNG thumbnails.
- **Prov-GigaPath** – vision foundation model embeddings at slide level
  (one vector per WSI) and tile level (one vector per 256 × 256 px
  tile).

### Browse the data catalogue

``` r
library(imageFeatureTCGA)
library(dplyr)

# Full catalogue (54,253 files)
getCatalog()

# HoVerNet files only
getCatalog("hovernet")

# ProvGigaPath files only
getCatalog("provgigapath")
```

Key columns: `pipeline`, `format`, `filename`, `tcga_barcode`,
`Case.ID`, `Project.ID`, `level`.

### Import HoVerNet data

HoVerNet outputs are imported as `SpatialExperiment` where columns are
nuclei and spatial coordinates are centroids in slide pixel space.

``` r
hspe <- getCatalog("hovernet") |>
    dplyr::filter(
        filename == paste(
            "TCGA-VG-A8LO-01A-01-DX1",
            "B39A4D64-82A1-4A04-8AB6-918F3058B83B",
            "json", "gz", sep = "."
        )
    ) |>
    getFileURLs() |>
    HoverNet(outClass = "SpatialExperiment") |>
    import()

hspe
colData(hspe)
```

### Import ProvGigaPath embeddings

``` r
# Slide-level: one 768-dim vector per WSI
getCatalog("provgigapath") |>
    dplyr::filter(
        filename == paste(
            "TCGA-VG-A8LO-01A-01-DX1",
            "B39A4D64-82A1-4A04-8AB6-918F3058B83B",
            "csv", "gz", sep = "."
        ),
        level == "slide_level"
    ) |>
    getFileURLs() |>
    ProvGiga() |>
    import()

# Tile-level: one vector per 256×256 px tile, with tile_x and tile_y coordinates
getCatalog("provgigapath") |>
    dplyr::filter(
        filename == paste(
            "TCGA-VG-A8LO-01A-01-DX1",
            "B39A4D64-82A1-4A04-8AB6-918F3058B83B",
            "csv", "gz", sep = "."
        ),
        level == "tile_level"
    ) |>
    getFileURLs() |>
    ProvGiga() |>
    import()
```

### Import multiple slides

``` r
# Import slide-level embeddings for multiple TCGA-GBM slides
pgl <- getCatalog("provgigapath") |>
    dplyr::filter(level == "slide_level", Project.ID == "TCGA-GBM") |>
    dplyr::slice(1:3) |>
    getFileURLs() |>
    ProvGigaList() |>
    import()

pgl
```

------------------------------------------------------------------------

## Step 3: Spatial Analysis with imageTCGAutils

`imageTCGAutils` provides utility functions for dimensionality reduction
and spatial autocorrelation on tile- and nucleus-level data imported by
`imageFeatureTCGA`. Key capabilities include:

- **PCA** on Prov-GigaPath tile embeddings
- **Global Moran’s I** and **Geary’s C** — spatial autocorrelation
  statistics
- **Local Moran’s I (LISA)** — identify spatially clustered tiles
- **Nucleus-to-tile matching** — assign HoVerNet nuclei to their
  containing Prov-GigaPath tile

See the [imageTCGAutils package
page](https://billila.github.io/imageTCGAWorkflow/articles/pkg-imageTCGAutils.md)
and the [package
documentation](https://github.com/waldronlab/imageTCGAutils) for
function-level usage.

------------------------------------------------------------------------

## Step 4: Visualization with HistoImagePlot

`HistoImagePlot` renders side-by-side panels of the tissue thumbnail
image and HoVerNet nucleus segmentation with coloured cell type labels.
It accepts HoVerNet outputs in JSON and H5AD formats and fetches the
corresponding thumbnail automatically via `BiocFileCache`.

See the [HistoImagePlot package
page](https://billila.github.io/imageTCGAWorkflow/articles/pkg-HistoImagePlot.md)
and the [package
documentation](https://github.com/waldronlab/HistoImagePlot) for
function-level usage.

------------------------------------------------------------------------

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
