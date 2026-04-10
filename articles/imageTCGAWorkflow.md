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
imageTCGA::run_app()
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
# Inspect the full data catalogue
cat <- imageFeatureTCGA::imageCatalogue()
head(cat)
```

The catalogue columns include `file_id`, `cases.submitter_id` (TCGA
barcode), `cancer_type`, `pipeline`, `data_type`, and the remote `url`.

``` r
# Filter to HoVerNet H5AD files for ovarian cancer
ov_h5ad <- cat |>
    dplyr::filter(cancer_type == "OV",
                  pipeline == "HoVerNet",
                  data_type == "h5ad")

nrow(ov_h5ad)
head(ov_h5ad[, c("cases.submitter_id", "pipeline", "data_type")])
```

### Import into a SpatialExperiment (HoVerNet)

HoVerNet outputs are imported as `SpatialExperiment` objects where each
row is a nucleus and spatial coordinates reflect its centroid on the
slide.

``` r
# Import HoVerNet H5AD data for a single slide
spe_hover <- imageFeatureTCGA::importHoVerNet(
    file_id  = ov_h5ad$file_id[1],
    format   = "h5ad"
)
spe_hover
```

``` r
# Import HoVerNet JSON (includes contours) for the same slide
spe_json <- imageFeatureTCGA::importHoVerNet(
    file_id  = ov_h5ad$file_id[1],
    format   = "json"
)
```

The resulting `SpatialExperiment` contains: - `assay("X")` – per-nucleus
feature matrix (intensity statistics, morphology) - `spatialCoords()` –
centroid (x, y) in slide pixel coordinates - `colData()` – cell type
label, nucleus area, and QC flags

### Import tile embeddings (Prov-GigaPath)

``` r
# Filter catalogue to tile-level Prov-GigaPath files for one OV slide
pg_tiles <- cat |>
    dplyr::filter(cases.submitter_id == ov_h5ad$cases.submitter_id[1],
                  pipeline == "ProvGigaPath",
                  data_type == "tile_embeddings")

spe_tiles <- imageFeatureTCGA::importProvGigaPath(
    file_id  = pg_tiles$file_id[1],
    level    = "tile"
)
spe_tiles
```

Each column is a 256 × 256 px tile; the 1,536-dimensional embedding
vector from Prov-GigaPath is stored in `assay("embeddings")`.

### Assemble a MultiAssayExperiment

Both assays can be linked to TCGA clinical metadata in a single
`MultiAssayExperiment`:

``` r
mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(
        HoVerNet    = spe_hover,
        GigaPath    = spe_tiles
    )
)
mae
```

------------------------------------------------------------------------

## Step 3: Spatial Analysis with imageTCGAutils

`imageTCGAutils` provides utility functions for dimensionality reduction
and spatial autocorrelation on the tile- and nucleus-level data.

### PCA on tile embeddings

``` r
# Run PCA on Prov-GigaPath tile embeddings and store result in reducedDim
spe_tiles <- imageTCGAutils::runTilePCA(spe_tiles, ncomponents = 30)

# Visualise first two PCs coloured by tissue region
imageTCGAutils::plotTileReducedDim(spe_tiles, dimred = "PCA",
                                    colour_by = "tissue_class")
```

### Spatial autocorrelation

Moran’s I and Geary’s C quantify whether a feature (e.g. proportion of
tumour nuclei in a tile) is spatially clustered, dispersed, or random.
Local Moran’s I (LISA) identifies individual tiles that drive global
clustering.

``` r
# Compute global Moran's I for the tumour-nucleus proportion per tile
mi <- imageTCGAutils::globalMoransI(
    spe   = spe_tiles,
    assay = "embeddings",
    feat  = "tumour_prop"
)
mi

# Local Moran's I (LISA) — returns a SpatialExperiment with LISA scores
spe_tiles <- imageTCGAutils::localMoransI(spe_tiles, feat = "tumour_prop")

# Geary's C
gc <- imageTCGAutils::globalGearysC(spe_tiles, feat = "tumour_prop")
gc
```

### Matching HoVerNet nuclei to Prov-GigaPath tiles

``` r
# Assign each nucleus (spe_hover) to its containing tile (spe_tiles)
spe_hover <- imageTCGAutils::matchNucleiToTiles(
    nuclei = spe_hover,
    tiles  = spe_tiles
)
# colData(spe_hover)$tile_id now contains the parent tile identifier
```

------------------------------------------------------------------------

## Step 4: Visualization with HistoImagePlot

`HistoImagePlot` renders side-by-side panels of the tissue thumbnail and
nucleus segmentation overlay, useful for quality control and for
illustrating findings.

### Side-by-side thumbnail + segmentation

``` r
# Import HoVerNet JSON for the slide (includes nucleus contours)
spe_json <- HistoImagePlot::importHoVerNetJSON(
    file_id = ov_h5ad$file_id[1]
)

# Plot thumbnail next to coloured segmentation map
HistoImagePlot::plotHistoImage(
    spe        = spe_json,
    colour_by  = "cell_type",
    point_size = 0.4,
    title      = paste("OV –", ov_h5ad$cases.submitter_id[1])
)
```

Supported colour palettes can be passed via `palette =`. The thumbnail
is fetched from the catalogue and cached automatically by
`BiocFileCache`.

### Overlay on a region of interest

``` r
# Restrict to a 2000 × 2000 px region starting at (5000, 3000)
HistoImagePlot::plotHistoImage(
    spe       = spe_json,
    colour_by = "cell_type",
    xlim      = c(5000, 7000),
    ylim      = c(3000, 5000)
)
```

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
