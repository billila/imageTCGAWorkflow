# Package: imageFeatureTCGA

## Overview

**GitHub:** <https://github.com/waldronlab/imageFeatureTCGA>

`imageFeatureTCGA` provides convenient access to histopathology-derived
features from TCGA through two complementary pipelines:

- **HoVerNet** — cell segmentation and classification at single-nucleus
  level
- **ProvGigaPath** — slide- and tile-level deep learning embeddings

The data catalogue contains **54,253 pre-computed files**: 33,177
HoVerNet files and 21,076 ProvGigaPath files. All files can be imported
directly into Bioconductor objects without requiring local image files
or GPU resources.

``` r
library(imageFeatureTCGA)
library(SummarizedExperiment)
library(dplyr)
```

## Data Catalogue

Use `getCatalog()` to download the full catalogue of available files:

``` r
getCatalog()
```

Filter by pipeline:

``` r
# HoVerNet files (33,177)
getCatalog("hovernet")

# ProvGigaPath files (21,076)
getCatalog("provgigapath")
```

Key columns in the catalogue: `pipeline`, `format`, `filename`,
`tcga_barcode`, `Case.ID`, `Project.ID`, `level`, `lat`, `lon`.

**HoVerNet formats:** `json`, `h5ad`, `thumb` (PNG thumbnail)

**ProvGigaPath formats:** `csv` (slide-level and tile-level,
distinguished by the `level` column: `"slide_level"` or `"tile_level"`)

## Importing HoVerNet Data

HoVerNet segmentation results are imported as `SpatialExperiment` or
`SpatialFeatureExperiment`. Columns represent individual nuclei.

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
#> class: SpatialExperiment
#> dim: 0 67081
#> metadata(1): type_map
#> assays(1): counts
#> colData names(10): cell_id x ... B sample_id
#> spatialCoords names(2): x y
```

Each nucleus has:

- `x`, `y` spatial coordinates (pixels relative to the slide)
- `type` — integer cell type code
- `type_prob` — classification probability
- `label` — cell type name (e.g. `"neopla"`, `"inflam"`, `"connec"`,
  `"necros"`)

``` r
colData(hspe)
```

## Importing ProvGigaPath Embeddings

### Slide-level embeddings

Each row is a slide; the embedding vector (768 dimensions, named
`V1`…`V768`) summarises the entire WSI. `ProvGiga()` for a single file
returns a tibble.

``` r
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
#> # A tibble: 1 × 771
#>   slideName  tumorType fileName     V1     V2     V3  ...
#>   <chr>      <chr>     <chr>     <dbl>  <dbl>  <dbl>
#> 1 TCGA-VG-…  <NA>      TCGA-VG… -0.355 0.584 -0.402 …
#> # ℹ 768 embedding columns (V1–V768) + 3 metadata columns
```

### Tile-level embeddings

Each row is a 256 × 256 px tile; columns include `tile_x`, `tile_y`, and
1,536 embedding dimensions (`0`…`1535`).

``` r
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
#> # A tibble: 211 × 1,543
#>   slide_name  tile_id tile_name tile_x tile_y    `0`    `1`  ...
#>   <chr>         <dbl> <chr>      <dbl>  <dbl>  <dbl>  <dbl>
#> 1 TCGA-VG-…        0 02155x_…   2155  18612  0.177  -1.42  …
#> # ℹ 201 more rows; 1,536 embedding columns (named `0`–`1535`) + tile_x, tile_y
```

## Importing Multiple Files

`ProvGigaList()` accepts multiple files at once. When all files are the
same level, `import()` returns a single `SummarizedExperiment`; for
mixed levels it returns a named list.

``` r
# Three slide-level files for TCGA-GBM → single SummarizedExperiment (768 × 3)
pgl <- getCatalog("provgigapath") |>
    dplyr::filter(level == "slide_level", Project.ID == "TCGA-GBM") |>
    dplyr::slice(1:3) |>
    getFileURLs() |>
    ProvGigaList() |>
    import()

pgl
#> class: SummarizedExperiment
#> dim: 768 3
#> metadata(5): slideName tumorType fileName patientIds sampleIds
#> assays(1): embeddings
#> colnames(3): TCGA-02-0001-01Z-00-DX1 TCGA-02-0001-01Z-00-DX2 TCGA-02-0001-01Z-00-DX3
```

``` r
# Both levels for the same slide → list with $slide_level and $tile_level
pgl_mixed <- getCatalog("provgigapath") |>
    dplyr::filter(
        filename == paste(
            "TCGA-VG-A8LO-01A-01-DX1",
            "B39A4D64-82A1-4A04-8AB6-918F3058B83B",
            "csv", "gz", sep = "."
        ),
        level %in% c("slide_level", "tile_level")
    ) |>
    getFileURLs() |>
    ProvGigaList() |>
    import()

pgl_mixed$slide_level
#> class: SummarizedExperiment
#> dim: 768 1
#> assays(1): embeddings
#> colnames(1): TCGA-VG-A8LO-01A-01-DX1

pgl_mixed$tile_level
#> class: SummarizedExperiment
#> dim: 1538 1
#> assays(1): tiles
#> colnames(1): TCGA-VG-A8LO-01A-01-DX1
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
#>  date     2026-04-13
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
