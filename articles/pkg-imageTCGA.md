# Package: imageTCGA

## Overview

**GitHub:** <https://github.com/billila/imageTCGA>

`imageTCGA` is a Shiny application for interactive exploration of the
TCGA diagnostic image database. It provides clinical and geographical
filtering of ~11,765 images from ~9,640 patients, and generates
ready-to-run `GenomicDataCommons` R code for downloading selected
images.

## Launch the app

``` r
library(imageTCGA)
imageTCGA::imageTCGA()
```

## Features

### Clinical filtering

The left panel allows filtering by:

- **Project** (cancer type, e.g. TCGA-OV, TCGA-BRCA)
- **Histological type** (morphology code)
- **Tumour grade** and **stage**

### Dot plots

For gynaecological tumours (BRCA, OV, UCS, UCEC), dot plots show the
distribution of images across clinical strata.

### Interactive maps

Images are annotated with the contributing institution. The map tab
displays image origins at three levels:

- Institution (centre name)
- US state
- Country

### Code generation

After applying filters, the **Generate code** panel produces a
`GenomicDataCommons` query that reproduces the current selection:

``` r
# Example of generated code:
library(GenomicDataCommons)
files() |>
    filter(
        ~ cases.project.project_id == "TCGA-OV" &
          data_type == "Slide Image" &
          experimental_strategy == "Diagnostic Slide"
    ) |>
    results(size = 100)
```

Copy this code into your R session, then pass the `file_id` column to
`gdcdata()` to download the SVS files.

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
