# Downloading TCGA Diagnostic Images

## Introduction

The TCGA diagnostic image database contains approximately **11,765
whole-slide images (WSI)** in SVS format from ~9,640 patients across 32
cancer types. Images are freely available from the [NIH Genomic Data
Commons (GDC)](https://portal.gdc.cancer.gov/).

Two complementary ways to find and download images are covered here:

1.  **`imageTCGA` Shiny app** — interactive point-and-click interface
    with clinical filters and code generation
2.  **`GenomicDataCommons` R package** — programmatic access for
    scripted or large-scale downloads

## Option 1: Shiny App (imageTCGA)

The `imageTCGA` app provides an interactive interface to the full TCGA
image catalogue. It is the recommended starting point for exploratory
data selection.

``` r

imageTCGA::imageTCGA()
```

### App features

**Filter panel** — narrow down images by:

- Cancer type (project ID, e.g. TCGA-OV, TCGA-BRCA)
- Histological type
- Tumour grade and stage
- Patient sex and age range

**Dot plots** — available for gynaecological tumours (BRCA, OV, UCS,
UCEC), showing image counts stratified by clinical variables.

**Interactive maps** — display the geographic origin of each image at
three levels: institution (centre), US state, and country.

**Code generation panel** — once you are satisfied with the filter
selection, click *Generate code*. The app outputs a ready-to-run
`GenomicDataCommons` query (see Option 2 below) that you can copy
directly into your R session.

## Option 2: GenomicDataCommons (programmatic)

The `GenomicDataCommons` Bioconductor package provides direct
programmatic access to the GDC REST API.

### Installation

``` r
BiocManager::install("GenomicDataCommons")
library(GenomicDataCommons)
```

### Query the image catalogue

``` r
# Authenticate (only needed for controlled-access data; images are open-access)
# gdc_token() # not required for diagnostic slides

# Build a query for diagnostic slides from ovarian cancer
q <- files() |>
    GenomicDataCommons::filter(
        ~ cases.project.project_id == "TCGA-OV" &
          data_type             == "Slide Image" &
          experimental_strategy == "Diagnostic Slide"
    )

# Preview results
res <- results(q, size = 10)
res[, c("file_id", "file_name", "file_size")]
```

### Filter by clinical metadata

``` r
# Add clinical filters: stage IV, high-grade serous
q_filtered <- files() |>
    GenomicDataCommons::filter(
        ~ cases.project.project_id               == "TCGA-OV" &
          data_type                              == "Slide Image" &
          experimental_strategy                  == "Diagnostic Slide" &
          cases.diagnoses.tumor_stage            == "stage iv" &
          cases.diagnoses.morphology             == "8441/3"   # serous cystadenocarcinoma
    )

res_filtered <- results(q_filtered, size = 100)
nrow(res_filtered)
```

### Check available metadata fields

``` r
# See all filterable fields for files
available_fields("files") |> grep("cases", x = _, value = TRUE) |> head(20)
```

### Download images

``` r
# Download a single image
file_ids <- results(q, size = 1)$file_id

gdcdata(
    uuids           = file_ids,
    destination_dir = "tcga_images/",
    progress        = TRUE
)
```

The downloaded SVS files can be very large (500 MB – 3 GB each). For
large-scale downloads, consider running on an HPC system and using the
[GDC Data Transfer
Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
(command-line client).

### GDC Data Transfer Tool (batch download)

For hundreds of images, create a manifest file from the GDC portal or
via R and use the CLI tool:

``` r
# Write a manifest
manifest(q_filtered) |> write.table("gdc_manifest.txt",
    sep = "\t", quote = FALSE, row.names = FALSE)
```

``` bash
# Then download with the GDC client
gdc-client download -m gdc_manifest.txt -d tcga_images/
```

### Metadata retrieval

Clinical and demographic metadata can be retrieved alongside file
metadata:

``` r
# Get associated case metadata
case_ids <- res$cases |> sapply(\(x) x[[1]])

cases() |>
    GenomicDataCommons::filter(~ case_id %in% case_ids) |>
    select(c("case_id",
             "diagnoses.tumor_stage",
             "diagnoses.morphology",
             "demographic.gender",
             "demographic.age_at_index")) |>
    results(size = length(case_ids))
```

## Image Format

TCGA diagnostic images are stored as **SVS** (Aperio ScanScope) files, a
multi-resolution pyramidal TIFF format. In R, raw SVS files are
typically not opened directly — instead, pre-computed features (HoVerNet
nuclei segmentation, Prov-GigaPath embeddings) are accessed via
`imageFeatureTCGA` without requiring the original images.

If you need to inspect SVS files directly, use the
[openslide-python](https://openslide.org/api/python/) library:

``` bash
# Install the system library and Python bindings
pip install openslide-bin openslide-python

# Inspect a slide in Python
python3 -c "
import openslide
slide = openslide.OpenSlide('tcga_images/TCGA-XX-XXXX.svs')
print('Levels:', slide.level_count)
print('Dimensions:', slide.dimensions)
print('MPP:', slide.properties.get('openslide.mpp-x'))
"
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
