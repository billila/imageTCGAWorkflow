# Downstream Analyses

## Introduction

Once histopathology features have been extracted and loaded into
Bioconductor objects, they can be integrated with molecular omics data
(genomics, transcriptomics, proteomics) for downstream multi-modal
analyses.

This vignette illustrates analytical patterns used in the imageTCGA
manuscript. Full reproducible code is available in the manuscript
repository: <https://github.com/billila/manuscript_imageTCGA/>

The image-derived features used are:

- **Slide-level Prov-GigaPath embeddings** (768-dim vector per patient)
- **Tile-level Prov-GigaPath embeddings** (768-dim vector per tile)
- **HoVerNet nuclei classification** (cell type proportions per
  slide/tile)
- **Spatial autocorrelation statistics** from `imageTCGAutils`

## Data preparation

``` r
library(imageFeatureTCGA)
library(dplyr)
library(SummarizedExperiment)

# Import slide-level embeddings for a cohort of OV slides
pgl <- getCatalog("provgigapath") |>
  dplyr::filter(level == "slide_level", Project.ID == "TCGA-OV") |>
  dplyr::slice(1:50) |>
  getFileURLs() |>
  ProvGigaList() |>
  import()

pgl
# SummarizedExperiment: 768 features × 50 slides
# assay "embeddings" contains the Prov-GigaPath vectors
```

## Dimensionality reduction

### PCA on slide embeddings

``` r
emb_mat <- t(assay(pgl, "embeddings"))  # slides × 768

pca_res <- prcomp(emb_mat, scale. = TRUE, center = TRUE)
summary(pca_res)$importance[, 1:10]

library(ggplot2)
pca_df <- data.frame(
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  sample = rownames(pca_res$x)
)

ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(title = "PCA — slide-level GigaPath embeddings (TCGA-OV)")
```

### UMAP

``` r
library(uwot)

umap_res <- umap(
  emb_mat,
  n_neighbors = 15,
  min_dist    = 0.1,
  metric      = "cosine"
)

umap_df <- data.frame(
  UMAP1  = umap_res[, 1],
  UMAP2  = umap_res[, 2],
  sample = rownames(emb_mat)
)

ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(title = "UMAP — slide-level GigaPath embeddings (TCGA-OV)")
```

## Clustering

``` r
km <- kmeans(pca_res$x[, 1:20], centers = 4, nstart = 25)
pca_df$cluster <- factor(km$cluster)

ggplot(pca_df, aes(PC1, PC2, colour = cluster)) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw()
```

## Multi-omics integration with MOFA+

**MOFA+** (Multi-Omics Factor Analysis, [Argelaguet et
al. 2020](https://doi.org/10.1186/s13059-020-02015-1)) identifies latent
factors that explain variance across multiple data modalities
simultaneously. Here image embeddings are used as one modality alongside
molecular omics data.

``` r
library(MOFA2)

# Prepare data as a named list of matrices (features × samples)
# Each matrix must have the same column names (sample IDs)
data_list <- list(
  image = assay(pgl, "embeddings"),  # 768 × n_slides
  rna   = rna_mat                    # genes × n_slides (your RNA matrix)
)

mofa <- create_mofa(data_list)

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15

train_opts <- get_default_training_options(mofa)
train_opts$seed <- 42

mofa <- prepare_mofa(
  mofa,
  model_options   = model_opts,
  training_options = train_opts
)

mofa <- run_mofa(mofa)

# Variance explained per factor per modality
plot_variance_explained(mofa, max_r2 = 15)
```

See the manuscript repository for the full MOFA+ analysis applied to
TCGA-OV: <https://github.com/billila/manuscript_imageTCGA/>

## Point Pattern Analysis (PPA)

Point pattern analysis treats nucleus locations as spatial point
processes and tests for clustering, regularity, or complete spatial
randomness. Applied to HoVerNet nuclei coordinates, it characterises the
spatial organisation of different cell types across the tissue section.

Key metrics used in the manuscript:

- **Ripley’s K / L function** — tests for clustering vs. regularity at
  multiple spatial scales
- **Cross-type K function** — spatial interaction between pairs of cell
  types (e.g. neoplastic and inflammatory)
- **Nearest-neighbour distance distribution** — distribution of
  distances to the nearest nucleus of the same or a different type

These analyses are performed with the `spatstat` R package on the
`(x, y)` coordinates from HoVerNet `SpatialExperiment` objects.

``` r
library(spatstat)
library(imageFeatureTCGA)

# Import HoVerNet data for one slide
hspe <- getCatalog("hovernet") |>
  dplyr::filter(format == "json") |>
  head(1) |>
  getFileURLs() |>
  HoverNet(outClass = "SpatialExperiment") |>
  import()

# Extract coordinates and cell types
coords <- as.data.frame(spatialCoords(hspe))
types  <- colData(hspe)$label

# Create a marked point pattern
win <- owin(
  xrange = range(coords$x),
  yrange = range(coords$y)
)
pp <- ppp(x = coords$x, y = coords$y, marks = factor(types), window = win)

# Ripley's L for neoplastic cells
L_neo <- Lest(subset(pp, marks == "neopla"), correction = "isotropic")
plot(L_neo, main = "Ripley's L — neoplastic nuclei")

# Cross-type K: neoplastic vs. inflammatory
K_cross <- Kcross(pp, i = "neopla", j = "inflam", correction = "isotropic")
plot(K_cross, main = "Cross-type K: neoplastic vs. inflammatory")
```

## Association with Clinical Outcomes

``` r
library(survival)
library(survminer)

# Merge cluster assignment with clinical data
surv_df <- data.frame(
  sample  = names(km$cluster),
  cluster = factor(km$cluster)
) |>
  dplyr::left_join(clinical_data, by = "sample")  # your clinical data frame

fit <- survfit(Surv(os_days, os_status) ~ cluster, data = surv_df)

ggsurvplot(
  fit, data = surv_df,
  pval       = TRUE,
  risk.table = TRUE,
  palette    = "Set1",
  title      = "Overall survival by image-derived cluster (TCGA-OV)"
)
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
