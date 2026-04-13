# Prov-GigaPath Feature Extraction

## Overview

**Prov-GigaPath** ([Xu et al.,
2024](https://doi.org/10.1038/s41586-024-07441-w)) is a whole-slide
pathology foundation model trained on over 1.3 billion pathology image
tiles from Providence Health & Services. It produces high-quality,
general-purpose embeddings for computational pathology.

We ran Prov-GigaPath on all available TCGA diagnostic WSIs and extracted
embeddings at two granularities:

| Level           | Unit                    | Embedding dim | Files in catalogue |
|-----------------|-------------------------|---------------|--------------------|
| **Tile level**  | 256 × 256 px tile       | 1,536         | ~14,000 CSV files  |
| **Slide level** | Entire WSI (aggregated) | 768           | ~7,000 CSV files   |

The resulting features are accessible without any GPU via
`imageFeatureTCGA` — see the [main workflow
vignette](https://billila.github.io/imageTCGAWorkflow/articles/imageTCGAWorkflow.md)
for import examples.

This document describes how the extraction was performed, for users who
wish to re-run or extend the pipeline.

## Model

Prov-GigaPath uses a two-stage architecture:

1.  **Tile encoder** — a DINOv2 ViT-g/14 model fine-tuned on 1.3 B
    pathology tiles. Each 256 × 256 px tile at 0.5 µm/px is encoded into
    a 1,536-dimensional vector.
2.  **Slide encoder** — a LongNet transformer that aggregates tile-level
    tokens into a single slide-level representation, preserving spatial
    context via dilated attention.

The model weights are distributed via Hugging Face
(`prov-gigapath/prov-gigapath`) under the [CC-BY-NC-4.0
license](https://creativecommons.org/licenses/by-nc/4.0/). Access
requires a Hugging Face account and accepting the model terms.

## Container

Inference was run inside an Apptainer/Singularity container to ensure
reproducibility across different compute environments.

### Definition file

``` singularity
Bootstrap: docker
From: mambaorg/micromamba:1.5.8

%files
    environment.yaml /environment.yaml
    . /app

%post
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

    apt-get update && apt-get install -y git && apt-get clean

    # Create the gigapath conda environment
    micromamba create -y -n gigapath -f /environment.yaml

    # Install the package and dependencies
    cd /app
    micromamba run -n gigapath pip install -e .
    micromamba run -n gigapath pip install "timm>=1.0.3"
    micromamba run -n gigapath pip install "openslide-bin"

    micromamba clean --all --yes

%environment
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH=$MAMBA_ROOT_PREFIX/envs/gigapath/bin:$PATH
```

> **Note:** The Hugging Face token (`HF_TOKEN`) is required to download
> the model weights and must be set as an environment variable at
> runtime — do **not** embed it in the container definition or commit it
> to version control. Set it at runtime:
> `apptainer run --env HF_TOKEN=<your_token> gigapath.sif`

### Build the container

``` bash
apptainer build gigapath.sif gigapath.def
```

## Running Tile-Level Extraction

Tiles are extracted at **0.5 µm/px** (20× equivalent magnification),
tiled into non-overlapping 256 × 256 px patches, and encoded by the
DINOv2 tile encoder.

``` python
import gigapath
import torch
from PIL import Image
import pandas as pd
import glob, os

# Load the tile encoder
tile_encoder = gigapath.load_tile_encoder("hf_hub:prov-gigapath/prov-gigapath")
tile_encoder.eval()

transform = gigapath.tile_transform(img_size=256)

def encode_tiles(svs_path, output_csv):
    tiles = gigapath.tile_slide(svs_path, mpp=0.5, tile_size=256)
    embeddings = []
    for tile_img, coords in tiles:
        t = transform(tile_img).unsqueeze(0)
        with torch.no_grad():
            emb = tile_encoder(t).squeeze().numpy()
        embeddings.append({"x": coords[0], "y": coords[1], **{f"d{i}": v for i, v in enumerate(emb)}})
    pd.DataFrame(embeddings).to_csv(output_csv, index=False)
```

Output: one CSV per slide, rows = tiles, columns = `x`, `y`,
`d0`…`d1535`.

## Running Slide-Level Extraction

The slide encoder aggregates all tile embeddings using LongNet
attention.

``` python
slide_encoder = gigapath.load_slide_encoder("hf_hub:prov-gigapath/prov-gigapath")
slide_encoder.eval()

def encode_slide(tile_csv, output_csv):
    df = pd.read_csv(tile_csv)
    coords = torch.tensor(df[["x", "y"]].values, dtype=torch.float32)
    tile_embs = torch.tensor(df[[f"d{i}" for i in range(1536)]].values,
                             dtype=torch.float32).unsqueeze(0)
    with torch.no_grad():
        slide_emb = slide_encoder(tile_embs, coords.unsqueeze(0)).squeeze().numpy()
    pd.DataFrame([slide_emb]).to_csv(output_csv, index=False)
```

Output: one CSV per slide, a single row of 1,536 values.

## Batch Processing on a Compute Cluster

For large-scale processing across all TCGA slides:

``` bash
#!/bin/bash
# Run on a SLURM cluster, one job per slide

apptainer run \
    --nv \
    --env HF_TOKEN="${HF_TOKEN}" \
    --bind /data/tcga:/data/tcga \
    --bind /data/output:/data/output \
    gigapath.sif \
    python extract_features.py \
        --input  /data/tcga/${FILE_ID}.svs \
        --output /data/output/${FILE_ID}_tiles.csv \
        --level  tile
```

## Accessing Pre-Computed Embeddings

If you do not have GPU resources, the pre-computed embeddings for all
TCGA slides are available via `imageFeatureTCGA`:

``` r
library(imageFeatureTCGA)
library(dplyr)

# Tile-level embeddings
getCatalog("provgigapath") |>
    dplyr::filter(level == "tile_level", Project.ID == "TCGA-OV") |>
    head(1) |>
    getFileURLs() |>
    ProvGiga() |>
    import()

# Slide-level embeddings
getCatalog("provgigapath") |>
    dplyr::filter(level == "slide_level", Project.ID == "TCGA-OV") |>
    head(1) |>
    getFileURLs() |>
    ProvGiga() |>
    import()
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
