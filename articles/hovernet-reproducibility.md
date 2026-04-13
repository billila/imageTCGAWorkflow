# Reproducing HoVerNet Nuclei Segmentation on TCGA Images

## Overview

This document describes how to reproduce the HoVerNet nuclei
segmentation results that underlie the pre-computed features available
in [imageFeatureTCGA](https://github.com/waldronlab/imageFeatureTCGA).

The complete setup instructions, scripts, and environment files are
maintained in the
**[hovernethelp](https://github.com/hpages/hovernethelp)** repository by
Hervé Pagès.

HoVer-Net ([Graham et al.,
2019](https://doi.org/10.1016/j.media.2019.101563)) performs
simultaneous nuclei **detection**, **segmentation**, and
**classification** on histopathology whole-slide images (WSI). The
pipeline was run on all TCGA diagnostic images and the outputs (JSON,
GeoJSON, H5AD, thumbnails) are accessible via `imageFeatureTCGA`.

## Hardware Requirements

Running HoVerNet from scratch requires GPU resources:

| Resource | Minimum                       |
|----------|-------------------------------|
| GPU      | NVIDIA A100 or L40S           |
| CPUs     | 16                            |
| RAM      | 60 GB                         |
| OS       | Ubuntu 24.04                  |
| Storage  | ~780 GB for full TCGA outputs |

Processing time is approximately **76–81 minutes per image** on a single
GPU node. A 7-node mini-cluster (described below) achieves the full TCGA
dataset in 15–25 hours.

> If you only need the pre-computed results, skip to the
> [imageFeatureTCGA workflow
> vignette](https://billila.github.io/imageTCGAWorkflow/articles/imageTCGAWorkflow.md)
> — no GPU is required to access the stored features.

## Environment Setup

All instructions below follow the [hovernethelp setup
guide](https://github.com/hpages/hovernethelp/blob/main/setup_hovernet_Ubuntu2404.txt).

### GPU drivers

For **NVIDIA A100**:

``` bash
sudo apt-get install nvidia-linux-grid-535
```

For **NVIDIA L40S**:

``` bash
sudo apt-get install nvidia-driver-550
```

Verify the installation:

``` bash
nvidia-smi
```

### Miniconda and HoVer-Net

``` bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda init bash && source ~/.bashrc

# Clone HoVer-Net
git clone https://github.com/vqdang/hover_net
cd hover_net

# Create conda environment from hovernethelp
wget https://raw.githubusercontent.com/hpages/hovernethelp/main/environment.yml
conda env create -f environment.yml
conda activate hovernet
```

The `environment.yml` pins key packages:

- Python 3.12.3
- PyTorch 2.5.1
- `opencv-python`, `scikit-image`, `scikit-learn`, `scipy`, `pandas`
- `openslide` for whole-slide image reading

### Pre-trained model

The pannuke_type model (TF2→PyTorch conversion) is downloaded with
`gdown`:

``` bash
pip install gdown
gdown <model_gdrive_id> -O pannuke_type_tf2pytorch.tar
tar -xf pannuke_type_tf2pytorch.tar
```

See the [hovernethelp setup
file](https://github.com/hpages/hovernethelp/blob/main/setup_hovernet_Ubuntu2404.txt)
for the exact model ID and path configuration.

## Downloading TCGA Images

Images are downloaded using the `GenomicDataCommons` Bioconductor
package. The `imageTCGA` Shiny app generates the query code
interactively (see the [main workflow
vignette](https://billila.github.io/imageTCGAWorkflow/articles/imageTCGAWorkflow.html#step-1-data-discovery-with-imagetcga)).

``` r
# Example: download a single OV diagnostic image
library(GenomicDataCommons)

q <- files() |>
    filter(~ cases.project.project_id == "TCGA-OV" &
               data_type == "Slide Image" &
               experimental_strategy == "Diagnostic Slide") |>
    results(size = 1)

gdcdata(q$file_id, destination_dir = "tcga_images/")
```

## Running HoVerNet Inference

### Single image

``` bash
python run_infer.py \
  --gpu='0' \
  --nr_types=6 \
  --type_info_path=type_info.json \
  --batch_size=12 \
  --model_mode=original \
  --model_path=pannuke_type_tf2pytorch/hovernet_fast_pannuke_type_tf2pytorch.tar \
  --nr_inference_workers=1 \
  --nr_post_proc_workers=3 \
  wsi \
  --input_dir=tcga_images/ \
  --output_dir=hovernet_output/ \
  --save_thumb \
  --save_mask
```

Key parameters:

| Parameter              | Value | Note                                                                                |
|------------------------|-------|-------------------------------------------------------------------------------------|
| `nr_types`             | 6     | Nuclei classes (neoplastic, inflammatory, connective, dead, epithelial, background) |
| `batch_size`           | 12    | Adjust to GPU memory                                                                |
| `nr_inference_workers` | 1     | GPU workers                                                                         |
| `nr_post_proc_workers` | 3     | CPU post-processing workers                                                         |

### Batch processing (multiple images)

The
[`infer_batch2.sh`](https://github.com/hpages/hovernethelp/blob/main/infer_batch2.sh)
script automates the full pipeline:

1.  Reads a manifest of TCGA image IDs to process
2.  Downloads each SVS file via `GenomicDataCommons` (with rsync retry,
    up to 10 attempts with 5-minute intervals)
3.  Runs HoVerNet inference and validates JSON output
4.  Pushes results to a central storage node (`hoverboss`)
5.  Updates the manifest with success/failure status and timestamps

``` bash
# Usage
bash infer_batch2.sh manifest.tsv
```

## Output Files

For each input SVS image, HoVerNet produces:

| File           | Format  | Content                                            |
|----------------|---------|----------------------------------------------------|
| `<id>.json`    | JSON    | Per-nucleus centroid, contour, type, probabilities |
| `<id>.geojson` | GeoJSON | Same data in GIS-compatible format                 |
| `<id>.h5ad`    | H5AD    | Feature matrix (intensity, morphology statistics)  |
| `<id>.png`     | PNG     | Tissue thumbnail with segmentation overlay         |

These files are exactly what `imageFeatureTCGA` imports via
`HoverNet() |> import()` — see the [imageFeatureTCGA package
page](https://billila.github.io/imageTCGAWorkflow/articles/pkg-imageFeatureTCGA.md)
for how to load them into a `SpatialExperiment`.

## Scaling with a Mini-Cluster

The [hovernethelp minicluster
setup](https://github.com/hpages/hovernethelp/tree/main/minicluster)
describes a 7-node architecture:

- **Worker nodes** (`hovernet1–4`, `hovernet6–7`, `kakapo1`): GPU nodes
  running inference
- **Central node** (`hoverboss`): 4 CPUs, 779 GB storage, collects
  results

Workers sync completed JSON/H5AD files to `hoverboss` via cron every 15
minutes. Each worker processes a non-overlapping range of the TCGA image
list. Total throughput: **100+ images/day**, completing the full TCGA
dataset in 15–25 hours.

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
