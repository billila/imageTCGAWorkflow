# imageTCGAWorkflow

An end-to-end Bioconductor workflow for histopathology image analysis
using TCGA diagnostic whole-slide images.

## Ecosystem packages

[TABLE]

Click here to explore the shiny app:
[imageTCGA](https://shiny.sph.cuny.edu/app/imageTCGA/)

## Overview

The TCGA image database contains ~11,765 diagnostic whole-slide images
(WSI) from ~9,640 patients. This workflow covers:

1.  **Data discovery** — browse and select images with the `imageTCGA`
    Shiny app
2.  **Data access** — download WSIs from GDC using `GenomicDataCommons`
3.  **Feature retrieval** — import pre-computed HoVerNet and
    Prov-GigaPath features via `imageFeatureTCGA`
4.  **Spatial analysis** — PCA, Moran’s I, LISA with `imageTCGAutils`
5.  **Visualization** — overlay cell segmentation on tissue thumbnails
    with `HistoImagePlot`
6.  **Downstream analyses** — multi-omics integration (MOFA+), point
    pattern analysis, survival

## Workshop Docker

docker run -e PASSWORD=bioc -p 8787:8787
ghcr.io/billila/imagetcgaworkflow:latest

This will start an RStudio Server instance accessible at
<http://localhost:8787> (username: rstudio, password: bioc)

\`\`\`

## Website

Full workflow documentation:
<https://billila.github.io/imageTCGAWorkflow>
