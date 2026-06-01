FROM bioconductor/bioconductor_docker:RELEASE_3_23

WORKDIR /home/rstudio

# System libraries for spatial packages (sf, spdep, sfdep, spatstat)
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3-pip \
        libicu-dev \
        libhdf5-dev \
        libgdal-dev \
        libudunits2-dev \
        libproj-dev \
        libgeos-dev \
    && pip3 install --break-system-packages tensorflow \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy workshop package source
COPY --chown=rstudio:rstudio . /home/rstudio/imageTCGAWorkflow

# Install the package + all Depends/Imports/Suggests from DESCRIPTION
RUN Rscript -e "devtools::install('imageTCGAWorkflow', dependencies = TRUE, \
    build_vignettes = FALSE)"

# Pre-download workshop data as rstudio user
USER rstudio

RUN Rscript -e ' \
    library(imageFeatureTCGA); library(dplyr); \
    sid <- "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D"; \
    getCatalog("hovernet") |> \
        filter(filename == paste0(sid, ".json.gz")) |> \
        getFileURLs() |> HoverNet(outClass = "SpatialExperiment") |> import()'

RUN Rscript -e ' \
    library(imageFeatureTCGA); library(dplyr); \
    sid <- "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D"; \
    getCatalog("hovernet") |> \
        filter(filename == paste0(sid, ".h5ad.gz")) |> \
        getFileURLs() |> HoverNet(outClass = "SpatialExperiment") |> import()'

RUN Rscript -e ' \
    library(imageFeatureTCGA); library(dplyr); \
    sid <- "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D"; \
    getCatalog("hovernet") |> \
        filter(filename == paste0(sid, ".h5ad.gz")) |> \
        getFileURLs() |> HoverNet(outClass = "SpatialFeatureExperiment") |> import()'

RUN Rscript -e ' \
    library(imageFeatureTCGA); library(dplyr); \
    sid <- "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D"; \
    getCatalog("provgigapath") |> \
        filter(filename == paste0(sid, ".csv.gz"), level == "slide_level") |> \
        getFileURLs() |> ProvGiga() |> import()'

RUN Rscript -e ' \
    library(imageFeatureTCGA); library(dplyr); \
    sid <- "TCGA-23-1021-01Z-00-DX1.F07C221B-D401-47A5-9519-10DE59CA1E9D"; \
    getCatalog("provgigapath") |> \
        filter(filename == paste0(sid, ".csv.gz"), level == "tile_level") |> \
        getFileURLs() |> ProvGiga() |> import()'

USER root

EXPOSE 8787
CMD ["/init"]
