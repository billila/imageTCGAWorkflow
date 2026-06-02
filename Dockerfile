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
    build_vignettes = TRUE)"
