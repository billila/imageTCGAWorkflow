FROM bioconductor/bioconductor_docker:RELEASE_3_23
WORKDIR /home/rstudio

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-setuptools \
        python3-dev \
        python3-pip \
        libicu-dev \
        libhdf5-dev \
        libgdal-dev \
        libudunits2-dev \
        libproj-dev \
        libgeos-dev \
        libcairo2-dev libfreetype6-dev libpng-dev \
        libtiff5-dev libjpeg-dev libxt-dev libharfbuzz-dev libfribidi-dev \
    && pip3 install --break-system-packages --ignore-installed mofapy2==0.7.1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY --chown=rstudio:rstudio . /home/rstudio/imageTCGAWorkflow

RUN Rscript -e "reticulate::use_python('/usr/bin/python3', force=TRUE)"

RUN Rscript -e "devtools::install('imageTCGAWorkflow', dependencies = TRUE, \
    build_vignettes = TRUE)"
