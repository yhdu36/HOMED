
FROM bioconductor/bioconductor_docker:RELEASE_3_22

WORKDIR /app

# Install HOMED dependencies from Bioconductor and CRAN
RUN Rscript -e "BiocManager::install( \
    c('limma', 'lumi', 'nlme', 'genefilter', \
      'matrixStats', 'EpiDISH', 'doParallel', \
      'foreach', 'ggplot2', 'dplyr', \
      'RColorBrewer'), \
    ask=FALSE, update=FALSE)"

# Copy HOMED R package source
COPY . /app/HOMED

# Install HOMED itself
RUN R CMD INSTALL /app/HOMED

# Verify the installation
RUN Rscript -e "library(HOMED); \
    cat('HOMED installed successfully\n')"

# Default to an interactive R session
CMD ["R"]

