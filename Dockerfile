FROM bioconductor/bioconductor_docker:RELEASE_3_17

ARG VERSION

LABEL org.opencontainers.image.source=https://github.com/js2264/HiCExperiment
LABEL org.opencontainers.image.documentation=https://js2264.github.io/HiCExperiment
LABEL org.opencontainers.image.authors="HiCExperiment authors"
LABEL org.opencontainers.image.description="Import Hi-C matrix file formats in R and performs common array operations on them"
LABEL org.opencontainers.image.licenses=MIT
LABEL org.opencontainers.image.version ${VERSION}

COPY . /app/
WORKDIR /app

RUN make install