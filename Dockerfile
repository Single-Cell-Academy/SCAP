# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:3.6.3

# system libraries of general use
## install debian packages
RUN apt-get update -qq && \
    apt-get upgrade -y && \
    apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libtiff-dev \
    libjpeg-dev \
    libudunits2-dev \
    libgdal-dev \
    gcc \
    libpq-dev -y

## update system libraries
RUN apt-get install -y git && \
    apt-get install python3 -y && \
    apt-get install python3-pip -y && \
    apt-get install python3-venv -y && \
    apt-get install python3-wheel -y && \
    pip3 install -U pip

# clone SCAP repo
RUN git clone https://github.com/Single-Cell-Academy/SCAP.git

WORKDIR "/SCAP"

RUN R -e 'renv::use_python()' && \
    R -e 'renv::restore()'

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('/SCAP/R', host = '0.0.0.0', port = 3838)"]
