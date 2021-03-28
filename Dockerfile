# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny-verse:3.6.3

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
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
    llvm \
    git \
    python-dev \
    python-pip \
    python3 \
    python3-dev\
    python3-dbg \
    python3-pip \
    python3-venv \
    python3-wheel \
    python3-setuptools \
    python3-llvmlite
    
# clone SCAP repo
COPY . ./SCAP

WORKDIR "/SCAP"

RUN python3 -m venv ./renv/python/virtualenvs/renv-python-3.7.3 && \
    ./renv/python/virtualenvs/renv-python-3.7.3/bin/pip3 install --upgrade pip && \
    ./renv/python/virtualenvs/renv-python-3.7.3/bin/pip3 install wheel setuptools &&\
    ./renv/python/virtualenvs/renv-python-3.7.3/bin/pip3 install -r requirements.txt &&\
    R -e 'renv::use_python(python = "./renv/python/virtualenvs/renv-python-3.7.3/bin/python3")' && \
    R -e 'renv::restore()'

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('/SCAP/R', host = '0.0.0.0', port = 3838)"]
