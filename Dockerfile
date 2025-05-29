FROM debian:bullseye-slim

ARG DEBIAN_FRONTEND="noninteractive"
ARG TZ="UTC"
ENV TZ="${TZ}"

WORKDIR /app

RUN apt-get update && apt-get install -y wget bash

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
RUN bash /tmp/miniconda.sh -b -p /app/miniconda
RUN rm /tmp/miniconda.sh

# Add Conda to PATH
ENV PATH=/app/miniconda/bin:$PATH

# Initialize Conda
RUN /app/miniconda/bin/conda init

# Configure conda channels
RUN /app/miniconda/bin/conda config --add channels bioconda && \
    /app/miniconda/bin/conda config --add channels conda-forge

# Copy and install environment specification
COPY conda_env_explicit.txt .

RUN /app/miniconda/bin/conda install --file conda_env_explicit.txt 

# Activate environment by default
SHELL ["/bin/bash", "-c"]

# Copy script into container
COPY screen_tradb.py /usr/local/bin/screen_tradb.py

# Make it executable and symlink to /usr/local/bin
RUN chmod +x /usr/local/bin/screen_tradb.py && \
    ln -s /usr/local/bin/screen_tradb.py /usr/local/bin/screen_tradb

