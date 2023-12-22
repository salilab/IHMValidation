FROM ubuntu:22.04
MAINTAINER Arthur Zalevsky <aozalevsky@gmail.com>

ARG ATSAS
ARG TZ=America/Los_Angeles
ENV TZ=${TZ}
ENV DEBIAN_FRONTEND=noninteractive

# Install common dependencies
RUN apt-get update -y && apt-get install -y curl git subversion build-essential

# Set up basic conda environment using conda-forge packages
RUN curl -L https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Mambaforge-23.3.1-1-Linux-x86_64.sh -o miniforge.sh && echo "dbe92c011a1315b9626e2f93a165892f4b89177145bc350b3859a483a3642a24  miniforge.sh" sha256sum -c --status && bash ./miniforge.sh -b -p /opt/conda && rm -f miniforge.sh && /opt/conda/bin/conda init -q bash

# Ensure that we always use a shell with conda
ENV PATH /opt/conda/bin:${PATH}
RUN echo "conda activate" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"] 

# Build molprobity from source code and put the binaries we need in /usr/bin/
ADD "https://api.github.com/repos/rlabduke/MolProbity/commits?per_page=1" latest_commit
RUN mkdir /opt/molprobity && cd /opt/molprobity && curl -O https://raw.githubusercontent.com/rlabduke/MolProbity/master/install_via_bootstrap.sh && bash ./install_via_bootstrap.sh 4
ENV PATH ${PATH}:/opt/molprobity/build/bin

# Install SAS tools: ATSAS and SASCIFTOOLS
COPY ${ATSAS} /.
RUN apt-get update -y && apt-get install -y /${ATSAS} && rm -f /${ATSAS}
RUN pip install git+https://git.embl.de/grp-svergun/sasciftools3.git@master#egg=sasciftools
# Install additional certificate for SASBDB
RUN curl https://cacerts.digicert.com/ThawteEVRSACAG2.crt.pem -o ThawteEVRSACAG2.crt --output-dir /usr/local/share/ca-certificates/ && update-ca-certificates
ENV REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt

# Install Python packages and other dependencies needed by IHMValidation
RUN apt-get update -y && apt-get install -y xvfb libasound2 libdbus-glib-1-2 libgtk-3-0 && mamba install -y pip requests numpy pandas scikit-learn matplotlib "bokeh<3.0" selenium pyvirtualdisplay firefox geckodriver jinja2 && pip install -U iqplot pdfkit ihm && curl -L -O https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6.1-2/wkhtmltox_0.12.6.1-2.jammy_amd64.deb && apt install -y ./wkhtmltox_0.12.6.1-2.jammy_amd64.deb && rm -rf wkhtmltox_0.12.6.1-2.jammy_amd64.deb
# Hardcode selenium path
ENV SE_MANAGER_PATH=/opt/conda/bin/selenium-manager

# Get IHMValidation
ADD "https://api.github.com/repos/salilab/IHMValidation/commits?per_page=1" latest_commit
RUN git clone --depth 1 https://github.com/salilab/IHMValidation.git /opt/IHMValidation

# Create input and output dirs
RUN mkdir -p /ihmv/input /ihmv/output /ihmv/cache /ihmv/databases
