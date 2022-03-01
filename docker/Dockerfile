FROM centos:7
MAINTAINER Ben Webb <ben@salilab.org>

# Install ATSAS and its dependencies
COPY ATSAS-3.0.5-2.el7.x86_64.rpm /
RUN yum -y localinstall /ATSAS-3.0.5-2.el7.x86_64.rpm && yum -y install libgfortran5 git && rm -f /ATSAS*.rpm

# Set up basic conda environment using conda-forge packages
RUN curl -L https://github.com/conda-forge/miniforge/releases/download/4.11.0-0/Miniforge3-4.11.0-0-Linux-x86_64.sh -o miniforge.sh && echo "4c24b38969ac413efa3a984290296f72578340d06004f2c7ba5efcbf828ec86b  miniforge.sh" sha256sum -c --status && bash ./miniforge.sh -b -p /root/miniforge && rm -f miniforge.sh && /root/miniforge/bin/conda init -q bash

# Install Python packages needed by IHMValidation and install extra RPMs
# needed to run firefox
RUN source /root/.bashrc && eval "$(conda shell.bash hook)" && conda install -y pandas flask ihm scikit-learn selenium bokeh requests firefox geckodriver && pip install python-decouple pdfkit && yum -y install gtk3 alsa-lib

# Get wkhtmltopdf directly from web site (there is a conda package, but it
# does not handle https: links)
RUN curl -LO https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6-1/wkhtmltox-0.12.6-1.centos7.x86_64.rpm && yum -y localinstall wkhtmltox-0.12.6-1.centos7.x86_64.rpm && rm -f wkhtmltox-0.12.6-1.centos7.x86_64.rpm

# Build molprobity from source code and put the binaries we need in /usr/bin/
RUN source /root/.bashrc && eval "$(conda shell.bash hook)" && yum -y install mesa-libGL-devel subversion ca-certificates git gcc-c++ && mkdir molprobity && cd molprobity && curl -O https://raw.githubusercontent.com/rlabduke/MolProbity/master/install_via_bootstrap.sh && bash install_via_bootstrap.sh 12 && ln -sf /molprobity/build/bin/molprobity.clashscore /usr/bin/ && ln -sf /molprobity/build/bin/molprobity.ramalyze /usr/bin/ && ln -sf /molprobity/build/bin/molprobity.rotalyze /usr/bin/ && ln -sf /molprobity/build/bin/molprobity.molprobity /usr/bin/

# Get IHMValidation and configure
RUN git clone https://github.com/salilab/IHMValidation.git && source /root/.bashrc && eval "$(conda shell.bash hook)" && yum -y install which && echo "ATSAS=\"`which datcmp`\"" > IHMValidation/.env && echo "Molprobity_ramalyze=\"`which molprobity.ramalyze`\"" >> IHMValidation/.env && echo "Molprobity_molprobity=\"`which molprobity.molprobity`\"" >> IHMValidation/.env && echo "Molprobity_clashscore=\"`which molprobity.clashscore`\"" >> IHMValidation/.env && echo "Molprobity_rotalyze=\"`which molprobity.rotalyze`\"" >> IHMValidation/.env && echo "wkhtmltopdf=\"`which wkhtmltopdf`\"" >> IHMValidation/.env
