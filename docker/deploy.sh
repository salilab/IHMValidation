#!/bin/bash

apt update -y
apt install -y wget git

# Download prebuilt singularity package for Ubuntu 20.04
wget https://salilab.org/~arthur/ihmv/packages/singularity_3.8.4-2_amd64.deb
# backup link
# wget https://vsb.fbb.msu.ru/share/aozalevsky/pdb_dev/ihmv/packages/singularity_3.8.4-2_amd64.deb

# Install singularity
apt install -yf ./singularity_3.8.4-2_amd64.deb

# Remove deb
rm ./singularity_3.8.4-2_amd64.deb

# Download precompiled singularity image
wget https://salilab.org/~arthur/ihmv/prebuilt_containers/ihmv_20250205.sif

# Clone code
git clone --branch dev_2.0 https://github.com/salilab/IHMValidation.git

# Create dirs
mkdir -p input output cache

# Test installation
singularity exec --pid --bind IHMValidation/:/opt/IHMValidation,input:/ihmv/input,output:/ihmv/output,cache:/ihmv/cache ihmv_20250205.sif /opt/IHMValidation/ihm_validation/ihm_validator.py --output-root /ihmv/output --cache-root /ihmv/cache --force -h
