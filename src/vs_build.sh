#!/bin/bash

# This is a build script for setting up the appropriate directories/files prior to running PEPSI
# To run:
#
#       bash vs_build.sh
#

# Set up working directories
src_dir=$(pwd)
data="../data"
cadd="$data/cadd"

# Download CADD v1.3 annotations
cd "$cadd"
wget "https://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3_inclAnno.tsv.gz"
wget "https://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3_inclAnno.tsv.gz.tbi"

# Download MaxEntScan
cd ".."
wget "http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz"
tar -xvf "fordownload.tar.gz"
rm "fordownload.tar.gz"
mv "fordownload" "MaxEntScan"

cd "$src_dir"
