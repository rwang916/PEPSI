#!/bin/bash

# This is a build script for setting up the appropriate directories/files prior to running PEPSI
# To run:
#
#	bash vs_build.sh
#

# Set up working directories
data="../data"
cadd="$data/cadd"

# Download CADD v1.3 annotations
wget "https://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3_inclAnno.tsv.gz" -P "$cadd"
wget "https://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3_inclAnno.tsv.gz.tbi" -P "$cadd"

# Download MaxEntScan
wget "http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz" -P "$data"
tar -xvf "$data/fordownload.tar.gz"
rm "$data/fordownload.tar.gz"
mv "$data/fordownload" "$data/MaxEntScan"

# Download SVMBP
wget "https://bitbucket.org/regulatorygenomicsupf/svm-bpfinder/get/727e2d8ea4ad.zip" -P "$data"
unzip "$data/727e2d8ea4ad.zip"
rm "$data/727e2d8ea4ad.zip"
mv "$data/regulatorygenomicsupf-svm-bpfinder-727e2d8ea4ad" "$data/SVMBP"
wget "http://regulatorygenomics.upf.edu/Software/SVM_BP/calculate_best_BP_per_intron.pl" -P "$data/SVMBP"
