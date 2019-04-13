# PEPSI
### Prediction of variant Effect on Percent Spliced In

https://github.com/rwang916/PEPSI

Authors: Robert Wang, Yaqiong Wang, Zhiqiang Hu

PEPSI is a tool created for predicting the impact of coding and noncoding variation on 
pre-mRNA splicing based on sequence conservation, RNA secondary structure, and regulatory 
sequence elements.

Please contact [Robert Wang](mailto:rwang916@berkeley.edu) for feedback, questions, or bugs.

## Table of Contents
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)

## Requirements

#### python
PEPSI requires python2 (>= 2.7.9)

#### R
PEPSI requires R (>= 3.2.2) and it is available on the [R website](https://www.r-project.org/).

Additional R packages needed are:
* withr
* ggplot2
* randomForest
* dplyr

#### parallel
GNU parallel is available on their [website](https://www.gnu.org/software/parallel/).
The latest version of GNU parallel should be downloaded and installed. Make sure that
the path to the 'parallel' executable is included in your 'PATH' environment variable.

#### tabix
tabix can be downloaded from their [website](https://github.com/samtools/tabix).
Follow the installation guidelines on their page. The latest version of tabix should
be downloaded and installed. Have the path to the 'tabix' executable in your 'PATH' 
environment variable.

#### RNAplfold
RNAplfold is part of the ViennaRNA suite, which can be downloaded from their 
[website](https://www.tbi.univie.ac.at/RNA/). Download the latest release of the 
ViennaRNA suite and follow the instructions for installation. Include path to RNAplfold 
executable in your 'PATH' environment variable. 

## Installation

PEPSI can be downloaded (cloned) using the git command.
```bash
git clone https://github.com/rwang916/PEPSI.git
```

Once downloaded, run the following commands to download necessary files for running PEPSI
```bash
cd PEPSI/src
bash vs_build.sh
```

## Usage

To run PEPSI, run:
```bash
bash vs_main.sh <path_to_training_set> <path_to_test_set> <number_of_threads> <path_to_training_set_sequences> <path_to_test_set_sequences>
```

Example usage:
```bash
bash vs_main.sh ../data/vexseq/HepG2_delta_PSI_CAGI_training.tsv ../data/vexseq/HepG2_delta_PSI_CAGI_test.tsv 20 
	../data/vexseq/vs.sequences.tsv ../data/vexseq/vs.sequences.tsv
```

