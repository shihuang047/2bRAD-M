# 2b-RAD-M
The computational pipeline for microbiome analysis on 2b-RAD data

## How it works
 
## Installation
 ### Requirements
 * Operating systems: Unix, OSX
 * Conda/Miniconda >= 3
 ### Speed and memory usage
 Most steps of the program are quite fast, require < 2Gb of RAM, and are compatible with multithreading. About 20 minutes are required for loading the 2bTag     database. For a typical gut metagenome, ~1-5 minutes are required for species profiling.
 ### Download the pipeline
 Clone the latest version from GitHub (recommended):  
`git clone https://github.com/shihuang047/2bRAD-M/`  
 This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added.
 ### Install 2bRAD-M pipeline in a conda environment 
 * Conda installation
   [Miniconda] (https://docs.conda.io/en/latest/miniconda.html) provides the conda environment and package manager, and is the recommended way to install 2bRAD-M. 
 * Create a conda environment for 2bRAD-M pipeline:
   After installing Miniconda and opening a new terminal, make sure you’re running the latest version of conda:
   
   `conda update conda`
   
   Once you have Miniconda installed, create a conda environment with the yml file at the directory: `tools/2bRAD-M-2020.11.24-conda.yml`.
   
   `conda env create -n 2bRAD-M-2020.11.24 --file tools/2bRAD-M-2020.11.24-conda.yml`
   
 * Activate the 2bRAD-M conda environment by running the following command:
 
   `conda activate 2bRAD-M-2020.11.24`
   
 * Perl dependency installation:
 
   `cpan Parallel::ForkManager`
   
 ### Fetch the reference database (required) and download the example data
 `gunzip tools/BuildstructureAndDownload.mk.gz`
 `make -f tools/BuildstructureAndDownload.mk`
 The path of 2bRAD-M pipeline can be assign by "Database_path=", for instance:
 `make -f tools/BuildstructureAndDownload.mk Database_path=./2B-RAD-M-ref_db/`
 
## 2bRAD-M pipeline tutorial
 * [Analyze the MOCK-MSA1002 community (sequenceing data)](docs/analyze_mock.md)
 perl bin/2bRADM_Pipline.pl -t 3 -l list -d 2B-RAD-M-ref_db -o output
 * [Analyze in silico mock community (synthetic shogun data)](docs/snp_diversity.md)
 perl bin/2bRADM_Pipline.pl -t 1 -l list -d 2B-RAD-M-ref_db -o output -s1 17 -s2 17 -gsc 60
 
## 2bRAD-M scripts for customized analyses 
 * [Extract 2b tags](docs/extract_2b.md)
 * [Build your own customized 2bTag database](docs/build_db.md)
 * [Species profiling for a single sample](doc/profile_single_sample.md)
 * [Merge species profiles for multiple samples](doc/profile_single_sample.md)
 
## Citation

