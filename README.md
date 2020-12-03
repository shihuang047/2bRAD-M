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
 ### Conda activation and 2bRAD-M pipeline installation
 * Installation of Conda: XXX (link)
 * Environment creation:
 `conda env create -n 2bRAD-M-2020.11.24 --file tools/2bRAD-M-2020.11.24-conda.yml`
 * Environment activation:
 `conda activate 2bRAD-M-2020.11.24`
 * Perl dependency installation:
 `cpan Parallel::ForkManager`
 2bRAD-M environment need to be activated every time you use it by "conda activate 2bRAD-M-2020.11.24"
 ### Reference database and example data downloading
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
 
## Citing

