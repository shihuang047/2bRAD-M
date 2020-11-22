# 2b-RAD-M
The computational pipeline for microbiome analysis on 2b-RAD data

## How it works
 
## Installation
 ### Requirements
 * Operating systems: Unix, OSX
 * Perl=5.26.2
 ### Speed and memory usage
 Most steps of the program are quite fast, require < 2Gb of RAM, and are compatible with multithreading. About 20 minutes are required for loading the 2bTag     database. For a typical gut metagenome, ~1-5 minutes are required for species profiling.
 ### Reference database
 * [Download reference full genome database from NCBI](docs/ref_db.md) 
 * [Download the prebuilt 2bTag database](docs/2bTag_db.md)
 ### Download the software
 Clone the latest version from GitHub (recommended):  
`git clone https://github.com/shihuang047/2bRAD-M/`  
 This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added.

 ### Dependencies and update your environmental variables
 * Install the latest version of PERL (5.26.2) via https://www.perl.org/get.html
 * Install PEAR via https://anaconda.org/bioconda/pear/0.9.6/download/osx-64/pear-0.9.6-h977ceac_6.tar.bz2
 
## 2bRAD-M pipeline tutorial
 * [Analyze the mock community (MOCK-MSA1002)](docs/analyze_mock.md)
 * [Analyze in silico mock community](docs/snp_diversity.md)
 
## 2bRAD-M scripts for customized analyses 
 * [Extract 2b tags](docs/extract_2b.md)
 * [Build your own custom 2bTag database](docs/build_db.md)
 * [Species profiling for a single sample](doc/profile_single_sample.md)
 * [Merge species profiles for multiple samples](doc/profile_single_sample.md)
 
## Citing

