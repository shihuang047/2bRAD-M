# 2bRAD-M
----------------------------
This repository provides the computational pipeline for microbiome analysis on 2b-RAD data presented in the paper below:

[Species-resolved sequencing of low-biomass microbiomes by 2bRAD-M](https://www.biorxiv.org/content/10.1101/2020.12.01.405647v1) by Zheng Sun, Shi Huang, Pengfei Zhu, Lam Tzehau, Helen Zhao, Jia Lv, Rongchao Zhang, Lisha Zhou, Qianya Niu, Xiuping Wang, Meng Zhang, Gongchao Jing, Zhenmin Bao, Jiquan Liu, Shi Wang, Jian Xu. BioRxiv, doi: https://doi.org/10.1101/2020.12.01.405647

## How it works
 The principle of 2bRAD-M on microbiome analyses on low-biomass samples: 
 
 (1) reliable enzyme-digested sequence tags can be derived that are specific to high-resolution taxa (e.g., species or strain) yet universally applicable for a broad range or all of bacterial, archaeal and fungal genomes; 
 
 (2) these taxa-specific, iso-length sequence tags can be evenly amplified and sequenced;
 
 (3) the tag sequences can be mapped to reference genomes to reconstruct faithfully the taxonomic composition.
 
 You can also find more details for the 2bRAD-M workflow below. 
 
 ![workflow](2bRAD-M_workflow.png)
 
 * The experimental workflow has two steps: 
 
   (1) BcgI (a commercially available Type IIB restriction enzymes) is used, as an example, to digest total genomic DNA extracted from microbiome samples. BcgI recognizes the sequence of CGA-N6-TGC in the genomic DNA and cleaves on both upstream (12-10 bp) and downstream (10-12 bp) of this signature, producing short and iso-length DNA (32bp without sticky ends) across all loci. 
   
   (2) These so-called “2bRAD fragments” are ligated to adaptors, amplified and then sequenced. 
   
 * The computational workflow. The foundation here is a unique 2bRAD tag database (“2b-Tag-DB”), which contains taxa-specific 2bRAD tags identified from all the sequenced bacteria, fungi and archaea genomes. Mapping the 2bRAD reads against 2b-Tag-DB thus identifies the presence of species in a sample. Subsequently, to estimate relative abundance of the identified taxa, the mean read coverage of all 2bRAD tags specific to each taxon is derived. To improve utilization rate of reads and classification accuracy, a secondary, sample-specific 2b-Tag-DB was dynamically derived from only those candidate taxa identified in a particular sample, which produces more species-specific 2bRAD tags than the original 2b-Tag-DB and results in more accurate modeling of relative abundance of taxa.

## Installation
 ### System requirements
 * Operating systems: Unix, OSX
 * Conda/Miniconda >= 3
 ### Speed and memory usage
 Most steps of the program are quite fast, require < 2Gb of RAM, and are compatible with multithreading. About 20 minutes are required for loading the 2bTag     database. For a typical gut metagenome, ~1-5 minutes are required for species profiling.
 ### Download the pipeline
 * Clone the latest version from GitHub (recommended):  
 
   `git clone https://github.com/shihuang047/2bRAD-M/`
   
   `cd 2bRAD-M`
   
    This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added.
 * Directly download without installing GitHub:
 
   `wget https://github.com/shihuang047/2bRAD-M/archive/master.zip`
   
   `unzip master.zip`
   
   `cd 2bRAD-M-master`
   
 ### Install 2bRAD-M pipeline in a conda environment 
 * Conda installation
   [Miniconda](https://docs.conda.io/en/latest/miniconda.html) provides the conda environment and package manager, and is the recommended way to install 2bRAD-M. 
 * Create a conda environment for 2bRAD-M pipeline:
   After installing Miniconda and opening a new terminal, make sure you’re running the latest version of conda:
   
   `conda update conda`
   
   Once you have Miniconda installed, create a conda environment with the yml file ``tools/2bRAD-M-2020.11.24-conda.yml``.
   
   `conda env create -n 2bRAD-M-2020.11.24 --file tools/2bRAD-M-2020.11.24-conda.yml`
   
 * Activate the 2bRAD-M conda environment by running the following command:
 
   `source activate 2bRAD-M-2020.11.24`

 ### Construct the reference 2B-Tag database (required)
 * Construct the 2bRAD-M species unique marker database (2B-Tag-DB)
   
   (1) Download the prebuilt 2b-Tag-DB from Figshare based on the NCBI Refseq (Oct., 2019).
   
   (2) Download the full genomes from NCBI Refseq for secondary 2b-Tag-DB construction for each sample.
   
   (3) Download the example datasets for pipeline tutorial
 
    `perl tools/Download_2bRADTagDB.pl your_database_path(default:./2B-RAD-M-ref_db/)`
    
    It usually can take around 30 mins to save all files in the `your_database_path`, but it still depends on your internet connenction speed and stability.
 
## 2bRAD-M pipeline tutorial
 * [Analyze the MOCK-MSA1002 community (sequencing data (SE 150bp): `MSA1002_R1.fq.gz`)](docs/analyze_mock.md)
 
  `perl bin/2bRADM_Pipline.pl -t 1 -l your_database_path/list_simulation -d your_database_path -o outdir -s1 5,13 -s2 5,13 -gsc 60`

 * [Analyze in silico mock community (synthetic sequencing data (SE 150bp): `simulate_50.fa.gz`)](docs/analyze_mock.md)
 
  `perl bin/2bRADM_Pipline.pl -t 3 -l your_database_path/list_mock -d your_database_path -o outdir`
 
## 2bRAD-M scripts for customized analyses 
 * [Extract 2b tags](docs/extract_2b.md)
 * [Build your own customized 2bTag database](docs/build_db.md)
 * [Species profiling for a single sample](doc/profile_single_sample.md)
 * [Merge species profiles for multiple samples](doc/profile_single_sample.md)
 
## Acknowledgement

   This work was funded by Grant 31800088 from National Natural Science Foundation and 2019M652501 from China Postdoctoral Science Foundation, and Taishan Scholar Fund of Shandong Province of China. 

