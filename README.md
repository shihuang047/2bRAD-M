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
 #### Dependencies
 All scripts in 2bRAD-M are written using Perl and recommended to run in a conda environment. In the Unix systems, or Mac OSX, the program should work properly as all required packages can be appropreiately download and installed in the conda environment. 
 #### Disk space
 Construction of a 2bRAD-M standard database (i.e., 2b-Tag-DB) requires approximately 10 GB of disk space. 
 #### Memory usage
 Running the standard pipeline requires < 30Gb of RAM, which is also compatible with multithreading. For example, the BcgI-derived (default) database size is 9.32 GB, and you will need more than that in RAM if you want to build the default database. In a test early on, the peak memory can reach up to 29GB.
 #### Speed 
About 20 minutes are required for loading the 2b-Tag-DB. For a typical gut metagenome, ~40 minutes are required for species profiling.
 
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
   
   Once you have Miniconda installed, create a conda environment with the yml file `tools/2bRAD-M-2020.11.24-conda.yml`.
   
   `conda env create -n 2bRAD-M-20201225 --file tools/2bRAD-M-20201225-conda.yml`
   
 * Activate the 2bRAD-M conda environment by running the following command:
 
   `source activate 2bRAD-M-20201225`

 ### Construct the reference 2B-Tag database (required)
 * Construct the 2bRAD-M species unique marker database (2B-Tag-DB)
   
   (1) Download the prebuilt 2b-Tag-DB from Figshare based on the NCBI Refseq (Oct., 2019).
   
   (2) Download the full genomes from NCBI Refseq for secondary 2b-Tag-DB construction for each sample.
   
   (3) Download the example datasets for pipeline tutorial
    
    `perl tools/Download_2bRADTagDB.pl $your_database_path`
    
    By default, `$your_database_path=./2B-RAD-M-ref_db/`
    
    It usually can take around 30 mins to save all files in the `your_database_path`, but it still depends on your internet connenction speed and stability.
 
## 2bRAD-M pipeline tutorial
* The 2bRAD-M analysis pipeline comprises a combination of 2bRAD-M scripts and optimized parameters for analyzing the 2bRAD or shotgun metagenomics sequencing data, which can output the most comprehensive output on each sample. The pipeline includes:

    (1) **The digital restriction digestion** It is required when input DNA sequences are longer than 31bp or 33bp (e.g., 150bp) or derived from the common shotgun sequencing protocols. If input DNA sequences were produced by the 2bRAD sequencing protocol this step will be skipped. 
       
    (2) **Qualitative analysis** Identify the microbes and estimate their abundances based on the 2bRAD (such as. BcgI derived) species-specific markers of a prebuilt and comprehensive 2b-Tag-DB based on the NCBI Refseq (Oct., 2019).
       
    (3) **Quantitative analysis** Estimate the microbial abundances more precisely based on the 2bRAD species-specific markers in a sample-specific 2b-Tag-DB that contains only condidate genomes that identified in the given biolgical sample. Given the sample-specific 2b-Tag-DB is more compact, it produces more species-specific 2bRAD markers than the original 2b-Tag-DB and results in more accurate modeling of relative abundance of taxa.
    	


```
DESCRIPTION
    We here provided a streamlined 2bRAD pipeline for analyzing microbial compositions from the 2bRAD/shotgun metagenomics data based on the species-specific 2bRAD markers.

PARAMETERS
  -t   <int>    The acceptable types of an input sequencing data file. The file path should be also listed in the sample list file (para -l)
                [1] generic genome data in a fasta format
                [2] shotgun metagenomic data in a fastq Format(either SE or PE platform is accepted)
                [3] 2bRAD data from a SE sequencing platform in a fastq format
                [4] 2bRAD data from a PE sequencing platform in a fastq format
  -l   <file>   The filepath of the sample list. Each line includes an input sample ID and the file path of corresponding DNA sequence data where each field should be separated by <tab>. A line in this file that begins with # will be ignored. Only four formats of a sample list file are accepted and should match with parameter -t: 
                [1] sample<tab>sample.fa(.gz)
                [2] sample<tab>shotgun.1.fq(.gz)(<tab>shotgun.2.fq.gz)
                [3] sample<tab>2bsingle.fq(.gz or 2bsingle.1.fq.gz)
                [4] sample1<tab>sample2<tab>sample3<tab>sample4<tab>sample5<tab>R1.fq(.gz)<tab>R2.fq(.gz)
  -d   <dir>    The working path of 2B-Tag-DB database
  -o   <dir>    The output directory (if not exists, it will be created automatically as 'outdir')
OPTIONS of Qualitative Analysis
  -p   <str>   If qualitative analysis applies or not [default: $qual] (yes or no)
  -s1  <str>   The enzyme site(s) for the qualitative analysis. One or more of site can be specified. (comma separated) [default: $site1]
               It also represents enzymatic digestion or data splitting, and combining qualitative analysis results(for quantitative analysis).
               [1]CspCI  [5]BcgI  [9]BplI     [13]CjePI  [17]AllEnzyme
               [2]AloI   [6]CjeI  [10]FalI    [14]Hin4I
               [3]BsaXI  [7]PpiI  [11]Bsp24I  [15]AlfI
               [4]BaeI   [8]PsrI  [12]HaeIV   [16]BslFI
  -t1  <str>   The taxonomic level for 2bRAD markers in the qualitative database, which should be one of the following: kingdom,phylum,class,order,family,genus,species,strain. [default: $level1]
OPTIONS of Quantitative Analysis
  -q   <str>   If the quantitative analysis applies or not [default: $quan] (yes or no)
  -gsc <int>   G score threshold for identifying the condidate microbes present in a sample in qualitative analysis, which also determines the membership of sample-specific 2B-Tag-DB database in the quantitative analysis step. [default: $g_score_threshold, it means >$g_score_threshold]
  -gcf <int>   The threshold of the 2bRAD tag number for the presence of a microbial genome (i.e., GCF) in the qualitative analysis, which also determines the membership of sample-specific 2B-Tag-DB database in the quantitative analysis step. [default: $GCF_threshold, it means >$GCF_threshold]
  -s2  <str>   The enzyme site for the quantitative analysis. (refer to -s1) [default: $site2, must be included in para -s1]
  -t2  <str>   The taxonomic level for 2bRAD markers in the quantitative database, which should be one of the following: kingdom,phylum,class,order,family,genus,species,strain. [default: $level2]
OPTIONS of CPU
  -c1  <int>   The number of CPUs used in the digital digestion step for multiple samples [default: $cpu1]
  -c2  <int>   The number of CPUs used for calculating the abundance for multiple samples [default: $cpu2] (each CPU needs about 15~65G of memory)
OPTIONS of Quality Control
  -qc  <str>   If quality control apply or not [default: $qc] (yes or no)
  -qcn <float> The maximum ratio of base \"N\" [default: $qc_n]
  -qcq <int>   The minimum quality score to keep [default: $qc_q]
  -qcp <int>   The minimum percentage of bases that must have [-qcq] quality [default: $qc_p]
  -qcb <int>   The quality values of base [default: $qc_b]
OPTIONS of Abundance Stat
  -ms  <str>   Mock Sample Name (separated by commas)
  -ncs <str>   Negative Control Sample Name (separated by commas)
  -h|help   Print this help information.
```
 
   * [Analyze in silico mock community (synthetic 2bRAD sequencing data: `simulate_50.fa.gz`)](docs/analyze_mock.md)
 
  `perl bin/2bRADM_Pipline.pl -t 1 -l $your_database_path/list_simulation -d $your_database_path -o outdir -s1 5,13 -s2 5,13 -gsc 60`

  * [Analyze the 2bRAD sequencing data of a mock microbial community: MSA1002 (`MSA1002_R1.fq.gz`)](docs/analyze_mock.md)
   [MSA1002](https://www.atcc.org/en/Global/Products/MSA-1002.aspx) comprises the genomic material from 20 microbial strains that are evenly mixed. We sequenced the DNA in this sample using our 2bRAD protocol for optimizing and testing the bioinformatic pipeline.
 
  `perl bin/2bRADM_Pipline.pl -t 3 -l $your_database_path/list_mock -d $your_database_path -o outdir`
 
## 2bRAD-M scripts for customized analyses 
 * [Extract 2b tags](docs/extract_2b.md)
 * [Build your own customized 2bTag database](docs/build_db.md)
 * [Species profiling for a single sample](doc/profile_single_sample.md)
 * [Merge species profiles for multiple samples](doc/profile_single_sample.md)
 
## Acknowledgement

   This work was funded by Grant 31800088 from National Natural Science Foundation and 2019M652501 from China Postdoctoral Science Foundation, and Taishan Scholar Fund of Shandong Province of China. 

