# `sm-SNIPER`

![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)

<p align="center">
  <img src="./images/sm-SNIPER.png" width=40% height=40% >
</p>

*Construction still in progress*

A Snakemake workflow for **Highly accurate Single Nucleotide polymorphisms calling and Implication of haplotypes in Probe-capture based long-Read Nanopore sequencing** from raw FastQ to phased VCF. It is especially designed to process amplicon sequencing data for non-model organisms such as *Plasmodium* spp.

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#dependencies)
  * [Features](#features)
  * [Usage](#usage)
  * [Installation](#installation)
  * [Configuration](#resultsration)
  * [Execution](#execution)
  * [Examples Dataset](#examples-dataset)
  * [Result](#result)
  * [Benchmarking](#benchmarking)
  * [Haplotype phasing](#Haplotype-phasing)
  * [Tips & FAQs](#tips)
  * [Limitation](#limitation)
  * [Reporting Issue](#reporting-issue)
  * [Link](#link)



## Authors
- [Myo T. Naung](https://github.com/myonaung)
- [Swapnil Tichkule](https://github.com/Swap90)
- [Andrew Guy](https://github.com/andrewguy)
- [Somesh Mehra](https://github.com/smehra34)
- [Somya Mehra](https://github.com/somyamehra)

## Software
This project is written based on the following software

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| BCFtools       | https://doi.org/10.1093/gigascience/giab008       |
| bedtools       | https://doi.org/10.1093/bioinformatics/btq033     |
| Freebayes      | https://doi.org/10.48550/arXiv.1207.3907          |
| GATK-4         | https://doi.org/10.1038/ng.806                    |
| longshot       | https://doi.org/10.1038/s41467-019-12493-y        |
| minimap2       | https://doi:10.1093/bioinformatics/btab705        |
| NanoSim        | https://doi.org/10.1093/gigascience/gix010        |
| PEPPER         | https://doi.org/10.5281/zenodo.5275510            |
| R-tidyverse    | https://doi.org/10.21105/joss.01686               |
| samtools       | https://doi.org/10.1093/bioinformatics/btp352     |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |
| VCFtools       | https://doi.org/10.1093/bioinformatics/btr330     |
| WhatsHap       | https://doi.org/10.1186/s13059-020-02158-1        |


## Dependencies 
The following software are required to install prior to running sm-SNIPER
* [Conda](https://docs.conda.io/en/latest/) >= 4.12.0
* [Singularity](https://github.com/apptainer/singularity) >= 3.8.5 
* [graphviz](https://graphviz.org/) (optional)

## Features
- **minimap2**: mapping to the reference genome
- **Quality control**: Read-depth and coverage calculation with samtools 
- **Variant calling with [longshot](https://github.com/pjedge/longshot)**: Support only SNVs and restricted to high quality variants from high coverage regions (at least 10X)
- **Variant calling with [PEPPER](https://github.com/kishwarshafin/pepper) which is run on singularity container**: Support SNVs and Indels including variants from low quality and low coverage regions
- **Merging SVNs supported by both callers**: highly confidents SNPs from high quality regions
- **Classification of variants (SVNs) by Support Vector Machine based on in-house reference database (it is genereated by sequencing of reference strains using [STAR-seq protocol](https://www.protocols.io/private/ACE2C16BC17D11EC94CE0A58A9FEAC02))**: rescues variants that failed the first filter using reference database
- **Final variant call**: final variant call 
- **Generation of only primary alignment BAM files**: to use for haplotype phasing
- **Haplotype phasing with [WhatsHap](https://github.com/whatshap/whatshap)**: Haplotype phasing to final variant callset  

**Note on in-house reference database**: The data includes true and false variants from the targeted amplicon panel (detected by variant calling methods - longshot or PEPPER) for known reference strains (3D7 and BB12) at 100X coverage. Two clinical samples are also sequenced at 200X coverage and variants detected by each variant method are cross-checked with variants from [Naung et al., 2022](https://github.com/myonaung/Naung-et-al-2021), and concensus variants are assumed as true variants from clinical samples. 

*Disclaimer*: These samples are not used in the following benchmarking to avoid overfitting.

<p align="center">
  <img src="./images/sm-SNIPER.svg" >
</p>

<h4 align="center">Figure 1: Framework for variation discovery with sm-SNIPER.</h4>


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

## Installation
1. Install snakemake, which requires conda & mamba, according to the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. Clone/download this repository (e.g. `git clone https://github.com/myonaung/sm-SNIPER.git`)
## Configuration
### Sample annotation specifications
* `Input (FastQ)` files must be annotated with sample name, and thus `(sample_name).fastq`.
* Based on the analyses, the following parameters in the `workflow/config/config.yaml` file and resource files in `workflow/resources/` are to be adjusted 

  - `reference` - name of the target reference genome along with index .fai file from workflow/resources/ folder (e.g. resources/2. PlasmoDB-46_Pfalciparum3D7_Genome.fasta)
  - `bed` - bed coordinate files for region of interest
  - `data` - file path to folder that contains fastq files (e.g. desktop/fastq)
  - `ont_chemistry` - the chemistry of flowcell used for sequencing (default is R9 flowcell that is `ont_r9_guppy5_sup`, other options include `ont_r10_q20` for R10 chemistry or `hifi` (for Hifi). 
  - `min_coverage`: minimum coverage used for variant calling
  - `max_coverage`: maximum coverage used for variant calling
  - `min_alt_frac`: specification of a potential SNV (or minor clones in the case of malaria multiclonal infection) to have at least this fraction of alternate allele observations
* Based on nature of data to be analysed, it is recommended to change `svm_training_longshot.txt` and `svm_training_pepper.txt` from `workflow/resources/` folder. However, in the absence of reference dataset, the `training dataset` from  `workflow/resources/` should suffice.
* Failed log from downstream `PEPPER` variant calling steps are to be ignored at the moment since they are not relying for the pipeline.


### Execution

#### 1. Install and activate conda environment
It is recommended to execute always from within top level of the pipeline directory (i.e `sm-SNIPER/`). Firstly, conda environment that includes all the core software has to be installed upon the first run of the workflow. It might take several minutes.

```
###envname can be replaced by any name

git clone https://github.com/myonaung/sm-SNIPER.git

conda env create --name envname --file=workflow/envs/default.yml
conda activate envname
```
#### 2. Download singularity image
Singularity image for PEPPER variant calling step has to be downloaded, and placed it under `workflow/envs` folder.

```
cd sm-SNIPER
singularity pull docker://kishwars/pepper_deepvariant:r0.7
mv pepper_deepvariant_r0.7.sif workflow/envs
```

#### 3. Build folder structure (will be replaced with a different feature soon)
Relevent folders to run sm-SNIPER is created using `init.sh`. The path to data (i.e. fastq files) has to be added to the `init.sh` file.

```
sh init.sh
```

#### 4. Execute a dry-run
Checking the pipeline with dry-run options. It is to print a summary of the DAG of jobs
```
cd workflow
snakemake -p -n
```
#### 5. Execute workflow local
Command for execution with two cores
```
cd workflow
snakemake -p -c2 -k
```

#### 6. Execute workflow on a cluster
#### 6a. Slurm system

We can invoke snakemake with the profile to automatically submit to the cluster with `sbatch` as follow:
```
snakemake --profile slurm/
```

Parameters can be customised by editing to `config.yaml` file from `workflow/slurm` folder. Details customisation can be found at this [repository](https://github.com/jdblischak/smk-simple-slurm).

## Examples Dataset
To ensure reproducibility of results and to make the pipeline easy-to-replicate, we provide all required reference data for the analysis on Zendodo: 
- [nanopore simulated BAM](https://zenodo.org/deposit/6571220)
- [nanopore amplicon-seq BAM file from 3D7 reference strain mock infection](https://zenodo.org/record/6622424#.YqA7-RNByu4)
- [nanopore ama1 BAM file from synthetic mocked infection of 3D7 and BB12 mixture](https://zenodo.org/record/6630467#.YqL1ehNBz_Q)

## Result

The final variant output can be found in the `sm-SNIPER/workflow/results/{SAMPLE}/out/final_vcf_SNIPER`.

## Benchmarking
Benchmarking of the reference strain shown in the following example was done based on the amplicons from Naung *et al.,* 2023. Genomic coordinates of these can be found in the [bed file](https://github.com/myonaung/sm-SNIPER/blob/main/workflow/resources/myo_panel_bed.bed).

### 1. Quality control

Information on read-depth and coverage can be found at `workflow/results/{SAMPLE}/out/QC`. 

* _amplicon_depth.txt represents per base-pair sequenced read-depth.
* _mapping_summary.txt represents mapping summary.
* _amplicon_total_depth.txt represents average read-depth for each amplicon specified in bed file.
* _mean_depth.txt represents mean read-depth of the amplicon panel.


### 2. Baseline error or false discovery rate (FDR)
To quantify baseline error rates for sm-SNIPER in the context of [STAR-seq](https://www.protocols.io/private/ACE2C16BC17D11EC94CE0A58A9FEAC02), we mapped raw amplicon sequencing data from probe-capture based 3D7 mocked infection (*with high human dna background*) to the publicly available *P. falciparum 3D7* reference genome (version 3). The sequencing was done an average of 250X coverage. The SNVs identification of sm-SNIPER with  popular variant calling methods are evaluated.

<!--- (FDR = expected (# false predictions/ # total predictions)) ---> 
<!--- Precision = TruePositives / (TruePositives + FalsePositives) ---> 

|Method           | coverage (X)|length (bps) |No. expected SNVs|No. observed SNVs   |
|:---------------:| :----------:| :---------: | :-------------: | :-----------------:|
|Freebayes v1.3.6 |     250     | 140629      |             0   |              75    |  
|BCFtools v1.15.1 |     250     | 140629      |             0   |              16    | 
|PEPPER           |     250     | 140629      |             0   |               8    | 
|Longshot         |     250     | 140629      |             0   |               1    | 
|sm-SNIPER        |     250     | 140629      |             0   |               0    | 
|sm-SNIPER        |     50      | 140629      |             0   |               0    | 
|sm-SNIPER        |     30      | 140629      |             0   |               0    |  
|sm-SNIPER        |     15      | 140629      |             0   |               0    |  
|sm-SNIPER        |      5      | 140629      |             0   |              10    |  

[Freebayes](https://github.com/freebayes/freebayes) was executed as follows:

```
freebayes -f reference.fasta -t bed.bed -p 2 -m 20 -q 7 --limit-coverage 3000 --min-coverage 10 -C 10 --throw-away-indel-obs bam.bam > out.vcf
```

[BCFtools](http://samtools.github.io/bcftools/bcftools.html) was executed as follows:

```
bcftools mpileup -d 3000 -Q 7 -Ou -I -f reference.fasta -R bed.bed bam.bam | bcftools call -mv -Ov -o out_raw.vcf 
bcftools filter  -sLowQual -e'%QUAL<20 & MQ < 10 | AC<=3 & DP<5 | %MAX(AC)/%MAX(DP) < 0.01' out_raw.vcf >  out_filtered.vcf
```
*Execution of longshot and PEPPER was similar to the ones used in sm-SNIPER.*

### 3. Capacity to identify true variants
We used synthetic mocked infection of 3D7 and BB12 mixture to evaluate true variant discovery with *sm-SNIPER*, and minor clone is BB12 strains. We used highly polymorphic and repeat-free *ama1* (PF3D7_1133400) sequenced at 200X coverage with STAR-seq as a benchmarking amplicon. 

|Mocked infection|Minor clone prop  | Gene ID         | bps   |No. expected SNVs| No. detected true SNVs|No. detected false SNVs| FDR| Precision| Recall | F1 score|
|:--------------:| :--------------:| :--------------: | :---: |:---------------:| :------------------:  | :--------------------:|:--:|:--------:|:------:|:-------:|
|mixed-1| 0.33            |   PF3D7_1133400  |   1869|         32     |            31          |           0           | 0  |    1.00  |    0.97    |  0.98   |
|mixed-2| < 0.1           |   PF3D7_1133400  |   1869|         32     |            14          |           0           | 0  |    1.00  |    0.43    |  0.60   |


### 4. Samples with lower coverage 
To understand the extent of baseline error and SNV detection at lower coverage, we subsampled the alignments file (mixed-1) to 30X, 15X, and 5X coverage with  seqtk v1.3. Variant calls were evaluated at each respective coverage.

```
#subset 10000 reads from fastq
seqtk sample sample.fastq 10000 > subsample.fastq
```

|Coverage  | Gene ID         | bps   |No. expected SNVs| No. detected true SNVs|No. detected false SNVs|  Precision| Recall     | F1 score|
| :--------------:| :--------------: | :---: |:---------------:| :------------------:  | :--------------------:|:--------:|:----------:|:-------:|
|  50X            |   PF3D7_1133400  |   1869|         32      |           30          |             0         |    1      |  0.94      |   0.97     |
|  30X            |   PF3D7_1133400  |   1869|         32      |           30          |             0         |    1      |  0.94      |   0.97     |
|  15X            |   PF3D7_1133400  |   1869|         32      |           30          |             0         |    1      |  0.94      |   0.97     |
|  5X             |   PF3D7_1133400  |   1869|         32      |           13          |             1         |    0.93   |  0.41      |   0.57     |



## Haplotype phasing  

Haplotype phasing method used in `sm-SNIPER` was based on `longshot` and [WhatsHap](https://github.com/whatshap/whatshap), and thus haplotype switch error may still present in the cases where there are more than 2 infected strains in a single sample. Therefore, for the best results in haplotype phasing, we recommend to use [SHAPEIT4](https://odelaneau.github.io/shapeit4/#documentation) which uses phase information (*value of phase tag (PS) block*) from WhatsHap to extract phase information along with our filtered bam file. See [here](https://github.com/myonaung/Mapping-within-host-antigenic-escape) for example workflow of *Plasmodium falciparum*.

See phased results on the example vcf files at [ama1_final_SNIPER.vcf](https://github.com/myonaung/sm-SNIPER/blob/main/workflow/results/ama1_final_SNIPER.vcf), where major clone is the reference strain and thus clear seperation of `0|1` is expected for all variants in genotype (GT)field.

To reconstruct the fasta files from a phased VCF, the `bcftools consensus` command can be used. 

```
bgzip phased.vcf
tabix phased.vcf.gz
bcftools consensus -H 1 -f reference.fasta phased.vcf.gz > haplotype1.fasta
bcftools consensus -H 2 -f reference.fasta phased.vcf.gz > haplotype2.fasta
```
## Tips
Here are some tips for troubleshooting & FAQs:
- always first perform a dry-run with option `-n`
- always run the pipeline with `-k` options to complete independent steps if an upstream step fails
- in case the pipeline crashes, manually cancel the pipeline as follow
```
snakemake --unlock 
snakemake -k -c8 --rerun-incomplete
```
- command for generating the directed acyclic graph (DAG) of all jobs with current resultsration (installation of [graphviz](https://graphviz.org/) will be required)
```
snakemake --dag --forceall | dot -Tsvg > workflow/dags/all_DAG.svg
```

## Limitation
Our sm-SNIPER pipeline cannot perform indel discovery and SNV calling is restricted to MOI = 2.


## Reporting Issue
Please create issue [here](https://github.com/myonaung/sm-SNIPER/issues/new) for any problem developed from using `sm-SNIPER` or to request a new features.

## Link
[Github page](https://myonaung.github.io/sm-SNIPER/)



