# `sm-SNIPER`

![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)


A Snakemake workflow for **Highly accurate Single Nucleotide polymorphisms calling and Implication of haplotypes in Probe-capture based long-Read Nanopore sequencing**

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#dependencies)
  * [Features](#features)
  * [Usage](#usage)
  * [Installation](#installation)
  * [Configuration](#configuration)
  * [Execution](#execution)
  * [Examples](#examples)
  * [Results](#results)
  * [Tips & FAQs](#tips)



## Authors
- [Myo T. Naung](https://github.com/myonaung)

## Software
This project is written based on the following software

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| BCFtools       | https://doi.org/10.1093/gigascience/giab008       |
| bedtools       | https://doi.org/10.1093/bioinformatics/btq033     |
| GATK-4         | https://doi.org/10.1038/ng.806                    |
| longshot       | https://doi.org/10.1038/s41467-019-12493-y        |
| minimap2       | https://doi:10.1093/bioinformatics/btab705        |
| pandas         | https://doi.org/10.5281/zenodo.3509134            |
| PEPPER         | https://doi.org/10.5281/zenodo.5275510            |
| R-tidyverse    | https://doi.org/10.21105/joss.01686               |
| samtools       | https://doi.org/10.1093/bioinformatics/btp352     |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |
| VCFtools       | https://doi.org/10.1093/bioinformatics/btr330     |


## Dependencies 
The following software are required to install prior to running sm-SNIPER
* Conda
* Singularity >= 3.8.5

## Features
- Mapping to the reference genome
- Quality control: Read-depth and coverage calculation
- Variant calling via [longshot](https://github.com/pjedge/longshot)
- Variant calling via [PEPPER](https://github.com/kishwarshafin/pepper)
- Merging SVNs supported by both callers
- Support Vector Machine
- Final variant call
- Generation of only primary alignment BAM files

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

## Installation
1. Install snakemake, which requires conda & mamba, according to the [documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. Clone/download this repository (e.g. git clone https://github.com/myonaung/sm-SNIPER.git)
## Configuration

### Execution

#### 1. Install and activate conda environment
It is recommended to execute always from within top level of the pipeline directory (i.e sm-SNIPER/)

#### 2. Execute a dry-run
Checking the pipeline with dry-run options. It is to print print a summary of the DAG of jobs
```
snakemake -p -n
```
#### 3. Execute workflow local
#### 4. Execute workflow on a cluster

## Examples
To ensure reproducibility of results and to make the pipeline easy-to-replicate, we provide all required reference data for the analysis on Zendodo: 

## Results
## Tips
Here are some tips for troubleshooting & FAQs:
- always first perform a dry-run with option -n
- always run the pipeline with -k options to complete independent steps if an upstream step fails