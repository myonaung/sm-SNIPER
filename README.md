# Snakemake workflow: `sm-SNIPER`

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

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

## Installation

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
## Results
## Tips