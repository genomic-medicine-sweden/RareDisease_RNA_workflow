[![CI test](https://github.com/genomic-medicine-sweden/RareDisease_RNA_workflow/actions/workflows/ci_test.yml/badge.svg?branch=main)](https://github.com/genomic-medicine-sweden/RareDisease_RNA_workflow/actions/workflows/ci_test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.3-brightgreen.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# RareDisease_RNA_workflow

nextflow main.nf --help

To test the pipeline run:

```Console
nextflow main.nf -profile test,<singularity/docker>
```

run a single sample:

```Console
nextflow main.nf -profile singularity --r1 read1.fq.gz --r2 read2.fq.gz --sample sampleID --output output_directory
```

run a single sample with multiple fastq files

```Console
nextflow main.nf -profile singularity --r1 "folder/*R1*.fq.gz" --r2 "folder/*R2*.fq.gz" --sample sampleID --output output_directory
```

NOTE: you need to add quotation marks around the search pattern

run all samples in a samplesheet:

```Console
nextflow main.nf -profile singularity --samplesheet sample.csv --output output_directory
```

The samplesheet is a comma-separated file with the following header:

```Console
id,r1,r2
```

## Setup
Modify the config file:

```
reference_dir : specify the folder with all your references

star_index : the star reference index folder

fasta : the reference fasta file

gtf : gene annotations in gtf format

strandedness : library strandedness <forward/reverse>, can be empty
```

The pipeline will automatically download and cache the singularity/docker images.

You can also install all dependencies, as listed in dependencies

## Dependencies
When using singularity/docker:

	nextflow
	singularity/docker

otherwise:

	nextflow
	bcftools
    BootstrapAnn (https://github.com/J35P312/BootstrapAnn)
	fastQC
	gatk
	gffcompare
	multiQC
	picard
	samtools
	STAR
	stringtie
