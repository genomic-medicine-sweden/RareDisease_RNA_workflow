[![CI test](https://github.com/genomic-medicine-sweden/RareDisease_RNA_workflow/actions/workflows/ci_test.yml/badge.svg?branch=main)](https://github.com/genomic-medicine-sweden/RareDisease_RNA_workflow/actions/workflows/ci_test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.3-brightgreen.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
# RareDisease_RNA_workflow

nextflow main.nf --help

run a single sample:

	nextflow main.nf -profile singularity --r1 read1.fq.gz --r2 read2.fq.gz --sample sampleID --output output_directory -c config.conf

run a single sample with multiple fastq files

	nextflow main.nf -profile singularity --r1 "folder/*R1*.fq.gz" --r2 "folder/*R2*.fq.gz" --sample sampleID --output output_directory -c config.conf

NOTE: you need to add quotation marks around the search pattern

run all samples in a samplesheet:

	nextflow main.nf -profile singularity --samplesheet sample.csv --output output_directory -c config.conf

the samplesheet is a comma-separated file with the following header:

	id,r1,r2

The id, r1 and r2 are mandatory, the vcf column may be left empty

# setup
Modify the config file:

    reference_dir : specify the folder with all your references

	star_index : the star reference index folder

	fasta : the reference fasta file

	gtf : gene annotations in gtf format

The pipeline will automatically download and cache the latest singularity image.

Alternatively you can download the singularity collection:

	singularity pull shub://J35P312/RareDisease_RNA_workflow

Or install all dependencies, as listed in dependencies

# dependencies
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
