# RareDisease_RNA_workflow

nextflow main.nf --help

run a single sample:

	nextflow main.nf -profile singularity --r1 read1.fq.gz --r2 --read2.fq.gz --sample sampleID --output output_directory -c config.conf

	optionally, a vcf file may be provided:
	
	nextflow main.nf -profile singularity --samplesheet sample.csv --output output_directory --vcf input.vcf -c config.conf


run all samples in a samplesheet:

	nextflow main.nf -profile singularity --samplesheet sample.csv --output output_directory -c config.conf

the samplesheet is a comma-separated file with the following header:

	sample,r1,r2,vcf

The sample, r1 and r2 are mandatory, the vcf column may be left empty	

# setup
Modify the config file:

    reference_dir : specify the folder with all your references 

	STAR_ref_dir : the star reference index folder

	ref :the reference fasta file (dict and fai file required)

The pipeline will automatically download and cache the latest singularity image. 

Alternatively you can download the singularity collection:

	singularity pull shub://J35P312/RareDisease_RNA_workflow

Or install all dependencies, as listed in dependencies

# dependencies
When using singularity:

	nextflow
	singularity

otherwise:

	nextflow
	samtools
	STAR
	gatk
	stringtie
	picard
	star-fusion
	fusioncatcher
	Arriba	
	multiQC
	fastQC
    BootstrapAnn (https://github.com/J35P312/BootstrapAnn)

