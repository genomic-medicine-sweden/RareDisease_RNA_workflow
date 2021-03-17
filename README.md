# RareDisease_RNA_workflow

nextflow main.nf --help

run a single sample:

	nextflow main.nf --r1 read1.fq.gz --r2 --read2.fq.gz --sample sampleID --output output_directory -c config.conf

run all samples in a samplesheet:

	nextflow main.nf --samplesheet sample.csv --output output_directory -c config.conf

# setup
Modify the config file:

	STAR_ref_dir : the path to the star reference index folder

	ref : the path to the reference fasta file (dict and fai file required)

Download the singularity collection:

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

