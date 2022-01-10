#!/usr/bin/dev nextflow

params.help=false
params.r1=false
params.r2=false
params.samplesheet=false

if(params.help){
    println "GMS-RNA workflow"
    println "Usage: nextflow main.nf --r1 read1.fq.gz --r2 --read2.fq.gz --sample sampleID --output output_directory -c config.conf"
    println "or provide a samplesheet"
    println "Usage: nextflow main.nf --samplesheet sample.csv --output output_directory -c config.conf"
    exit 0

}else if(params.r1 && params.r2 && params.sample){

     reads_align=Channel.of( [params.sample, params.r1, params.r2] )

}else if (params.samplesheet){
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
        .set {reads_align}

}else{
    println "Nope."
    exit 0
}

process STAR_Aln{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
	    tuple val(sample),val(r1),val(r2) from reads_align
           
    output:
        tuple val(sample), file("${sample}.Aligned.sortedByCoord.out.bam"), file("${sample}.Aligned.sortedByCoord.out.bam.bai") into STAR_output
        tuple val(sample), file("${sample}.ReadsPerGene.out.tab")
        tuple val(sample), file('*Log.out'), file('*Log.final.out'), file('*Log.progress.out') into star_multiqc

    """
    zcat ${r1} > R1.fastq
    zcat ${r2} > R2.fastq

    STAR --genomeDir ${params.STAR_ref_dir} \\
         --readFilesIn R1.fastq R2.fastq\\
         --twopassMode Basic \\
         --outReadsUnmapped None \\
         --runThreadN ${task.cpus} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMattrRGline ID:${sample} PL:${params.platform} SM:${sample} \\
         --outFileNamePrefix ${sample}. \\
         --quantMode GeneCounts \\
         --outSAMstrandField intronMotif

    samtools index ${sample}.Aligned.sortedByCoord.out.bam 
    """

}

STAR_output.into {gatk_split_input; metrics_input}

process picard_collectrnaseqmetrics{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(bam), path(bai) from metrics_input
    
    output:
        tuple val(sample), path("${sample}_rna_metrics.txt") into metric_multiqc

    """
    picard CollectRnaSeqMetrics \\
        --STRAND_SPECIFICITY ${params.strandedness} \\
        --REF_FLAT ${params.annotation_refflat} \\
        --INPUT ${bam} \\
        --OUTPUT ${sample}_rna_metrics.txt \\
    """

}

process gatk_split{

    input:
        tuple val(sample) , file(bam), file(bai) from gatk_split_input

    output:
        tuple val(sample), file("${sample}.RG.split.Aligned.sortedByCoord.out.bam"), file("${sample}.RG.split.Aligned.sortedByCoord.out.bai") into gatk_split_output

    """
    gatk SplitNCigarReads -R ${params.ref} -I ${bam} -O ${sample}.RG.split.Aligned.sortedByCoord.out.bam --create-output-bam-index
    """

}

process GATK_ASE{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
	        tuple val(sample) , file(bam),file(bai) from gatk_split_output

    output:
        tuple val(sample), file("${sample}.ase.vcf"), file("${sample}.GATKASE.csv") into gatk_hc

//  TODO: Add bcftools to container
    script:

        """
        gatk HaplotypeCaller -R ${params.ref} -I ${bam} -stand-call-conf 10 -O ${sample}.vcf --minimum-mapping-quality 10
        bcftools view --genotype het --max-alleles --min-alleles 2 --types snps -O z -o {sample}.vcf.gz ${sample}.vcf
        tabix -p vcf ${sample}.vcf.gz
        gatk ASEReadCounter -R ${params.ref} -O ${sample}.GATKASE.csv -I ${bam} -V ${sample}.vcf.gz
        python ${params.bootstrapann} --vcf ${sample}.vcf.gz --ase ${sample}.GATKASE.csv > ${sample}.ase.vcf
        """

}

// Combine metric output files to one channel
multiqc_input = star_multiqc.join(metric_multiqc).collect{it[1..-1]}

process multiqc{
    publishDir "${params.output}", mode: 'copy', overwrite: true
    
    input: 
        //tuple val(sample), file(picard_metrics) from metric_multiqc
        path(qc_files) from multiqc_input

    output:
        path "*multiqc_report.html"
        path "*_data"
    """
    multiqc .
    """
}

process utilities_DROP{
    // this process will generate 1) gene counts from STAR output 2) sample_annotations.tsv for DROP 3) config.yml for DROP

    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(tab) from star_gene_counts
    
    output:
        file("external_geneCounts.tsv") into external_geneCounts

    """
        generate_gene_counts.py --sample $sample --star $star_gene_counts --output external_geneCounts.tsv
        generate_drop_sample_annot.py --samples --counts --output
        generate_drop_config.py
    """

}
// instructions 1: process *ReadsPerGene.out.tab from STAR and produce a gene_counts.tsv.gz with the following headers:
// geneID, sample1, sample2, ..

// instructions 2: construct sample_annotation.tsv with the following headers:
// if gene counts: col_1=rna_id, ..., col_5=drop_group, col_11=absolute path, col_12=gene_annotation
// if for AS+MAE: RNA_ID,RNA_BAM_FILE,DNA_VCF_FILE,DNA_ID,DROP_GROUP,PAIRED_END,COUNT_MODE,COUNT_OVERLAPS,STRAND,HPO_TERMS,GENE_COUNTS_FILE,GENE_ANNOTATION,GENOME

// instructions 3: construct config.yml

// instructions 4: use 2 + 3 files as input for running DROP
// snakemake will need to be in a docker that inherits
