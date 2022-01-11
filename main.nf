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

     ch_reads = Channel.of( [params.sample, params.r1, params.r2] )

}else if (params.samplesheet){
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
        .set { ch_reads }

}else{
    println "Nope."
    exit 0
}

ch_reads.into { ch_reads_align; ch_reads_qc; ch_sampel_id }
ch_multiqc_input = ch_sampel_id.map{ it.first() }

process fastqc{

    input:
	    tuple val(sample), val(r1), val(r2) from ch_reads_qc
    
    output:
        tuple val(sample), file('*1_fastqc.zip'), file('*2_fastqc.zip') into fastqc_multiqc

    script:

    def R1 = r1.getName()
    def R2 = r2.getName()

    """
    ln -s ${r1} ${R1} 
    ln -s ${r2} ${R2} 
    fastqc --threads ${task.cpus} ${R1} ${R2}
    """

}
ch_multiqc_input = ch_multiqc_input.join(fastqc_multiqc)

process STAR_Aln{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
	    tuple val(sample), val(r1), val(r2) from ch_reads_align
           
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
ch_multiqc_input = ch_multiqc_input.join(star_multiqc)

STAR_output.into {gatk_split_input; metrics_input; stringtie_input}

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
ch_multiqc_input = ch_multiqc_input.join(metric_multiqc)

process stringtie{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(bam),file(bai) from stringtie_input

    output:
        tuple val(sample), file("${sample}_stringtie.gtf") into stringtie_output

    script:

    """
    stringtie ${bam} -p ${task.cpus} ${params.stranded} -G ${params.gtf} > ${sample}_stringtie.gtf
    """

}

process gffcompare{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(stringtie_gtf) from stringtie_output

    output:
        tuple val(sample), file("${sample}_stringtie.annotated.gtf") into gffcompare_output
        tuple val(sample), file("${sample}_stringtie.stats") into gffcompare_multiqc

    script:
    """
    gffcompare -r ${params.gtf} -o ${sample}_stringtie ${stringtie_gtf}
    """
}
ch_multiqc_input = ch_multiqc_input.join(gffcompare_multiqc)

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
ch_multiqc_input = ch_multiqc_input.collect{it[1..-1]}
process multiqc{
    publishDir "${params.output}", mode: 'copy', overwrite: true
    
    input: 
        //tuple val(sample), file(picard_metrics) from metric_multiqc
        path(qc_files) from ch_multiqc_input

    output:
        path "*multiqc_report.html"
        path "*_data"
    """
    multiqc .
    """
}
