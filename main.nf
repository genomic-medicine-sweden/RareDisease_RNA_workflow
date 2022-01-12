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

     reads_align=Channel.of( [params.sample, file(params.r1), file(params.r2)] )

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

process stringtie{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample) , file(bam),file(bai)  from stringtie_input

    output:
        tuple val(sample), file("${sample}.stringtie.gtf") into stringtie_output

    script:

    """
    stringtie ${bam} -p ${task.cpus} ${params.stranded} -G ${params.gtf} > ${sample}.stringtie.gtf
    """

}

process gffcompare{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(stringtie_gtf) from stringtie_output

    output:
        tuple val(sample), file("${sample}.stringtie.stats"), file("${sample}.stringtie.annotated.gtf") into gffcompare_output

    script:
    """
    gffcompare -r ${params.gtf} -o ${sample}.stringtie ${stringtie_gtf}
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
