#!/usr/bin/dev nextflow

params.help=false
params.r1=false
params.vcf=false
params.r2=false
params.samplesheet=false
params.platform="ILLUMINA"

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
        .map{ row-> tuple(row.id, row.read1, row.read2 ) }
        .into {reads_align}

}else{
    println "Nope."
    exit 0
}
print(params.r1)
process STAR_Aln{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
	        tuple val(sample),val(r1),val(r2) from reads_align
           
    output:
        tuple val(sample), file("${sample}.Aligned.sortedByCoord.out.bam"), file("${sample}.Aligned.sortedByCoord.out.bam.bai"), file("${sample}.ReadsPerGene.out.tab") into STAR_output

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

process gatk_split{

    input:
        tuple val(sample) , file(bam), file(bai), file(counts) from STAR_output

    output:
        tuple val(sample), file("${sample}.RG.split.Aligned.sortedByCoord.out.bam"), file("${sample}.RG.split.Aligned.sortedByCoord.out.bai") into gatk_split_output

    """
    gatk SplitNCigarReads -R ${params.ref} -I ${bam} -O ${sample}.RG.split.Aligned.sortedByCoord.out.bam --create-output-bam-index
    """
}

process GATK_ASE{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
	        tuple val(sample) , file(bam),file(bai)  from gatk_split_output


    output:
        tuple val(sample), file("${sample}.ase.vcf"), file("${sample}.GATKASE.csv") into gatk_hc

    script:

        """
        gatk HaplotypeCaller -R ${params.ref} -I ${bam} -stand-call-conf 10 -O ${sample}.vcf --minimum-mapping-quality 10
        bgzip ${sample}.vcf
        tabix -p vcf ${sample}.vcf.gz
        gatk ASEReadCounter  -R ${params.ref} -O ${sample}.GATKASE.csv -I ${bam} -V ${sample}.vcf.gz
        python ${params.bootstrapann} --vcf ${sample}.vcf.gz --ase ${sample}.GATKASE.csv > ${sample}.ase.vcf
        """


}
