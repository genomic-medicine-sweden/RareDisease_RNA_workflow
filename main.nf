#!/usr/bin/dev nextflow

if(params.help){
    println "GMS-RNA workflow"
    println "Usage: nextflow main.nf --r1 read1.fq.gz --r2 --read2.fq.gz --sample sampleID --output output_directory -c config.conf"
    println "or provide a samplesheet"
    println "Usage: nextflow main.nf --samplesheet sample.csv --output output_directory -c config.conf"
    exit 0

}else if(params.r1 && params.r2 && params.sample){

     reads_align=Channel( tuple(params.sample, file(params.r1), file(params.r2)) )

}else if (params.samplesheet){
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
        .into {reads_align}

}else{
    println "Nope."
    exit 0
}

process STAR_Aln{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        set val(sample) , file(r1), file(r2) from reads_align
           
    output:
        file "${sample}.Chimeric.out.junction" into junctions
        file "${sample}.Aligned.sortedByCoord.out.bam" into star_bam
        file "${sample}.Aligned.sortedByCoord.out.bam.bai" into bai
        file "${sample}.ReadsPerGene.out.tab" into geneCounts

    """
    STAR --genomeDir ${params.STAR_ref_dir} \\
         --readFilesIn ${r1} ${r2}  \\
         --twopassMode Basic \\
         --outReadsUnmapped None \\
         --runThreadN ${task.cpus} \\
         --limitBAMsortRAM ${params.STAR_bam_sort_ram} \\
         --outSAMtype BAM SortedByCoordinate \\ 
         --outFileNamePrefix ${sample}. \\
         --quantMode GeneCounts \\
         --outSAMstrandField intronMotif \\
         --readFilesCommand gunzip -c

    samtools index ${sample}.Aligned.sortedByCoord.out.bam
    """

}

