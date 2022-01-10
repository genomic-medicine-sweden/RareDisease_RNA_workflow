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

     reads_align=Channel.of( [params.sample, file(params.r1), file(params.r2), file(params.vcf)] )

}else if (params.samplesheet){
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ row-> tuple(row.id, file(row.read1), file(row.read2),file(row.vcf) ) }
        .into {reads_align}

}else{
    println "Nope."
    exit 0
}


process STAR_Aln{
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "*.{bam,bai,ReadsPerGene.out.tab}" 

    input:
        tuple val(sample), file(r1), file(r2), file(vcf) from reads_align

           
    output:
        tuple val(sample), file("${sample}.Aligned.sortedByCoord.out.bam"), file("${sample}.Aligned.sortedByCoord.out.bam.bai"), file("${sample}.ReadsPerGene.out.tab"), file(vcf) into STAR_output
        tuple val(sample), file("${sample}.Signal.UniqueMultiple.str1.out.wig") into sj_wig_input
        tuple val(sample), file("${sample}.SJ.out.tab") into sj_bed_input

    """
    STAR --genomeDir ${params.STAR_ref_dir} \\
         --readFilesIn ${r1} ${r2}  \\
         --twopassMode Basic \\
         --outReadsUnmapped None \\
         --runThreadN ${task.cpus} \\
         --limitBAMsortRAM ${params.STAR_bam_sort_ram} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMattrRGline ID:${sample} PL:${params.platform} SM:${sample} \\
         --outFileNamePrefix ${sample}. \\
         --quantMode GeneCounts \\
         --outSAMstrandField intronMotif \\
         --outWigNorm None \\
         --outWigType wiggle \\
         --outWigStrand Unstranded \\
         --readFilesCommand gunzip -c

    samtools index ${sample}.Aligned.sortedByCoord.out.bam
    """

}

process gatk_split{
   
    input:
        tuple val(sample) , file(bam), file(bai), file(counts), file(vcf) from STAR_output

    output:
        tuple val(sample), file("${sample}.RG.split.Aligned.sortedByCoord.out.bam"), file("${sample}.RG.split.Aligned.sortedByCoord.out.bai"), file(vcf) into gatk_split_output

    """
    gatk SplitNCigarReads -R ${params.ref} -I ${bam} -O ${sample}.RG.split.Aligned.sortedByCoord.out.bam --create-output-bam-index
    """

}

process sj_wig{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(wig) from sj_wig_input
    
    output: 
        tuple val(sample), file("${sample}.bigwig")

    """
        cut -f1,2 ${params.ref}.fai > contig_file_size.tsv &&
        wigToBigWig ${wig} contig_file_size.tsv ${sample}.bigwig
    """
    
}

process sj_bed{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(sj) from sj_bed_input
    
    output: 
        tuple val(sample), file("${sample}_sj.bed.gz"), file("${sample}_sj.bed.gz.tbi")

    """
        star_sj_tab_to_bed.pl ${sj} | bgzip --stdout > ${sample}_sj.bed.gz
        tabix -p bed ${sample}_sj.bed.gz       
    """
    
}

process GATK_ASE{
   
   publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample) , file(bam),file(bai) , file(vcf) from gatk_split_output

    output:
        tuple val(sample), file("${sample}.ase.vcf"), file("${sample}.GATKASE.csv") into gatk_hc

    script:

    //the input is uncompressed vcf   
    if ( vcf.exists() && !(vcf.getExtension() == "gz" ) )

        """
        bgzip ${vcf}
        tabix -p vcf ${vcf}.gz
        gatk ASEReadCounter  -R ${params.ref} -O ${sample}.GATKASE.csv -I ${bam} -V ${vcf}.gz
        python ${params.bootstrapann} --vcf ${vcf}.gz --ase ${sample}.GATKASE.csv > ${sample}.ase.vcf
        """
    //the input is gzipped vcf
    else if (vcf.exists())

        """
        tabix -p vcf ${vcf}
        gatk ASEReadCounter  -R ${params.ref} -O ${sample}.GATKASE.csv -I ${bam} -V ${vcf}
        python ${params.bootstrapann} --vcf ${vcf}.gz --ase ${sample}.GATKASE.csv > ${sample}.ase.vcf
        """
    // no input vcf, run SNV calling on RNA-seq data
    else

        """
        gatk HaplotypeCaller -R ${params.ref} -I ${bam} -stand-call-conf 10 -O ${sample}.vcf --minimum-mapping-quality 10
        bgzip ${sample}.vcf
        tabix -p vcf ${sample}.vcf.gz
        gatk ASEReadCounter  -R ${params.ref} -O ${sample}.GATKASE.csv -I ${bam} -V ${sample}.vcf.gz
        python ${params.bootstrapann} --vcf ${sample}.vcf.gz --ase ${sample}.GATKASE.csv > ${sample}.ase.vcf
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