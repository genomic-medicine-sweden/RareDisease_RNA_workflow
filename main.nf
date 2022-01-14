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
        .set { ch_reads }

}else{
    println "Nope."
    exit 0
}

ch_fasta = Channel.value(params.fasta)
ch_reads.into { ch_cat_reads; ch_sampel_id }
ch_multiqc_input = ch_sampel_id.map{ it.first() }

if (params.star_index.endsWith('.tar.gz')) {

    process untar_star_index{

        input:
            file(star_index_compressed) from file("${params.star_index}")

        output:
            file("${star_index_decompressed}") into ch_star_index

        script:

        star_index_decompressed = star_index_compressed.getSimpleName()

        """
        tar -xzvf ${star_index_compressed} --no-same-owner
        """

    }
} else {

    ch_star_index = Channel.value(params.star_index)
}

if (params.gtf.endsWith('.gz')) {

    process gunzip_gtf{

        input:
            file(gtf_compressed) from file("${params.gtf}")

        output:
            path('*.gtf') into ch_gtf

        """
        gunzip -f ${gtf_compressed}
        """
    }
} else {

    ch_gtf = Channel.value(params.gtf)
}

if (!file(params.fasta + '.fai').exists()) {

    process index_fasta{

        input:
            path(fasta) from ch_fasta

        output:
            path('*.fai') into ch_fai

        """
        samtools faidx "${fasta}"
        """
    }
}
else{

    ch_fai = Channel.value(params.fasta + '.fai')
}

fasta_dict = file(params.fasta).getSimpleName() + '.dict'
if (!file(fasta_dict).exists()) {

    process generate_fasta_dict{

        input:
            path(fasta) from ch_fasta

        output:
            path('*.dict') into ch_dict

        """
        gatk CreateSequenceDictionary -R "${fasta}"
        """
    }
}
else{

    ch_dict= Channel.value(fasta_dict)
}


process cat_fastq{

    input:
        tuple val(sample), file(r1), file(r2) from ch_cat_reads

    output:
        tuple val(sample), file("${sample}_1.fastq"), file("${sample}_2.fastq") into ch_fastq

    """
    zcat ${r1} > ${sample}_1.fastq
    zcat ${r2} > ${sample}_2.fastq

    """
}
ch_fastq.into { ch_reads_align; ch_reads_qc }

process fastqc{

    input:
	    tuple val(sample), val(r1), val(r2) from ch_reads_qc

    output:
        tuple val(sample), file("${r1.simpleName}*_fastqc.zip"), file("${r2.simpleName}*_fastqc.zip") into fastqc_multiqc

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

if (!params.annotation_refflat || file(params.annotation_refflat).isEmpty()) {

    process gtf2refflat{

        input:
            path(gtf) from ch_gtf

        output:
            path('*.refflat') into ch_refflat

        script:

        def genepred = gtf.getSimpleName() + '.genepred'
        def refflat = gtf.getSimpleName() + '.refflat'

        """
        gtfToGenePred -genePredExt -geneNameAsName2 ${gtf} ${genepred}
        paste ${genepred} ${genepred} | cut -f12,16-25 > ${refflat}
        """

    }
} else {

    ch_refflat = Channel.value(params.annotation_refflat)
}


process STAR_Aln{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
	    tuple val(sample), val(r1), val(r2) from ch_reads_align
        path(star_index) from ch_star_index

    output:
        tuple val(sample), file("${sample}.Aligned.sortedByCoord.out.bam") into STAR_output
        tuple val(sample), file("${sample}.ReadsPerGene.out.tab")
        tuple val(sample), file('*Log.out'), file('*Log.final.out'), file('*Log.progress.out') into star_multiqc

    """

    STAR --genomeDir ${star_index} \\
         --readFilesIn ${r1} ${r2} \\
         --twopassMode Basic \\
         --outReadsUnmapped None \\
         --runThreadN ${task.cpus} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMattrRGline ID:${sample} PL:${params.platform} SM:${sample} \\
         --outFileNamePrefix ${sample}. \\
         --quantMode GeneCounts \\
         --outSAMstrandField intronMotif

    """
}
ch_multiqc_input = ch_multiqc_input.join(star_multiqc)

process index_bam{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file("${sample}.Aligned.sortedByCoord.out.bam") from STAR_output

    output:
        tuple val(sample), file("${sample}.Aligned.sortedByCoord.out.bam"), file("${sample}.Aligned.sortedByCoord.out.bam.bai") into ch_mapped_reads

    """
    samtools index ${sample}.Aligned.sortedByCoord.out.bam
    """
}
ch_mapped_reads.into {metrics_input; gatk_split_input; stringtie_input}

process picard_collectrnaseqmetrics{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(bam), path(bai) from metrics_input
        path(refflat) from ch_refflat

    output:
        tuple val(sample), path("${sample}_rna_metrics.txt") into metric_multiqc

    """
    picard CollectRnaSeqMetrics \\
        --STRAND_SPECIFICITY ${params.strandedness} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam} \\
        --OUTPUT ${sample}_rna_metrics.txt \\
    """
}
ch_multiqc_input = ch_multiqc_input.join(metric_multiqc)

process stringtie{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(bam),file(bai) from stringtie_input
        path(gtf) from ch_gtf

    output:
        tuple val(sample), file("${sample}_stringtie.gtf") into stringtie_output

    script:

    """
    stringtie ${bam} -p ${task.cpus} ${params.stranded} -G ${gtf} > ${sample}_stringtie.gtf
    """

}

process gffcompare{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(stringtie_gtf) from stringtie_output
        path(gtf) from ch_gtf

    output:
        tuple val(sample), file("${sample}_stringtie.annotated.gtf") into gffcompare_output
        tuple val(sample), file("${sample}_stringtie.stats") into gffcompare_multiqc

    script:
    """
    gffcompare -r ${gtf} -o ${sample}_stringtie ${stringtie_gtf}
    """
}
ch_multiqc_input = ch_multiqc_input.join(gffcompare_multiqc)

process gatk_split{

    input:
        tuple val(sample) , file(bam), file(bai) from gatk_split_input
        path(fasta) from ch_fasta
        path(fai) from ch_fai
        path(dict) from ch_dict

    output:
        tuple val(sample), file("${sample}.RG.split.Aligned.sortedByCoord.out.bam"), file("${sample}.RG.split.Aligned.sortedByCoord.out.bai") into ch_gatk_split_output

    """
    gatk SplitNCigarReads -R ${fasta} -I ${bam} -O ${sample}.RG.split.Aligned.sortedByCoord.out.bam --create-output-bam-index
    """

}
ch_gatk_split_output.into { ch_gatk_split_hk; ch_gatk_split_ase }


process gatk_haplotypecaller{

    input:
	    tuple val(sample), file(bam),file(bai) from ch_gatk_split_hk
        path(fasta) from ch_fasta
        path(fai) from ch_fai
        path(dict) from ch_dict

    output:
        tuple val(sample), file("${sample}.vcf") into ch_haplotypecaller

    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -stand-call-conf 10 -O ${sample}.vcf --minimum-mapping-quality 10
    """
}

process bcftools_prep_vcf{

    input:
        tuple val(sample), file(vcf) from ch_haplotypecaller

    output:
        tuple val(sample), file("${sample}_biallelic.vcf.gz"), file("${sample}_biallelic.vcf.gz.tbi") into ch_prep_vcf

    """
    bcftools view --genotype het --max-alleles 2 --min-alleles 2 --types snps -O z -o ${sample}_biallelic.vcf.gz ${sample}.vcf
    bcftools index --tbi ${sample}_biallelic.vcf.gz
    """
}
ch_prep_vcf.into { ch_prep_vcf_ase; ch_prep_vcf_bootstrapann }


process gatk_asereadcounter{

    input:
        tuple val(sample), file(vcf), file(tbi) from ch_prep_vcf_ase
	    tuple val(sample) , file(bam),file(bai) from ch_gatk_split_ase
        path(fasta) from ch_fasta
        path(fai) from ch_fai
        path(dict) from ch_dict

    output:
        tuple val(sample), file("${sample}_ase.csv") into ch_asereadcounter

    """
    gatk ASEReadCounter -R ${fasta} -O ${sample}_ase.csv -I ${bam} -V ${vcf}
    """
}

process bootstrapann{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), file(vcf), file(tbi) from ch_prep_vcf_bootstrapann
        tuple val(sample), file(csv) from ch_asereadcounter

    output:
        tuple val(sample), file("${sample}_ase.vcf") into ch_ase_vcf

    """
    BootstrapAnn.py --vcf ${vcf} --ase ${csv} > ${sample}_ase.vcf
    """
}

// Combine metric output files to one channel
ch_multiqc_input = ch_multiqc_input.collect{it[1..-1]}
process multiqc{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        path(qc_files) from ch_multiqc_input

    output:
        path "*multiqc_report.html"
        path "*_data"
    """
    multiqc .
    """
}
