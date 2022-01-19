#!/usr/bin/dev nextflow

nextflow.enable.dsl=2
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

     ch_reads = Channel.of( [params.sample, file(params.r1), file(params.r2)] )


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

// Initiate reference channels
ch_fasta = file(params.fasta)
ch_star_index = file(params.star_index)
ch_gtf = file(params.gtf)
ch_fai = file(params.fasta + '.fai')
ch_dict = file(params.fasta).getParent()/file(params.fasta).getSimpleName() + '.dict'
fasta_dict = file(params.fasta).getParent()/file(params.fasta).getSimpleName() + '.dict'
ch_refflat = params.annotation_refflat
    ? Channel.fromPath( params.annotation_refflat, checkIfExists: true )
    : Channel.empty()

ch_multiqc_input = ch_reads.map{ it.first() }

// Setup tempdir - can be overridden in config
params.tmpdir = "${workflow.workDir}/run_tmp/${workflow.runName}"
file(params.tmpdir).mkdir()

process untar_star_index{

    input:
        path star_index

    output:
        path "${star_index_decompressed}", emit: star_index

    script:

    star_index_decompressed = star_index.getSimpleName()

    """
    tar -xzvf ${star_index} --no-same-owner
    """
}

process gunzip_gtf{

    input:
        path gtf

    output:
        path('*.gtf') , emit: gtf

    """
    gunzip -f ${gtf}
    """
}

process index_fasta{

    input:
        path fasta

    output:
        path('*.fai'), emit: fai

    """
    samtools faidx ${fasta}
    """
}

process build_fasta_dict{

    input:
        path fasta

    output:
        path('*.dict'), emit: dict

    """
    gatk CreateSequenceDictionary -R ${fasta}
    """
}

process gtf2refflat{

    input:
        path gtf

    output:
        path('*.refflat'), emit: refflat

    script:

    def genepred = gtf.getSimpleName() + '.genepred'
    def refflat = gtf.getSimpleName() + '.refflat'

    """
    gtfToGenePred -genePredExt -geneNameAsName2 ${gtf} ${genepred}
    paste ${genepred} ${genepred} | cut -f12,16-25 > ${refflat}
    """
}

process cat_fastq{

    input:
        tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), file("${sample}_1.fastq"), file("${sample}_2.fastq") , emit: fastq

    """
    zcat ${r1} > ${sample}_1.fastq
    zcat ${r2} > ${sample}_2.fastq

    """
}

process fastqc{

    input:
	    tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), path ("*1_fastqc.zip"), path ("*2_fastqc.zip"), emit: zip

    script:

    def R1 = r1.getName()
    def R2 = r2.getName()

    """
    [ ! -f ${R1} ] && ln -s ${r1} ${R1}
    [ ! -f ${R2} ] && ln -s ${r2} ${R2}
    fastqc --threads ${task.cpus} ${R1} ${R2}
    """
}

process STAR_Aln{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
	    tuple val(sample), path(r1), path(r2)
        path star_index

    output:
        tuple val(sample), file("${sample}.Aligned.sortedByCoord.out.bam") , emit: bam
        tuple val(sample), file("${sample}.ReadsPerGene.out.tab") , emit : counts
        tuple val(sample), file('*Log.out'), file('*Log.final.out'), file('*Log.progress.out') , emit: star_multiqc

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

process index_bam{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path('*.bai'), emit: bai

    """
    samtools index ${bam}
    """
}

process picard_collectrnaseqmetrics{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(bam), path(bai)
        path(refflat)

    output:
        tuple val(sample), path("${sample}_rna_metrics.txt")

    script:
    """
    picard CollectRnaSeqMetrics \\
        --TMP_DIR ${params.tmpdir} \\
        --STRAND_SPECIFICITY ${params.strandedness} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam} \\
        --OUTPUT ${sample}_rna_metrics.txt \\
    """
}

process stringtie{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(bam), path(bai)
        path gtf

    output:
        tuple val(sample), path("${sample}_stringtie.gtf"), emit: gtf

    script:

    """
    stringtie ${bam} -p ${task.cpus} ${params.stranded} -G ${gtf} > ${sample}_stringtie.gtf
    """
}

process gffcompare{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(stringtie_gtf)
        path gtf

    output:
        tuple val(sample), file("${sample}_stringtie.annotated.gtf"), emit: gtf
        tuple val(sample), file("${sample}_stringtie.stats"), emit: multiqc

    script:
    """
    gffcompare -r ${gtf} -o ${sample}_stringtie ${stringtie_gtf}
    """
}

process gatk_split{

    input:
        tuple val(sample), path(bam), path(bai)
        path fasta
        path fai
        path dict

    output:
        tuple val(sample), path("${sample}.RG.split.Aligned.sortedByCoord.out.bam"), path("${sample}.RG.split.Aligned.sortedByCoord.out.bai")

    """
    gatk SplitNCigarReads --tmp-dir ${params.tmpdir} -R ${fasta} -I ${bam} -O ${sample}.RG.split.Aligned.sortedByCoord.out.bam --create-output-bam-index
    """

}

process gatk_haplotypecaller{

    input:
	    tuple val(sample), path(bam), file(bai)
        path fasta
        path fai
        path dict

    output:
        tuple val(sample), path("${sample}.vcf")

    script:
        file(params.tmpdir).mkdir()

    """
    gatk HaplotypeCaller --tmp-dir ${params.tmpdir} -R ${fasta} -I ${bam} -stand-call-conf 10 -O ${sample}.vcf --minimum-mapping-quality 10
    """
}

process bcftools_compress_and_index{

    input:
        tuple val(sample), path(vcf)

    output:
        tuple val(sample), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    """
    bcftools view --output-type z --threads ${task.cpus} ${vcf} > ${vcf}.gz
    bcftools index --tbi ${vcf}.gz
    """
}

process bcftools_prep_vcf{

    input:
        tuple val(sample), path(vcf)

    output:
        tuple val(sample), path("${sample}_biallelic.vcf.gz"), path("${sample}_biallelic.vcf.gz.tbi")

    """
    bcftools view --genotype het --max-alleles 2 --min-alleles 2 --types snps -O z -o ${sample}_biallelic.vcf.gz ${sample}.vcf
    bcftools index --tbi ${sample}_biallelic.vcf.gz
    """
}


process gatk_asereadcounter{

    input:
        tuple val(sample), path(vcf), path(tbi)
	    tuple val(sample), path(bam), path(bai)
        path fasta
        path fai
        path dict

    output:
        tuple val(sample), file("${sample}_ase.csv")

    script:
        file(params.tmpdir).mkdir()

    """
    gatk ASEReadCounter --tmp-dir ${params.tmpdir} -R ${fasta} -O ${sample}_ase.csv -I ${bam} -V ${vcf}
    """
}

process bootstrapann{

    input:
        tuple val(sample), path(vcf), path(tbi)
        tuple val(sample), path(csv)

    output:
        tuple val(sample), path("${sample}_ase.vcf")

    """
    BootstrapAnn.py --vcf ${vcf} --ase ${csv} > ${sample}_ase.vcf
    """
}

process recompress_and_index_vcf{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(vcf)

    output:
        tuple val(sample), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    """
    bcftools view --output-type z --threads ${task.cpus} ${vcf} > ${vcf}.gz
    bcftools index --tbi ${vcf}.gz
    """
}

process multiqc{
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        path(qc_files)

    output:
        path "*multiqc_report.html"
        path "*_data"

    """
    multiqc .
    """
}

workflow {

    main:

    // preprocess references
    if ( ch_star_index.getExtension() == 'gz' ) {

        untar_star_index(ch_star_index)
        .set{ ch_star_index }
    }
    if ( ch_gtf.getExtension() == 'gz' ) {

        gunzip_gtf(ch_gtf)
        .set{ ch_gtf }
    }
    if ( ch_fai.isEmpty() ) {

        index_fasta(ch_fasta)
        .set{ ch_fai}
    }
    if ( ch_dict.isEmpty() ) {

        build_fasta_dict(ch_fasta)
        .set{ ch_dict }
    }
    if (ch_refflat.toList().isEmpty()) {

        gtf2refflat(ch_gtf)
        .set{ ch_refflat }
    }

    // main workflow
    cat_fastq(ch_reads)

    fastqc(cat_fastq.out)
    ch_multiqc_input = ch_multiqc_input.join(fastqc.out.zip)

    STAR_Aln(cat_fastq.out, ch_star_index)
    ch_multiqc_input = ch_multiqc_input.join(STAR_Aln.out.star_multiqc)

    index_bam(STAR_Aln.out.bam)
    ch_indexed_bam = STAR_Aln.out.bam.join(index_bam.out)

    picard_collectrnaseqmetrics(ch_indexed_bam, ch_refflat)
    ch_multiqc_input = ch_multiqc_input.join(picard_collectrnaseqmetrics.out)

    stringtie(ch_indexed_bam, ch_gtf)

    gffcompare(stringtie.out, ch_gtf)
    ch_multiqc_input = ch_multiqc_input.join(gffcompare.out.multiqc)

    gatk_split(ch_indexed_bam, ch_fasta, ch_fai, ch_dict)

    gatk_haplotypecaller(gatk_split.out, ch_fasta, ch_fai, ch_dict)

    bcftools_compress_and_index(gatk_haplotypecaller.out)

    bcftools_prep_vcf(gatk_haplotypecaller.out)

    gatk_asereadcounter(bcftools_prep_vcf.out, ch_indexed_bam, ch_fasta, ch_fai, ch_dict)

    bootstrapann(bcftools_compress_and_index.out, gatk_asereadcounter.out)

    recompress_and_index_vcf(bootstrapann.out)

    // Combine metric output files to one channel for multiqc
    ch_multiqc_input = ch_multiqc_input.collect{it[1..-1]}
    multiqc(ch_multiqc_input)
}