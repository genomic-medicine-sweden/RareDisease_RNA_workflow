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

// Initiate references
ch_fasta = file(params.fasta)
ch_star_index = file(params.star_index)
ch_gtf = file(params.gtf)
ch_fai = file(params.fasta + '.fai')
ch_dict = file("${ch_fasta.Parent}/${ch_fasta.SimpleName}.dict")
ch_refflat = params.annotation_refflat ? file(params.annotation_refflat) : file("${ch_gtf.Parent}/${ch_gtf.SimpleName}.refflat")
ch_rrna_intervals = params.rrna_intervals ? file(params.rrna_intervals) : params.rrna_intervals
ch_downsample_regions = params.downsample_regions ? file(params.downsample_regions) : params.downsample_regions
ch_vep_cache = file(params.vep_cache)
ch_reference_cnts = params.reference_count_file ? file(params.reference_count_file) : []

// Setup tempdir - can be overridden in config
params.tmpdir = "${workflow.workDir}/run_tmp/${workflow.runName}"
result = file(params.tmpdir).mkdirs()
if ( !result ) {
    exit 1, "Cannot create tmpdir: $params.tmpdir"
}

// Import processes
include {
    bcftools_compress_and_index as recompress_and_index_vcf;
    bcftools_compress_and_index;
    bcftools_prep_vcf;
    bcftools_variantcall;
    bootstrapann;
    build_fasta_dict;
    build_rrna_intervallist;
    cat_fastq;
    drop_aberrant_expression;
    fastp;
    filter_bam;
    gatk_asereadcounter;
    gatk_haplotypecaller;
    gatk_split;
    generate_SA4drop;
    generate_gene_counts4drop;
    get_rrna_transcripts;
    gffcompare;
    gtf2refflat;
    gunzip_gtf;
    index_bam;
    index_fasta;
    multiqc;
    picard_collectrnaseqmetrics;
    STAR_Aln;
    stringtie;
    untar_star_index;
    vep;
} from './modules/main'

workflow {

    main:

    // Preprocess references
    ch_star_index = untar_star_index(ch_star_index).ifEmpty(ch_star_index)
    ch_gtf = gunzip_gtf(ch_gtf).ifEmpty(ch_gtf)
    ch_fai = ch_fai.isEmpty() ? index_fasta(ch_fasta).collect() : ch_fai
    ch_dict = ch_dict.isEmpty() ? build_fasta_dict(ch_fasta).collect() : ch_dict
    ch_refflat = ch_refflat.isEmpty() ? gtf2refflat(ch_gtf).collect() : ch_refflat
    ch_rrna_intervals = ch_rrna_intervals ?: build_rrna_intervallist(ch_dict, get_rrna_transcripts(ch_gtf).collect()).ifEmpty([])

    // Create channel for multiqc
    ch_multiqc_input = ch_reads.map{ it.first() }

    cat_fastq(ch_reads)

    // Trimming
    fastp(cat_fastq.out)

    // Alignment
    STAR_Aln(fastp.out.reads, ch_star_index)
    index_bam(STAR_Aln.out.bam)
    ch_indexed_bam = STAR_Aln.out.bam.join(index_bam.out)

    // QC
    picard_collectrnaseqmetrics(ch_indexed_bam, ch_refflat, ch_rrna_intervals)

    // Assemble transcripts
    stringtie(ch_indexed_bam, ch_gtf)
    gffcompare(stringtie.out.gtf, ch_gtf)

    // Ready files for DROP tool
    STAR_Aln.out.counts.collect{ sample, cnt_file -> sample }.set{ ch_sample_collect }
    STAR_Aln.out.counts.collect{ sample, cnt_file -> cnt_file }.set{ ch_cnts_collect }

    // Drop
    generate_gene_counts4drop(ch_cnts_collect, ch_sample_collect, ch_gtf, ch_reference_cnts)
    generate_SA4drop( generate_gene_counts4drop.out.processed_gene_counts, params.gtf, [] )
    drop_aberrant_expression(
       generate_SA4drop.out.sample_annotation_drop,
       generate_gene_counts4drop.out.processed_gene_counts,
       ch_reference_cnts,
       ch_fasta,
       ch_gtf
    )

    // ASE subworkflow
    ch_indexed_bam = ch_downsample_regions ? filter_bam(ch_indexed_bam, ch_downsample_regions) : ch_indexed_bam
    gatk_split(ch_indexed_bam, ch_fasta, ch_fai, ch_dict)

    if (params.variantcaller == "gatk") {
        ch_vcf = gatk_haplotypecaller(gatk_split.out, ch_fasta, ch_fai, ch_dict)
    }
    else if (params.variantcaller == "bcftools") {
        ch_vcf = bcftools_variantcall(gatk_split.out, ch_fasta, ch_fai)
    }
    else {
        exit 1, 'Please provide a valid variantcaller [gatk/bcftools].'
    }
    bcftools_compress_and_index(ch_vcf)
    bcftools_prep_vcf(ch_vcf)
    gatk_asereadcounter(bcftools_prep_vcf.out, ch_indexed_bam, ch_fasta, ch_fai, ch_dict, ch_gtf)
    bootstrapann(bcftools_compress_and_index.out, gatk_asereadcounter.out)
    vep(bootstrapann.out, ch_fasta, ch_fai, ch_vep_cache)
    /* TODO: pass the bootstrapann vcf file to recompress and index in the case when
        we don't run VEP. Also the multiqc part will need to be fixed
    */
    recompress_and_index_vcf(vep.out.vcf)

    // Aggregate log files for MultiQC
    ch_multiqc_input = ch_multiqc_input.join(fastp.out.json)
    ch_multiqc_input = ch_multiqc_input.join(STAR_Aln.out.star_multiqc)
    ch_multiqc_input = ch_multiqc_input.join(picard_collectrnaseqmetrics.out)
    ch_multiqc_input = ch_multiqc_input.join(gffcompare.out.multiqc)
    ch_multiqc_input = ch_multiqc_input.join(vep.out.html)
    ch_multiqc_input = ch_multiqc_input.collect{it[1..-1]}
    multiqc(ch_multiqc_input)
}
