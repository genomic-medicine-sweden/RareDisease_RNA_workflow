
// Processes used in main workflow

process untar_star_index{

    input:
        path star_index

    output:
        path "${star_index_decompressed}", emit: star_index

    when:
        star_index.getExtension() == 'gz'

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

    when:
        gtf.getExtension() == 'gz'

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

process get_rrna_transcripts{

    input:
        path gtf

    output:
        path('rrna.bed'), emit: rrna_bed

    """
    $baseDir/bin/get_rrna_transcripts ${gtf} > rrna.gtf
    $baseDir/bin/gtf2bed rrna.gtf > rrna.bed
    """
}

process build_rrna_intervallist{

    input:
        path fasta_dict
        path bed

    output:
        path('rrna.interval_list'), emit: rrna_interval_list

    when:
        bed.size() > 1

    """
    gatk BedToIntervalList -INPUT ${bed} -SEQUENCE_DICTIONARY ${fasta_dict} -OUTPUT rrna.interval_list
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
        tuple val(sample), file("${sample}_1.fastq*"), file("${sample}_2.fastq*") , emit: fastq

    script:
    def extension = ( r1.getExtension() == "gz" ) ? ".fastq.gz" : ".fastq"
    def read_1 = sample + "_1" + extension
    def read_2 = sample + "_2" + extension

    """
    cat ${r1} > ${read_1}
    cat ${r2} > ${read_2}
    """
}

process fastp {

    input:
        tuple val(sample), file(r1), file(r2)

    output:
        tuple val(sample), path("${sample}_1.fastp.fastq.gz"), path("${sample}_2.fastp.fastq.gz"), emit: reads
        tuple val(sample), path('*.json')          , emit: json
        tuple val(sample), path('*.html')          , emit: html
        tuple val(sample), path('*.log')           , emit: log

    script:

    """
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${sample}_1.fastp.fastq.gz \\
        --out2 ${sample}_2.fastp.fastq.gz \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html \\
        --detect_adapter_for_pe \\
        --correction \\
        --overrepresentation_analysis \\
        --thread ${task.cpus} \\
        2> ${sample}.fastp.log
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

    input:
	    tuple val(sample), path(r1), path(r2)
        path star_index

    output:
        tuple val(sample), file("${sample}.bam") , emit: bam
        tuple val(sample), file("${sample}.ReadsPerGene.out.tab") , emit : counts
        tuple val(sample), file('*Log.out'), file('*Log.final.out'), file('*Log.progress.out') , emit: star_multiqc

    script:

    def read_cmd = (r1.getExtension() == "gz") ? "--readFilesCommand gunzip -c" : ""

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
         --outSAMstrandField intronMotif \\
         --peOverlapNbasesMin 10 \\
         --peOverlapMMp 0.1 \\
         --chimSegmentMin 12 \\
         --chimJunctionOverhangMin 12 \\
         --chimOutType WithinBAM \\
         $read_cmd

    mv ${sample}.Aligned.sortedByCoord.out.bam ${sample}.bam
    """
}

process index_bam{

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path('*.bai'), emit: bai

    """
    samtools index ${bam}
    """
}

process generate_gene_counts4drop{

	input:
		path counts
        val samples
        path gtf
        path reference_count_file

	output:
		path('processed_geneCounts.tsv'), emit: processed_gene_counts

    when:
        task.ext.when == null || task.ext.when

    script:
        def ref_counts = reference_count_file ? "--ref_count_file $reference_count_file" : ""

	"""
	generate_gene_counts.py \\
		--star $counts \\
		--sample $samples \\
		--strandedness $params.strandedness \\
        $ref_counts \\
		--output processed_geneCounts.tsv \\
        --gtf $gtf \\
	"""
}

process generate_SA4drop{

	input:
		path processed_gene_counts
        path gtf
        path reference_count_file

	output:
		path('sample_annotation.tsv'), emit: sample_annotation_drop

    when:
        task.ext.when == null || task.ext.when

    script:
        def ref_counts = reference_count_file ? "--ref_count_file $reference_count_file" : ""

	"""
	generate_drop_sample_annot.py \\
		--count_file $processed_gene_counts \\
        --gtf \$(basename $gtf) \\
        $ref_counts \\
		--output sample_annotation.tsv \\
	"""
}

process drop_aberrant_expression{

    beforeScript 'TMPDIR=\$PWD'

    input:
        path sample_annotation
        path gene_counts
        path reference_count_file
        path fasta
        path gtf

    output:
        path('config.yaml'), emit: config_drop
        path('output'), emit: drop_ae_out

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    drop init
    generate_drop_config.py \\
        --genome_fasta \$(basename $fasta) \\
        --gtf \$(basename $gtf) \\
        --output config.yaml
    snakemake aberrantExpression --cores ${task.cpus}
    """
}

process picard_collectrnaseqmetrics{

    input:
        tuple val(sample), path(bam), path(bai)
        path(refflat)
        path(rrna_intervals)

    output:
        tuple val(sample), path("${sample}_rna_metrics.txt")

    script:
    def strandedness = ''
    if (params.strandedness == 'forward') {
        strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    } else if (params.strandedness == 'reverse') {
        strandedness = '--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND'
    }
    def rrna = rrna_intervals == [] ? '' : "--RIBOSOMAL_INTERVALS ${rrna_intervals}"

    """
    picard CollectRnaSeqMetrics \\
        --TMP_DIR ${params.tmpdir} \\
        ${strandedness} \\
        ${rrna} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam} \\
        --OUTPUT ${sample}_rna_metrics.txt \\
    """
}

process stringtie{

    input:
        tuple val(sample), path(bam), path(bai)
        path gtf

    output:
        tuple val(sample), path("${sample}_stringtie.gtf"), emit: gtf
	    tuple val(sample), path("${sample}_stringtie.tab"), emit: tab

    script:
    def strandedness = ''
    if (params.strandedness == 'forward') {
        strandedness = '--fr'
    } else if (params.strandedness == 'reverse') {
        strandedness = '--rf'
    }

    """
    stringtie ${bam} -p ${task.cpus} ${strandedness} -G ${gtf} -A ${sample}_stringtie.tab > ${sample}_stringtie.gtf
    """
}

process gffcompare{

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

process filter_bam {

    input:
        tuple val(sample), path(bam), path(bai)
        path regions

    output:
        tuple val(sample), path("*_filtered.bam"), path("*_filtered.bam.bai")

    script:
    def seed_frac = Math.floor(Math.random() * (100 - 1) + 1) + 0.001

    """
    samtools view --threads ${task.cpus} --bam --unoutput non_select.bam --target-file ${regions} ${bam} | \\
    samtools view -s ${seed_frac} --threads ${task.cpus} --bam --output select.bam
    samtools merge --threads ${task.cpus} -o ${sample}_filtered.bam non_select.bam select.bam
    samtools index ${sample}_filtered.bam
    """
}

process gatk_split{

    input:
        tuple val(sample), path(bam), path(bai)
        path fasta
        path fai
        path dict

    output:
        tuple val(sample), path("${bam.baseName}_split.bam"), path("${bam.baseName}_split.bai")

    """
    gatk SplitNCigarReads --tmp-dir ${params.tmpdir} -R ${fasta} -I ${bam} -O ${bam.baseName}_split.bam --create-output-bam-index
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

process bcftools_variantcall{

    input:
        tuple val(sample), path(bam), path(bai)
        path fasta
        path fai

    output:
        tuple val(sample), path("${sample}.vcf")

    script:

    """
    bcftools mpileup --fasta-ref ${fasta} --output-type u  --max-depth 20000 ${bam} | \\
    bcftools call --pval-threshold 0.01  -mv --output-type v --threads ${task.cpus} --output ${sample}.vcf
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
    bcftools view --threads ${task.cpus} --genotype het --max-alleles 2 --min-alleles 2 --types snps -O z -o ${sample}_biallelic.vcf.gz ${sample}.vcf
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

process vep{

    input:
        tuple val(sample), path(vcf)
        path fasta
        path fai
        path cache

    output:
        tuple val(sample), path("${sample}_ann.vcf"), emit: vcf
        tuple val(sample), path("*_ann.html"), emit: html

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''

    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample}_ann.vcf \\
        --fork ${task.cpus} \\
        --dir_cache ${cache} \\
        --cache \\
        --vcf \\
        --stats_file ${sample}_ann.html \\
        --merged \\
        --compress_output bgzip \\
        --offline \\
        ${args}
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
