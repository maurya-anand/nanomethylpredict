#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    input_samples = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sampleid, file(row.cram), file(row.index)) }

    reference = file(params.reference)
    reference_fai = file(params.reference_fai)
    region_of_interest = file(params.regions_bed)
    EXTRACT_REGION(
        input_samples,
        reference,
        reference_fai,
        region_of_interest,
    )
    VARIANT_CALL(
        EXTRACT_REGION.out.bam,
        reference,
        reference_fai,
    )
    GET_CG_LOCI(
        VARIANT_CALL.out.vcf,
        reference,
        reference_fai,
    )
    chromosomes = GET_CG_LOCI.out.locifile
        .splitText()
        .map { row -> row.split('\t')[0] }
        .unique()
    bam_loci_chr = VARIANT_CALL.out.bam
        .combine(GET_CG_LOCI.out.locifile)
        .combine(chromosomes)
    PREDICT_METHYLATION(
        bam_loci_chr.map { sampleid, bam, bai, _locifile, chr -> tuple(sampleid, bam, bai, chr) },
        reference,
        reference_fai,
        GET_CG_LOCI.out.locifile,
    )
    merged_preds = PREDICT_METHYLATION.out.pred.groupTuple()
    MERGE_METHYLATION(merged_preds)
}

process EXTRACT_REGION {
    publishDir "${params.outdir}/${sampleid}/alignment", mode: 'copy'
    input:
    tuple val(sampleid), path(cram), path(cram_index)
    path reference
    path reference_fai
    path region_of_interest
    output:
    tuple val(sampleid), path("${sampleid}.region.sorted.bam"), path("${sampleid}.region.sorted.bam.bai"), emit: bam
    script:
    def input_flags = cram.toString().endsWith('.cram') ? "-T ${reference} ${cram}" : "${cram}"
    """
    total_threads=\$(nproc)
    available=\$((total_threads - 2))
    samtools_threads=\$(( available * 2 / 10 ))
    [ \$samtools_threads -lt 1 ] && samtools_threads=1
    samtools view -@ \$samtools_threads -b -L ${region_of_interest} ${input_flags} > ${sampleid}.region.bam
    samtools sort -@ \$samtools_threads -o ${sampleid}.region.sorted.bam ${sampleid}.region.bam
    samtools index -@ \$samtools_threads ${sampleid}.region.sorted.bam
    """
}

process VARIANT_CALL {
    publishDir "${params.outdir}/${sampleid}/variant_call", mode: 'copy'
    input:
    tuple val(sampleid), path(bam), path(bai)
    path reference
    path reference_fai
    output:
    tuple val(sampleid), path("${sampleid}.aligned.sorted.haplotagged.bam"), path("${sampleid}.aligned.sorted.haplotagged.bam.bai"), emit: bam
    tuple val(sampleid), path("${sampleid}.vcf.gz"), path("${sampleid}.vcf.gz.tbi"), emit: vcf
    path ("${sampleid}.visual_report.html"), emit: report
    path ("${sampleid}.pepper.margin.deepvariant.log"), emit: main_log
    path ("logs/*.log"), emit: logs
    script:
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    export TMPDIR=\$PWD/tmp
    export TEMP=\$TMPDIR
    export TMP=\$TMPDIR
    export PARALLEL_TMPDIR=\$TMPDIR
    mkdir -p \$TMPDIR
    run_pepper_margin_deepvariant call_variant \
    -b ${sampleid}.aligned.sorted.bam \
    -f ${reference} \
    -o . \
    --gvcf \
    --sample_name ${sampleid} \
    -p ${sampleid} \
    -t \${threads} \
    --dv_sort_by_haplotypes true \
    --keep_intermediate_bam_files true \
    --ont_r10_q20 >> ${sampleid}.pepper.margin.deepvariant.log 2>&1
    mv intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam \
    ${sampleid}.aligned.sorted.haplotagged.bam
    mv intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam.bai \
    ${sampleid}.aligned.sorted.haplotagged.bam.bai
    """
}

process GET_CG_LOCI {
    publishDir "${params.outdir}/${sampleid}/ref", mode: 'copy'
    input:
    tuple val(sampleid), path(vcf), path(idx)
    path reference
    path reference_fai
    output:
    tuple val(sampleid), path("${sampleid}.prediction.locifile.tsv"), emit: locifile
    script:
    """
    awk 'BEGIN{OFS="\t"} /^>/ {if (seq) {gsub("N", "", seq); for (i=1; i<=length(seq)-1; i++) {if (substr(seq, i, 2) == "CG") print chr, i, "5mC", "+"}} chr=substr(\$0, 2); seq=""} /^[^>]/ {seq = seq \$0} END {if (seq) {gsub("N", "", seq); for (i=1; i<=length(seq)-1; i++) {if (substr(seq, i, 2) == "CG") print chr, i, "5mC", "+"}}}' "${reference}" > "${sampleid}.generic.nfl.locifile.tsv"
    awk 'BEGIN{OFS="\t"} {print \$1, \$2-1, \$2, ".", "0", \$4}' "${sampleid}.generic.nfl.locifile.tsv" > "${sampleid}.generic_cpgs.bed"
    bcftools view -f PASS ${vcf} | bcftools query -f '%CHROM\t%POS0\t%POS\n' > ${sampleid}.variants.pass.bed
    bedtools subtract -a "${sampleid}.generic_cpgs.bed" -b "${sampleid}.variants.pass.bed" > "${sampleid}.sample_specific_cpgs.bed"
    awk 'BEGIN{OFS="\t"} {print \$1, \$2+1, "5mC", \$6}' "${sampleid}.sample_specific_cpgs.bed" > "${sampleid}.prediction.locifile.tsv"
    """
}

process PREDICT_METHYLATION {
    publishDir "${params.outdir}/${sampleid}/methylation_pred/subset_${chr}", mode: 'copy'
    input:
    tuple val(sampleid), path(bam), path(bai), val(chr)
    path reference
    path reference_fai
    path locifile
    output:
    tuple val(sampleid), path("${sampleid}_${chr}_methylation_pred.bed"), emit: pred
    script:
    """
    mkdir -p ${sampleid}_prepdata_${chr}
    nfl prepdata -f -p -c ${chr} ${bam} ${reference} ${locifile} ${sampleid}_prepdata_${chr}
    nfl predict /app/NanoFreeLunch.jl-0.28.0/model/${params.nfl_pred_model} ${sampleid}_prepdata_${chr}/forward/Xdata ${sampleid}_${chr}_forward.bed
    nfl predict /app/NanoFreeLunch.jl-0.28.0/model/${params.nfl_pred_model} ${sampleid}_prepdata_${chr}/backward/Xdata ${sampleid}_${chr}_backward.bed
    cat ${sampleid}_${chr}_forward.bed ${sampleid}_${chr}_backward.bed | sort -k1,1 -k2,2n > ${sampleid}_${chr}_methylation_pred.bed
    """
}

process MERGE_METHYLATION {
    publishDir "${params.outdir}/${sampleid}/methylation_pred", mode: 'copy'
    input:
    tuple val(sampleid), path(beds)
    output:
    path ("${sampleid}_methylation_pred.bed"), emit: merged_bed
    script:
    """
    cat ${beds} | sort -k1,1 -k2,2n > ${sampleid}_methylation_pred.bed
    """
}
