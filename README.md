# NanoMethylPredict

> [!WARNING]
>
> This pipeline is under active development. It is not ready for production use.

A Nextflow pipeline for predicting DNA methylation from Oxford Nanopore Technologies (ONT) sequencing data **without requiring native methylation calls**. NanoMethylPredict integrates variant calling with PEPPER-Margin-DeepVariant and methylation prediction using [NanoFreeLunch.jl](https://gitee.com/zhixingfeng/NanoFreeLunch.jl) to get sample-specific methylation predictions from standard ONT basecalled data.

> [!NOTE]
>
> **NanoMethylPredict leverages NanoFreeLunch to enable methylation prediction from ONT data that lacks native methylation calls.** This is particularly valuable for:
>
> - **Legacy datasets**: Analyze older ONT data sequenced and basecalled without methylation-aware basecalling.
> - **Re-analysis**: Extract methylation information from existing BAM/CRAM files without re-basecalling.
>
> **Article** : Zhixing Feng, Chenxi Zhang, Shuo Jin, Jiale Niu, Huijuan Feng, Quantitative detection of DNA methylation from nanopore sequencing data without raw signals, GigaScience, Volume 14, 2025, giaf113, [https://doi.org/10.1093/gigascience/giaf113](https://doi.org/10.1093/gigascience/giaf113)

## Features

- **No methylation calls required** - Works with standard ONT basecalled data containing base quality values (Guppy, Dorado).
- **Variant-aware methylation prediction** - excludes CpG sites affected by variants.

## Requirements

- **Nextflow** version 22.10.0 or higher
- **Container runtime**: Docker, Singularity, or Apptainer

## Quick Start

### 1. Prepare your sample sheet (CSV)

Create a CSV file with the following columns:

```csv
sampleid,cram,index
sample1,/path/to/sample1.cram,/path/to/sample1.cram.crai
sample2,/path/to/sample2.bam,/path/to/sample2.bam.bai
```

### 2. Prepare target regions (BED)

Create a BED file defining regions of interest:

```tsv
chr1    1000000    2000000
chr2    5000000    6000000
```

### 3. Configure parameters in `nextflow.config`

```groovy
params {
    sample_sheet = "/path/to/samples.csv"
    reference = "/path/to/reference.fasta"
    reference_fai = "/path/to/reference.fasta.fai"
    regions_bed = "/path/to/regions.bed"
    nfl_pred_model = "r10_dorado0.model"
    outdir = "/path/to/output"
}
```

> [!NOTE]
>
> Check the [NanoFreeLunch.jl documentation](https://gitee.com/zhixingfeng/NanoFreeLunch.jl) for generating models from scratch. Pre-trained methylation models are available in the [model directory](https://gitee.com/zhixingfeng/NanoFreeLunch.jl/tree/main/model)
> Example model names: `r10_dorado0.model`, `r9_guppy3_or_4.model`, `r9_guppy5_or_6.model`.

### 4. Run the pipeline

**Using Docker:**

```bash
nextflow run main.nf -profile docker
```

## Workflow Overview

The pipeline consists of the following main steps:

1. **Region Extraction** (`EXTRACT_REGION`):

   - Extracts reads overlapping target regions from input CRAM/BAM files.
   - Handles both CRAM (with reference) and BAM formats.
   - Output: Sorted BAM with index for target regions.

2. **Variant Calling and Phasing** (`VARIANT_CALL`):

   - Calls variants and haplotags reads using PEPPER-Margin-DeepVariant.
   - Phases variants and assigns haplotype tags to reads.
   - Optimized for ONT R10 chemistry with Q20+ basecalling.
   - Output: Haplotagged BAM, VCF with variants, and QC reports.

3. **CpG Loci Identification** (`GET_CG_LOCI`):

   - Identifies all CpG dinucleotides in the reference genome.
   - Filters out CpG sites overlapping PASS variants to avoid confounding methylation predictions.
   - Creates sample-specific loci file accounting for genetic variants.
   - Output: Tab-separated loci file for methylation prediction.

4. **Methylation Prediction** (`PREDICT_METHYLATION`):

   - Filters reads by quality (MAPQ ≥ 1, removes unmapped reads).
   - Processes each chromosome in parallel for efficiency.
   - Runs NanoFreeLunch prediction for both forward and reverse strands.
   - Uses deep learning models trained on specific ONT chemistry and basecaller versions.
   - Output: Per-chromosome BED files with methylation probabilities.

5. **Merge Predictions** (`MERGE_METHYLATION`):
   - Combines per-chromosome predictions into a single genome-wide file.
   - Sorts by chromosome and position for downstream analysis.
   - Output: Final merged BED file with all methylation predictions.

### Output Files

- **alignment/**: Region-extracted BAM files
- **variant_call/**:
  - Phased, haplotype-tagged BAM
  - Variant calls (VCF)
  - QC reports and logs
- **ref/**: Sample-specific CpG loci file (variant-filtered)
- **methylation_pred/**:
  - Per-chromosome methylation predictions (BED format)
  - Merged genome-wide methylation predictions

## Output Directory Structure

All outputs are organized **per sample**:

```text
outdir/
└── {sampleid}/
    ├── alignment/
    │   ├── {sampleid}.region.sorted.bam
    │   └── {sampleid}.region.sorted.bam.bai
    ├── variant_call/
    │   ├── {sampleid}.aligned.sorted.haplotagged.bam
    │   ├── {sampleid}.aligned.sorted.haplotagged.bam.bai
    │   ├── {sampleid}.vcf.gz
    │   ├── {sampleid}.vcf.gz.tbi
    │   ├── {sampleid}.visual_report.html
    │   ├── {sampleid}.pepper.margin.deepvariant.log
    │   └── logs/
    ├── ref/
    │   └── {sampleid}.prediction.locifile.tsv
    ├── methylation_pred/
    │   ├── subset_{chr1}/
    │   │   └── {sampleid}_{chr1}_methylation_pred.bed
    │   ├── subset_{chr2}/
    │   │   └── {sampleid}_{chr2}_methylation_pred.bed
    │   └── {sampleid}_methylation_pred.bed (merged)
```

## Components

| Component                                                            | Version |
| -------------------------------------------------------------------- | ------- |
| [bedtools](https://bedtools.readthedocs.io/)                         | 2.30.0  |
| [NanoFreeLunch](https://gitee.com/zhixingfeng/NanoFreeLunch.jl)      | 0.28.0  |
| [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) | r0.8    |
| [samtools](http://www.htslib.org/)                                   | 1.13    |

## References

- Shafin, K., Pesout, T., Chang, PC. et al. Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads. Nat Methods 18, 1322–1332 (2021). [https://doi.org/10.1038/s41592-021-01299-w](https://doi.org/10.1038/s41592-021-01299-w)
- Zhixing Feng, Chenxi Zhang, Shuo Jin, Jiale Niu, Huijuan Feng, Quantitative detection of DNA methylation from nanopore sequencing data without raw signals, GigaScience, Volume 14, 2025, giaf113, [https://doi.org/10.1093/gigascience/giaf113](https://doi.org/10.1093/gigascience/giaf113)

## License

This pipeline is licensed under the **GNU General Public License v3.0** due to its use of GPL-licensed components.

- **NanoFreeLunch** (GPL v3)
- **PEPPER-Margin-DeepVariant** (MIT)

Users must comply with the respective licenses of these tools when using this pipeline.
