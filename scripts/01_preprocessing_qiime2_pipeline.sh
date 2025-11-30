#!/bin/bash
# =============================================================================
# 16S rRNA Gene Sequencing Analysis Pipeline for Illumina MiSeq Paired-End Data
# =============================================================================
#
# Study: Salivary estrogens and oral microbiota in aging women
# Platform: Illumina MiSeq
# Target region: 16S rRNA gene V3-V4
# Database: Human Oral Microbiome Database (HOMD) v15.23
#
# Pipeline steps:
#   1. Quality control with FastQC
#   2. Adapter trimming and quality filtering with Cutadapt
#   3. Import into QIIME2
#   4. Denoising with DADA2 (paired-end)
#   5. Taxonomic classification with HOMD
#   6. Phylogenetic tree construction
#   7. Diversity analysis
#   8. Export results for R analysis
#
# Requirements:
#   - FastQC v0.11.9+
#   - Cutadapt v4.4+
#   - QIIME2 v2023.5+
#   - HOMD database (v15.23)
#
# =============================================================================

set -e  # Exit on error

# =============================================================================
# CONFIGURATION - Modify these parameters for your analysis
# =============================================================================

# Input/Output directories
INPUT_DIR="raw_sequences"           # Directory with raw FASTQ files
TRIMMED_DIR="trimmed_sequences"     # Output directory for trimmed reads
QIIME_DIR="qiime2_output"           # Output directory for QIIME2 results
EXPORT_DIR="exported_results"       # Final exported results for R

# Sample manifest file (tab-separated: sample-id, forward-filepath, reverse-filepath)
MANIFEST="manifest.tsv"

# Metadata file
METADATA="metadata.tsv"

# Reference database paths (HOMD v15.23)
HOMD_SEQS="homd_v15.23_seqs.qza"
HOMD_TAX="homd_v15.23_tax.qza"

# Cutadapt parameters (based on quality analysis)
TRIM_FRONT_F=17                     # Trim from 5' end of forward reads
TRIM_FRONT_R=0                      # Trim from 5' end of reverse reads
TRIM_TAIL_F=12                      # Trim from 3' end of forward reads
TRIM_TAIL_R=95                      # Trim from 3' end of reverse reads
MIN_LENGTH=100                      # Minimum read length after trimming
QUALITY_CUTOFF=20                   # Quality score cutoff
MAX_EE=0.1                          # Maximum expected errors

# DADA2 parameters (based on quality profiles)
TRUNC_LEN_F=288                     # Truncation length for forward reads
TRUNC_LEN_R=205                     # Truncation length for reverse reads
TRIM_LEFT_F=0                       # Additional trim from left (forward)
TRIM_LEFT_R=0                       # Additional trim from left (reverse)
MAX_EE_F=2                          # Max expected errors (forward)
MAX_EE_R=2                          # Max expected errors (reverse)
TRUNC_Q=2                           # Truncate at first base with quality <= Q

# Diversity analysis parameters
SAMPLING_DEPTH=5000                 # Rarefaction depth for diversity metrics

# Number of threads
THREADS=8

# =============================================================================
# STEP 0: Setup and directory creation
# =============================================================================

echo "=== Step 0: Setting up directories ==="

mkdir -p ${TRIMMED_DIR}
mkdir -p ${QIIME_DIR}/{denoising,taxonomy,phylogeny,diversity,exports}
mkdir -p ${EXPORT_DIR}

# =============================================================================
# STEP 1: Quality Control with FastQC
# =============================================================================

echo "=== Step 1: Quality control with FastQC ==="

mkdir -p fastqc_reports/raw
mkdir -p fastqc_reports/trimmed

# Run FastQC on raw reads
fastqc -t ${THREADS} -o fastqc_reports/raw ${INPUT_DIR}/*.fastq

# Generate MultiQC report
multiqc fastqc_reports/raw -o fastqc_reports/raw

echo "Review FastQC reports in fastqc_reports/raw/ before proceeding"
echo "Adjust TRIM and TRUNC parameters based on quality profiles if needed"

# =============================================================================
# STEP 2: Adapter Trimming and Quality Filtering with Cutadapt
# =============================================================================

echo "=== Step 2: Trimming with Cutadapt ==="

# Process each sample (paired-end)
for forward_file in ${INPUT_DIR}/*_1.fastq; do
    # Extract sample name
    filename=$(basename "$forward_file" _1.fastq)
    reverse_file="${INPUT_DIR}/${filename}_2.fastq"

    echo "Processing sample: ${filename}"

    # Run Cutadapt
    cutadapt \
        --cores ${THREADS} \
        -q ${QUALITY_CUTOFF} \
        -u ${TRIM_FRONT_F} \
        -u -${TRIM_TAIL_F} \
        -U ${TRIM_FRONT_R} \
        -U -${TRIM_TAIL_R} \
        --minimum-length ${MIN_LENGTH} \
        --max-expected-errors ${MAX_EE} \
        --pair-filter=any \
        -o "${TRIMMED_DIR}/${filename}_1.trim.fastq" \
        -p "${TRIMMED_DIR}/${filename}_2.trim.fastq" \
        "${forward_file}" \
        "${reverse_file}"
done

# Run FastQC on trimmed reads
fastqc -t ${THREADS} -o fastqc_reports/trimmed ${TRIMMED_DIR}/*.trim.fastq
multiqc fastqc_reports/trimmed -o fastqc_reports/trimmed

echo "Trimming complete. Review FastQC reports in fastqc_reports/trimmed/"

# =============================================================================
# STEP 3: Import into QIIME2
# =============================================================================

echo "=== Step 3: Importing data into QIIME2 ==="

# Activate QIIME2 environment (uncomment if needed)
# conda activate qiime2-2023.5

# Import paired-end reads using manifest file
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path ${MANIFEST} \
    --output-path ${QIIME_DIR}/paired-end-demux.qza \
    --input-format PairedEndFastqManifestPhred33V2

# Generate visualization of demultiplexed sequences
qiime demux summarize \
    --i-data ${QIIME_DIR}/paired-end-demux.qza \
    --o-visualization ${QIIME_DIR}/paired-end-demux.qzv

echo "Import complete. View paired-end-demux.qzv to verify quality profiles"

# =============================================================================
# STEP 4: Denoising with DADA2 (Paired-End)
# =============================================================================

echo "=== Step 4: Denoising with DADA2 ==="

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ${QIIME_DIR}/paired-end-demux.qza \
    --p-trunc-len-f ${TRUNC_LEN_F} \
    --p-trunc-len-r ${TRUNC_LEN_R} \
    --p-trim-left-f ${TRIM_LEFT_F} \
    --p-trim-left-r ${TRIM_LEFT_R} \
    --p-max-ee-f ${MAX_EE_F} \
    --p-max-ee-r ${MAX_EE_R} \
    --p-trunc-q ${TRUNC_Q} \
    --p-n-threads ${THREADS} \
    --o-table ${QIIME_DIR}/denoising/dada2-table.qza \
    --o-representative-sequences ${QIIME_DIR}/denoising/dada2-rep-seqs.qza \
    --o-denoising-stats ${QIIME_DIR}/denoising/dada2-stats.qza

# Generate visualizations
qiime feature-table summarize \
    --i-table ${QIIME_DIR}/denoising/dada2-table.qza \
    --o-visualization ${QIIME_DIR}/denoising/dada2-table-summary.qzv \
    --m-sample-metadata-file ${METADATA}

qiime feature-table tabulate-seqs \
    --i-data ${QIIME_DIR}/denoising/dada2-rep-seqs.qza \
    --o-visualization ${QIIME_DIR}/denoising/dada2-rep-seqs.qzv

qiime metadata tabulate \
    --m-input-file ${QIIME_DIR}/denoising/dada2-stats.qza \
    --o-visualization ${QIIME_DIR}/denoising/dada2-stats.qzv

echo "Denoising complete. Check dada2-stats.qzv for read retention statistics"

# =============================================================================
# STEP 5: Taxonomic Classification with HOMD
# =============================================================================

echo "=== Step 5: Taxonomic classification with HOMD ==="

# Classify sequences using VSEARCH consensus classifier
qiime feature-classifier classify-consensus-vsearch \
    --i-query ${QIIME_DIR}/denoising/dada2-rep-seqs.qza \
    --i-reference-reads ${HOMD_SEQS} \
    --i-reference-taxonomy ${HOMD_TAX} \
    --p-perc-identity 0.99 \
    --p-threads ${THREADS} \
    --o-classification ${QIIME_DIR}/taxonomy/taxonomy-homd.qza \
    --o-search-results ${QIIME_DIR}/taxonomy/vsearch-results.qza

# Generate taxonomy visualization
qiime metadata tabulate \
    --m-input-file ${QIIME_DIR}/taxonomy/taxonomy-homd.qza \
    --o-visualization ${QIIME_DIR}/taxonomy/taxonomy-homd.qzv

# Generate taxa bar plots
qiime taxa barplot \
    --i-table ${QIIME_DIR}/denoising/dada2-table.qza \
    --i-taxonomy ${QIIME_DIR}/taxonomy/taxonomy-homd.qza \
    --m-metadata-file ${METADATA} \
    --o-visualization ${QIIME_DIR}/taxonomy/taxa-barplot.qzv

echo "Taxonomic classification complete"

# =============================================================================
# STEP 6: Phylogenetic Tree Construction
# =============================================================================

echo "=== Step 6: Building phylogenetic tree ==="

# Align sequences with MAFFT
qiime alignment mafft \
    --i-sequences ${QIIME_DIR}/denoising/dada2-rep-seqs.qza \
    --o-alignment ${QIIME_DIR}/phylogeny/aligned-rep-seqs.qza \
    --p-n-threads ${THREADS}

# Mask highly variable positions
qiime alignment mask \
    --i-alignment ${QIIME_DIR}/phylogeny/aligned-rep-seqs.qza \
    --o-masked-alignment ${QIIME_DIR}/phylogeny/masked-aligned-rep-seqs.qza

# Build tree with FastTree
qiime phylogeny fasttree \
    --i-alignment ${QIIME_DIR}/phylogeny/masked-aligned-rep-seqs.qza \
    --o-tree ${QIIME_DIR}/phylogeny/unrooted-tree.qza \
    --p-n-threads ${THREADS}

# Root tree at midpoint
qiime phylogeny midpoint-root \
    --i-tree ${QIIME_DIR}/phylogeny/unrooted-tree.qza \
    --o-rooted-tree ${QIIME_DIR}/phylogeny/rooted-tree.qza

echo "Phylogenetic tree construction complete"

# =============================================================================
# STEP 7: Diversity Analysis
# =============================================================================

echo "=== Step 7: Diversity analysis ==="

# Generate alpha rarefaction curves
qiime diversity alpha-rarefaction \
    --i-table ${QIIME_DIR}/denoising/dada2-table.qza \
    --i-phylogeny ${QIIME_DIR}/phylogeny/rooted-tree.qza \
    --p-max-depth ${SAMPLING_DEPTH} \
    --m-metadata-file ${METADATA} \
    --o-visualization ${QIIME_DIR}/diversity/alpha-rarefaction.qzv

# Core metrics (phylogenetic)
qiime diversity core-metrics-phylogenetic \
    --i-table ${QIIME_DIR}/denoising/dada2-table.qza \
    --i-phylogeny ${QIIME_DIR}/phylogeny/rooted-tree.qza \
    --p-sampling-depth ${SAMPLING_DEPTH} \
    --m-metadata-file ${METADATA} \
    --p-n-jobs-or-threads ${THREADS} \
    --output-dir ${QIIME_DIR}/diversity/core-metrics

# Additional alpha diversity metrics
qiime diversity alpha \
    --i-table ${QIIME_DIR}/denoising/dada2-table.qza \
    --p-metric 'simpson' \
    --o-alpha-diversity ${QIIME_DIR}/diversity/simpson-vector.qza

qiime diversity alpha \
    --i-table ${QIIME_DIR}/denoising/dada2-table.qza \
    --p-metric 'dominance' \
    --o-alpha-diversity ${QIIME_DIR}/diversity/dominance-vector.qza

echo "Diversity analysis complete"

# =============================================================================
# STEP 8: Export Results for R Analysis
# =============================================================================

echo "=== Step 8: Exporting results ==="

# Export feature table
qiime tools export \
    --input-path ${QIIME_DIR}/denoising/dada2-table.qza \
    --output-path ${EXPORT_DIR}/feature-table

# Convert to TSV
biom convert \
    -i ${EXPORT_DIR}/feature-table/feature-table.biom \
    -o ${EXPORT_DIR}/feature-table/feature-table.tsv \
    --to-tsv

# Export taxonomy
qiime tools export \
    --input-path ${QIIME_DIR}/taxonomy/taxonomy-homd.qza \
    --output-path ${EXPORT_DIR}/taxonomy

# Export phylogenetic tree
qiime tools export \
    --input-path ${QIIME_DIR}/phylogeny/rooted-tree.qza \
    --output-path ${EXPORT_DIR}/phylogeny

# Export alpha diversity metrics
for metric in faith_pd shannon evenness observed_features simpson dominance; do
    if [ -f "${QIIME_DIR}/diversity/core-metrics/${metric}_vector.qza" ]; then
        qiime tools export \
            --input-path ${QIIME_DIR}/diversity/core-metrics/${metric}_vector.qza \
            --output-path ${EXPORT_DIR}/alpha-diversity/${metric}
    fi
done

# Export Simpson (if generated separately)
if [ -f "${QIIME_DIR}/diversity/simpson-vector.qza" ]; then
    qiime tools export \
        --input-path ${QIIME_DIR}/diversity/simpson-vector.qza \
        --output-path ${EXPORT_DIR}/alpha-diversity/simpson
fi

# Export beta diversity matrices
for metric in bray_curtis jaccard weighted_unifrac unweighted_unifrac; do
    if [ -f "${QIIME_DIR}/diversity/core-metrics/${metric}_distance_matrix.qza" ]; then
        qiime tools export \
            --input-path ${QIIME_DIR}/diversity/core-metrics/${metric}_distance_matrix.qza \
            --output-path ${EXPORT_DIR}/beta-diversity/${metric}
    fi
done

echo "=== Pipeline complete ==="
echo "Results exported to: ${EXPORT_DIR}/"
echo ""
echo "Output files for R analysis:"
echo "  - Feature table: ${EXPORT_DIR}/feature-table/feature-table.tsv"
echo "  - Taxonomy: ${EXPORT_DIR}/taxonomy/taxonomy.tsv"
echo "  - Phylogenetic tree: ${EXPORT_DIR}/phylogeny/tree.nwk"
echo "  - Alpha diversity: ${EXPORT_DIR}/alpha-diversity/"
echo "  - Beta diversity: ${EXPORT_DIR}/beta-diversity/"
