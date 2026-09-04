#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# ITS amplicon sequencing processing pipeline
#
# Workflow:
#   FASTQ
#     -> fastp quality control
#     -> Cutadapt primer trimming
#     -> QIIME 2 import
#     -> DADA2 denoising
#     -> UNITE taxonomy assignment
#     -> fungal sequence filtering
#     -> phylogenetic tree construction
#     -> rarefaction and diversity analysis
#
# Reference database:
#   UNITE v10 dynamic (2024-04-04), with entries containing
#   "Incertae sedis" removed before classifier training
# ============================================================


# -------------------- parameters --------------------

# Input directory.
# Each sample has its own subdirectory containing:
#   *.R1.fq.gz
#   *.R2.fq.gz
RAW_DIR="raw_data/ITS"

METADATA="metadata.tsv"

# UNITE classifier
ITS_CLASSIFIER="database/ITS-classifier_ver10_dynamic_all_04.04.2024_dev_no_Incertae_sedis.qza"

# Output directory
OUT_DIR="results/ITS"

# Threads
FASTP_THREADS=54
CUTADAPT_THREADS=54
DADA2_THREADS=68
CLASSIFIER_THREADS=10
TREE_THREADS=64

# Rarefaction depth used in the final analysis
SAMPLING_DEPTH=44680


# -------------------- Output directories --------------------

QC_DIR="${OUT_DIR}/01_quality_control"
QIIME_DIR="${OUT_DIR}/02_qiime2"

mkdir -p "${QC_DIR}" "${QIIME_DIR}"


# ============================================================
# 1. Quality control and primer trimming
# ============================================================

echo "[1/6] Quality control and primer trimming"

for sample_dir in "${RAW_DIR}"/*; do
    sample=$(basename "${sample_dir}")
    R1=("${sample_dir}"/*.R1.fq.gz)
    R2=("${sample_dir}"/*.R2.fq.gz)

    # Quality control
    fastp -w "${FASTP_THREADS}" -i "${R1[0]}" -I "${R2[0]}" -o "${QC_DIR}/${sample}.R1.fq.gz" -O "${QC_DIR}/${sample}.R2.fq.gz"

    # Remove ITS primers from both read directions
    cutadapt \
        --cores "${CUTADAPT_THREADS}" \
        -a '^CTTGGTCATTTAGAGGAAGTAA...GCATCGATGAAGAACGCAGC' \
        -A '^GCTGCGTTCTTCATCGATGC...TTACTTCCTCTAAATGACCAAG' \
        -o "${QC_DIR}/${sample}.R1.cutadapt.fq.gz" \
        -p "${QC_DIR}/${sample}.R2.cutadapt.fq.gz" \
        --untrimmed-output "${QC_DIR}/${sample}.R1.untrimmed.fq.gz" \
        --untrimmed-paired-output "${QC_DIR}/${sample}.R2.untrimmed.fq.gz" \
        "${QC_DIR}/${sample}.R1.fq.gz" \
        "${QC_DIR}/${sample}.R2.fq.gz"

    # Remove intermediate files
    rm \
        "${QC_DIR}/${sample}.R1.fq.gz" \
        "${QC_DIR}/${sample}.R2.fq.gz"
done

rm -f fastp.html fastp.json


# ============================================================
# 2. Create QIIME 2 manifest
# ============================================================

echo "[2/6] Creating QIIME 2 manifest"

MANIFEST="${QIIME_DIR}/ITS_manifest.tsv"
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${MANIFEST}"

for forward_file in "${QC_DIR}"/*.R1.cutadapt.fq.gz; do
    sample=$(basename "${forward_file}" ".R1.cutadapt.fq.gz")
    reverse_file="${QC_DIR}/${sample}.R2.cutadapt.fq.gz"

    forward_path=$(realpath "${forward_file}")
    reverse_path=$(realpath "${reverse_file}")

    echo -e "${sample}\t${forward_path}\t${reverse_path}" >> "${MANIFEST}"
done


# ============================================================
# 3. Import reads and DADA2 denoising
# ============================================================

echo "[3/6] DADA2 denoising"

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${MANIFEST}" \
    --output-path "${QIIME_DIR}/ITS.paired-end-demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${QIIME_DIR}/ITS.paired-end-demux.qza" \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads "${DADA2_THREADS}" \
    --o-table "${QIIME_DIR}/ITS.table.qza" \
    --o-representative-sequences "${QIIME_DIR}/ITS.rep-seqs.qza" \
    --o-denoising-stats "${QIIME_DIR}/ITS.denoising-stats.qza"


# ============================================================
# 4. UNITE taxonomy assignment and fungal filtering
# ============================================================

echo "[4/6] UNITE taxonomy assignment and fungal filtering"

qiime feature-classifier classify-sklearn \
    --p-n-jobs "${CLASSIFIER_THREADS}" \
    --i-classifier "${ITS_CLASSIFIER}" \
    --i-reads "${QIIME_DIR}/ITS.rep-seqs.qza" \
    --o-classification "${QIIME_DIR}/ITS.taxonomy-euk-no_Incertae_sedis.qza"

# Retain fungal ASVs only
qiime taxa filter-table \
    --i-table "${QIIME_DIR}/ITS.table.qza" \
    --i-taxonomy "${QIIME_DIR}/ITS.taxonomy-euk-no_Incertae_sedis.qza" \
    --p-include 'k__Fungi' \
    --o-filtered-table "${QIIME_DIR}/ITS.table-fungi-no_Incertae_sedis.qza"

qiime feature-table filter-seqs \
    --i-data "${QIIME_DIR}/ITS.rep-seqs.qza" \
    --i-table "${QIIME_DIR}/ITS.table-fungi-no_Incertae_sedis.qza" \
    --o-filtered-data "${QIIME_DIR}/ITS.rep-seqs-fungi-no_Incertae_sedis.qza"


# ============================================================
# 5. Phylogenetic tree
# ============================================================

echo "[5/6] Constructing phylogenetic tree"

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "${QIIME_DIR}/ITS.rep-seqs-fungi-no_Incertae_sedis.qza" \
    --p-n-threads "${TREE_THREADS}" \
    --o-alignment "${QIIME_DIR}/ITS.aligned-rep-seqs-fungi.qza" \
    --o-masked-alignment "${QIIME_DIR}/ITS.masked-aligned-rep-seqs-fungi.qza" \
    --o-tree "${QIIME_DIR}/ITS.unrooted-tree-fungi.qza" \
    --o-rooted-tree "${QIIME_DIR}/ITS.rooted-tree-fungi.qza"


# ============================================================
# 6. Rarefaction and diversity analysis
# ============================================================

echo "[6/6] Rarefaction and diversity analysis"

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${QIIME_DIR}/ITS.rooted-tree-fungi.qza" \
    --i-table "${QIIME_DIR}/ITS.table-fungi-no_Incertae_sedis.qza" \
    --p-sampling-depth "${SAMPLING_DEPTH}" \
    --m-metadata-file "${METADATA}" \
    --output-dir "${QIIME_DIR}/core-metrics-results"


# Save the rarefied ASV table used for downstream analyses
cp \
    "${QIIME_DIR}/core-metrics-results/rarefied_table.qza" \
    "${QIIME_DIR}/ITS.rarefied-table-fungi-no_Incertae_sedis.qza"


# Filter representative sequences according to the rarefied ASV table
qiime feature-table filter-seqs \
    --i-data "${QIIME_DIR}/ITS.rep-seqs-fungi-no_Incertae_sedis.qza" \
    --i-table "${QIIME_DIR}/ITS.rarefied-table-fungi-no_Incertae_sedis.qza" \
    --o-filtered-data "${QIIME_DIR}/ITS.rarefied-rep-seqs-fungi-no_Incertae_sedis.qza"


echo "ITS amplicon processing completed successfully."
