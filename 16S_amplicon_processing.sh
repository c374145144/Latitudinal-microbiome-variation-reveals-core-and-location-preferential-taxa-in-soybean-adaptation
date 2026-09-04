#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 16S rRNA amplicon sequencing processing pipeline
#
# Workflow:
#   FASTQ
#     -> fastp quality control
#     -> Cutadapt primer trimming
#     -> PEAR paired-end read merging
#     -> QIIME 2 import
#     -> DADA2 denoising
#     -> Greengenes2 mapping and taxonomy assignment
#     -> phylogenetic tree construction
#     -> rarefaction and diversity analysis
#
# Reference database:
#   Greengenes2 2022.10
# ============================================================


# -------------------- parameters --------------------

# Input directory.
# Each sample have its own subdirectory containing:
#   *.R1.fq.gz
#   *.R2.fq.gz
RAW_DIR="raw_data/16S"

METADATA="metadata.tsv"

# Greengenes2 reference files
GG2_BACKBONE="database/2022.10.backbone.full-length.fna.qza"
GG2_TAXONOMY="database/2022.10.taxonomy.asv.nwk.qza"

# Output directory
OUT_DIR="results/16S"

# Threads
FASTP_THREADS=8
CUTADAPT_THREADS=32
PEAR_THREADS=32
DADA2_THREADS=54
TREE_THREADS=56
GG2_THREADS=64

# Rarefaction depth used in the final analysis
SAMPLING_DEPTH=25585


# -------------------- Output directories --------------------

QC_DIR="${OUT_DIR}/01_quality_control"
PEAR_DIR="${OUT_DIR}/02_merged_reads"
QIIME_DIR="${OUT_DIR}/03_qiime2"

mkdir -p "${QC_DIR}" "${PEAR_DIR}" "${QIIME_DIR}"


# ============================================================
# 1. Quality control and primer trimming
# ============================================================

echo "[1/7] Quality control and primer trimming"

for sample_dir in "${RAW_DIR}"/*; do
    sample=$(basename "${sample_dir}")
    R1=("${sample_dir}"/*.R1.fq.gz)
    R2=("${sample_dir}"/*.R2.fq.gz)

    # Quality control
    fastp -w "${FASTP_THREADS}" -i "${R1[0]}" -I "${R2[0]}" -o "${QC_DIR}/${sample}.R1.fq.gz" -O "${QC_DIR}/${sample}.R2.fq.gz"

    # Remove 16S primers
    cutadapt \
        -j "${CUTADAPT_THREADS}" \
        -g AACMGGATTAGATACCCKG \
        -G ACGTCATCCCCACCTTCC \
        -m 1 \
        -o "${QC_DIR}/${sample}_qc.R1.fq.gz" \
        -p "${QC_DIR}/${sample}_qc.R2.fq.gz" \
        --too-short-output "${QC_DIR}/${sample}_short.R1.fq.gz" \
        --too-short-paired-output "${QC_DIR}/${sample}_short.R2.fq.gz" \
        --untrimmed-output "${QC_DIR}/${sample}_untrimmed.R1.fq.gz" \
        --untrimmed-paired-output "${QC_DIR}/${sample}_untrimmed.R2.fq.gz" \
        "${QC_DIR}/${sample}.R1.fq.gz" \
        "${QC_DIR}/${sample}.R2.fq.gz"

    # Remove intermediate files
    rm \
        "${QC_DIR}/${sample}.R1.fq.gz" \
        "${QC_DIR}/${sample}.R2.fq.gz"
done

rm -f fastp.html fastp.json


# ============================================================
# 2. Merge paired-end reads using PEAR
# ============================================================

echo "[2/7] Merging paired-end reads"

for forward_file in "${QC_DIR}"/*_qc.R1.fq.gz; do
    sample=$(basename "${forward_file}" "_qc.R1.fq.gz")
    reverse_file="${QC_DIR}/${sample}_qc.R2.fq.gz"
    pear -j "${PEAR_THREADS}" -f "${forward_file}" -r "${reverse_file}" -o "${PEAR_DIR}/${sample}"
done


# ============================================================
# 3. Create QIIME 2 manifest
# ============================================================

echo "[3/7] Creating QIIME 2 manifest"

MANIFEST="${QIIME_DIR}/16S_manifest.tsv"
echo -e "sample-id\tabsolute-filepath\tdirection" > "${MANIFEST}"

for file in "${PEAR_DIR}"/*.assembled.fastq; do
    sample=$(basename "${file}" ".assembled.fastq")
    filepath=$(realpath "${file}")
    echo -e "${sample}\t${filepath}\tforward" >> "${MANIFEST}"
done


# ============================================================
# 4. Import reads and DADA2 denoising
# ============================================================

echo "[4/7] DADA2 denoising"

qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "${MANIFEST}" \
    --output-path "${QIIME_DIR}/16S.assembled-demux.qza" \
    --input-format SingleEndFastqManifestPhred33V2

qiime dada2 denoise-single \
    --p-n-threads "${DADA2_THREADS}" \
    --p-trunc-len 0 \
    --p-trim-left 0 \
    --i-demultiplexed-seqs "${QIIME_DIR}/16S.assembled-demux.qza" \
    --o-table "${QIIME_DIR}/16S.table.qza" \
    --o-representative-sequences "${QIIME_DIR}/16S.rep-seqs.qza" \
    --o-denoising-stats "${QIIME_DIR}/16S.denoising-stats.qza"

# ============================================================
# 5. Greengenes2 taxonomy assignment
# ============================================================

echo "[5/7] Greengenes2 taxonomy assignment"

qiime greengenes2 non-v4-16s \
    --i-table "${QIIME_DIR}/16S.table.qza" \
    --i-sequences "${QIIME_DIR}/16S.rep-seqs.qza" \
    --i-backbone "${GG2_BACKBONE}" \
    --o-mapped-table "${QIIME_DIR}/16S.gg2.table.qza" \
    --o-representatives "${QIIME_DIR}/16S.gg2.rep-seqs.qza" \
    --p-threads "${GG2_THREADS}"

qiime greengenes2 taxonomy-from-table \
    --i-reference-taxonomy "${GG2_TAXONOMY}" \
    --i-table "${QIIME_DIR}/16S.gg2.table.qza" \
    --o-classification "${QIIME_DIR}/16S.gg2.taxonomy.qza"


# ============================================================
# 6. Phylogenetic tree
# ============================================================

echo "[6/7] Constructing phylogenetic tree"

qiime phylogeny align-to-tree-mafft-fasttree \
    --p-n-threads "${TREE_THREADS}" \
    --p-parttree \
    --i-sequences "${QIIME_DIR}/16S.gg2.rep-seqs.qza" \
    --o-alignment "${QIIME_DIR}/16S.gg2.aligned-rep-seqs.qza" \
    --o-masked-alignment "${QIIME_DIR}/16S.gg2.masked-aligned-rep-seqs.qza" \
    --o-tree "${QIIME_DIR}/16S.gg2.unrooted-tree.qza" \
    --o-rooted-tree "${QIIME_DIR}/16S.gg2.rooted-tree.qza"


# ============================================================
# 7. Rarefaction and diversity analysis
# ============================================================

echo "[7/7] Rarefaction and diversity analysis"

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${QIIME_DIR}/16S.gg2.rooted-tree.qza" \
    --i-table "${QIIME_DIR}/16S.gg2.table.qza" \
    --p-sampling-depth "${SAMPLING_DEPTH}" \
    --m-metadata-file "${METADATA}" \
    --output-dir "${QIIME_DIR}/core-metrics-results"


# Save the rarefied ASV table used for downstream analyses
cp \
    "${QIIME_DIR}/core-metrics-results/rarefied_table.qza" \
    "${QIIME_DIR}/16S.gg2.rarefied-table.qza"


# Filter representative sequences according to the rarefied ASV table
qiime feature-table filter-seqs \
    --i-data "${QIIME_DIR}/16S.gg2.rep-seqs.qza" \
    --i-table "${QIIME_DIR}/16S.gg2.rarefied-table.qza" \
    --o-filtered-data "${QIIME_DIR}/16S.gg2.rarefied-rep-seqs.qza"


echo "16S amplicon processing completed successfully."