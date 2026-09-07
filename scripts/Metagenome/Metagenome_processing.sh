#!/bin/bash

# ============================================================
# Metagenomic assembly, annotation, and coverage analysis
# ============================================================

# 1. Metagenome co-assembly
megahit -t 72 --presets meta-large \
  -1 [sample1_1.fastq],[sample2_1.fastq],[sample3_1.fastq],... \
  -2 [sample1_2.fastq],[sample2_2.fastq],[sample3_2.fastq],... \
  -o [output_dir]


# 2. Structural annotation
prokka --cpus 18 --metagenome --outdir [output_dir] [assembly.fasta]

# 3. KEGG Orthology (KO) annotation
kofam_scan/bin/exec_annotation \
  -o [output] \
  -p [kofam_db/profiles] \
  -k [kofam_db/ko_list] \
  [annotated_proteins.fasta]

# 4. Read sketching
fairy sketch \
  -1 [sample_1.fastq] \
  -2 [sample_2.fastq] \
  -d [sample_sketches]


# 5. Contig coverage estimation
fairy coverage \
  -t 36 \
  [sample_sketches/*.bcsp] \
  [assembly.fasta] \
  -o [assembly.coverage.tsv]
