# ============================================================
# 03_network_analysis.R
#
# 1. Overall network
# 2. Location × compartment networks
# 3. Zi/Pi node roles
#
# Network thresholds:
# Spearman; filter_thres = 0.0005; P < 0.01; |rho| > 0.6
# Module detection: cluster_fast_greedy
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(file2meco)
library(microeco)
library(meconetcomp)
library(magrittr)
library(igraph)

MARKERS <- c("16S", "ITS")

marker_files <- list(
  `16S` = list(
    domain = "Bacteria",
    abundance = "results/16S/03_qiime2/16S.gg2.rarefied-table.qza",
    repseq = "results/16S/03_qiime2/16S.gg2.rarefied-rep-seqs.qza",
    taxonomy = "results/16S/03_qiime2/16S.gg2.taxonomy.qza",
    tree = "results/16S/03_qiime2/16S.gg2.rooted-tree.qza"
  ),
  ITS = list(
    domain = "Fungi",
    abundance = "results/ITS/02_qiime2/ITS.rarefied-table-fungi-no_Incertae_sedis.qza",
    repseq = "results/ITS/02_qiime2/ITS.rarefied-rep-seqs-fungi-no_Incertae_sedis.qza",
    taxonomy = "results/ITS/02_qiime2/ITS.taxonomy-euk-no_Incertae_sedis.qza",
    tree = "results/ITS/02_qiime2/ITS.rooted-tree-fungi.qza"
  )
)

metadata_file <- "metadata.tsv"
location_levels <- c("He", "Ha", "BJ", "WH", "GZ")
latitude_levels <- c("50.26° N", "45.70° N", "40.17° N", "30.58° N", "23.16° N")
cultivar_levels <- c("H43", "H84", "JD17", "TL1", "HC6")
compartment_levels <- c("B", "R", "Ro", "L")

load_amplicon_data <- function(marker) {
  f <- marker_files[[marker]]
  test1 <- qiime2meco(f$abundance, sample_table = metadata_file, taxonomy_table = f$taxonomy,
                      phylo_tree = f$tree, rep_fasta = f$repseq, auto_tidy = TRUE)
  test1$sample_table$compartment <- factor(test1$sample_table$compartment, levels = compartment_levels)
  test1$sample_table$location <- factor(test1$sample_table$location, levels = location_levels)
  test1$sample_table$variety <- factor(test1$sample_table$variety, levels = cultivar_levels)
  test1$sample_table$lat_label <- factor(test1$sample_table$lat_label, levels = latitude_levels)
  if (!"loc-com" %in% colnames(test1$sample_table))
    test1$sample_table$`loc-com` <- paste(test1$sample_table$location, test1$sample_table$compartment, sep = "-")
  return(test1)
}

for (marker in MARKERS) {

  cat("\n====================", marker, "====================\n")
  test1 <- load_amplicon_data(marker)
  outdir <- file.path("results", "downstream", marker, "03_network")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # ==========================================================
  # 1. Overall network
  # ==========================================================

  t1 <- trans_network$new(dataset = test1, cor_method = "spearman", filter_thres = 0.0005)
  t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
  t1$cal_module(method = "cluster_fast_greedy")
  t1$cal_network_attr()
  t1$res_network_attr
  t1$save_network(filepath = file.path(outdir, paste0(marker, "_overall_network.gexf")))


  # ==========================================================
  # 2. Location × compartment networks
  # ==========================================================

  loccom_attr <- data.frame()

  for (loccom in unique(test1$sample_table$`loc-com`)) {
    sub_group <- clone(test1)
    sub_group$sample_table <- subset(sub_group$sample_table, `loc-com` == loccom)
    sub_group$tidy_dataset()

    t1 <- trans_network$new(dataset = sub_group, cor_method = "spearman", filter_thres = 0.0005)
    t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
    t1$cal_module(method = "cluster_fast_greedy")
    t1$cal_network_attr()

    tmp <- as.data.frame(t1$res_network_attr)
    tmp$Combination <- loccom
    loccom_attr <- rbind(loccom_attr, tmp)

    t1$save_network(filepath = file.path(outdir, paste0(marker, "_", loccom, "_network.gexf")))
  }

  write.table(loccom_attr, file.path(outdir, paste0(marker, "_network_attributes_by_location_compartment.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # ==========================================================
  # 3. Zi/Pi node roles
  # ==========================================================

  amp_network <- list()

  for (loccom in unique(test1$sample_table$`loc-com`)) {
    tmp <- clone(test1)
    tmp$sample_table <- subset(tmp$sample_table, `loc-com` == loccom)
    tmp$tidy_dataset()

    tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
    tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
    amp_network[[loccom]] <- tmp
  }

  amp_network %<>% cal_module(undirected_method = "cluster_fast_greedy")
  amp_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

  node_role_df <- data.frame()

  for (loccom in names(amp_network)) {
    tmp <- data.frame(Combination = loccom, amp_network[[loccom]]$res_node_table, check.names = FALSE)
    node_role_df <- rbind(node_role_df, tmp)
  }
  
  write.table(node_role_df, file.path(outdir, paste0(marker, "_node_roles_Zi_Pi.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

}
