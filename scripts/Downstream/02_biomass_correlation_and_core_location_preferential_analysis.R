# ============================================================
# 02_core_preferential_and_biomass_analysis.R
#
# 16S / ITS manuscript-specific analyses:
# 1. Biomass-associated genera and context dependence
# 2. Correlation strength and abundance-rank percentile
# 3. Geographic distribution of genera
# 4. Core genera and their abundance profiles
# 5. Core vs non-core abundance rank
# 6. Location-preferential genera and overlap
# 7. LEfSe enrichment at home sites
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(file2meco)
library(microeco)
library(aplot)
library(UpSetR)
library(grid)

MARKERS <- c("16S", "ITS")

marker_files <- list(
  `16S` = list(
    abundance = "results/16S/03_qiime2/16S.gg2.rarefied-table.qza",
    repseq = "results/16S/03_qiime2/16S.gg2.rarefied-rep-seqs.qza",
    taxonomy = "results/16S/03_qiime2/16S.gg2.taxonomy.qza",
    tree = "results/16S/03_qiime2/16S.gg2.rooted-tree.qza"
  ),
  ITS = list(
    abundance = "results/ITS/02_qiime2/ITS.rarefied-table-fungi-no_Incertae_sedis.qza",
    repseq = "results/ITS/02_qiime2/ITS.rarefied-rep-seqs-fungi-no_Incertae_sedis.qza",
    taxonomy = "results/ITS/02_qiime2/ITS.taxonomy-euk-no_Incertae_sedis.qza",
    tree = "results/ITS/02_qiime2/ITS.rooted-tree-fungi.qza"
  )
)

metadata_file <- "metadata.tsv"
biomass_file <- "biomass.tsv"

location_levels <- c("He", "Ha", "BJ", "WH", "GZ")
latitude_levels <- c("50.26° N", "45.70° N", "40.17° N", "30.58° N", "23.16° N")
cultivar_levels <- c("H43", "H84", "JD17", "TL1", "HC6")
compartment_levels <- c("B", "R", "Ro", "L")
lat_map <- c(He = "50.26° N", Ha = "45.70° N", BJ = "40.17° N", WH = "30.58° N", GZ = "23.16° N")

load_amplicon_data <- function(marker) {
  f <- marker_files[[marker]]
  test1 <- qiime2meco(f$abundance, sample_table = metadata_file, taxonomy_table = f$taxonomy,
                      phylo_tree = f$tree, rep_fasta = f$repseq, auto_tidy = TRUE)
  test1$sample_table$compartment <- factor(test1$sample_table$compartment, levels = compartment_levels)
  test1$sample_table$location <- factor(test1$sample_table$location, levels = location_levels)
  test1$sample_table$variety <- factor(test1$sample_table$variety, levels = cultivar_levels)
  test1$sample_table$lat_label <- factor(test1$sample_table$lat_label, levels = latitude_levels)
  test1$sample_table$`loc-com` <- paste(test1$sample_table$location, test1$sample_table$compartment, sep = "-")
  return(test1)
}

env_data <- read.delim(biomass_file, header = TRUE, sep = "\t", check.names = FALSE)
env_data <- env_data %>% column_to_rownames(var = "sample-id")


for (marker in MARKERS) {

  cat("\n====================", marker, "====================\n")
  test1 <- load_amplicon_data(marker)
  domain <- ifelse(marker == "16S", "Bacteria", "Fungi")

  outdir <- file.path("results", "downstream", marker, "02_core_preferential_biomass")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


  # ==========================================================
  # 1. Biomass-associated genera and context dependence
  # ==========================================================

  all_pos_sig_res <- data.frame()

  for (com in c("R", "Ro", "L")) {
    for (loc in location_levels) {

      sub_group <- clone(test1)
      sub_group$sample_table <- subset(sub_group$sample_table, compartment == com & location == loc)
      sub_group$tidy_dataset()

      m1 <- sub_group$merge_taxa(taxa = "Genus")
      m1$tax_table <- m1$tax_table[m1$tax_table$Genus != "g__", ]
      m1$tax_table <- m1$tax_table %>% filter(!grepl("g__$", Genus) & !grepl("unidentified", Genus))
      m1$tidy_dataset()
      m1$cal_abund()
      rownames(m1$taxa_abund$Genus) <- gsub("^.*\\|g__(.*)", "\\1", rownames(m1$taxa_abund$Genus))

      t1 <- trans_env$new(dataset = m1, add_data = env_data[, c("Aboveground_biomass", "Belowground_biomass")])
      t1$cal_cor(use_data = "Genus", p_adjust_method = "none")

      pos_sig_res <- t1$res_cor %>% filter(Significance != "", Correlation > 0)

      t2 <- trans_abund$new(dataset = m1, taxrank = "Genus", ntaxa = 10)
      result <- t2$data_abund %>%
        group_by(Taxonomy) %>%
        summarise(mean_abund = mean(all_mean_abund, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_abund)) %>%
        mutate(Rank = rank(-mean_abund, ties.method = "average"), N = n(),
               Abundance_percentile = 100 * (N - Rank) / (N - 1)) %>%
        select(-N)

      pos_sig_res_genus_abund <- result %>% filter(Taxonomy %in% pos_sig_res$Taxa)

      if (nrow(pos_sig_res) > 0) {
        pos_sig_res$Taxa <- sub("^.*\\|g__(.*)", "\\1", pos_sig_res$Taxa)
        pos_sig_res$Compartment <- com
        pos_sig_res$Location <- loc
        pos_sig_res$Latitude <- lat_map[loc]
        pos_sig_res <- pos_sig_res %>%
          left_join(pos_sig_res_genus_abund, by = c("Taxa" = "Taxonomy")) %>%
          arrange(desc(Correlation))
        all_pos_sig_res <- rbind(all_pos_sig_res, pos_sig_res)
      }
    }
  }

  write.table(all_pos_sig_res, file.path(outdir, paste0(marker, "_biomass_positive_correlated_genera.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  # UpSet analysis across the 15 location × compartment combinations
  all_pos_sig_res$Combo <- paste(all_pos_sig_res$Compartment, all_pos_sig_res$Latitude, sep = ": ")

  full_combos <- c(
    "L: 50.26° N", "Ro: 50.26° N", "R: 50.26° N",
    "L: 45.70° N", "Ro: 45.70° N", "R: 45.70° N",
    "L: 40.17° N", "Ro: 40.17° N", "R: 40.17° N",
    "L: 30.58° N", "Ro: 30.58° N", "R: 30.58° N",
    "L: 23.16° N", "Ro: 23.16° N", "R: 23.16° N"
  )

  df_agb <- all_pos_sig_res[all_pos_sig_res$Env == "Aboveground_biomass", c("Taxa", "Combo")]
  df_bgb <- all_pos_sig_res[all_pos_sig_res$Env == "Belowground_biomass", c("Taxa", "Combo")]

  make_matrix <- function(df_env) {
    taxa <- sort(unique(df_env$Taxa))
    mat <- matrix(0L, nrow = length(taxa), ncol = length(full_combos), dimnames = list(taxa, full_combos))
    if (nrow(df_env) > 0) {
      for (i in seq_len(nrow(df_env))) {
        taxon <- df_env$Taxa[i]
        combo <- df_env$Combo[i]
        if (combo %in% full_combos) mat[taxon, combo] <- 1L
      }
    }
    as.data.frame(mat)
  }

  mat_agb <- make_matrix(df_agb)
  mat_bgb <- make_matrix(df_bgb)

  col_L <- "#4C78A8"; col_Ro <- "#54A24B"; col_R <- "#B279A2"

  for (trait in c("AGB", "BGB")) {

    mat_upset <- if (trait == "AGB") mat_agb else mat_bgb
    trait_title <- if (trait == "AGB") "Aboveground biomass" else "Belowground biomass"

    sets_use <- names(sort(colSums(mat_upset), decreasing = TRUE))
    sets_cols <- ifelse(grepl("^L: ", sets_use), col_L,
                        ifelse(grepl("^Ro: ", sets_use), col_Ro, col_R))

    pdf(file.path(outdir, paste0(marker, "_biomass_", trait, "_upset.pdf")), width = 9.33, height = 6)
    upset(mat_upset, sets = sets_use, nsets = length(sets_use), keep.order = TRUE,
          nintersects = 60, order.by = "freq", decreasing = TRUE, mb.ratio = c(0.56, 0.44),
          main.bar.color = "#E45756", sets.bar.color = sets_cols, matrix.color = "grey25",
          point.size = 2.5, line.size = 1.2,
          text.scale = c(2.3, 1.8, 2, 1.8, 1.6, 1.6),
          mainbar.y.label = "Intersection Size", sets.x.label = "Set Size")
    grid.text(trait_title, x = 0.67, y = 0.91, gp = gpar(fontsize = 24, fontface = "bold"))
    dev.off()
  }


  # ==========================================================
  # 2. Correlation strength and abundance-rank percentile
  # ==========================================================

  all_pos_sig_res$Env1 <- gsub("_", " ", all_pos_sig_res$Env)

  p1 <- ggplot(all_pos_sig_res, aes(x = Correlation, y = Abundance_percentile, shape = Env1, color = Env1)) +
    geom_point(size = 2, alpha = 0.9) +
    labs(x = "Correlation coefficient", y = "Abundance percentile (%)") +
    scale_shape_manual(values = c("Aboveground biomass" = 15, "Belowground biomass" = 17)) +
    scale_color_manual(values = c("Aboveground biomass" = "#fccccb", "Belowground biomass" = "#8CCDBF")) +
    theme_bw(base_size = 16) +
    theme(panel.border = element_blank(), axis.line = element_line(color = "black"),
          panel.grid = element_blank())

  p2 <- ggplot(all_pos_sig_res, aes(x = Env1, y = Correlation, fill = Env1)) +
    geom_boxplot(width = 0.8) +
    stat_compare_means(method = "wilcox", label = "p.signif", label.x = 1.5, size = 5) +
    scale_fill_manual(values = c("Aboveground biomass" = "#fccccb", "Belowground biomass" = "#8CCDBF")) +
    coord_flip() + theme_bw(base_size = 16) +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid = element_blank())

  p3 <- ggplot(all_pos_sig_res, aes(x = Env1, y = Abundance_percentile, fill = Env1)) +
    geom_boxplot(width = 0.8) +
    stat_compare_means(method = "wilcox", label = "p.signif", label.x = 1.28, size = 5) +
    scale_fill_manual(values = c("Aboveground biomass" = "#fccccb", "Belowground biomass" = "#8CCDBF")) +
    theme_bw(base_size = 16) +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), panel.grid = element_blank())

  g <- p1 %>% insert_top(p2, height = 0.2) %>% insert_right(p3, width = 0.2)

  ggsave(file.path(outdir, paste0(marker, "_biomass_correlation_vs_abundance_percentile.pdf")), g, width = 9.6, height = 6.5)

  # ==========================================================
  # 3. Geographic distribution of genera
  # ==========================================================

  unique_loccom <- unique(test1$sample_table$`loc-com`)
  details_list <- list()
  abundance_list <- list()

  for (loccom in unique_loccom) {

    sub_group <- clone(test1)
    sub_group$sample_table <- subset(sub_group$sample_table, `loc-com` == loccom)
    sub_group$tidy_dataset()

    m1 <- sub_group$merge_taxa(taxa = "Genus")
    m1$tax_table <- m1$tax_table[m1$tax_table$Genus != "g__", ]
    m1$tidy_dataset()

    m2 <- m1$merge_samples("variety")
    m3 <- trans_venn$new(dataset = m2, ratio = NULL)
    m3$data_details <- m3$data_details[m3$data_details$`H43&H84&HC6&JD17&TL1` != "", ,]
    details_list[[loccom]] <- m3$data_details$`H43&H84&HC6&JD17&TL1`

    all_variety_shared_genus <- m2$otu_table[
      m2$otu_table$H43 != 0 & m2$otu_table$H84 != 0 & m2$otu_table$HC6 != 0 & m2$otu_table$JD17 != 0 & m2$otu_table$TL1 != 0, ]
    abundance_list[[loccom]] <- all_variety_shared_genus
  }

  all_count_stacked_data <- data.frame()
  all_abund_stacked_data <- data.frame()

  for (compartment in compartment_levels) {

    loc_comp_details <- details_list[paste(location_levels, compartment, sep = "-")]

    shared_elements_all <- Reduce(intersect, loc_comp_details)

    unique_elements <- lapply(names(loc_comp_details), function(name) {
      current <- loc_comp_details[[name]]
      others <- unlist(loc_comp_details[names(loc_comp_details) != name])
      setdiff(current, others)
    })
    names(unique_elements) <- names(loc_comp_details)

    combn_indices <- unlist(lapply(2:4, function(n) {
      combn(names(loc_comp_details), n, simplify = FALSE)
    }), recursive = FALSE)

    shared_elements <- lapply(combn_indices, function(combo) {
      shared_genus <- Reduce(intersect, loc_comp_details[combo])
      remaining <- setdiff(names(loc_comp_details), combo)
      setdiff(shared_genus, unlist(loc_comp_details[remaining]))
    })
    names(shared_elements) <- sapply(combn_indices, paste, collapse = "&")

    # Genus counts
    stacked_count <- data.frame(location = character(), Type = character(), Genus_count = integer())

    stacked_count <- rbind(stacked_count, data.frame(
      location = gsub(paste0("-", compartment), "", names(loc_comp_details)),
      Type = "Shared Across All Locations",
      Genus_count = length(shared_elements_all)
    ))

    for (combo in names(shared_elements)) {
      stacked_count <- rbind(stacked_count, data.frame(
        location = gsub(paste0("-", compartment), "", combo),
        Type = "Shared Across Some Locations",
        Genus_count = length(shared_elements[[combo]])
      ))
    }

    for (name in names(unique_elements)) {
      stacked_count <- rbind(stacked_count, data.frame(
        location = gsub(paste0("-", compartment), "", name),
        Type = "Location-specific",
        Genus_count = length(unique_elements[[name]])
      ))
    }

    stacked_count <- stacked_count %>%
      separate_rows(location, sep = "&") %>%
      group_by(location, Type) %>%
      summarise(Genus_count = sum(Genus_count), .groups = "drop") %>%
      mutate(Latitude = lat_map[location], Compartment = compartment)

    all_count_stacked_data <- rbind(all_count_stacked_data, stacked_count)

    # Cumulative genus abundance
    loc_comp_abundance <- abundance_list[paste(location_levels, compartment, sep = "-")]

    shared_abundance <- lapply(loc_comp_abundance, function(df) df[shared_elements_all, , drop = FALSE])
    shared_abundance_sum <- sum(do.call(rbind, shared_abundance))
    shared_abundance_df <- data.frame(Combination = "He&Ha&BJ&WH&GZ", Abundance_Sum = shared_abundance_sum, row.names = NULL)

    unique_abundance <- lapply(names(unique_elements), function(name) {
      abundance_data <- abundance_list[[name]]
      abundance_data[rownames(abundance_data) %in% unique_elements[[name]], ]})
    names(unique_abundance) <- names(unique_elements)
    unique_abundance_df <- data.frame(Combination = names(unique_abundance), Abundance_Sum = sapply(unique_abundance, sum))

    shared_combinations <- lapply(combn_indices, function(combo) {
      unique_shared <- shared_elements[[paste(combo, collapse = "&")]]
      combo_abundance <- lapply(combo, function(loc) abundance_list[[loc]])
      combined_abundance <- do.call(rbind, lapply(combo_abundance, function(df) {
        df[unique_shared, , drop = FALSE]
      }))
      sum(combined_abundance)
    })
    shared_combinations_df <- data.frame(Combination = sapply(combn_indices, paste, collapse = "&"), Abundance_Sum = unlist(shared_combinations))

    final_abund <- rbind(unique_abundance_df, shared_combinations_df, shared_abundance_df)
    final_abund$Combination <- gsub(paste0("-", compartment), "", final_abund$Combination)

    specific_locations <- final_abund[final_abund$Combination %in% location_levels, ]
    partial_shared <- final_abund[grepl("&", final_abund$Combination) & !grepl("He&Ha&BJ&WH&GZ", final_abund$Combination), ]

    merged_data <- data.frame()
    for (loc in location_levels) {
      Abundance <- sum(partial_shared[grepl(loc, partial_shared$Combination), "Abundance_Sum"], na.rm = TRUE)
      merged_data <- rbind(merged_data, data.frame(Combination = paste(loc, "Partially Shared"), Abundance_Sum = Abundance))
    }

    all_shared <- final_abund[grepl("He&Ha&BJ&WH&GZ", final_abund$Combination), ]
    final_summary <- rbind(specific_locations, merged_data, all_shared)

    stacked_abund <- data.frame(
      location = rep(location_levels, each = 3),
      Type = rep(c("Location-specific", "Shared Across Some Locations", "Shared Across All Locations"), times = 5),
      Abundance = c(
        final_summary[1, "Abundance_Sum"], final_summary[6, "Abundance_Sum"], final_summary[11, "Abundance_Sum"],
        final_summary[2, "Abundance_Sum"], final_summary[7, "Abundance_Sum"], final_summary[11, "Abundance_Sum"],
        final_summary[3, "Abundance_Sum"], final_summary[8, "Abundance_Sum"], final_summary[11, "Abundance_Sum"],
        final_summary[4, "Abundance_Sum"], final_summary[9, "Abundance_Sum"], final_summary[11, "Abundance_Sum"],
        final_summary[5, "Abundance_Sum"], final_summary[10, "Abundance_Sum"], final_summary[11, "Abundance_Sum"]
      )
    ) %>% mutate(Latitude = lat_map[location], Compartment = compartment)

    all_abund_stacked_data <- rbind(all_abund_stacked_data, stacked_abund)
  }

  all_count_stacked_data$Type <- factor(all_count_stacked_data$Type, levels = c("Location-specific", "Shared Across Some Locations", "Shared Across All Locations"))
  all_count_stacked_data$Latitude <- factor(all_count_stacked_data$Latitude, levels = latitude_levels)
  all_count_stacked_data$Compartment <- factor(all_count_stacked_data$Compartment, levels = compartment_levels)

  all_abund_stacked_data$Type <- factor(all_abund_stacked_data$Type, levels = c("Location-specific", "Shared Across Some Locations", "Shared Across All Locations"))
  all_abund_stacked_data$Latitude <- factor(all_abund_stacked_data$Latitude, levels = latitude_levels)
  all_abund_stacked_data$Compartment <- factor(all_abund_stacked_data$Compartment, levels = compartment_levels)

  distribution_colors <- c("Location-specific" = "lightcoral",
                           "Shared Across Some Locations" = "lightgreen",
                           "Shared Across All Locations" = "lightblue")

  p_count <- ggplot(all_count_stacked_data, aes(x = Latitude, y = Genus_count, fill = Type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Compartment, nrow = 1) +
    labs(x = NULL, y = "Number of genera", fill = "Type") +
    scale_fill_manual(values = distribution_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))

  p_abund <- ggplot(all_abund_stacked_data, aes(x = Latitude, y = Abundance, fill = Type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Compartment, nrow = 1) +
    labs(x = NULL, y = "Abundance of genera", fill = "Type") +
    scale_fill_manual(values = distribution_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))

  g <- ggarrange(p_count, p_abund, nrow = 2, common.legend = TRUE, legend = "top")
  ggsave(file.path(outdir, paste0(marker, "_genus_geographic_distribution.pdf")),
         g, width = 15, height = 8)

  # Mean abundance per genus for the three distribution categories
  final_df <- left_join(all_abund_stacked_data, all_count_stacked_data, by = c("location", "Type", "Latitude", "Compartment")) %>%
    mutate(Mean_abundance_per_genus = Abundance / Genus_count) %>%
    arrange(desc(Mean_abundance_per_genus))

  plot_df <- final_df %>%
    mutate(Type = factor(Type, levels = c("Location-specific", "Shared Across Some Locations", "Shared Across All Locations")),
           Compartment = factor(Compartment, levels = compartment_levels))

  my_comparisons <- list(
    c("Location-specific", "Shared Across Some Locations"),
    c("Location-specific", "Shared Across All Locations"),
    c("Shared Across Some Locations", "Shared Across All Locations")
  )

  p <- ggplot(plot_df, aes(x = Type, y = Mean_abundance_per_genus, fill = Type)) +
    geom_boxplot() +
    facet_wrap(~ Compartment, nrow = 1, scales = "free_y") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
    labs(x = NULL, y = "Mean abundance per genus", fill = "Type") +
    scale_fill_manual(values = distribution_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(outdir, paste0(marker, "_mean_abundance_by_distribution_category.pdf")), p, width = 12, height = 8)

  # ==========================================================
  # 4. Core genera and their abundance profiles
  # ==========================================================

  test1$cal_abund()
  genus_table <- test1$taxa_abund$Genus

  mat <- as.matrix(genus_table)
  tax <- rownames(mat)
  genus <- ifelse(grepl("\\|g__[^|]+$", tax), sub(".*\\|g__", "", tax), "Unidentified")
  genus_table_merged <- rowsum(mat, group = genus, reorder = TRUE)

  type <- sub("-\\d+$", "", colnames(genus_table_merged))
  sum_by_type <- rowsum(t(genus_table_merged), group = type, reorder = TRUE)
  counts <- table(type)
  counts <- as.numeric(counts[rownames(sum_by_type)])
  genus_mean_by_type <- t(sweep(sum_by_type, 1, counts, "/"))

  core_flag <- rowSums(genus_mean_by_type > 0) == ncol(genus_mean_by_type)
  core_genera <- rownames(genus_mean_by_type)[core_flag]
  core_df <- data.frame(core_genus = core_genera)

  write.table(core_df, file.path(outdir, paste0(marker, "_core_genera.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  # Relative abundance of core genera
  t1 <- trans_abund$new(dataset = test1, taxrank = "Genus", ntaxa = nrow(test1$taxa_abund$Genus), groupmean = "group", group_morestats = TRUE)
  t1$data_abund$Taxonomy <- sub("^.*\\|g__", "", t1$data_abund$Taxonomy)
  t1$data_abund <- t1$data_abund[t1$data_abund$Taxonomy %in% core_genera, ]

  split_names <- strsplit(as.character(t1$data_abund$Sample), "-")
  t1$data_abund$location <- sapply(split_names, function(x) x[1])
  t1$data_abund$variety <- sapply(split_names, function(x) x[2])
  t1$data_abund$compartment <- sapply(split_names, function(x) x[3])
  t1$data_abund$lat_label <- lat_map[t1$data_abund$location]

  t1$data_abund$compartment <- factor(t1$data_abund$compartment, levels = compartment_levels)
  t1$data_abund$location <- factor(t1$data_abund$location, levels = location_levels)
  t1$data_abund$variety <- factor(t1$data_abund$variety, levels = cultivar_levels)
  t1$data_abund$lat_label <- factor(t1$data_abund$lat_label, levels = latitude_levels)

  sample_order <- c()
  for (loc in location_levels) {
    for (var in cultivar_levels) sample_order <- c(sample_order, paste(loc, var, compartment_levels, sep = "-"))
  }

  core_colors <- c("#8cd0c3", "#3d9f3c", "#7087e4", "#cb8d8b", "#b2446b", "#0c54a5",
                   "#faf5b5", "#bcb9d8", "#f18072", "#e31a1c", "#80b1d2",
                   "#d7d7d5", "#fcc3b4", "#f49513", "#ba7fb5", "#b3d46b")

  p <- t1$plot_bar(facet = c("compartment", "lat_label"), xtext_angle = 45, xtext_size = 7, barwidth = 1, order_x = sample_order) +
    scale_fill_manual(values = core_colors) +
    theme_bw() + theme(panel.grid = element_blank())

  ggsave(file.path(outdir, paste0(marker, "_core_genera_abundance.pdf")), p, width = 18, height = 5.5, limitsize = FALSE)


  # ==========================================================
  # 5. Core vs non-core abundance rank
  # ==========================================================

  genus_mean_rank <- genus_mean_by_type[rownames(genus_mean_by_type) != "Unidentified", , drop = FALSE]
  core_in_mat <- intersect(setdiff(core_genera, "Unidentified"), rownames(genus_mean_rank))

  rank_mat <- apply(genus_mean_rank, 2, function(x) {
    x[x == 0] <- NA
    r <- rank(-x, ties.method = "average", na.last = "keep")
    N <- sum(!is.na(r))
    100 * (N - r) / (N - 1)
  })

  plot_rank_df <- data.frame(genus = rep(rownames(rank_mat), times = ncol(rank_mat)), type = rep(colnames(rank_mat), each = nrow(rank_mat)), rank = as.vector(rank_mat) )
  plot_rank_df$core_flag <- ifelse(plot_rank_df$genus %in% core_in_mat, "Core Genus", "Non-core Genus")
  plot_rank_df <- plot_rank_df %>% filter(!is.na(rank))

  wilcox_result <- wilcox.test(rank ~ core_flag, data = plot_rank_df)
  p_label <- paste0("Wilcoxon P = ", formatC(wilcox_result$p.value, format = "e", digits = 2))

  p <- ggplot(plot_rank_df, aes(x = core_flag, y = rank, fill = core_flag)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black") +
    scale_fill_manual(values = c("Core Genus" = "#97cd99", "Non-core Genus" = "#fe7a7e")) +
    labs(x = NULL, y = "Relative abundance percentile (%)") +
    annotate("text", x = 1.5, y = 105, label = p_label, size = 5) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())

  ggsave(file.path(outdir, paste0(marker, "_core_noncore_abundance_rank.pdf")), p, width = 4, height = 4.5)

  # ==========================================================
  # 6. Location-preferential genera and overlap
  # ==========================================================

  home_map <- data.frame(location = c("He", "Ha", "BJ", "WH", "GZ"), variety = c("H43", "H84", "JD17", "TL1", "HC6"), stringsAsFactors = FALSE)

  col_parts <- do.call(rbind, strsplit(colnames(genus_mean_by_type), "-", fixed = TRUE))
  col_parts <- as.data.frame(col_parts, stringsAsFactors = FALSE)
  colnames(col_parts) <- c("location", "variety", "compartment")
  col_parts$colname <- colnames(genus_mean_by_type)

  res_list <- list()
  k <- 1

  for (i in seq_len(nrow(home_map))) {

    loc <- home_map$location[i]
    var <- home_map$variety[i]

    for (comp in compartment_levels) {

      home_col <- paste(loc, var, comp, sep = "-")
      home_mean <- genus_mean_by_type[, home_col]
      other_cols <- col_parts$colname[col_parts$variety == var & col_parts$compartment == comp & col_parts$location != loc]

      other_mat <- genus_mean_by_type[, other_cols, drop = FALSE]
      max_other <- apply(other_mat, 1, max, na.rm = TRUE)
      keep <- (home_mean > max_other) & (home_mean > 0)

      if (any(keep)) {
        tmp <- data.frame(
          genus = rownames(genus_mean_by_type)[keep],
          home_loc = loc,
          variety = var,
          compartment = comp,
          home_mean = home_mean[keep],
          stringsAsFactors = FALSE
        )
        res_list[[k]] <- tmp
        k <- k + 1
      }
    }
  }

  preferred_df <- if (length(res_list) == 0) data.frame() else do.call(rbind, res_list)
  preferred_df <- preferred_df[order(preferred_df$home_loc, preferred_df$compartment, -preferred_df$home_mean), ]
  preferred_df <- preferred_df[preferred_df$genus != "Unidentified", ]

  write.table(preferred_df, file.path(outdir, paste0(marker, "_location_preferential_genera.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)


  # UpSet plot of location-preferential genera
  df <- preferred_df[preferred_df$compartment != "B", ]
  df$Comb <- paste(df$variety, df$compartment, sep = "_")

  full_preferred_combos <- c(
    "H43_L", "H43_Ro", "H43_R",
    "H84_L", "H84_Ro", "H84_R",
    "JD17_L", "JD17_Ro", "JD17_R",
    "TL1_L", "TL1_Ro", "TL1_R",
    "HC6_L", "HC6_Ro", "HC6_R"
  )

  taxa <- unique(df$genus)
  mat_pref <- matrix(0, nrow = length(taxa), ncol = length(full_preferred_combos), dimnames = list(taxa, full_preferred_combos))

  for (i in seq_len(nrow(df))) mat_pref[df$genus[i], df$Comb[i]] <- 1
  mat_pref <- as.data.frame(mat_pref)

  sets_order <- c(
    "H43_L", "H84_L", "JD17_L", "TL1_L", "HC6_L",
    "H43_Ro", "H84_Ro", "JD17_Ro", "TL1_Ro", "HC6_Ro",
    "H43_R", "H84_R", "JD17_R", "TL1_R", "HC6_R"
  )

  sets_cols <- ifelse(grepl("_L$", sets_order), col_L, ifelse(grepl("_Ro$", sets_order), col_Ro, col_R))

  pdf(file.path(outdir, paste0(marker, "_location_preferential_upset.pdf")), width = 10, height = 6.67)
  upset(mat_pref, sets = sets_order, nsets = length(sets_order), keep.order = TRUE,
        nintersects = 60, order.by = "freq", decreasing = TRUE, mb.ratio = c(0.56, 0.44),
        main.bar.color = "#E45756", sets.bar.color = sets_cols, matrix.color = "grey25",
        point.size = 2.5, line.size = 1.2,
        text.scale = c(2.3, 1.8, 2, 1.8, 1.6, 1.25),
        mainbar.y.label = "Intersection Size", sets.x.label = "Set Size")
  grid.text(domain, x = 0.65, y = 0.92, gp = gpar(fontsize = 22, fontface = "bold"))
  dev.off()

  # ==========================================================
  # 7. LEfSe enrichment at home sites
  # ==========================================================

  final_lefse_loc_pref_df <- data.frame()

  for (com in compartment_levels) {

    sub_group <- clone(test1)
    sub_group$sample_table <- subset(sub_group$sample_table, compartment == com)
    sub_group$tidy_dataset()

    for (main_cultivars in c("He-H43", "Ha-H84", "BJ-JD17", "WH-TL1", "GZ-HC6")) {

      m1 <- clone(sub_group)
      split_names <- strsplit(main_cultivars, "-")
      loc <- sapply(split_names, function(x) x[1])
      var <- sapply(split_names, function(x) x[2])
      lat <- lat_map[loc]

      m1$sample_table <- subset(m1$sample_table, variety == var)
      m1$tidy_dataset()
      m1$cal_abund()

      t1 <- trans_diff$new(dataset = m1, method = "lefse", taxa_level = "Genus", group = "lat_label")
      t1$res_diff$Taxa <- gsub("^.*\\|g__(.*)", "\\1", t1$res_diff$Taxa)
      rownames(t1$res_diff) <- gsub("^.*\\|g__(.*)", "\\1", rownames(t1$res_diff))

      preferred_microbes <- preferred_df %>% filter(home_loc == loc, variety == var, compartment == com)

      lefse_loc_pref <- t1$res_diff %>%
        filter(Group == lat, Taxa %in% preferred_microbes$genus) %>%
        mutate(Compartment = com, Cultivar = var, Domain = domain, `Home site` = loc) %>%
        select(Domain, Compartment, Cultivar, `Home site`, Taxa, Method, LDA, P.unadj, P.adj, Significance)

      final_lefse_loc_pref_df <- rbind(final_lefse_loc_pref_df, lefse_loc_pref)
    }
  }
  
  write.table(final_lefse_loc_pref_df, file.path(outdir, paste0(marker, "_LEfSe_home_site_enriched_genera.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
}
