# ============================================================
# 01_amplicon_community_analysis.R
#
# 16S / ITS community analyses:
# 1. Alpha diversity
# 2. Linear mixed-effects model
# 3. Beta diversity
# 4. Cultivar effects within location × compartment
# 5. Taxonomic composition
# 6. Distance-decay
# ============================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(file2meco)
library(microeco)
library(ggpubr)
library(lme4)
library(lmerTest)
library(car)
library(broom)
library(reshape2)

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
  return(test1)
}

for (marker in MARKERS) {

  cat("\n====================", marker, "====================\n")
  test1 <- load_amplicon_data(marker)
  outdir <- file.path("results", "downstream", marker, "01_community")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # ==========================================================
  # 1. Alpha diversity
  # ==========================================================

  test1$cal_alphadiv(measures = "Shannon")

  alpha_test <- list()
  for (group_var in c("variety", "lat_label", "compartment")) {
    t1 <- trans_alpha$new(dataset = test1, group = group_var)

    t1$cal_diff(method = "KW")
    total_p <- t1$res_diff$P.adj
    alpha_test[[paste0(group_var, "_KW")]] <- t1$res_diff

    t1$cal_diff(method = "KW_dunn")
    alpha_test[[paste0(group_var, "_Dunn")]] <- t1$res_diff

    p <- t1$plot_alpha(measure = "Shannon", color_values = c("#fccccb", "#bdb5e1", "#b0d992", "#f9d580", "#8ccDBF"), xtext_angle = 0, add_sig_text_size = 6.5) +
      ggtitle(paste0(group_var, ": p = ", formatC(total_p, format = "e", digits = 2))) + theme_bw() +
      theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
      
    ggsave(file.path(outdir, paste0(marker, "_Shannon_", group_var, ".pdf")), p, width = 5, height = 4)
  }
  alpha_test

  # ==========================================================
  # 2. Linear mixed-effects model for Shannon diversity
  # ==========================================================

  data <- test1$alpha_diversity
  split_data <- do.call(rbind, strsplit(as.character(rownames(data)), "-"))
  data$location <- split_data[, 1]
  data$variety <- split_data[, 2]
  data$compartment <- split_data[, 3]
  data$group <- paste(data$location, data$variety, data$compartment, sep = "-")

  data$location <- factor(data$location, levels = location_levels)
  data$compartment <- factor(data$compartment, levels = compartment_levels)
  data$variety <- factor(data$variety, levels = cultivar_levels)  

  control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  model <- lmerTest::lmer(Shannon ~ variety * location * compartment + (1 | group), data = data, REML = TRUE, control = control)
  anova_result <- Anova(model, test.statistic = "F", multivariate = TRUE, p.adjust.methods = "Holm")
  anova_result

  # ==========================================================
  # 3. Beta diversity: Bray-Curtis / PCoA / PERMANOVA / ANOSIM
  # ==========================================================

  test1$cal_betadiv(unifrac = TRUE)
  t1 <- trans_beta$new(dataset = test1, measure = "bray")
  t1$cal_ordination(method = "PCoA")

  p <- t1$plot_ordination(plot_color = "lat_label", plot_shape = "variety",
                          plot_type = "point", point_alpha = 1, point_size = 3) +
    geom_point(aes(fill = compartment), size = 3) +
    stat_ellipse(aes(group = compartment), size = 1, linetype = 2, show.legend = FALSE) +
    theme_bw() + theme(panel.grid = element_blank())

  ggsave(file.path(outdir, paste0(marker, "_PCoA.pdf")), p, width = 8, height = 6)

  beta_test <- list()
  for (group_var in c("location", "variety", "compartment")) {
    t1$cal_manova(manova_all = TRUE, group = group_var)
    beta_test[[paste0(group_var, "_PERMANOVA")]] <- t1$res_manova

    t1$cal_anosim(paired = FALSE, group = group_var)
    beta_test[[paste0(group_var, "_ANOSIM")]] <- t1$res_anosim
  }
  beta_test

  # ==========================================================
  # 4. Cultivar effects within location × compartment
  # ==========================================================

  test1$sample_table$Combination <- paste(test1$sample_table$compartment, test1$sample_table$lat_label, sep = ": ")
  unique_loccom <- unique(test1$sample_table$Combination)

  manova_all_result <- data.frame()
  manova_pair_result <- data.frame()
  plot_list <- list()

  for (loccom in unique_loccom) {
    sub_group <- clone(test1)
    sub_group$sample_table <- subset(sub_group$sample_table, Combination == loccom)
    sub_group$tidy_dataset()
    sub_group$cal_betadiv(unifrac = TRUE)

    t1 <- trans_beta$new(dataset = sub_group, group = "variety", measure = "bray")

    t1$cal_manova(manova_all = TRUE)
    manova_all_result <- rbind(manova_all_result, cbind(Combination = loccom, t1$res_manova[1, ]))

    t1$cal_manova(manova_all = FALSE, p_adjust_method = "none")
    manova_pair_result <- rbind(manova_pair_result, cbind(Combination = loccom, t1$res_manova))
  }

  if ("Significance" %in% colnames(manova_all_result)) manova_all_result <- select(manova_all_result, -Significance)
  manova_all_result$P.adj <- p.adjust(manova_all_result$`Pr(>F)`, method = "fdr")

  if ("p.adjusted" %in% colnames(manova_pair_result)) manova_pair_result <- select(manova_pair_result, -p.adjusted)
  if ("Significance" %in% colnames(manova_pair_result)) manova_pair_result <- select(manova_pair_result, -Significance)
  manova_pair_result$P.adj <- p.adjust(manova_pair_result$p.value, method = "fdr")

  write.table(manova_all_result, file.path(outdir, paste0(marker, "_cultivar_effect_location_compartment.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(manova_pair_result, file.path(outdir, paste0(marker, "_pair_cultivar_effect_location_compartment.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  
  for (loccom in unique_loccom) {
    sub_group <- clone(test1)
    sub_group$sample_table <- subset(sub_group$sample_table, Combination == loccom)
    colnames(sub_group$sample_table)[colnames(sub_group$sample_table) == "variety"] <- "Cultivar"
    sub_group$tidy_dataset()
    sub_group$cal_betadiv(unifrac = TRUE)
    t1 <- trans_beta$new(dataset = sub_group, group = "Cultivar", measure = "bray")
    t1$cal_ordination(method  = "PCoA")
    loccom_row <- manova_all_result[manova_all_result$Combination  == loccom, ]
    R2_value<-round(loccom_row$R2,digits = 3)
    p_adj<-round(loccom_row$P.adj,digits = 3)
    p1 <- t1$plot_ordination(plot_color = "Cultivar",
                            plot_type = c("point", "ellipse"),point_alpha = 1,point_size = 2,
                            color_values = c("#fccccb","#bdb5e1","#b0d992","#f9d580","#8ccDBF")) +
      theme_bw() + theme(panel.grid = element_blank()) +
      ggtitle(paste0(loccom," R2=",R2_value, ", p=",p_adj)) +
      theme(legend.text = element_text(size = 18),
            legend.title = element_text(size = 20),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 13),
            plot.title = element_text(size = 14,face = "bold"))
    plot_list[[loccom]] <- p1
    }

  new_order <- c("B: 50.26° N", "R: 50.26° N", "Ro: 50.26° N", "L: 50.26° N",
                "B: 45.70° N", "R: 45.70° N", "Ro: 45.70° N", "L: 45.70° N",
                "B: 40.17° N", "R: 40.17° N", "Ro: 40.17° N", "L: 40.17° N",
                "B: 30.58° N", "R: 30.58° N", "Ro: 30.58° N", "L: 30.58° N",
                "B: 23.16° N", "R: 23.16° N", "Ro: 23.16° N", "L: 23.16° N" )
  sorted_plot_list <- plot_list[new_order]
  g<-ggarrange(plotlist = sorted_plot_list, common.legend = TRUE, legend="top",ncol = 4,nrow = 5)
  ggsave(file.path(outdir, paste0(marker, "_sub_cultivar_beta_div", ".pdf")), g, width = 14.5, height = 16)
  # ==========================================================
  # 5. Taxonomic composition
  # ==========================================================

  test1$cal_abund()

  # Phylum
  t1 <- trans_abund$new(dataset = test1, taxrank = "Phylum", ntaxa = 10, groupmean = "group", group_morestats = TRUE)
  split_names <- strsplit(as.character(t1$data_abund$Sample), "-")
  t1$data_abund$location <- sapply(split_names, function(x) x[1])
  t1$data_abund$variety <- sapply(split_names, function(x) x[2])
  t1$data_abund$compartment <- sapply(split_names, function(x) x[3])
  t1$data_abund$lat_label <- c(He = "50.26° N", Ha = "45.70° N", BJ = "40.17° N",
                               WH = "30.58° N", GZ = "23.16° N")[t1$data_abund$location]

  p <- t1$plot_bar(facet = c("compartment", "lat_label"), xtext_angle = 45, barwidth = 1) + theme_bw()
  ggsave(file.path(outdir, paste0(marker, "_phylum_composition.pdf")), p, width = 20, height = 6, limitsize = FALSE)

  # Genus
  t1 <- trans_abund$new(dataset = test1, taxrank = "Genus", ntaxa = 15, groupmean = "group", group_morestats = TRUE)
  split_names <- strsplit(as.character(t1$data_abund$Sample), "-")
  t1$data_abund$location <- sapply(split_names, function(x) x[1])
  t1$data_abund$variety <- sapply(split_names, function(x) x[2])
  t1$data_abund$compartment <- sapply(split_names, function(x) x[3])
  t1$data_abund$lat_label <- c(He = "50.26° N", Ha = "45.70° N", BJ = "40.17° N",
                               WH = "30.58° N", GZ = "23.16° N")[t1$data_abund$location]

  p <- t1$plot_bar(facet = c("compartment", "lat_label"), xtext_angle = 45, barwidth = 1) + theme_bw()
  ggsave(file.path(outdir, paste0(marker, "_genus_composition.pdf")), p, width = 20, height = 6, limitsize = FALSE)

  # ==========================================================
  # 6. Bray-Curtis similarity vs latitude difference
  # ==========================================================

loc_lat_mapping <- data.frame(location = c("He", "Ha", "BJ", "WH", "GZ"), latitude = c("50.2551", "45.6993", "40.1741", "30.5816", "23.1645"))
col <- c("#FB9A99", "#bdb5e1", "#b0d992", "#f9d580", "#8ccDBF")

for (com in c("B", "R", "Ro", "L")) {

  # ----------------------------------------------------------
  # Samples from the five cultivars within one compartment
  # ----------------------------------------------------------

  unique_varcom <- paste(c("H43", "H84", "JD17", "TL1", "HC6"), com, sep = "-")
  df_all <- data.frame()


  # Calculate Bray-Curtis similarity
  for (varcom in unique_varcom) {
    sub_group <- clone(test1)
    sub_group$sample_table <- subset(sub_group$sample_table, `var-com` == varcom)
    sub_group$tidy_dataset()
    sub_group$cal_betadiv(unifrac = TRUE)

    bray_sim <- 1 - as.matrix(sub_group$beta_diversity$bray)
    diag(bray_sim) <- 0
    bray_sim[upper.tri(bray_sim)] <- 0
    bray_sim <- reshape2::melt(bray_sim)
    bray_sim <- subset(bray_sim, value != 0)
    bray_sim$loc1 <- substr(bray_sim$Var1, 1, 2)
    bray_sim$loc2 <- substr(bray_sim$Var2, 1, 2)

    df <- bray_sim %>%
      filter(loc1 != loc2) %>%
      left_join(loc_lat_mapping, by = c("loc1" = "location")) %>%
      left_join(loc_lat_mapping, by = c("loc2" = "location"), suffix = c("_1", "_2")) %>%
      mutate(lat_diff = abs(as.numeric(latitude_1) - as.numeric(latitude_2)), variety = strsplit(varcom, "-")[[1]][1]) %>%
      select(variety, Var1, Var2, value, lat_diff)

    df_all <- rbind(df_all, df)
  }


  # ----------------------------------------------------------
  # Export Bray-Curtis similarity and latitude differences
  # ----------------------------------------------------------

  final_df_all <- df_all
  final_df_all$Compartment <- com
  final_df_all <- final_df_all %>% rename(Cultivar = variety, Value = value, delta_Latitude = lat_diff)

  final_df_all$Location1 <- sub("-.*", "", final_df_all$Var1)
  final_df_all$Location2 <- sub("-.*", "", final_df_all$Var2)
  final_df_all$Latitude1 <- loc_lat_mapping$latitude[match(final_df_all$Location1, loc_lat_mapping$location)]
  final_df_all$Latitude2 <- loc_lat_mapping$latitude[match(final_df_all$Location2, loc_lat_mapping$location)]
  final_df_all <- final_df_all[, c("Compartment", "Cultivar", "Var1", "Var2", "Value", "Latitude1", "Latitude2", "delta_Latitude")]

  write.table(final_df_all, file.path(outdir, paste0(marker, "_", com, "_distance_decay_pairs.tsv")), sep = "\t",row.names = FALSE, quote = FALSE)

  # ----------------------------------------------------------
  # Plot Bray-Curtis similarity against latitude difference
  # ----------------------------------------------------------

  df_all$variety <- factor(df_all$variety, levels = cultivar_levels)
  colnames(df_all)[colnames(df_all) == "variety"] <- "Variety"

  p <- ggplot(df_all, aes(x = lat_diff, y = value, color = Variety)) +
    geom_point(size = 2, shape = 1) +
    stat_smooth(method = lm, aes(color = Variety, fill = Variety)) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    labs(x = "ΔLatitude", y = "Bray-Curtis similarity", color = "Cultivar", fill = "Cultivar") +
    scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 30, 35)) +
    theme_classic() + ggtitle(paste0("Com: ", com))

  ggsave(file.path(outdir, paste0(marker, "_", com, "_bray_latitude.pdf")), plot = p, width = 5, height = 4)


  # ----------------------------------------------------------
  # Linear regression equations for each cultivar
  # ----------------------------------------------------------

  fit_list <- df_all %>% group_by(Variety) %>% do(broom::tidy(lm(value ~ lat_diff, data = .)))

  equations <- fit_list %>%
    group_by(Variety) %>%
    summarize(
      equation = paste0(
        "y=",
        round(estimate[term == "(Intercept)"], 6),
        ifelse(estimate[term == "lat_diff"] >= 0, "+", "-"),
        abs(round(estimate[term == "lat_diff"], 10)),
        "x"
      )
    )

  write.table(equations, file.path(outdir, paste0(marker, "_", com, "_equations_latitude.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
