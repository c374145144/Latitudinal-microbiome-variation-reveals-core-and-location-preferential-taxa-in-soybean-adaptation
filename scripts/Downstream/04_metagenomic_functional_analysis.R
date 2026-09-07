# ============================================================
# 04_metagenomic_functional_analysis.R
#
# Shotgun metagenomic downstream analyses based on KO abundance:
# 1. Functional alpha diversity
# 2. Dominant KO composition
# 3. Functional beta diversity
# 4. Location-preferential KOs and overlap
# 5. Nitrogen, phosphorus and potassium cycling functions
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(microeco)
library(UpSetR)
library(pheatmap)

ko_abundance_file <- "Final_Full_KO_Abundance_Matrix_p1.tsv"
metadata_file <- "meta_metadata.tsv"
ko_annotation_file <- "tax_table.tsv"
npk_ko_file <- "NPK_cycle_wang.txt"

location_levels <- c("He", "Ha", "BJ", "WH", "GZ")
latitude_levels <- c("50.26° N", "45.70° N", "40.17° N", "30.58° N", "23.16° N")
cultivar_levels <- c("H43", "H84", "JD17", "TL1", "HC6")
compartment_levels <- c("B", "R")
lat_map <- c(He = "50.26° N", Ha = "45.70° N", BJ = "40.17° N", WH = "30.58° N", GZ = "23.16° N")
var2loc <- c(H43 = "He", H84 = "Ha", JD17 = "BJ", TL1 = "WH", HC6 = "GZ")

outdir <- file.path("results", "downstream", "metagenomics", "04_functional_analysis")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

load_metagenome_data <- function() {
  ko_tab <- read.table(ko_abundance_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  sample_tab <- read.table(metadata_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  tax_tab <- read.table(ko_annotation_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, quote = "")

  test1 <- microtable$new(otu_table = ko_tab, sample_table = sample_tab, tax_table = tax_tab)
  test1$sample_table$compartment <- factor(test1$sample_table$compartment, levels = compartment_levels)
  test1$sample_table$location <- factor(test1$sample_table$location, levels = location_levels)
  test1$sample_table$variety <- factor(test1$sample_table$variety, levels = cultivar_levels)
  test1$sample_table$lat_label <- factor(test1$sample_table$lat_label, levels = latitude_levels)
  return(test1)
}

test1 <- load_metagenome_data()


# ==========================================================
# 1. Functional alpha diversity
# ==========================================================
test1$cal_alphadiv(measures = "Shannon")

# Compartment
t_alpha_comp <- trans_alpha$new(dataset = test1, group = "compartment")
t_alpha_comp$cal_diff(method = "KW")
kw_comp <- t_alpha_comp$res_diff

p_comp <- t_alpha_comp$plot_alpha(measure = "Shannon", color_values = c("#fccccb", "#8ccDBF"), xtext_angle = 0, add_sig = FALSE) +
  ggtitle(paste0("Compartment: p = ", formatC(kw_comp$P.adj[1], format = "e", digits = 2))) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())
ggsave(file.path(outdir, "KO_functional_alpha_diversity_compartment.pdf"), p_comp, width = 5, height = 4)

# Location
t_alpha_loc <- trans_alpha$new(dataset = test1, group = "lat_label")
t_alpha_loc$cal_diff(method = "KW")
kw_loc <- t_alpha_loc$res_diff

t_alpha_loc$cal_diff(method = "KW_dunn")

p_loc <- t_alpha_loc$plot_alpha(measure = "Shannon", color_values = c("#fccccb", "#bdb5e1", "#b0d992", "#f9d580", "#8ccDBF"), xtext_angle = 0, add_sig_text_size = 6.5) +
  ggtitle(paste0("Latitude: p = ", formatC(kw_loc$P.adj[1], format = "e", digits = 2))) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())
ggsave(file.path(outdir, "KO_functional_alpha_diversity_location.pdf"), p_loc, width = 5, height = 4)

# Cultivar
t_alpha_var <- trans_alpha$new(dataset = test1, group = "variety")
t_alpha_var$cal_diff(method = "KW")
kw_var <- t_alpha_var$res_diff

t_alpha_var$cal_diff(method = "KW_dunn")
p_var <- t_alpha_var$plot_alpha(measure = "Shannon", color_values = c("#fccccb", "#bdb5e1", "#b0d992", "#f9d580", "#8ccDBF"), xtext_angle = 0, add_sig_text_size = 6.5) +
  ggtitle(paste0("Cultivar: p = ", formatC(kw_var$P.adj[1], format = "e", digits = 2))) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())
ggsave(file.path(outdir, "KO_functional_alpha_diversity_cultivar.pdf"), p_var, width = 5, height = 4)


# ==========================================================
# 2. Dominant KO composition
# ==========================================================

test1$cal_abund(rel = FALSE)

t_abund <- trans_abund$new(dataset = test1, taxrank = "func_descrip", ntaxa = 15, use_percentage = FALSE, groupmean = "group", group_morestats = TRUE)

split_names <- strsplit(as.character(t_abund$data_abund$Sample), "-")
t_abund$data_abund$location <- sapply(split_names, function(x) x[1])
t_abund$data_abund$variety <- sapply(split_names, function(x) x[2])
t_abund$data_abund$compartment <- sapply(split_names, function(x) x[3])
t_abund$data_abund$lat_label <- lat_map[t_abund$data_abund$location]

t_abund$data_abund$compartment <- factor(t_abund$data_abund$compartment, levels = compartment_levels)
t_abund$data_abund$location <- factor(t_abund$data_abund$location, levels = location_levels)
t_abund$data_abund$variety <- factor(t_abund$data_abund$variety, levels = cultivar_levels)
t_abund$data_abund$lat_label <- factor(t_abund$data_abund$lat_label, levels = latitude_levels)

sample_order <- c()
for (loc in location_levels) {
  for (var in cultivar_levels) sample_order <- c(sample_order, paste(loc, var, compartment_levels, sep = "-"))
}

ko_colors <- c("#8cd0c3", "#3d9f3c", "#7087e4", "#cb8d8b", "#b2446b", "#0c54a5", "#faf5b5", "#bcb9d8",
               "#f18072", "#e31a1c", "#d7d7d5", "#fcc3b4", "#f49513", "#ba7fb5", "#b3d46b")

p_abund <- t_abund$plot_bar(facet = c("compartment", "lat_label"), bar_full = FALSE, xtext_angle = 45, xtext_size = 7, barwidth = 1, order_x = sample_order) +
  scale_fill_manual(values = ko_colors) +
  ylab("Abundance (TPM)") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 8.5), axis.text.y = element_text(size = 15))

write.table(t_abund$data_abund, file.path(outdir, "top15_KO_abundance.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
ggsave(file.path(outdir, "top15_KO_abundance.pdf"), p_abund, width = 14.5, height = 6)


# ==========================================================
# 3. Functional beta diversity
# ==========================================================

test1$cal_betadiv(unifrac = FALSE)

t_beta <- trans_beta$new(dataset = test1, measure = "bray")
t_beta$cal_ordination(method = "PCoA")

p_beta <- t_beta$plot_ordination(plot_color = "lat_label", plot_shape = "variety", plot_type = "point",
                                 point_alpha = 1, point_size = 3,
                                 color_values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
                                 shape_values = c(21, 22, 23, 24, 25)) +
  geom_point(aes(color = lat_label, shape = variety, fill = compartment), size = 4, stroke = 1.5) +
  scale_fill_manual(values = c(B = "#fccccb", R = "#8ccDBF")) +
  stat_ellipse(aes(group = compartment), level = 0.90, linewidth = 1, linetype = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")) +
  labs(color = "Latitude", fill = "Compartment", shape = "Cultivar") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), legend.position = "right")

ggsave(file.path(outdir, "KO_BrayCurtis_PCoA.pdf"), p_beta, width = 10, height = 7)

permanova_all <- data.frame()
anosim_all <- data.frame()

for (group_var in c("location", "variety", "compartment")) {
  t_beta$cal_manova(manova_all = TRUE, group = group_var)
  tmp <- as.data.frame(t_beta$res_manova)
  tmp$Factor <- group_var
  permanova_all <- rbind(permanova_all, tmp)

  t_beta$cal_anosim(paired = FALSE, group = group_var)
  tmp <- as.data.frame(t_beta$res_anosim)
  tmp$Factor <- group_var
  anosim_all <- rbind(anosim_all, tmp)
}

write.table(permanova_all, file.path(outdir, "KO_BrayCurtis_PERMANOVA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(anosim_all, file.path(outdir, "KO_BrayCurtis_ANOSIM.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# ==========================================================
# 4. Location-preferential KOs and overlap
# ==========================================================

varcom_levels <- c("H43-B", "H43-R", "H84-B", "H84-R", "JD17-B", "JD17-R", "TL1-B", "TL1-R", "HC6-B", "HC6-R")
all_loc_pref <- data.frame()

for (varcom in varcom_levels) {

  sub_group <- clone(test1)
  sub_group$sample_table <- subset(sub_group$sample_table, `var-com` == varcom)
  sub_group$tidy_dataset()
  sub_group$sample_table$location <- factor(sub_group$sample_table$location, levels = location_levels)
  sub_group$cal_abund(rel = FALSE)
  var <- strsplit(varcom, "-", fixed = TRUE)[[1]][1]
  com <- strsplit(varcom, "-", fixed = TRUE)[[1]][2]
  home_loc <- as.character(var2loc[var])

  t_pref_abund <- trans_abund$new(dataset = sub_group, taxrank = "KO_gene", ntaxa = Inf, use_percentage = FALSE, groupmean = "location", group_morestats = TRUE)

  highest_abund_loc_ko <- t_pref_abund$data_abund %>%
    group_by(Taxonomy) %>%
    filter(Abundance[Sample == home_loc] > max(Abundance[Sample != home_loc])) %>%
    ungroup() %>%
    distinct(Taxonomy) %>%
    pull(Taxonomy)

  sub_group1 <- clone(sub_group)
  sub_group1$tax_table <- sub_group1$tax_table %>% filter(KO_gene %in% highest_abund_loc_ko)
  sub_group1$tidy_dataset()
  sub_group1$cal_abund(rel = FALSE)

  t_pref_diff <- trans_diff$new(dataset = sub_group1, method = "lefse", group = "location", taxa_level = "func_descrip")

  loc_pref_results <- t_pref_diff$res_diff %>%
    filter(Group == home_loc) %>%
    mutate(Compartment = com, Cultivar = var, `Home site` = home_loc) %>%
    select(Compartment, Cultivar, `Home site`, Taxa, Method, LDA, P.unadj, P.adj, Significance)

  all_loc_pref <- rbind(all_loc_pref, loc_pref_results)
}

all_loc_pref_final <- all_loc_pref %>%
  separate(Taxa, into = c("KO info", "Function description"), sep = "\\|", extra = "merge", fill = "right") %>%
  separate(`KO info`, into = c("KO ID", "Gene symbol"), sep = ":", extra = "merge", fill = "right") %>%
  select(Compartment, Cultivar, `Home site`, `KO ID`, `Gene symbol`, `Function description`, Method, LDA, P.unadj, P.adj, Significance)

write.table(all_loc_pref_final, file.path(outdir, "location_preferential_KOs.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# UpSet plot
all_loc_pref_final$Comb <- paste(all_loc_pref_final$Cultivar, all_loc_pref_final$Compartment, sep = "_")
full_combos <- c("H43_B", "H43_R", "H84_B", "H84_R", "JD17_B", "JD17_R", "TL1_B", "TL1_R", "HC6_B", "HC6_R")

df_use <- all_loc_pref_final[, c("KO ID", "Comb")]
taxa <- unique(df_use$`KO ID`)
mat_pref <- matrix(0, nrow = length(taxa), ncol = length(full_combos), dimnames = list(taxa, full_combos))

for (i in seq_len(nrow(df_use))) mat_pref[df_use$`KO ID`[i], df_use$Comb[i]] <- 1
mat_pref <- as.data.frame(mat_pref)

sets_order <- c("H43_B", "H84_B", "JD17_B", "TL1_B", "HC6_B", "H43_R", "H84_R", "JD17_R", "TL1_R", "HC6_R")
sets_cols <- ifelse(grepl("_B$", sets_order), "#FFDD93", "#B279A2")

pdf(file.path(outdir, "location_preferential_KO_upset.pdf"), width = 8.33, height = 6.67)
upset(mat_pref, sets = sets_order, nsets = length(sets_order), keep.order = TRUE,
      nintersects = 30, order.by = "freq", decreasing = TRUE, mb.ratio = c(0.52, 0.48),
      main.bar.color = "skyblue", sets.bar.color = sets_cols, matrix.color = "grey25",
      point.size = 2.5, line.size = 1.2,
      text.scale = c(2.3, 1.8, 2, 1.8, 1.6, 1.25),
      mainbar.y.label = "Intersection Size", sets.x.label = "Set Size")
dev.off()

# ==========================================================
# 5. Nitrogen, phosphorus and potassium cycling functions
# ==========================================================

npk_ko <- read.table(npk_ko_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

npk_test <- clone(test1)
npk_test$tax_table <- npk_test$tax_table[rownames(npk_test$tax_table) %in% npk_ko$KO, ]
npk_test$tidy_dataset()
npk_test$cal_abund(rel = FALSE)

func_levels <- c("N mineralization", "N fixation", "N reduction", "Nitrification", "Denitrification",
                 "Nitrogen limitation", "Organic P mineralization", "Inorganic P solubilization",
                 "P transportation", "P starvation", "K transport")

func_col <- c(
  "N mineralization" = "#7b3294",
  "N fixation" = "#c2a5cf",
  "N reduction" = "#fdae61",
  "Nitrification" = "#1b9e77",
  "Denitrification" = "#d95f02",
  "Nitrogen limitation" = "#7570b3",
  "Organic P mineralization" = "#e7298a",
  "Inorganic P solubilization" = "#66a61e",
  "P transportation" = "#e6ab02",
  "P starvation" = "#a6761d",
  "K transport" = "#1f78b4"
)

lefse_col <- c(He = "#E41A1C", Ha = "#377EB8", BJ = "#4DAF4A", WH = "#984EA3", GZ = "#FF7F00", NS = "grey85")
hm_col <- colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100)

for (com in compartment_levels) {

  sub_group <- clone(npk_test)
  sub_group$sample_table <- subset(sub_group$sample_table, compartment == com)
  sub_group$tidy_dataset()
  sub_group$cal_abund(rel = FALSE)

  otu <- as.matrix(sub_group$taxa_abund$KO_gene)
  meta <- sub_group$sample_table

  stopifnot(all(colnames(otu) %in% rownames(meta)))
  meta <- meta[colnames(otu), , drop = FALSE]

  loc_samples <- split(rownames(meta), meta$location)
  mean_mat <- sapply(loc_samples, function(samps) rowMeans(otu[, samps, drop = FALSE], na.rm = TRUE))
  mean_mat <- as.matrix(mean_mat)
  mean_mat <- mean_mat[, location_levels, drop = FALSE]

  ko_id <- sub(":.*", "", rownames(mean_mat))
  func <- npk_ko$Function[match(ko_id, npk_ko$KO)]
  func <- factor(func, levels = func_levels)

  ord <- order(func, ko_id)
  mean_mat <- mean_mat[ord, , drop = FALSE]

  anno_row <- data.frame(Function = func[ord])
  rownames(anno_row) <- rownames(mean_mat)

  mat_plot <- t(scale(t(mean_mat)))
  mat_plot[is.na(mat_plot)] <- 0

  t_npk_diff <- trans_diff$new(dataset = sub_group, method = "lefse", group = "location", taxa_level = "KO_gene")
  lefse_results <- t_npk_diff$res_diff[, c("Taxa", "Group")]
  lefse_map <- setNames(as.character(lefse_results$Group), lefse_results$Taxa)

  anno_row$LEfSe <- lefse_map[rownames(mean_mat)]
  anno_row$LEfSe[is.na(anno_row$LEfSe)] <- "NS"

  pref_df <- all_loc_pref %>% filter(Compartment == com) %>% transmute(Group = `Home site`, Taxa = sub("\\|.*", "", Taxa)) %>% distinct()

  tri_mat <- matrix("", nrow = nrow(mean_mat), ncol = ncol(mean_mat), dimnames = dimnames(mean_mat))

  for (i in seq_len(nrow(pref_df))) {
    ko <- pref_df$Taxa[i]
    loc <- pref_df$Group[i]
    if (ko %in% rownames(tri_mat) & loc %in% colnames(tri_mat)) tri_mat[ko, loc] <- "★"
  }

  ann_colors <- list(Function = func_col, LEfSe = lefse_col)

  g <- pheatmap(mat_plot,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_row = anno_row,
                annotation_colors = ann_colors,
                color = hm_col,
                display_numbers = tri_mat,
                number_color = "black",
                fontsize_number = 10,
                fontsize_row = 9,
                fontsize_col = 10,
                cellwidth = 15,
                cellheight = 10)

  write.table(mean_mat, file.path(outdir, paste0(com, "_NPK_mean_TPM_by_location.tsv")), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
  ggsave(file.path(outdir, paste0(com, "_NPK_cycle_heatmap.pdf")), g, width = 5.5, height = 10)
}
