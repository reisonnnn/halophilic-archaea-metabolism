library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Install eulerr for proportional Venn diagrams
if (!require("eulerr", quietly = TRUE)) {
  install.packages("eulerr")
  library(eulerr)
}

# Install pheatmap for heatmaps
if (!require("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
  library(pheatmap)
}

project_root <- getwd()

cat("=========================================\n")
cat("COMPARATIVE STRAIN ANALYSIS - UPDATED\n")
cat("With Proportional Venn Diagrams!\n")
cat("=========================================\n\n")

# Strain information
strain_info <- data.frame(
  strain_id = c("hs", "hm", "ha"),
  full_name = c("Halobacterium salinarum", "Halomicrobium kobeituziens", "Haloarcula sp."),
  kbtz_id = c("KBTZ01", "KBTZ05", "KBTZ06"),
  stringsAsFactors = FALSE
)

cat("Analyzing strains:\n")
print(strain_info)
cat("\n")

all_annotations <- read.csv(
  file.path(project_root, "data/processed/metabolic_pathways/all_eggnog_annotations.csv"),
  stringsAsFactors = FALSE
)

cat("1. ANNOTATION COVERAGE\n")
cat("----------------------\n\n")

annotation_stats <- all_annotations %>%
  group_by(strain_id) %>%
  summarise(
    total_proteins = n(),
    annotated = sum(!is.na(Description) & Description != "" & Description != "-"),
    with_kegg_pathway = sum(!is.na(KEGG_Pathway) & KEGG_Pathway != "" & KEGG_Pathway != "-"),
    with_ko = sum(!is.na(KEGG_ko) & KEGG_ko != "" & KEGG_ko != "-"),
    with_cog = sum(!is.na(COG_category) & COG_category != "" & COG_category != "-"),
    unannotated = sum(is.na(Description) | Description == "" | Description == "-")
  ) %>%
  mutate(
    pct_annotated = round(annotated / total_proteins * 100, 1),
    pct_with_pathway = round(with_kegg_pathway / total_proteins * 100, 1),
    pct_with_ko = round(with_ko / total_proteins * 100, 1)
  ) %>%
  left_join(strain_info, by = "strain_id")

cat("Annotation Statistics:\n")
print(annotation_stats)
cat("\n")

write.csv(annotation_stats, 
          file.path(project_root, "results/tables/annotation_coverage.csv"),
          row.names = FALSE)

# Create plot with proper strain labels
p1 <- annotation_stats %>%
  select(strain_id, full_name, annotated, unannotated) %>%
  pivot_longer(cols = c(annotated, unannotated), 
               names_to = "status", values_to = "count") %>%
  ggplot(aes(x = strain_id, y = count, fill = status)) +
  geom_col(position = "stack") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), 
            color = "white", fontface = "bold", size = 5) +
  scale_fill_manual(values = c("annotated" = "#2E7D32", "unannotated" = "#D32F2F"),
                    labels = c("Annotated", "Unannotated")) +
  labs(title = "Annotation Coverage by Strain",
       subtitle = paste0("HS: ", annotation_stats$pct_annotated[1], "% | ",
                         "HM: ", annotation_stats$pct_annotated[2], "% | ",
                         "HA: ", annotation_stats$pct_annotated[3], "%"),
       x = "Strain", y = "Number of Proteins", fill = "Status") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(project_root, "results/figures/annotation_coverage.png"),
       p1, width = 10, height = 6, dpi = 300)

cat("✓ Annotation coverage plot saved\n\n")

cat("2. PROPORTIONAL VENN DIAGRAM (Euler Diagram)\n")
cat("---------------------------------------------\n\n")

# Get orthologs for each strain
hs_orthologs <- all_annotations %>%
  filter(strain_id == "hs", 
         !is.na(eggNOG_OGs) & eggNOG_OGs != "" & eggNOG_OGs != "-") %>%
  pull(eggNOG_OGs) %>%
  unique()

hm_orthologs <- all_annotations %>%
  filter(strain_id == "hm",
         !is.na(eggNOG_OGs) & eggNOG_OGs != "" & eggNOG_OGs != "-") %>%
  pull(eggNOG_OGs) %>%
  unique()

ha_orthologs <- all_annotations %>%
  filter(strain_id == "ha",
         !is.na(eggNOG_OGs) & eggNOG_OGs != "" & eggNOG_OGs != "-") %>%
  pull(eggNOG_OGs) %>%
  unique()

cat("Unique orthologs per strain:\n")
cat(sprintf("HS (Halobacterium): %d\n", length(hs_orthologs)))
cat(sprintf("HM (Halomicrobium): %d\n", length(hm_orthologs)))
cat(sprintf("HA (Haloarcula): %d\n", length(ha_orthologs)))
cat("\n")

# Calculate overlaps for Euler diagram (proportional)
euler_fit <- euler(list(
  "HS" = hs_orthologs,
  "HM" = hm_orthologs,
  "HA" = ha_orthologs
))

# Create proportional Euler diagram
png(file.path(project_root, "results/figures/venn_ortholog_overlap_proportional.png"),
    width = 3000, height = 3000, res = 300)

plot(euler_fit,
     fills = list(fill = c("#FF6B6B", "#95E1D3", "#A8E6CF"), alpha = 0.6),
     quantities = list(fontsize = 14, fontface = "bold"),
     labels = list(fontsize = 16, fontface = "bold"),
     edges = list(col = "white", lwd = 2),
     main = "Metabolic Gene Overlaps",
     cex.main = 1.5)

dev.off()

cat("✓ Proportional Euler diagram saved\n\n")

# Calculate overlap statistics
overlap_stats <- list(
  HS_only = length(setdiff(setdiff(hs_orthologs, hm_orthologs), ha_orthologs)),
  HM_only = length(setdiff(setdiff(hm_orthologs, hs_orthologs), ha_orthologs)),
  HA_only = length(setdiff(setdiff(ha_orthologs, hs_orthologs), hm_orthologs)),
  HS_HM = length(setdiff(intersect(hs_orthologs, hm_orthologs), ha_orthologs)),
  HS_HA = length(setdiff(intersect(hs_orthologs, ha_orthologs), hm_orthologs)),
  HM_HA = length(setdiff(intersect(hm_orthologs, ha_orthologs), hs_orthologs)),
  All_three = length(Reduce(intersect, list(hs_orthologs, hm_orthologs, ha_orthologs)))
)

cat("Overlap Statistics:\n")
cat(sprintf("Core (all 3 strains): %d genes\n", overlap_stats$All_three))
cat(sprintf("HS unique: %d genes\n", overlap_stats$HS_only))
cat(sprintf("HM unique: %d genes\n", overlap_stats$HM_only))
cat(sprintf("HA unique: %d genes\n", overlap_stats$HA_only))
cat(sprintf("HS-HM shared: %d genes\n", overlap_stats$HS_HM))
cat(sprintf("HS-HA shared: %d genes\n", overlap_stats$HS_HA))
cat(sprintf("HM-HA shared: %d genes\n", overlap_stats$HM_HA))
cat("\n")

# Save overlap stats
overlap_df <- data.frame(
  Category = names(overlap_stats),
  Gene_Count = unlist(overlap_stats)
)
write.csv(overlap_df,
          file.path(project_root, "results/tables/ortholog_overlap_stats.csv"),
          row.names = FALSE)

cat("3. GLYCEROL vs GLUCOSE METABOLISM\n")
cat("---------------------------------\n\n")

glycerol_metabolism <- all_annotations %>%
  filter(grepl("ko00561", KEGG_Pathway, ignore.case = TRUE)) %>%
  select(strain_id, X.query, KEGG_ko, KEGG_Pathway, EC, Description)

cat("Glycerol metabolism genes by strain:\n")
print(table(glycerol_metabolism$strain_id))
cat("\n")

glucose_metabolism <- all_annotations %>%
  filter(grepl("ko00010", KEGG_Pathway, ignore.case = TRUE)) %>%
  select(strain_id, X.query, KEGG_ko, KEGG_Pathway, EC, Description)

cat("Glucose/Glycolysis metabolism genes by strain:\n")
print(table(glucose_metabolism$strain_id))
cat("\n")

metabolism_comparison <- data.frame(
  strain_id = c("hs", "hm", "ha"),
  glycerol_genes = c(
    sum(glycerol_metabolism$strain_id == "hs"),
    sum(glycerol_metabolism$strain_id == "hm"),
    sum(glycerol_metabolism$strain_id == "ha")
  ),
  glucose_genes = c(
    sum(glucose_metabolism$strain_id == "hs"),
    sum(glucose_metabolism$strain_id == "hm"),
    sum(glucose_metabolism$strain_id == "ha")
  )
) %>%
  mutate(
    gly_glu_ratio = round(glycerol_genes / glucose_genes, 2)
  )

cat("Glycerol vs Glucose gene comparison:\n")
print(metabolism_comparison)
cat("\n")

p2 <- metabolism_comparison %>%
  pivot_longer(cols = c(glycerol_genes, glucose_genes),
               names_to = "pathway", values_to = "count") %>%
  ggplot(aes(x = strain_id, y = count, fill = pathway)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7),
            vjust = -0.5, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("glycerol_genes" = "#FF6B6B", 
                               "glucose_genes" = "#4ECDC4"),
                    labels = c("Glucose/Glycolysis", "Glycerol")) +
  labs(title = "Glycerol vs Glucose Metabolism Genes",
       subtitle = "Comparison across three halophilic archaea strains",
       x = "Strain", y = "Number of Genes", fill = "Pathway") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(project_root, "results/figures/glycerol_glucose_comparison.png"),
       p2, width = 10, height = 6, dpi = 300)

cat("✓ Metabolism comparison plot saved\n\n")

write.csv(glycerol_metabolism,
          file.path(project_root, "results/tables/glycerol_metabolism_detailed.csv"),
          row.names = FALSE)

write.csv(glucose_metabolism,
          file.path(project_root, "results/tables/glucose_metabolism_detailed.csv"),
          row.names = FALSE)

write.csv(metabolism_comparison,
          file.path(project_root, "results/tables/metabolism_comparison_summary.csv"),
          row.names = FALSE)

cat("=========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=========================================\n\n")

cat("Files saved:\n")
cat("- results/figures/annotation_coverage.png\n")
cat("- results/figures/venn_ortholog_overlap_proportional.png (NEW!)\n")
cat("- results/figures/glycerol_glucose_comparison.png\n")
cat("- results/tables/annotation_coverage.csv\n")
cat("- results/tables/ortholog_overlap_stats.csv (NEW!)\n")
cat("- results/tables/glycerol_metabolism_detailed.csv\n")
cat("- results/tables/glucose_metabolism_detailed.csv\n")
cat("- results/tables/metabolism_comparison_summary.csv\n\n")

cat("Key Findings:\n")
cat("1. Annotation coverage ranges from", min(annotation_stats$pct_annotated), 
    "% to", max(annotation_stats$pct_annotated), "%\n")
cat("2. Core ortholog group (all 3 strains):", overlap_stats$All_three, "genes\n")
cat("3. Glycerol metabolism genes:", sum(metabolism_comparison$glycerol_genes), "total\n")
cat("4. Glucose metabolism genes:", sum(metabolism_comparison$glucose_genes), "total\n\n")
