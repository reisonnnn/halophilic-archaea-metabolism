library(tidyverse)
library(pheatmap)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
}

cat("=========================================\n")
cat("AMINO ACID BIOSYNTHESIS - GHOSTKOALA\n")
cat("Based on KEGG map01230\n")
cat("=========================================\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

load_ko <- function(filepath) {
  read.delim(filepath, header = FALSE, 
             col.names = c("gene_id", "KO"),
             stringsAsFactors = FALSE) %>%
    filter(KO != "") %>%
    pull(KO) %>%
    unique()
}

hs_kos <- load_ko("annotation_outputs/GhostKOALA/HS_GhostKOALA.txt")
hm_kos <- load_ko("annotation_outputs/GhostKOALA/HM_GhostKOALA.txt")
ha_kos <- load_ko("annotation_outputs/GhostKOALA/HA_GhostKOALA.txt")

cat(sprintf("HS: %d KOs\n", length(hs_kos)))
cat(sprintf("HM: %d KOs\n", length(hm_kos)))
cat(sprintf("HA: %d KOs\n\n", length(ha_kos)))

# ==============================================================================
# AMINO ACID FAMILIES (following KEGG map01230 logic)
# ==============================================================================

aa_families <- list(
  
  "ASPARTATE FAMILY" = list(
    
    "Aspartate" = list(
      steps = list(
        "AspC" = c("K00812", "K00813")
      )
    ),
    
    "Asparagine" = list(
      steps = list(
        "AspC" = c("K00812", "K00813"),
        "AsnA_or_AsnB" = c("K01914", "K01953")
      )
    ),
    
    "Lysine (DAP)" = list(
      steps = list(
        "AspC" = c("K00812", "K00813"),
        "LysC_AK" = c("K00928"),
        "Asd" = c("K00133"),
        "DapA" = c("K00674"),
        "DapB" = c("K00215"),
        "DapD" = c("K00674"),
        "DapE_or_DapF" = c("K01439", "K01778"),
        "LysA" = c("K01586")
      )
    ),
    
    "Threonine" = list(
      steps = list(
        "AspC" = c("K00812", "K00813"),
        "ThrA_AK" = c("K00928"),
        "Asd" = c("K00133"),
        "ThrA_HDH" = c("K00133"),
        "ThrB" = c("K00133"),
        "ThrC" = c("K01733")
      )
    ),
    
    "Isoleucine" = list(
      steps = list(
        "AspC" = c("K00812", "K00813"),
        "ThrA_AK" = c("K00928"),
        "Asd" = c("K00133"),
        "ThrA_HDH" = c("K00133"),
        "ThrB" = c("K00133"),
        "ThrC" = c("K01733"),
        "IlvA" = c("K01754"),
        "IlvB_or_IlvI" = c("K01652", "K01653"),
        "IlvC" = c("K00053"),
        "IlvD" = c("K01687"),
        "IlvE" = c("K00826")
      )
    ),
    
    "Methionine" = list(
      steps = list(
        "AspC" = c("K00812", "K00813"),
        "ThrA_AK" = c("K00928"),
        "Asd" = c("K00133"),
        "MetA" = c("K00651"),
        "MetB_or_MetY" = c("K00548", "K17217"),
        "MetC" = c("K01760"),
        "MetE_or_MetH" = c("K00549", "K00548")
      )
    )
  ),
  
  "GLUTAMATE FAMILY" = list(
    
    "Glutamate" = list(
      steps = list(
        "GltB_or_GltS" = c("K00260", "K00261", "K00262", "K15371")
      )
    ),
    
    "Glutamine" = list(
      steps = list(
        "GltB_or_GltS" = c("K00260", "K00261", "K00262", "K15371"),
        "GlnA" = c("K01915")
      )
    ),
    
    "Proline" = list(
      steps = list(
        "GltB_or_GltS" = c("K00260", "K00261", "K00262", "K15371"),
        "ProB" = c("K00931"),
        "ProA" = c("K00147"),
        "ProC" = c("K00286")
      )
    ),
    
    "Arginine" = list(
      steps = list(
        "GltB_or_GltS" = c("K00260", "K00261", "K00262", "K15371"),
        "ArgB_or_ArgJ" = c("K00145", "K00620"),
        "ArgC" = c("K00821"),
        "ArgD" = c("K00821"),
        "ArgE" = c("K00611"),
        "ArgG_or_ArgH" = c("K01940", "K01755")
      )
    )
  ),
  
  "PYRUVATE FAMILY" = list(
    
    "Alanine" = list(
      steps = list(
        "AlaT" = c("K00259", "K14260")
      )
    ),
    
    "Valine" = list(
      steps = list(
        "IlvB_or_IlvI" = c("K01652", "K01653"),
        "IlvC" = c("K00053"),
        "IlvD" = c("K01687"),
        "IlvE" = c("K00826")
      )
    ),
    
    "Leucine" = list(
      steps = list(
        "IlvB_or_IlvI" = c("K01652", "K01653"),
        "IlvC" = c("K00053"),
        "IlvD" = c("K01687"),
        "LeuA" = c("K01649"),
        "LeuC" = c("K01703"),
        "LeuD" = c("K01704"),
        "LeuB" = c("K00052"),
        "IlvE" = c("K00826")
      )
    )
  ),
  
  "3PG FAMILY" = list(
    
    "Serine" = list(
      steps = list(
        "SerA" = c("K00058"),
        "SerC" = c("K00831")
      )
    ),
    
    "Glycine" = list(
      steps = list(
        "SerA" = c("K00058"),
        "SerC" = c("K00831"),
        "GlyA" = c("K00600")
      )
    ),
    
    "Cysteine" = list(
      steps = list(
        "SerA" = c("K00058"),
        "SerC" = c("K00831"),
        "CysE" = c("K00640"),
        "CysK_or_CysM" = c("K01738", "K12339")
      )
    )
  ),
  
  "AROMATIC FAMILY" = list(
    
    "Phenylalanine" = list(
      steps = list(
        "AroG_or_AroF" = c("K00766", "K03856"),
        "Chorismate_pathway" = c("K01626", "K00800", "K01735", "K01736"),
        "PheA" = c("K14170"),
        "TyrB_or_PheA" = c("K00812", "K00813", "K14170")
      )
    ),
    
    "Tyrosine" = list(
      steps = list(
        "AroG_or_AroF" = c("K00766", "K03856"),
        "Chorismate_pathway" = c("K01626", "K00800", "K01735", "K01736"),
        "TyrA" = c("K00220"),
        "TyrB" = c("K00812", "K00813")
      )
    ),
    
    "Tryptophan" = list(
      steps = list(
        "AroG_or_AroF" = c("K00766", "K03856"),
        "Chorismate_pathway" = c("K01626", "K00800", "K01735", "K01736"),
        "TrpE_or_TrpG" = c("K01657", "K01658"),
        "TrpD_or_TrpG" = c("K00766", "K01658"),
        "TrpC_or_TrpF" = c("K01609", "K01696"),
        "TrpA" = c("K01695"),
        "TrpB" = c("K01696")
      )
    )
  ),
  
  "HISTIDINE" = list(
    
    "Histidine" = list(
      steps = list(
        "HisG" = c("K00765", "K16841"),
        "HisI" = c("K01814"),
        "HisA" = c("K02502"),
        "HisH" = c("K00817"),
        "HisF" = c("K02500"),
        "HisB" = c("K00766"),
        "HisC" = c("K00817"),
        "HisD" = c("K00013")
      )
    )
  )
)

# ==============================================================================
# CALCULATE COMPLETENESS
# ==============================================================================

calc_completeness <- function(pathway_steps, strain_kos) {
  steps_present <- 0
  total_steps <- length(pathway_steps)
  
  for (step_kos in pathway_steps) {
    if (any(step_kos %in% strain_kos)) {
      steps_present <- steps_present + 1
    }
  }
  
  (steps_present / total_steps) * 100
}

results <- data.frame()

cat("Amino acid biosynthesis by family:\n")
cat("=====================================\n\n")

for (family_name in names(aa_families)) {
  cat(family_name, "\n")
  cat(strrep("-", nchar(family_name)), "\n")
  
  family <- aa_families[[family_name]]
  
  for (aa_name in names(family)) {
    pw <- family[[aa_name]]
    
    hs_pct <- calc_completeness(pw$steps, hs_kos)
    hm_pct <- calc_completeness(pw$steps, hm_kos)
    ha_pct <- calc_completeness(pw$steps, ha_kos)
    
    results <- rbind(results, data.frame(
      Family = family_name,
      Amino_acid = aa_name,
      Strain = c("HS", "HM", "HA"),
      Completeness = c(hs_pct, hm_pct, ha_pct)
    ))
    
    cat(sprintf("  %-20s HS: %3.0f%%  HM: %3.0f%%  HA: %3.0f%%\n",
                aa_name, hs_pct, hm_pct, ha_pct))
  }
  cat("\n")
}

results$Strain <- factor(results$Strain, levels = c("HS", "HM", "HA"))

# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

cat("Creating figures...\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Heatmap
heatmap_data <- results %>%
  pivot_wider(names_from = Strain, values_from = Completeness) %>%
  unite("AA", Family, Amino_acid, sep = " - ") %>%
  column_to_rownames("AA") %>%
  as.matrix()

heatmap_data <- heatmap_data[, c("HS", "HM", "HA")]

png("results/figures/aa_biosynthesis_heatmap.png",
    width = 2800, height = 3500, res = 300)

pheatmap(heatmap_data,
         color = colorRampPalette(c("#FFFFFF", "#FFF59D", "#FFD54F",
                                    "#FFA726", "#FF6F00", "#D84315"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 9,
         fontsize_col = 13,
         border_color = "grey60",
         cellwidth = 60,
         cellheight = 22,
         main = "Amino Acid Biosynthesis Completeness (%)",
         labels_col = c("HS\n(Halobacterium)", "HM\n(Halomicrobium)", 
                       "HA\n(Haloarcula)"),
         display_numbers = TRUE,
         number_format = "%.0f",
         fontsize_number = 8,
         breaks = seq(0, 100, length.out = 101))

dev.off()
cat("  Heatmap saved\n")

# Bar plot
p <- ggplot(results, aes(x = Amino_acid, y = Completeness, fill = Strain)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", Completeness)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2) +
  geom_hline(yintercept = 100, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  scale_fill_manual(values = c("HS" = "#FF6B6B", 
                               "HM" = "#4ECDC4", 
                               "HA" = "#95E1D3")) +
  facet_wrap(~Family, scales = "free_x", ncol = 2) +
  labs(title = "Amino Acid Biosynthesis by Family",
       x = "", y = "Completeness (%)", fill = "") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        strip.text = element_text(face = "bold")) +
  ylim(0, 110)

ggsave("results/figures/aa_biosynthesis_barplot.png", p,
       width = 16, height = 12, dpi = 300)
cat("  Bar plot saved\n")

# Comprehensive 20 amino acid comparison
aa_order <- c("Alanine", "Arginine", "Asparagine", "Aspartate",
              "Cysteine", "Glutamate", "Glutamine", "Glycine",
              "Histidine", "Isoleucine", "Leucine", "Lysine (DAP)",
              "Methionine", "Phenylalanine", "Proline", "Serine",
              "Threonine", "Tryptophan", "Tyrosine", "Valine")

p_comprehensive <- results %>%
  mutate(Amino_acid = factor(Amino_acid, levels = aa_order)) %>%
  ggplot(aes(x = Amino_acid, y = Completeness, fill = Strain)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", Completeness)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.5, angle = 0) +
  geom_hline(yintercept = 100, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  geom_hline(yintercept = 50, linetype = "dotted", 
             color = "orange", alpha = 0.3) +
  scale_fill_manual(values = c("HS" = "#FF6B6B", 
                               "HM" = "#4ECDC4", 
                               "HA" = "#95E1D3"),
                    labels = c("HS (Halobacterium)", 
                              "HM (Halomicrobium)", 
                              "HA (Haloarcula)")) +
  labs(title = "Amino Acid Biosynthesis Capacity Across Three Halophilic Archaea",
       subtitle = "Completeness of biosynthetic pathways for all 20 amino acids",
       x = "", 
       y = "Pathway Completeness (%)", 
       fill = "") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank()) +
  ylim(0, 115)

ggsave("results/figures/aa_biosynthesis_comprehensive.png", 
       p_comprehensive,
       width = 16, height = 8, dpi = 300)
cat("  Comprehensive 20-AA plot saved\n")

# Summary by family
family_summary <- results %>%
  group_by(Family, Strain) %>%
  summarise(
    Mean_completeness = mean(Completeness),
    Complete_AAs = sum(Completeness == 100),
    Total_AAs = n(),
    .groups = 'drop'
  )

write.csv(family_summary, 
          "results/tables/aa_biosynthesis_family_summary.csv",
          row.names = FALSE)

write.csv(results, 
          "results/tables/aa_biosynthesis_completeness.csv",
          row.names = FALSE)
cat("  Tables saved\n")

cat("\n=========================================\n")
cat("FAMILY SUMMARY\n")
cat("=========================================\n\n")
print(family_summary, n = Inf)

cat("\n=========================================\n")
cat("COMPLETE\n")
cat("=========================================\n")