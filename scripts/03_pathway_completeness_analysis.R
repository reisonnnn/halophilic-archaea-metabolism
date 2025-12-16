library(tidyverse)
library(pheatmap)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
}

cat("=========================================\n")
cat("PATHWAY COMPLETENESS - GHOSTKOALA\n")
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
# PATHWAY DEFINITIONS (all 10 steps for glycolysis)
# ==============================================================================

pathways <- list(
  
  "Classical glycolysis" = list(
    steps = list(
      "1_Hexokinase" = c("K00844"),
      "2_GPI" = c("K01810"),
      "3_PFK" = c("K00850"),
      "4_Aldolase" = c("K01623", "K01624"),
      "5_TPI" = c("K01803"),
      "6_GAPDH" = c("K00134", "K00150"),
      "7_PGK" = c("K00927"),
      "8_PGAM" = c("K01834", "K01689"),
      "9_Enolase" = c("K15633"),
      "10_PK" = c("K00873")
    )
  ),
  
  "Archaeal glycolysis" = list(
    steps = list(
      "1_Glucokinase" = c("K00845"),
      "2_GPI" = c("K01810"),
      "3_ADP-PFK" = c("K00918"),
      "4_Aldolase" = c("K01622", "K01623", "K01624"),
      "5_TPI" = c("K01803"),
      "6_GAPDH" = c("K00134", "K00150"),
      "7_PGK" = c("K00927"),
      "8_PGAM" = c("K01689"),
      "9_Enolase" = c("K15633"),
      "10_PK" = c("K00873")
    )
  ),
  
  "Entner-Doudoroff" = list(
    steps = list(
      "1_GPI" = c("K01810"),
      "2_G6PDH" = c("K00036"),
      "3_6PGL" = c("K01057"),
      "4_KDPG_aldolase" = c("K01625"),
      "5_TPI" = c("K01803"),
      "6_GAPDH" = c("K00134", "K00150"),
      "7_PGK" = c("K00927"),
      "8_PGAM" = c("K01834", "K01689"),
      "9_Enolase" = c("K15633"),
      "10_PK" = c("K00873")
    )
  ),
  
  "Gluconeogenesis" = list(
    steps = list(
      "1_PEPCK" = c("K01596"),
      "2_PGAM" = c("K01689", "K01834"),
      "3_PGK" = c("K00927"),
      "4_GAPDH" = c("K00134", "K00150"),
      "5_Aldolase" = c("K01623", "K01624"),
      "6_PFK" = c("K00850"),
      "7_GPI" = c("K01810"),
      "8_Hexokinase" = c("K00844")
    )
  ),
  
  "Glycerol metabolism" = list(
    steps = list(
      "1_GlpK" = c("K00864"),
      "2_GlpD" = c("K00111"),
      "3_TPI" = c("K01803"),
      "4_GAPDH" = c("K00134", "K00150"),
      "5_PGK" = c("K00927"),
      "6_PGAM" = c("K01689", "K01834"),
      "7_Enolase" = c("K15633"),
      "8_PK" = c("K00873")
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

cat("Pathway completeness:\n")
for (name in names(pathways)) {
  pw <- pathways[[name]]
  
  hs_pct <- calc_completeness(pw$steps, hs_kos)
  hm_pct <- calc_completeness(pw$steps, hm_kos)
  ha_pct <- calc_completeness(pw$steps, ha_kos)
  
  results <- rbind(results, data.frame(
    Pathway = rep(name, 3),
    Strain = c("HS", "HM", "HA"),
    Completeness = c(hs_pct, hm_pct, ha_pct)
  ))
  
  cat(sprintf("%-25s HS: %3.0f%%  HM: %3.0f%%  HA: %3.0f%%\n",
              paste0(name, ":"), hs_pct, hm_pct, ha_pct))
}

results$Strain <- factor(results$Strain, levels = c("HS", "HM", "HA"))

# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

cat("\nCreating figures...\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Heatmap
heatmap_data <- results %>%
  pivot_wider(names_from = Strain, values_from = Completeness) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

heatmap_data <- heatmap_data[, c("HS", "HM", "HA")]

png("results/figures/ghostkoala_pathway_heatmap.png",
    width = 2400, height = 2000, res = 300)

pheatmap(heatmap_data,
         color = colorRampPalette(c("#FFFFFF", "#FFF59D", "#FFD54F",
                                    "#FFA726", "#FF6F00", "#D84315"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 11,
         fontsize_col = 13,
         border_color = "grey60",
         cellwidth = 60,
         cellheight = 35,
         main = "Pathway Completeness (%)",
         labels_col = c("HS\n(Halobacterium)", "HM\n(Halomicrobium)", 
                       "HA\n(Haloarcula)"),
         display_numbers = TRUE,
         number_format = "%.0f",
         fontsize_number = 10,
         breaks = seq(0, 100, length.out = 101))

dev.off()
cat("  Heatmap saved\n")

# Bar plot
p <- ggplot(results, aes(x = Pathway, y = Completeness, fill = Strain)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", Completeness)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3) +
  geom_hline(yintercept = 100, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  scale_fill_manual(values = c("HS" = "#FF6B6B", 
                               "HM" = "#4ECDC4", 
                               "HA" = "#95E1D3")) +
  labs(title = "Pathway Completeness",
       x = "", y = "Completeness (%)", fill = "") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  ylim(0, 110)

ggsave("results/figures/ghostkoala_pathway_barplot.png", p,
       width = 12, height = 7, dpi = 300)
cat("  Bar plot saved\n")

# Save table
write.csv(results, "results/tables/ghostkoala_pathway_completeness.csv",
          row.names = FALSE)
cat("  Table saved\n")

cat("\n=========================================\n")
cat("COMPLETE\n")
cat("=========================================\n")