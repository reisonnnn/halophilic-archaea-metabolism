
# Load required libraries
library(tidyverse)
library(readxl)
library(VennDiagram)
library(ggplot2)
library(gridExtra)

# Set working directory (adjust as needed)
# setwd("path/to/halophilic-archaea-metabolism")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading annotation files...\n")

# Load GhostKOALA output
ghostkoala <- read.table(
  "annotation_outputs/GhostKOALA/HS_GhostKOALA.txt",
  sep = "\t",
  header = FALSE,
  col.names = c("protein_id", "ko_id"),
  fill = TRUE,
  stringsAsFactors = FALSE
)

# Load eggNOG output
eggnog <- read_excel(
  "annotation_outputs/eggnog/hs.emapper.annotations.xlsx",
  skip = 2  # Skip the header comment lines
)

# Check column names and find the KO column
cat("eggNOG columns found:\n")
cat(paste(names(eggnog), collapse = ", "), "\n\n")

cat("Data loaded successfully!\n")
cat(sprintf("GhostKOALA: %d proteins\n", nrow(ghostkoala)))
cat(sprintf("eggNOG: %d proteins\n", nrow(eggnog)))

# ============================================================================
# 2. EXTRACT AND CLEAN KO IDs
# ============================================================================

cat("\nExtracting KO identifiers...\n")

# Process GhostKOALA data
ghostkoala_clean <- ghostkoala %>%
  filter(!is.na(ko_id) & ko_id != "") %>%
  mutate(ko_id = str_trim(ko_id))

# Process eggNOG data - extract KO column and clean
# Find the correct column names (they might have special characters)
query_col <- names(eggnog)[1]  # First column is the query
ko_col <- names(eggnog)[grepl("KEGG_ko", names(eggnog), ignore.case = TRUE)][1]  # Find KEGG_ko column

cat(sprintf("Using query column: %s\n", query_col))
cat(sprintf("Using KO column: %s\n", ko_col))

eggnog_clean <- eggnog %>%
  filter(!is.na(.data[[ko_col]]) & .data[[ko_col]] != "-") %>%
  select(query = 1, ko_raw = all_of(ko_col)) %>%
  # Split multiple KO IDs (some entries have multiple KOs)
  separate_rows(ko_raw, sep = ",") %>%
  mutate(
    ko_id = str_remove(ko_raw, "ko:"),  # Remove "ko:" prefix
    ko_id = str_trim(ko_id)
  ) %>%
  filter(ko_id != "")

# Get unique KO sets for each tool
ghostkoala_kos <- unique(ghostkoala_clean$ko_id)
eggnog_kos <- unique(eggnog_clean$ko_id)

cat(sprintf("\nGhostKOALA unique KOs: %d\n", length(ghostkoala_kos)))
cat(sprintf("eggNOG unique KOs: %d\n", length(eggnog_kos)))

# ============================================================================
# 3. COMPARE KO ASSIGNMENTS
# ============================================================================

cat("\nComparing KO assignments...\n")

# Find overlap and unique KOs
ko_overlap <- intersect(ghostkoala_kos, eggnog_kos)
ko_ghostkoala_only <- setdiff(ghostkoala_kos, eggnog_kos)
ko_eggnog_only <- setdiff(eggnog_kos, ghostkoala_kos)

# Calculate statistics
total_unique_kos <- length(union(ghostkoala_kos, eggnog_kos))
overlap_percent <- (length(ko_overlap) / total_unique_kos) * 100

cat(sprintf("\n--- COMPARISON RESULTS ---\n"))
cat(sprintf("Total unique KOs (union): %d\n", total_unique_kos))
cat(sprintf("Overlapping KOs: %d (%.1f%%)\n", length(ko_overlap), overlap_percent))
cat(sprintf("GhostKOALA only: %d (%.1f%%)\n", 
            length(ko_ghostkoala_only), 
            (length(ko_ghostkoala_only)/total_unique_kos)*100))
cat(sprintf("eggNOG only: %d (%.1f%%)\n", 
            length(ko_eggnog_only), 
            (length(ko_eggnog_only)/total_unique_kos)*100))

# ============================================================================
# 4. CHECK HOUSEKEEPING GENES AS CONTROL
# ============================================================================

cat("\n--- CHECKING HOUSEKEEPING GENES (CONTROL) ---\n")

# Define core housekeeping genes that should be present
housekeeping_kos <- list(
  # DNA replication
  DNA_replication = c("K02314", "K02315", "K02316", "K02337"),
  # Ribosomal proteins
  Ribosomal = c("K02945", "K02946", "K02948", "K02950", "K02952"),
  # ATP synthase
  ATP_synthase = c("K02111", "K02112", "K02113", "K02114", "K02115"),
  # Translation factors
  Translation = c("K02355", "K02358", "K02992", "K02994")
)

# Check presence in each annotation
housekeeping_results <- data.frame(
  Category = character(),
  KO = character(),
  GhostKOALA = logical(),
  eggNOG = logical(),
  Both = logical(),
  stringsAsFactors = FALSE
)

for (category in names(housekeeping_kos)) {
  for (ko in housekeeping_kos[[category]]) {
    in_ghost <- ko %in% ghostkoala_kos
    in_egg <- ko %in% eggnog_kos
    housekeeping_results <- rbind(
      housekeeping_results,
      data.frame(
        Category = category,
        KO = ko,
        GhostKOALA = in_ghost,
        eggNOG = in_egg,
        Both = in_ghost & in_egg
      )
    )
  }
}

# Print housekeeping results
cat("\nHousekeeping Gene Detection:\n")
print(table(housekeeping_results$Category, housekeeping_results$Both))

cat("\nMissing housekeeping genes:\n")
missing <- housekeeping_results %>% filter(!Both)
if (nrow(missing) > 0) {
  print(missing)
} else {
  cat("All checked housekeeping genes found in both annotations!\n")
}

# ============================================================================
# 5. CREATE VISUALIZATIONS
# ============================================================================

cat("\nGenerating visualizations...\n")

# Create output directories if they don't exist
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

# Venn Diagram
venn.plot <- venn.diagram(
  x = list(
    GhostKOALA = ghostkoala_kos,
    eggNOG = eggnog_kos
  ),
  category.names = c("GhostKOALA", "eggNOG"),
  filename = NULL,
  fill = c("#440154ff", "#21908dff"),
  alpha = 0.5,
  cex = 2,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.dist = c(0.05, 0.05),
  cat.pos = c(-20, 20)
)

# Save Venn diagram
png("results/figures/hs_ko_venn_diagram.png", width = 800, height = 800)
grid.draw(venn.plot)
dev.off()

# Bar plot comparison
comparison_data <- data.frame(
  Category = c("Total KOs", "Overlap", "GhostKOALA only", "eggNOG only"),
  Count = c(
    total_unique_kos,
    length(ko_overlap),
    length(ko_ghostkoala_only),
    length(ko_eggnog_only)
  )
)

p1 <- ggplot(comparison_data, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "KO Assignment Comparison: GhostKOALA vs eggNOG",
    subtitle = "Halobacterium salinarum (HS)",
    x = "",
    y = "Number of KOs"
  )

ggsave("results/figures/hs_ko_comparison_barplot.png", p1, width = 10, height = 6, dpi = 300)

# ============================================================================
# 6. SAVE RESULTS TO TABLES
# ============================================================================

cat("\nSaving results to tables...\n")

# Summary statistics
summary_stats <- data.frame(
  Metric = c(
    "GhostKOALA unique KOs",
    "eggNOG unique KOs",
    "Total unique KOs (union)",
    "Overlapping KOs",
    "Overlap percentage",
    "GhostKOALA only KOs",
    "eggNOG only KOs",
    "GhostKOALA proteins annotated",
    "eggNOG proteins annotated"
  ),
  Value = c(
    length(ghostkoala_kos),
    length(eggnog_kos),
    total_unique_kos,
    length(ko_overlap),
    round(overlap_percent, 2),
    length(ko_ghostkoala_only),
    length(ko_eggnog_only),
    nrow(ghostkoala_clean),
    nrow(eggnog_clean)
  )
)

write.csv(summary_stats, 
          "results/tables/hs_annotation_comparison_summary.csv", 
          row.names = FALSE)

# Save KO lists
write.csv(data.frame(KO = sort(ko_overlap)), 
          "results/tables/hs_ko_overlap.csv", 
          row.names = FALSE)

write.csv(data.frame(KO = sort(ko_ghostkoala_only)), 
          "results/tables/hs_ko_ghostkoala_only.csv", 
          row.names = FALSE)

write.csv(data.frame(KO = sort(ko_eggnog_only)), 
          "results/tables/hs_ko_eggnog_only.csv", 
          row.names = FALSE)

# Save housekeeping results
write.csv(housekeeping_results, 
          "results/tables/hs_housekeeping_genes.csv", 
          row.names = FALSE)

# ============================================================================
# 7. FINAL SUMMARY
# ============================================================================

cat("\n============================================================================\n")
cat("Done\n")
cat("============================================================================\n")
cat("\nOutput files created:\n")
cat("  Figures:\n")
cat("    - results/figures/hs_ko_venn_diagram.png\n")
cat("    - results/figures/hs_ko_comparison_barplot.png\n")
cat("  Tables:\n")
cat("    - results/tables/hs_annotation_comparison_summary.csv\n")
cat("    - results/tables/hs_ko_overlap.csv\n")
cat("    - results/tables/hs_ko_ghostkoala_only.csv\n")
cat("    - results/tables/hs_ko_eggnog_only.csv\n")
cat("    - results/tables/hs_housekeeping_genes.csv\n")
cat("\n============================================================================\n")