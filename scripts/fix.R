# ==============================================================================
# FIX: GhostKOALA Venn Diagram
# Run this after setting working directory to your project root
# ==============================================================================

# Set working directory - EDIT THIS PATH if needed
# setwd("C:/Users/reison/Desktop/12.12/halophilic-archaea-metabolism")

library(eulerr)

# Define the counts from your analysis
# (These match your output: Core=755, HS-HM=22, HS-HA=41, HM-HA=133, HS=71, HM=67, HA=129)

# Method 1: Using named vector with intersection notation
venn_counts <- c(
  "HS" = 71,
  "HM" = 67,
  "HA" = 129,
  "HS&HM" = 22,
  "HS&HA" = 41,
  "HM&HA" = 133,
  "HS&HM&HA" = 755
)

# Create euler fit
euler_fit <- euler(venn_counts)

# Print to check
print(euler_fit)

# Create the plot object
venn_plot <- plot(
  euler_fit,
  fills = list(fill = c("#FF6B6B", "#4ECDC4", "#95E1D3"), alpha = 0.6),
  labels = list(fontsize = 14, fontface = "bold"),
  quantities = list(fontsize = 12, fontface = "bold"),
  edges = list(col = "white", lwd = 2),
  main = list(label = "GhostKOALA KO Distribution", fontsize = 16, fontface = "bold")
)

# Save to PNG - KEY: must use print() inside png device!
png("results/figures/ghostkoala_venn_diagram.png",
    width = 2400, height = 2400, res = 300)
print(venn_plot)  # THIS IS THE FIX - must explicitly print!
dev.off()

cat("✓ Venn diagram saved to results/figures/ghostkoala_venn_diagram.png\n")

# Also save a PDF version (often looks better)
pdf("results/figures/ghostkoala_venn_diagram.pdf",
    width = 8, height = 8)
print(venn_plot)
dev.off()

cat("✓ PDF version saved to results/figures/ghostkoala_venn_diagram.pdf\n")

# Display in RStudio viewer
print(venn_plot)

# ==============================================================================
# BONUS: Create a more detailed version with strain labels
# ==============================================================================

venn_plot_detailed <- plot(
  euler_fit,
  fills = list(fill = c("#FF6B6B", "#4ECDC4", "#95E1D3"), alpha = 0.6),
  labels = list(
    labels = c("HS\n(Halobacterium)", "HM\n(Halomicrobium)", "HA\n(Haloarcula)"),
    fontsize = 12, 
    fontface = "bold"
  ),
  quantities = list(fontsize = 11, fontface = "bold"),
  edges = list(col = "white", lwd = 2),
  main = list(
    label = "GhostKOALA KO Distribution Across Strains", 
    fontsize = 14, 
    fontface = "bold"
  )
)

png("results/figures/ghostkoala_venn_detailed.png",
    width = 2800, height = 2400, res = 300)
print(venn_plot_detailed)
dev.off()

cat("✓ Detailed version saved to results/figures/ghostkoala_venn_detailed.png\n")

# ==============================================================================
# Summary statistics
# ==============================================================================

cat("\n")
cat("===========================================\n")
cat("GhostKOALA KO Distribution Summary\n")
cat("===========================================\n")
cat(sprintf("Total unique KOs: %d\n", sum(venn_counts)))
cat(sprintf("Core (all 3 strains): %d (%.1f%%)\n", 
            venn_counts["HS&HM&HA"], 
            100 * venn_counts["HS&HM&HA"] / sum(venn_counts)))
cat(sprintf("HS unique: %d\n", venn_counts["HS"]))
cat(sprintf("HM unique: %d\n", venn_counts["HM"]))
cat(sprintf("HA unique: %d\n", venn_counts["HA"]))
cat(sprintf("HS-HM shared (not HA): %d\n", venn_counts["HS&HM"]))
cat(sprintf("HS-HA shared (not HM): %d\n", venn_counts["HS&HA"]))
cat(sprintf("HM-HA shared (not HS): %d\n", venn_counts["HM&HA"]))
cat("===========================================\n")