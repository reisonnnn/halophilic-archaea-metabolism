# ==============================================================================
# 02 - GhostKOALA Comparison Analysis
# ==============================================================================

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
}

library(tidyverse)
library(eulerr)

cat("=========================================\n")
cat("GHOSTKOALA COMPARISON ANALYSIS\n")
cat("=========================================\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

ghostkoala_dir <- "annotation_outputs/GhostKOALA"

load_ghostkoala <- function(filepath, strain_name) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  
  data <- read.delim(filepath, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("gene_id", "KO"))
  
  total_proteins <- nrow(data)
  
  data_annotated <- data %>% filter(!is.na(KO) & KO != "")
  data_annotated$strain <- strain_name
  
  cat(sprintf("  %s: %d/%d proteins with KO (%.1f%%)\n", 
              strain_name, nrow(data_annotated), total_proteins,
              100 * nrow(data_annotated) / total_proteins))
  
  return(list(data = data_annotated, 
              total = total_proteins, 
              annotated = nrow(data_annotated)))
}

cat("Loading GhostKOALA annotations...\n")
hs_result <- load_ghostkoala(file.path(ghostkoala_dir, "HS_GhostKOALA.txt"), "HS")
hm_result <- load_ghostkoala(file.path(ghostkoala_dir, "HM_GhostKOALA.txt"), "HM")
ha_result <- load_ghostkoala(file.path(ghostkoala_dir, "HA_GhostKOALA.txt"), "HA")

hs_data <- hs_result$data
hm_data <- hm_result$data
ha_data <- ha_result$data

all_data <- bind_rows(hs_data, hm_data, ha_data)
cat(sprintf("\nTotal KO assignments: %d\n\n", nrow(all_data)))

# ==============================================================================
# 2. EXTRACT UNIQUE KOs
# ==============================================================================

cat("Extracting unique KO numbers...\n")

hs_kos <- unique(hs_data$KO)
hm_kos <- unique(hm_data$KO)
ha_kos <- unique(ha_data$KO)

cat(sprintf("  HS: %d unique KOs\n", length(hs_kos)))
cat(sprintf("  HM: %d unique KOs\n", length(hm_kos)))
cat(sprintf("  HA: %d unique KOs\n", length(ha_kos)))

# ==============================================================================
# 3. CALCULATE OVERLAPS
# ==============================================================================

cat("\nCalculating overlaps...\n")

core_kos <- Reduce(intersect, list(hs_kos, hm_kos, ha_kos))
hs_hm_shared <- setdiff(intersect(hs_kos, hm_kos), core_kos)
hs_ha_shared <- setdiff(intersect(hs_kos, ha_kos), core_kos)
hm_ha_shared <- setdiff(intersect(hm_kos, ha_kos), core_kos)
hs_unique <- setdiff(hs_kos, union(hm_kos, ha_kos))
hm_unique <- setdiff(hm_kos, union(hs_kos, ha_kos))
ha_unique <- setdiff(ha_kos, union(hs_kos, hm_kos))

cat("\nKO distribution:\n")
cat(sprintf("  Core (all 3): %d\n", length(core_kos)))
cat(sprintf("  HS-HM shared: %d\n", length(hs_hm_shared)))
cat(sprintf("  HS-HA shared: %d\n", length(hs_ha_shared)))
cat(sprintf("  HM-HA shared: %d\n", length(hm_ha_shared)))
cat(sprintf("  HS unique: %d\n", length(hs_unique)))
cat(sprintf("  HM unique: %d\n", length(hm_unique)))
cat(sprintf("  HA unique: %d\n", length(ha_unique)))

# ==============================================================================
# 4. FETCH KEGG DESCRIPTIONS (WITH CACHING)
# ==============================================================================

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

cache_file <- "results/tables/kegg_ko_descriptions_cache.csv"
all_kos <- unique(c(hs_kos, hm_kos, ha_kos))

if (file.exists(cache_file)) {
  cat("\nLoading cached KEGG descriptions...\n")
  kegg_info <- read.csv(cache_file, stringsAsFactors = FALSE)
  
  # Check for new KOs not in cache
  missing_kos <- setdiff(all_kos, kegg_info$KO)
  
  if (length(missing_kos) > 0) {
    cat(sprintf("  Fetching %d new KOs not in cache...\n", length(missing_kos)))
  } else {
    cat(sprintf("  Cache complete: %d KOs loaded\n", nrow(kegg_info)))
  }
} else {
  cat("\nNo cache found. Fetching KEGG descriptions (this only happens once)...\n")
  missing_kos <- all_kos
  kegg_info <- data.frame()
}

# Fetch missing KOs
if (length(missing_kos) > 0) {
  cat(sprintf("Total to fetch: %d\n", length(missing_kos)))
  
  fetch_kegg_description <- function(ko_id) {
    url <- paste0("https://rest.kegg.jp/get/", ko_id)
    tryCatch({
      response <- readLines(url, warn = FALSE)
      name_line <- grep("^NAME", response, value = TRUE)
      def_line <- grep("^DEFINITION", response, value = TRUE)
      name <- ifelse(length(name_line) > 0, gsub("^NAME\\s+", "", name_line[1]), NA)
      definition <- ifelse(length(def_line) > 0, gsub("^DEFINITION\\s+", "", def_line[1]), NA)
      return(data.frame(KO = ko_id, name = name, definition = definition, stringsAsFactors = FALSE))
    }, error = function(e) {
      return(data.frame(KO = ko_id, name = NA, definition = NA, stringsAsFactors = FALSE))
    })
  }
  
  new_info <- data.frame()
  batch_size <- 50
  
  for (i in seq(1, length(missing_kos), by = batch_size)) {
    batch_end <- min(i + batch_size - 1, length(missing_kos))
    cat(sprintf("  Processing %d - %d of %d...\n", i, batch_end, length(missing_kos)))
    
    for (ko in missing_kos[i:batch_end]) {
      new_info <- bind_rows(new_info, fetch_kegg_description(ko))
      Sys.sleep(0.1)
    }
  }
  
  kegg_info <- bind_rows(kegg_info, new_info)
  write.csv(kegg_info, cache_file, row.names = FALSE)
  cat(sprintf("✓ Cache updated: %d KOs saved to %s\n", nrow(kegg_info), cache_file))
}

# ==============================================================================
# 5. CREATE OUTPUT TABLES
# ==============================================================================

cat("\nCreating output tables...\n")

create_annotated_table <- function(ko_list, category_name) {
  if (length(ko_list) == 0) {
    return(data.frame(KO = character(), name = character(), 
                      definition = character(), category = character()))
  }
  data.frame(KO = ko_list, stringsAsFactors = FALSE) %>%
    left_join(kegg_info, by = "KO") %>%
    mutate(category = category_name) %>%
    arrange(KO)
}

write.csv(create_annotated_table(core_kos, "Core"), 
          "results/tables/ghostkoala_core_all_three.csv", row.names = FALSE)
write.csv(create_annotated_table(hs_hm_shared, "HS-HM"), 
          "results/tables/ghostkoala_shared_HS_HM.csv", row.names = FALSE)
write.csv(create_annotated_table(hs_ha_shared, "HS-HA"), 
          "results/tables/ghostkoala_shared_HS_HA.csv", row.names = FALSE)
write.csv(create_annotated_table(hm_ha_shared, "HM-HA"), 
          "results/tables/ghostkoala_shared_HM_HA.csv", row.names = FALSE)
write.csv(create_annotated_table(hs_unique, "HS unique"), 
          "results/tables/ghostkoala_unique_HS.csv", row.names = FALSE)
write.csv(create_annotated_table(hm_unique, "HM unique"), 
          "results/tables/ghostkoala_unique_HM.csv", row.names = FALSE)
write.csv(create_annotated_table(ha_unique, "HA unique"), 
          "results/tables/ghostkoala_unique_HA.csv", row.names = FALSE)

cat("✓ Tables saved\n")

# ==============================================================================
# 6. VENN DIAGRAM
# ==============================================================================

cat("\nCreating Venn diagram...\n")

venn_counts <- c(
  "HS" = length(hs_unique),
  "HM" = length(hm_unique),
  "HA" = length(ha_unique),
  "HS&HM" = length(hs_hm_shared),
  "HS&HA" = length(hs_ha_shared),
  "HM&HA" = length(hm_ha_shared),
  "HS&HM&HA" = length(core_kos)
)

euler_fit <- euler(venn_counts)

venn_plot <- plot(
  euler_fit,
  fills = list(fill = c("#FF6B6B", "#4ECDC4", "#95E1D3"), alpha = 0.6),
  labels = list(fontsize = 14, fontface = "bold"),
  quantities = list(fontsize = 12, fontface = "bold"),
  edges = list(col = "white", lwd = 2)
)

png("results/figures/ghostkoala_venn_KO.png", width = 2400, height = 2400, res = 300)
print(venn_plot)
dev.off()

cat("✓ Venn diagram saved\n")

# ==============================================================================
# 7. ANNOTATION COVERAGE FIGURE
# ==============================================================================

cat("\nCreating annotation coverage figure...\n")

coverage_data <- data.frame(
  Strain = factor(c("HS", "HM", "HA"), levels = c("HS", "HM", "HA")),
  Total = c(hs_result$total, hm_result$total, ha_result$total),
  Annotated = c(hs_result$annotated, hm_result$annotated, ha_result$annotated)
) %>%
  mutate(
    Unannotated = Total - Annotated,
    Pct = round(100 * Annotated / Total, 1)
  )

coverage_long <- coverage_data %>%
  pivot_longer(cols = c(Annotated, Unannotated),
               names_to = "Status", values_to = "Count") %>%
  mutate(Status = factor(Status, levels = c("Unannotated", "Annotated")))

p_coverage <- ggplot(coverage_long, aes(x = Strain, y = Count, fill = Status)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.5) +
  geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 5) +
  geom_text(data = coverage_data,
            aes(x = Strain, y = Total + 150, label = paste0(Pct, "%")),
            inherit.aes = FALSE, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Unannotated" = "#E57373", "Annotated" = "#2E7D32"),
                    name = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "GhostKOALA Annotation Coverage",
       x = "Strain", y = "Number of Proteins") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

ggsave("results/figures/ghostkoala_annotation_coverage.png",
       p_coverage, width = 8, height = 6, dpi = 300)

cat("✓ Annotation coverage figure saved\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n=========================================\n")
cat("COMPLETE\n")
cat("=========================================\n\n")

cat("Tables:\n")
cat("  results/tables/ghostkoala_core_all_three.csv\n")
cat("  results/tables/ghostkoala_shared_*.csv\n")
cat("  results/tables/ghostkoala_unique_*.csv\n")
cat("  results/tables/kegg_ko_descriptions_cache.csv (reused on next run)\n\n")

cat("Figures:\n")
cat("  results/figures/ghostkoala_venn_KO.png\n")
cat("  results/figures/ghostkoala_annotation_coverage.png\n\n")

cat("Summary:\n")
cat(sprintf("  Total unique KOs: %d\n", sum(venn_counts)))
cat(sprintf("  Core (all 3): %d (%.1f%%)\n", length(core_kos), 
            100 * length(core_kos) / sum(venn_counts)))
cat(sprintf("  HS unique: %d | HM unique: %d | HA unique: %d\n",
            length(hs_unique), length(hm_unique), length(ha_unique)))
cat("=========================================\n")