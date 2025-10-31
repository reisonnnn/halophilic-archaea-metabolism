# ==============================================================================
# ARCHAEAL GLUCOSE & GLYCEROL METABOLISM - COMPREHENSIVE PATHWAY ANALYSIS
# Accounts for alternative archaeal pathways and pathway interconnections
# ==============================================================================

library(tidyverse)
library(pheatmap)

project_root <- getwd()

cat("=========================================\n")
cat("ARCHAEAL GLUCOSE METABOLISM ANALYSIS\n")
cat("=========================================\n\n")

# Load data
all_annotations <- read.csv(
  file.path(project_root, "data/processed/metabolic_pathways/all_eggnog_annotations.csv"),
  stringsAsFactors = FALSE
)

# CRITICAL FIX: Extract KO numbers properly
# KO numbers in the data have "ko:" prefix and may be comma-separated
# Example: "ko:K00844,ko:K01810" or "ko:K03413,ko02020" (mixed with pathways)
# We need to extract ONLY K##### format (K followed by exactly 5 digits)

cat("Preprocessing KO numbers...\n")

# Function to extract clean KO numbers from a string
extract_ko_numbers <- function(ko_string) {
  if(is.na(ko_string) || ko_string == "") {
    return(NA)
  }
  # Extract all matches of K followed by exactly 5 digits
  matches <- str_extract_all(ko_string, "K\\d{5}")[[1]]
  if(length(matches) == 0) {
    return(NA)
  }
  return(matches)
}

# Apply to all rows
all_annotations$KO_extracted <- lapply(all_annotations$KEGG_ko, extract_ko_numbers)

# Create long format: one row per KO per gene
all_annotations_long <- all_annotations %>%
  unnest(KO_extracted) %>%
  rename(KO_clean = KO_extracted) %>%
  filter(!is.na(KO_clean))

cat(sprintf("  Original rows: %d\n", nrow(all_annotations)))
cat(sprintf("  Expanded rows (one per KO): %d\n", nrow(all_annotations_long)))
cat(sprintf("  Unique KO numbers found: %d\n", 
            length(unique(all_annotations_long$KO_clean))))

# Show breakdown by strain
cat("\nKO numbers per strain:\n")
for(strain in c("hs", "hm", "ha")) {
  n_kos <- all_annotations_long %>%
    filter(strain_id == strain) %>%
    pull(KO_clean) %>%
    unique() %>%
    length()
  cat(sprintf("  %s: %d unique KO numbers\n", strain, n_kos))
}

# Sample of extracted KOs
cat("\nSample extracted KO numbers:\n")
sample_kos <- head(unique(all_annotations_long$KO_clean), 30)
cat(paste(sample_kos, collapse = ", "))
cat("\n\n")

# Verify key enzymes are found
test_enzymes <- c("K00844", "K00864", "K00134", "K01803", "K00927")
cat("Testing for key glycolysis/glycerol enzymes:\n")
for(ko in test_enzymes) {
  found <- all_annotations_long %>%
    filter(KO_clean == ko) %>%
    group_by(strain_id) %>%
    summarise(n = n(), .groups = 'drop')
  
  if(nrow(found) > 0) {
    cat(sprintf("  %s: Found in %d strains - ", ko, nrow(found)))
    cat(paste(found$strain_id, collapse = ", "), "\n")
  } else {
    cat(sprintf("  %s: NOT FOUND\n", ko))
  }
}
cat("\n")

# ==============================================================================
# DEFINE ALL POSSIBLE GLUCOSE METABOLISM PATHWAYS
# ==============================================================================

# 1. CLASSICAL EMP GLYCOLYSIS (Embden-Meyerhof-Parnas)
emp_classical <- data.frame(
  KO = c("K00844", "K01810", "K00850", "K01623", "K01624", 
         "K00134", "K00150", "K00927", "K01834", "K15633", "K15634",
         "K00873", "K12406"),
  Enzyme = c("HK - Hexokinase", "GPI - Glucose-6-P isomerase", 
             "PFK - Phosphofructokinase", "FBA class I", "FBA class II",
             "GAPDH (NAD+)", "GAPDH (NADP+)", "PGK - Phosphoglycerate kinase", 
             "PGAM", "ENO", "ENO (alternative)",
             "PK - Pyruvate kinase", "PK (alternative)"),
  Step = c("1", "2", "3", "4", "4", "6", "6", "7", "8", "9", "9", "10", "10"),
  Type = "Classical_EMP",
  stringsAsFactors = FALSE
)

# 2. ARCHAEAL MODIFIED EMP (uses different enzymes!)
emp_archaeal <- data.frame(
  KO = c("K00845", "K01810", "K00918", "K01622", "K01623", "K01624",
         "K00134", "K10705", "K00927", "K01689", "K15633", "K15634",
         "K00873", "K12406"),
  Enzyme = c("Glucokinase (archaeal)", "GPI", "ADP-dependent PFK (archaeal!)", 
             "FBA (archaeal)", "FBA class I", "FBA class II",
             "GAPDH", "GAPOR (archaeal)", "PGK", 
             "PGAM (archaeal 2,3-BPG dependent)", "ENO", "ENO (alternative)",
             "PK", "PK (alternative)"),
  Step = c("1", "2", "3", "4", "4", "4", "6", "6", "7", "8", "9", "9", "10", "10"),
  Type = "Archaeal_EMP",
  stringsAsFactors = FALSE
)

# 3. ENTNER-DOUDOROFF PATHWAY (common in archaea)
ed_pathway <- data.frame(
  KO = c("K00844", "K00845", "K01810", "K00036", "K01057", "K01625",
         "K00134", "K00150", "K00927", "K01834", "K15633",
         "K00873", "K12406"),
  Enzyme = c("Hexokinase", "Glucokinase", "GPI", 
             "G6P dehydrogenase", "6-phosphogluconolactonase",
             "KDG-6-P aldolase (ED key enzyme!)",
             "GAPDH", "GAPDH (NADP+)", "PGK", "PGAM", "ENO",
             "PK", "PK (alternative)"),
  Step = c("1", "1", "2", "3", "4", "5", "6", "6", "7", "8", "9", "10", "10"),
  Type = "ED_pathway",
  stringsAsFactors = FALSE
)

# 4. SEMI-PHOSPHORYLATED ED PATHWAY (archaeal variant)
ed_semi <- data.frame(
  KO = c("K00036", "K01057", "K01690", "K00865", "K00134", "K00927",
         "K01834", "K15633", "K00873", "K12406"),
  Enzyme = c("G6P dehydrogenase", "6-phosphogluconolactonase",
             "KDG aldolase (non-phosphorylated!)", "KDGK - KDG kinase",
             "GAPDH", "PGK", "PGAM", "ENO", "PK", "PK (alternative)"),
  Step = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "9"),
  Type = "Semi_ED",
  stringsAsFactors = FALSE
)

# 5. BRANCHED ED PATHWAY
ed_branched <- data.frame(
  KO = c("K00036", "K01057", "K00065", "K01625", "K01690",
         "K00134", "K00927", "K01834", "K15633", "K00873"),
  Enzyme = c("G6P dehydrogenase", "6-phosphogluconolactonase",
             "Gluconate dehydrogenase", "KDG-6-P aldolase",
             "KDG aldolase", "GAPDH", "PGK", "PGAM", "ENO", "PK"),
  Step = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
  Type = "Branched_ED",
  stringsAsFactors = FALSE
)

# 6. GLUCONEOGENESIS (reverse pathway - important!)
gluconeogenesis <- data.frame(
  KO = c("K01596", "K01610", "K01834", "K00927", "K00134", "K01624",
         "K00850", "K01810", "K00844"),
  Enzyme = c("PEP carboxykinase", "PEP carboxylase", "PGAM", "PGK", "GAPDH",
             "FBA", "PFK", "GPI", "Hexokinase"),
  Step = c("1", "1", "2", "3", "4", "5", "6", "7", "8"),
  Type = "Gluconeogenesis",
  stringsAsFactors = FALSE
)

# 7. GLYCEROL METABOLISM (INTERTWINED WITH GLYCOLYSIS!)
glycerol_pathway <- data.frame(
  KO = c("K00864", "K00111", "K00112", "K00865", "K00858", "K01803"),
  Enzyme = c("GlpK - Glycerol kinase (Glycerol -> Glycerol-3-P)",
             "GlpD - G3P dehydrogenase NAD+ (G3P -> DHAP)",
             "GlpD - G3P dehydrogenase (aerobic)",
             "GlpD (alternative)",
             "DHAK - Dihydroxyacetone kinase (DHA -> DHAP)",
             "TPI - Triose-phosphate isomerase (DHAP <-> G3P)"),
  Step = c("1", "2", "2", "2", "alt", "interchange"),
  Type = "Glycerol_to_DHAP",
  Note = c("Entry point", "Feeds into glycolysis at DHAP!", 
           "Feeds into glycolysis at DHAP!", "Feeds into glycolysis at DHAP!",
           "Alternative entry", "Connects DHAP and G3P (glycolysis step 5)"),
  stringsAsFactors = FALSE
)

# 8. LOWER GLYCOLYSIS (COMMON TO ALL - from G3P/DHAP to Pyruvate)
lower_glycolysis <- data.frame(
  KO = c("K01803", "K00134", "K00150", "K10705", "K00927", "K01834", 
         "K01689", "K15633", "K15634", "K00873", "K12406"),
  Enzyme = c("TPI - Triose-phosphate isomerase",
             "GAPDH (NAD+)", "GAPDH (NADP+)", "GAPOR (archaeal)",
             "PGK", "PGAM", "PGAM (archaeal)", 
             "ENO", "ENO (alternative)", "PK", "PK (alternative)"),
  Step = c("5", "6", "6", "6", "7", "8", "8", "9", "9", "10", "10"),
  Type = "Lower_Glycolysis",
  Note = "Common pathway from DHAP/G3P to Pyruvate",
  stringsAsFactors = FALSE
)

# ==============================================================================
# FUNCTION TO CHECK PATHWAY PRESENCE
# ==============================================================================

check_pathway <- function(annotations_long, pathway_df, pathway_name, strain_name) {
  
  # Get all KO numbers for this strain (use the cleaned KO_clean column)
  strain_kos <- annotations_long %>%
    filter(strain_id == strain_name) %>%
    pull(KO_clean) %>%
    unique()
  
  pathway_check <- pathway_df %>%
    mutate(
      Present = KO %in% strain_kos,
      Strain = strain_name
    )
  
  # Calculate presence by step (handling alternative enzymes)
  unique_steps <- unique(pathway_check$Step)
  steps_covered <- 0
  
  for(step in unique_steps) {
    step_enzymes <- pathway_check %>% filter(Step == step)
    if(any(step_enzymes$Present)) {
      steps_covered <- steps_covered + 1
    }
  }
  
  n_present <- sum(pathway_check$Present)
  n_total <- nrow(pathway_df)
  n_steps <- length(unique_steps)
  
  completeness <- (steps_covered / n_steps) * 100
  
  return(list(
    details = pathway_check,
    n_present = n_present,
    n_total = n_total,
    steps_covered = steps_covered,
    n_steps = n_steps,
    completeness = completeness
  ))
}

# ==============================================================================
# ANALYZE ALL PATHWAYS FOR ALL STRAINS
# ==============================================================================

pathways_list <- list(
  "Classical EMP\n(Glucose → Pyruvate)" = emp_classical,
  "Archaeal EMP\n(Modified EMP pathway)" = emp_archaeal,
  "Entner-Doudoroff\n(Glucose → Pyruvate, alternative)" = ed_pathway,
  "Semi-phosphorylated ED\n(Archaeal ED variant)" = ed_semi,
  "Branched ED\n(ED pathway variant)" = ed_branched,
  "Gluconeogenesis\n(Pyruvate → Glucose)" = gluconeogenesis,
  "Glycerol metabolism\n(Glycerol → DHAP)" = glycerol_pathway,
  "Lower glycolysis\n(G3P/DHAP → Pyruvate)" = lower_glycolysis
)

strains <- c("hs", "hm", "ha")
strain_labels <- c("HS (Halobacterium)", "HM (Halomicrobium)", "HA (Haloarcula)")

results_summary <- data.frame()
all_details <- list()

for(pathway_name in names(pathways_list)) {
  
  cat("\n")
  cat("==================================================\n")
  cat(toupper(pathway_name), "\n")
  cat("==================================================\n\n")
  
  pathway_df <- pathways_list[[pathway_name]]
  
  for(i in 1:length(strains)) {
    
    result <- check_pathway(all_annotations_long, pathway_df, pathway_name, strains[i])
    
    cat(sprintf("%s:\n", strain_labels[i]))
    cat(sprintf("  Enzymes present: %d/%d\n", result$n_present, result$n_total))
    cat(sprintf("  Steps covered: %d/%d\n", result$steps_covered, result$n_steps))
    cat(sprintf("  Completeness: %.1f%%\n", result$completeness))
    
    # Show present enzymes
    present_enzymes <- result$details %>% filter(Present)
    if(nrow(present_enzymes) > 0) {
      cat("  Present enzymes:\n")
      for(j in 1:nrow(present_enzymes)) {
        cat(sprintf("    %s: %s\n", present_enzymes$KO[j], present_enzymes$Enzyme[j]))
      }
    }
    
    # Show missing enzymes
    missing_enzymes <- result$details %>% filter(!Present)
    if(nrow(missing_enzymes) > 0) {
      cat("  Missing enzymes:\n")
      for(j in 1:nrow(missing_enzymes)) {
        cat(sprintf("    %s: %s\n", missing_enzymes$KO[j], missing_enzymes$Enzyme[j]))
      }
    }
    cat("\n")
    
    # Store summary
    results_summary <- rbind(results_summary, data.frame(
      Strain = strains[i],
      Strain_label = strain_labels[i],
      Pathway = pathway_name,
      Enzymes_present = result$n_present,
      Total_enzymes = result$n_total,
      Steps_covered = result$steps_covered,
      Total_steps = result$n_steps,
      Completeness = result$completeness,
      stringsAsFactors = FALSE
    ))
    
    # Store details
    result$details$Pathway <- pathway_name
    all_details[[paste(pathway_name, strains[i], sep = "_")]] <- result$details
  }
}

# ==============================================================================
# CREATE COMPREHENSIVE VISUALIZATION
# ==============================================================================

cat("\n")
cat("==================================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("==================================================\n\n")

# 1. Completeness heatmap for all pathways
completeness_matrix <- results_summary %>%
  select(Pathway, Strain, Completeness) %>%
  pivot_wider(names_from = Strain, values_from = Completeness) %>%
  as.data.frame()

rownames(completeness_matrix) <- completeness_matrix$Pathway
completeness_matrix <- completeness_matrix %>% select(-Pathway)
completeness_matrix <- as.matrix(completeness_matrix)
completeness_matrix <- completeness_matrix[, c("hs", "hm", "ha")]

# Order pathways logically
pathway_order <- c("Classical EMP\n(Glucose → Pyruvate)", 
                   "Archaeal EMP\n(Modified EMP pathway)", 
                   "Entner-Doudoroff\n(Glucose → Pyruvate, alternative)", 
                   "Semi-phosphorylated ED\n(Archaeal ED variant)", 
                   "Branched ED\n(ED pathway variant)", 
                   "Gluconeogenesis\n(Pyruvate → Glucose)", 
                   "Glycerol metabolism\n(Glycerol → DHAP)", 
                   "Lower glycolysis\n(G3P/DHAP → Pyruvate)")
completeness_matrix <- completeness_matrix[pathway_order, ]

png(file.path(project_root, "results/figures/glucose_pathways_completeness_heatmap.png"),
    width = 2400, height = 2000, res = 300)

pheatmap(completeness_matrix,
         color = colorRampPalette(c("#FFFFFF", "#FFF59D", "#FFD54F", "#FFA726", "#FF6F00", "#D84315", "#B71C1C"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 11,
         fontsize_col = 13,
         border_color = "grey60",
         cellwidth = 60,
         cellheight = 35,
         main = "Glucose Metabolism Pathway Completeness (%)",
         labels_col = c("HS (Halobacterium)", "HM (Halomicrobium)", "HA (Haloarcula)"),
         display_numbers = TRUE,
         number_format = "%.0f",
         fontsize_number = 10,
         number_color = "black",
         breaks = seq(0, 100, length.out = 101),
         legend_breaks = c(0, 25, 50, 75, 100),
         legend_labels = c("0%", "25%", "50%", "75%", "100%"))

dev.off()

cat("✓ Completeness heatmap saved\n")

# 2. Bar plot comparison
p_pathways <- ggplot(results_summary, 
                     aes(x = Pathway, y = Completeness, fill = Strain_label)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", Completeness)), 
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.8) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "orange", alpha = 0.3) +
  scale_fill_manual(values = c("HS (Halobacterium)" = "#FF6B6B", 
                               "HM (Halomicrobium)" = "#95E1D3", 
                               "HA (Haloarcula)" = "#A8E6CF")) +
  labs(title = "Glucose Metabolism Pathway Completeness",
       subtitle = "Multiple alternative pathways analyzed",
       x = "", y = "Pathway Completeness (%)", fill = "Strain") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom") +
  ylim(0, 110)

ggsave(file.path(project_root, "results/figures/glucose_pathways_completeness_barplot.png"),
       p_pathways, width = 12, height = 7, dpi = 300)

cat("✓ Bar plot saved\n")

# 3. Create detailed enzyme presence heatmap for key pathways
key_pathways <- c("Classical EMP\n(Glucose → Pyruvate)", 
                  "Archaeal EMP\n(Modified EMP pathway)", 
                  "Entner-Doudoroff\n(Glucose → Pyruvate, alternative)", 
                  "Glycerol metabolism\n(Glycerol → DHAP)")

for(pw in key_pathways) {
  
  pw_details <- bind_rows(
    all_details[[paste(pw, "hs", sep = "_")]],
    all_details[[paste(pw, "hm", sep = "_")]],
    all_details[[paste(pw, "ha", sep = "_")]]
  )
  
  enzyme_matrix <- pw_details %>%
    mutate(
      Enzyme_label = paste0(Step, ": ", KO, " - ", Enzyme),
      Present_num = as.numeric(Present)
    ) %>%
    select(Strain, Enzyme_label, Present_num) %>%
    pivot_wider(names_from = Strain, values_from = Present_num, values_fill = 0) %>%
    as.data.frame()
  
  rownames(enzyme_matrix) <- enzyme_matrix$Enzyme_label
  enzyme_matrix <- enzyme_matrix %>% select(-Enzyme_label)
  enzyme_matrix <- as.matrix(enzyme_matrix)
  enzyme_matrix <- enzyme_matrix[, c("hs", "hm", "ha")]
  
  # Check if there's variation in the data
  if(length(unique(as.vector(enzyme_matrix))) > 1) {
    
    png(file.path(project_root, 
                  paste0("results/figures/enzyme_presence_", 
                         tolower(gsub(" ", "_", pw)), ".png")),
        width = 2200, height = 800 + nrow(enzyme_matrix)*40, res = 300)
    
    pheatmap(enzyme_matrix,
             color = c("#FFEBEE", "#1B5E20"),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 9,
             fontsize_col = 12,
             border_color = "grey50",
             cellwidth = 55,
             cellheight = 25,
             main = paste0(pw, " - Enzyme Presence"),
             labels_col = c("HS (Halobacterium)", "HM (Halomicrobium)", "HA (Haloarcula)"),
             legend_breaks = c(0, 1),
             legend_labels = c("Absent", "Present"))
    
    dev.off()
    
    cat(sprintf("✓ %s enzyme presence heatmap saved\n", pw))
  } else {
    # All values are the same - create a simple visualization instead
    all_value <- unique(as.vector(enzyme_matrix))[1]
    status <- ifelse(all_value == 1, "ALL PRESENT", "ALL ABSENT")
    
    cat(sprintf("⚠ %s: %s - skipping heatmap (no variation)\n", pw, status))
  }
}

# ==============================================================================
# PATHWAY INTERCONNECTION ANALYSIS
# ==============================================================================

cat("\n")
cat("==================================================\n")
cat("PATHWAY INTERCONNECTIONS\n")
cat("==================================================\n\n")

# Check if glycerol pathway connects properly to lower glycolysis
cat("GLYCEROL -> GLYCOLYSIS CONNECTION:\n")
cat("Glycerol metabolism produces DHAP, which feeds into lower glycolysis\n")
cat("Key enzyme: K01803 (TPI - Triose-phosphate isomerase) converts DHAP <-> G3P\n\n")

for(i in 1:length(strains)) {
  strain_kos <- all_annotations_long %>%
    filter(strain_id == strains[i]) %>%
    pull(KO_clean) %>%
    unique()
  
  has_glycerol_kinase <- "K00864" %in% strain_kos
  has_g3p_dehydrogenase <- any(c("K00111", "K00112", "K00865") %in% strain_kos)
  has_tpi <- "K01803" %in% strain_kos
  has_gapdh <- any(c("K00134", "K00150", "K10705") %in% strain_kos)
  
  cat(sprintf("%s:\n", strain_labels[i]))
  cat(sprintf("  Glycerol kinase (K00864): %s\n", ifelse(has_glycerol_kinase, "PRESENT", "ABSENT")))
  cat(sprintf("  G3P dehydrogenase: %s\n", ifelse(has_g3p_dehydrogenase, "PRESENT", "ABSENT")))
  cat(sprintf("  TPI (DHAP <-> G3P): %s\n", ifelse(has_tpi, "PRESENT", "ABSENT")))
  cat(sprintf("  GAPDH (G3P -> pyruvate): %s\n", ifelse(has_gapdh, "PRESENT", "ABSENT")))
  
  if(has_glycerol_kinase && has_g3p_dehydrogenase && has_tpi && has_gapdh) {
    cat("  ✓ COMPLETE GLYCEROL -> PYRUVATE PATHWAY\n")
  } else {
    cat("  ✗ INCOMPLETE PATHWAY\n")
  }
  cat("\n")
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

write.csv(results_summary,
          file.path(project_root, "results/tables/glucose_pathways_completeness_summary.csv"),
          row.names = FALSE)

all_details_df <- bind_rows(all_details)
write.csv(all_details_df,
          file.path(project_root, "results/tables/glucose_pathways_enzyme_details.csv"),
          row.names = FALSE)

cat("\n")
cat("=========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=========================================\n\n")

cat("Files created:\n")
cat("1. results/figures/glucose_pathways_completeness_heatmap.png\n")
cat("2. results/figures/glucose_pathways_completeness_barplot.png\n")
cat("3. results/figures/enzyme_presence_classical_emp.png\n")
cat("4. results/figures/enzyme_presence_archaeal_emp.png\n")
cat("5. results/figures/enzyme_presence_entner-doudoroff.png\n")
cat("6. results/figures/enzyme_presence_glycerol_metabolism.png\n")
cat("7. results/tables/glucose_pathways_completeness_summary.csv\n")
cat("8. results/tables/glucose_pathways_enzyme_details.csv\n\n")

cat("=========================================\n")
cat("PATHWAY COMPLETENESS SUMMARY\n")
cat("=========================================\n\n")

# Show completeness by pathway
summary_table <- results_summary %>%
  group_by(Pathway) %>%
  summarise(
    HS = Completeness[Strain == "hs"],
    HM = Completeness[Strain == "hm"],
    HA = Completeness[Strain == "ha"],
    Average = mean(Completeness),
    .groups = 'drop'
  ) %>%
  arrange(desc(Average))

print(summary_table)
cat("\n")