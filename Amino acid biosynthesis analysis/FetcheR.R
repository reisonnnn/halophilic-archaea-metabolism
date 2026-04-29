# Amino Acid Biosynthesis Analysis with KEGG Module Logic Parsing
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, readr, stringr, ggplot2, KEGGREST)

# Load and validate data
if(!file.exists("annotated_pathway_genes.csv")) stop("Input file not found!")
data <- read_csv("annotated_pathway_genes.csv", show_col_types = FALSE)
if(nrow(data) == 0) stop("Input file is empty!")
if(!"KO" %in% colnames(data)) stop("Input file must contain 'KO' column!")
data <- data %>% mutate(KO = str_trim(KO)) %>% filter(!is.na(KO) & KO != "")
present_KOs <- unique(data$KO)
cat(sprintf("Found %d unique KO numbers\n\n", length(present_KOs)))

# Define amino acid modules
aa_modules <- list(
  Serine = list(modules = c("M00020"), names = "Serine biosynthesis"),
  Threonine = list(modules = c("M00018"), names = "Threonine biosynthesis"),
  Cysteine = list(
    modules = c("M00021", "M00338", "M00609"),
    names = c("Cysteine biosynthesis, serine pathway", "Cysteine biosynthesis, homocysteine pathway", "Cysteine biosynthesis, methionine pathway")
  ),
  Methionine = list(modules = c("M00017"), names = "Methionine biosynthesis"),
  Valine = list(modules = c("M00019"), names = "Valine biosynthesis"),
  Isoleucine = list(
    modules = c("M00019", "M00535", "M00570"),
    names = c("Isoleucine biosynthesis, shared with valine", "Isoleucine biosynthesis, pyruvate pathway", "Isoleucine biosynthesis, threonine pathway")
  ),
  Leucine = list(modules = c("M00432"), names = "Leucine biosynthesis"),
  Lysine = list(
    modules = c("M00016", "M00525", "M00526", "M00527", "M00030", "M00433", "M00031"),
    names = c("Lysine biosynthesis, succinyl-DAP pathway", "Lysine biosynthesis, acetyl-DAP pathway", "Lysine biosynthesis, DAP dehydrogenase pathway", "Lysine biosynthesis, DAP aminotransferase pathway", "Lysine biosynthesis, AAA pathway", "Lysine biosynthesis", "Lysine biosynthesis, mediated by LysW")
  ),
  Ornithine = list(
    modules = c("M00028", "M00763"),
    names = c("Ornithine biosynthesis, glutamate pathway", "Ornithine biosynthesis, mediated by LysW")
  ),
  Arginine = list(
    modules = c("M00844", "M00845"),
    names = c("Arginine biosynthesis, ornithine pathway", "Arginine biosynthesis, acetylcitrulline pathway")
  ),
  Proline = list(modules = c("M00015"), names = "Proline biosynthesis"),
  Histidine = list(modules = c("M00026"), names = "Histidine biosynthesis"),
  Tryptophan = list(modules = c("M00023"), names = "Tryptophan biosynthesis"),
  Phenylalanine = list(
    modules = c("M00024", "M00910"),
    names = c("Phenylalanine biosynthesis, phenylpyruvate pathway", "Phenylalanine biosynthesis, arogenate pathway")
  ),
  Tyrosine = list(
    modules = c("M00025", "M00040"),
    names = c("Tyrosine biosynthesis, HPP pathway", "Tyrosine biosynthesis, arogenate pathway")
  ),
  Alanine = list(modules = c(), names = "Alanine biosynthesis (from pyruvate)", key_enzymes = c("K00812", "K00813", "K14260")),
  Aspartate = list(modules = c(), names = "Aspartate biosynthesis (from oxaloacetate)", key_enzymes = c("K00812", "K00813", "K10206", "K00278")),
  Asparagine = list(modules = c(), names = "Asparagine biosynthesis (from aspartate)", key_enzymes = c("K01914", "K01915", "K01953")),
  Glutamate = list(modules = c(), names = "Glutamate biosynthesis (from 2-oxoglutarate)", key_enzymes = c("K00260", "K00261", "K00262", "K15371", "K00265", "K00266")),
  Glutamine = list(modules = c(), names = "Glutamine biosynthesis (from glutamate)", key_enzymes = c("K01915", "K01914")),
  Glycine = list(modules = c(), names = "Glycine biosynthesis (from serine or other precursors)", key_enzymes = c("K00600", "K00281", "K00282", "K00283", "K00605"))
)

# Function to evaluate KEGG module completeness
evaluate_module <- function(module_id, present_kos) {
  tryCatch({
    module_info <- keggGet(module_id)[[1]]
    definition <- module_info$DEFINITION
    
    if(is.null(definition)) {
      return(list(is_complete = FALSE, all_kos = character(0), present_kos = character(0), 
                  missing_kos = character(0), satisfied_steps = 0, total_steps = 0))
    }
    
    # Extract all KOs
    all_kos <- unique(str_extract_all(definition, "K\\d{5}")[[1]])
    present <- all_kos[all_kos %in% present_kos]
    missing <- all_kos[!all_kos %in% present_kos]
    
    # Parse steps - remove optional components
    def_clean <- gsub("-K\\d{5}", "", definition)
    def_clean <- gsub("--", "", def_clean)
    steps <- str_split(def_clean, "\\s+")[[1]]
    steps <- steps[steps != ""]
    
    total_steps <- length(steps)
    satisfied_steps <- 0
    
    # Evaluate each step (considering alternatives within each step)
    for(step in steps) {
      step_kos <- str_extract_all(step, "K\\d{5}")[[1]]
      if(length(step_kos) == 0) next
      # If ANY alternative enzyme is present, step is satisfied
      if(any(step_kos %in% present_kos)) satisfied_steps <- satisfied_steps + 1
    }
    
    is_complete <- (satisfied_steps == total_steps) && (total_steps > 0)
    
    return(list(
      is_complete = is_complete,
      all_kos = all_kos,
      present_kos = present,
      missing_kos = missing,
      satisfied_steps = satisfied_steps,
      total_steps = total_steps
    ))
  }, error = function(e) {
    return(list(is_complete = FALSE, all_kos = character(0), present_kos = character(0),
                missing_kos = character(0), satisfied_steps = 0, total_steps = 0))
  })
}

# Fetch and evaluate all modules
cat("Fetching KEGG module definitions...\n")
all_modules <- unique(unlist(lapply(aa_modules, function(x) x$modules)))
all_modules <- all_modules[all_modules != ""]

module_results <- list()
for(m in all_modules) {
  cat(sprintf("  %s...\n", m))
  module_results[[m]] <- evaluate_module(m, present_KOs)
  Sys.sleep(0.5)
}

# Evaluate each amino acid
results <- lapply(names(aa_modules), function(aa) {
  modules <- aa_modules[[aa]]$modules
  module_names <- aa_modules[[aa]]$names
  key_enzymes <- aa_modules[[aa]]$key_enzymes
  
  # Case 1: Amino acids with key enzymes (no modules)
  if(length(modules) == 0 && !is.null(key_enzymes)) {
    present <- key_enzymes[key_enzymes %in% present_KOs]
    missing <- key_enzymes[!key_enzymes %in% present_KOs]
    is_complete <- length(present) > 0  # ANY ONE enzyme is sufficient
    
    status <- ifelse(is_complete, "Complete", "No Pathway")
    
    return(tibble(
      Amino_Acid = aa,
      Module_ID = "Key enzyme",
      Pathway_Name = module_names,
      Total_Steps = 1,
      Satisfied_Steps = ifelse(is_complete, 1, 0),
      Total_KOs = length(key_enzymes),
      Present_KOs_count = length(present),
      Status = status,
      Present_KOs = paste(present, collapse = ", "),
      Missing_KOs = paste(missing, collapse = ", ")
    ))
  }
  
  # Case 2: No modules and no key enzymes
  if(length(modules) == 0) {
    return(tibble(
      Amino_Acid = aa,
      Module_ID = "None",
      Pathway_Name = module_names,
      Total_Steps = 0,
      Satisfied_Steps = 0,
      Total_KOs = 0,
      Present_KOs_count = 0,
      Status = "No Pathway",
      Present_KOs = "",
      Missing_KOs = ""
    ))
  }
  
  # Case 3: Regular KEGG modules
  lapply(seq_along(modules), function(i) {
    m <- modules[i]
    eval_result <- module_results[[m]]
    
    # Determine status based on presence of enzymes
    status <- case_when(
      eval_result$is_complete ~ "Complete",
      eval_result$satisfied_steps == 0 ~ "No Pathway",
      TRUE ~ "Incomplete"
    )
    
    tibble(
      Amino_Acid = aa,
      Module_ID = m,
      Pathway_Name = module_names[i],
      Total_Steps = eval_result$total_steps,
      Satisfied_Steps = eval_result$satisfied_steps,
      Total_KOs = length(eval_result$all_kos),
      Present_KOs_count = length(eval_result$present_kos),
      Status = status,
      Present_KOs = paste(eval_result$present_kos, collapse = ", "),
      Missing_KOs = paste(eval_result$missing_kos, collapse = ", ")
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# Summarize by amino acid (best pathway for each)
aa_summary <- results %>%
  group_by(Amino_Acid) %>%
  # First calculate completeness for each pathway
  mutate(
    Pathway_Completeness = ifelse(Total_Steps > 0, (Satisfied_Steps / Total_Steps) * 100, 0)
  ) %>%
  # Then select the best pathway (highest completeness)
  slice_max(order_by = Pathway_Completeness, n = 1, with_ties = FALSE) %>%
  summarise(
    Best_Module = first(Module_ID),
    Best_Pathway = first(Pathway_Name),
    Alternative_Pathways = n_distinct(Module_ID),
    Satisfied_Steps = first(Satisfied_Steps),
    Total_Steps = first(Total_Steps),
    Best_Steps = paste(first(Satisfied_Steps), "/", first(Total_Steps), sep = ""),
    Status = first(Status),
    Completeness_Pct = first(Pathway_Completeness),
    .groups = "drop"
  ) %>%
  arrange(Amino_Acid)

# Display results
cat("\n", rep("=", 80), "\n", sep = "")
cat("AMINO ACID BIOSYNTHESIS ANALYSIS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("SUMMARY (20 Amino Acids):\n")
print(aa_summary %>% select(Amino_Acid, Status, Completeness_Pct, Best_Steps, Best_Pathway), n = 20)

cat("\n\nSTATUS COUNTS:\n")
status_counts <- aa_summary %>% count(Status)
for(i in 1:nrow(status_counts)) {
  cat(sprintf("  %s: %d amino acids\n", status_counts$Status[i], status_counts$n[i]))
}

# Prototrophy assessment
complete_count <- sum(aa_summary$Status == "Complete")
if(complete_count == 20) {
  cat("\n✓ PROTOTROPH: Can synthesize all 20 amino acids!\n")
} else {
  essential <- aa_summary %>% filter(Status != "Complete") %>% pull(Amino_Acid)
  cat(sprintf("\n✗ AUXOTROPH: Requires %d amino acids from environment:\n", length(essential)))
  cat(paste("  -", essential, collapse = "\n"), "\n")
}

# Save detailed results
write_csv(aa_summary, "amino_acid_summary.csv")
write_csv(results, "amino_acid_module_details.csv")

# ============================================================================
# VISUALIZATIONS
# ============================================================================

# Prepare data for plotting - sort by completeness percentage
plot_data <- aa_summary %>%
  arrange(Completeness_Pct) %>%
  mutate(
    # Create display label with pathway info
    AA_Label = case_when(
      Amino_Acid == "Arginine" ~ paste0(Amino_Acid, " (from Ornithine)"),
      Amino_Acid == "Ornithine" ~ paste0(Amino_Acid, " (precursor)"),
      TRUE ~ Amino_Acid
    ),
    AA_Label = factor(AA_Label, levels = AA_Label),
    Status = factor(Status, levels = c("Complete", "Incomplete", "No Pathway"))
  )

# Plot 1: Bar plot sorted by completeness percentage
p1 <- ggplot(plot_data, 
             aes(x = AA_Label, y = Completeness_Pct, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "#27ae60", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(Completeness_Pct, 0), "%")),
            hjust = -0.1, size = 3, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Complete" = "#27ae60",
                               "Incomplete" = "#f39c12",
                               "No Pathway" = "white"),
                    drop = FALSE) +
  labs(title = "Amino Acid Biosynthesis Capability",
       subtitle = "Ranked by completeness | Green = complete; Orange = incomplete; Blank = no pathway",
       x = "Amino Acid", y = "Pathway Completeness (%)", fill = "Status") +
  scale_y_continuous(limits = c(0, 110),
                     breaks = seq(0, 100, 25),
                     labels = paste0(seq(0, 100, 25), "%"),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.y = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("amino_acid_biosynthesis_barplot.tiff", p1, width = 10, height = 8, dpi = 300, compression = "lzw")


# Display plots
print(p1)

cat("\n", rep("=", 80), "\n", sep = "")
cat("✅ Files saved:\n")
cat("   - amino_acid_summary.csv\n")
cat("   - amino_acid_module_details.csv\n")
cat("   - amino_acid_biosynthesis_barplot.png\n")
cat("   - amino_acid_biosynthesis_barplot.tiff\n")
cat(rep("=", 80), "\n", sep = "")