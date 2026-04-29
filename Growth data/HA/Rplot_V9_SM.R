library(tidyverse)
library(rstatix)
library(ggpubr)
library(car)
library(gridExtra)
library(readxl)

# File input for Excel files
file_path <- file.choose()
file_name <- basename(file_path)
output_dir <- tools::file_path_sans_ext(file_name)
if (!dir.exists(output_dir)) dir.create(output_dir)

# Read data from Excel
growth <- tryCatch(
  {
    read_excel(file_path, col_types = "guess")
  },
  error = function(e) {
    # Try alternative reading method if read_excel fails
    tryCatch(
      read.xlsx(file_path, check.names = FALSE),
      error = function(e2) stop("Error reading Excel file: ", e$message, "\nAlso tried alternative method: ", e2$message)
    )
  }
)

# Prepare data
names(growth)[1] <- "Time"
growth$Time <- as.numeric(gsub("\\D", "", growth$Time))

# Transform to long format with proper grouping
growth_long <- growth %>%
  pivot_longer(-Time, names_to = "Condition", values_to = "OD660") %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  group_by(Condition, Time) %>%
  mutate(
    Sample_ID = row_number(),
    BioRep = ceiling(Sample_ID / 3),
    TechRep = ((Sample_ID - 1) %% 3) + 1
  ) %>%
  ungroup()

# Average technical replicates
growth_long_avg <- growth_long %>%
  group_by(Time, Condition, BioRep) %>%
  summarise(OD660_avg = mean(OD660, na.rm = TRUE), .groups = "drop")

growth_long_avg <- growth_long_avg %>%
  mutate(
    Strain = factor(sub(" .*", "", Condition)),
    Media = factor(sub("^\\S+\\s+", "", Condition))
  )

# Get unique levels in order of appearance
strain_levels <- unique(growth_long_avg$Strain)
media_levels <- unique(growth_long_avg$Media)

# Summary statistics for all timepoints
growth_summary <- growth_long_avg %>%
  group_by(Time, Condition) %>%
  summarise(
    mean_OD = mean(OD660_avg),
    sd_OD = sd(OD660_avg),
    .groups = "drop"
  ) %>%
  mutate(
    Time = factor(Time),
    Strain = factor(sub(" .*", "", Condition), levels = strain_levels),
    Media = factor(sub("^\\S+\\s+", "", Condition), levels = media_levels)
  )

# Plot: All timepoints
p_all <- ggplot(growth_summary, aes(x = Media, y = mean_OD, fill = Time)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(mean_OD - sd_OD, 0), ymax = mean_OD + sd_OD),
    width = 0.2, position = position_dodge(0.8)
  ) +
  facet_wrap(~ Strain, scales = "free_x") +
  scale_fill_grey(start = 0.7, end = 0.2) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Archaea Growth (Mean ± SD)",
    y = "OD660", x = NULL, fill = "Time"
  )

print(p_all)

# ==================== STATISTICAL ANALYSIS (72h only) ====================

# Filter 72h data for statistics
data_72 <- growth_long_avg %>%
  filter(Time == 72)

# Check if we have multiple levels for analysis
n_strains <- n_distinct(data_72$Strain)
n_media <- n_distinct(data_72$Media)

cat(sprintf("\n=== DATA CHECK (72h) ===\n"))
cat(sprintf("Number of strains: %d\n", n_strains))
cat(sprintf("Number of media: %d\n", n_media))

if (n_strains < 1 || n_media < 2) {
  stop("Insufficient data for analysis. Need at least 1 strain and 2 media conditions.")
}

# Assumption testing
cat("\n=== NORMALITY TEST (72h) ===\n")

# Build appropriate model based on data structure
if (n_strains > 1) {
  aov_model <- aov(OD660_avg ~ Strain * Media, data = data_72)
} else {
  aov_model <- aov(OD660_avg ~ Media, data = data_72)
}

shapiro_result <- shapiro.test(residuals(aov_model))
print(shapiro_result)

is_normal <- shapiro_result$p.value >= 0.05

cat(sprintf("\nData distribution: %s (p = %.4f)\n",
            ifelse(is_normal, "NORMAL", "NON-NORMAL"),
            shapiro_result$p.value))

# QQ plot
p_qq <- ggqqplot(residuals(aov_model), title = "QQ Plot - 72h Model Residuals") +
  theme_minimal()
print(p_qq)

# Levene's test
cat("\n=== HOMOGENEITY OF VARIANCE (72h) ===\n")
levene_result <- tryCatch(
  leveneTest(OD660_avg ~ Condition, data = data_72),
  error = function(e) {
    warning("Levene's test failed: ", e$message)
    NULL
  }
)
print(levene_result)

# Choose statistical test based on normality
cat(sprintf("\n=== USING %s TESTS ===\n",
            ifelse(is_normal, "PARAMETRIC", "NON-PARAMETRIC")))

if (is_normal) {
  # Parametric: Welch ANOVA + Games-Howell
  test_results <- data_72 %>%
    group_by(Strain) %>%
    filter(n() >= 3, n_distinct(Condition) > 1) %>%
    group_modify(~ {
      test <- oneway.test(OD660_avg ~ Condition, data = ., var.equal = FALSE)
      data.frame(test_name = "Welch ANOVA", test_p = test$p.value)
    })
  
  posthoc_results <- data_72 %>%
    group_by(Strain) %>%
    filter(n_distinct(Condition) > 1, n() >= 3) %>%
    group_modify(~ tryCatch(
      games_howell_test(., OD660_avg ~ Condition),
      error = function(e) data.frame()
    ))
  
} else {
  # Non-parametric: Kruskal-Wallis + Dunn's test
  test_results <- data_72 %>%
    group_by(Strain) %>%
    filter(n() >= 3, n_distinct(Condition) > 1) %>%
    group_modify(~ {
      test <- kruskal.test(OD660_avg ~ Condition, data = .)
      data.frame(test_name = "Kruskal-Wallis", test_p = test$p.value)
    })
  
  posthoc_results <- data_72 %>%
    group_by(Strain) %>%
    filter(n_distinct(Condition) > 1, n() >= 3) %>%
    group_modify(~ tryCatch(
      dunn_test(., OD660_avg ~ Condition, p.adjust.method = "holm"),
      error = function(e) data.frame()
    ))
}

cat("\n=== MAIN TEST RESULTS ===\n")
print(test_results)

cat("\n=== POST-HOC RESULTS ===\n")
print(posthoc_results)

# Save results
write.csv(test_results, file.path(output_dir, "main_test_results.csv"), row.names = FALSE)

if (nrow(posthoc_results) > 0 && "p.adj" %in% names(posthoc_results)) {
  posthoc_results <- posthoc_results %>%
    mutate(
      Significance = case_when(
        p.adj <= 0.001 ~ "***",
        p.adj <= 0.01 ~ "**",
        p.adj <= 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  write.csv(posthoc_results, file.path(output_dir, "posthoc_results.csv"), row.names = FALSE)
}

# ==================== VISUALIZATION WITH 0h AND 72h SIDE-BY-SIDE ====================

# Always create the 0h/72h comparison plot, even without significance
# Prepare data with both 0h and 72h
data_0h_72h <- growth_long_avg %>%
  filter(Time %in% c(0, 72)) %>%
  mutate(
    Strain = factor(sub(" .*", "", Condition), levels = strain_levels),
    Media = factor(sub("^\\S+\\s+", "", Condition), levels = media_levels),
    Time_label = paste0(Time, "h")
  )

# Summary for plotting
data_summary <- data_0h_72h %>%
  group_by(Time_label, Strain, Condition, Media) %>%
  summarise(
    mean_OD = mean(OD660_avg, na.rm = TRUE),
    sd_OD = sd(OD660_avg, na.rm = TRUE),
    .groups = "drop"
  )

# Check if we have significant results to add brackets
has_significant_results <- (exists("posthoc_results") && 
                              nrow(posthoc_results) > 0 && 
                              "p.adj" %in% names(posthoc_results) && 
                              any(!is.na(posthoc_results$p.adj)) &&
                              any(posthoc_results$p.adj <= 0.05, na.rm = TRUE))

if (has_significant_results) {
  # Prepare significance table (72h only)
  res_table <- posthoc_results %>%
    mutate(
      P = p.adj,
      signif = case_when(
        P <= 0.001 ~ "***",
        P <= 0.01 ~ "**",
        P <= 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      label = ifelse(P < 0.001,
                     paste0("p=", formatC(P, format = "e", digits = 2)),
                     paste0("p=", round(P, 3))),
      group1_simple = factor(sub("^\\S+\\s+", "", group1), levels = media_levels),
      group2_simple = factor(sub("^\\S+\\s+", "", group2), levels = media_levels)
    ) %>%
    filter(P <= 0.05) %>%
    group_by(Strain) %>%
    mutate(y.position = max(data_0h_72h$OD660_avg, na.rm = TRUE) + row_number() * 0.15) %>%
    ungroup()
  
  # Base plot function
  base_plot <- function(show_xaxis = TRUE) {
    list(
      geom_col(position = position_dodge(0.8), width = 0.7, alpha = 1),
      geom_errorbar(
        aes(ymin = pmax(mean_OD - sd_OD, 0), ymax = mean_OD + sd_OD),
        width = 0.25, linewidth = 0.5, position = position_dodge(0.8)
      ),
      facet_wrap(~ Strain, scales = "free_x"),
      scale_fill_manual(values = c("0h" = "grey80", "72h" = "grey30")),
      theme_minimal(base_size = 14),
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = margin(10, 10, 5, 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5)
      ),
      labs(
        title = "Comparison of media formulations (0h vs 72h, Mean ± SD)",
        subtitle = "n=3 biological replicates; significance shown for 72h comparisons only",
        y = "OD660", x = NULL, fill = "Timepoint"
      )
    )
  }
  
  # Helper function to add significance brackets for 72h data
  add_significance <- function(p, res_tbl, label_type = "label") {
    # Adjust y.position for grouped bars
    res_tbl_adjusted <- res_tbl %>%
      mutate(
        xmin = as.numeric(group1_simple) + 0.2,  # Offset for 72h bar
        xmax = as.numeric(group2_simple) + 0.2
      )
    
    p + stat_pvalue_manual(
      res_tbl_adjusted, 
      label = label_type, 
      tip.length = 0.01, 
      y.position = "y.position",
      bracket.size = 0.5, 
      xmin = "xmin", 
      xmax = "xmax",
      hide.ns = TRUE, 
      label.size = if(label_type == "label") 2.5 else 3
    )
  }
  
  # Plot with p-values
  p_bar <- ggplot(data_summary, aes(x = Media, y = mean_OD, fill = Time_label)) +
    base_plot(TRUE)
  p_bar <- add_significance(p_bar, res_table, "label")
  
  # Plot with asterisks
  p_bar_ast <- ggplot(data_summary, aes(x = Media, y = mean_OD, fill = Time_label)) +
    base_plot(FALSE)
  p_bar_ast <- add_significance(p_bar_ast, res_table, "signif")
  
} else {
  # Create plots without significance brackets
  cat("\nNo significant results found for 72h comparisons.\n")
  
  # Base plot function (without significance addition)
  base_plot <- list(
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 1),
    geom_errorbar(
      aes(ymin = pmax(mean_OD - sd_OD, 0), ymax = mean_OD + sd_OD),
      width = 0.25, linewidth = 0.5, position = position_dodge(0.8)
    ),
    facet_wrap(~ Strain, scales = "free_x"),
    scale_fill_manual(values = c("0h" = "grey80", "72h" = "grey30")),
    theme_minimal(base_size = 14),
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(10, 10, 5, 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold", size = 14),
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5)
    ),
    labs(
      title = "Comparison of media formulations (0h vs 72h, Mean ± SD)",
      subtitle = "n=3 biological replicates; no significant differences at 72h (p > 0.05)",
      y = "OD660", x = NULL, fill = "Timepoint"
    )
  )
  
  # Create plots without significance
  p_bar <- ggplot(data_summary, aes(x = Media, y = mean_OD, fill = Time_label)) +
    base_plot
  
  p_bar_ast <- ggplot(data_summary, aes(x = Media, y = mean_OD, fill = Time_label)) +
    base_plot
}

# Media composition table
table_grob <- tableGrob(
  data.frame(
    Component = c("Salts solution", "Amino acids", "Phosphate", "Trace elements",
                  "Vitamins", "Glycerol", "Glucose", "Yeast extract", "Peptone"),
    `YEP` = c("+", "-", "-", "-", "-", "+", "-", "+", "+"),
    `YEP +TE` = c("+", "-", "-", "+", "-", "+", "-", "+", "+"),
    `SM` = c("+", "+", "+", "+", "+", "+", "+", "-", "-"),
    check.names = FALSE
  ),
  rows = NULL,
  theme = ttheme_minimal(
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, fontsize = 11),
      bg_params = list(fill = c("white", "grey95"))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 11),
      bg_params = list(fill = "grey90")
    )
  )
)

# Combined plot
combined_plot <- grid.arrange(
  p_bar_ast, table_grob,
  ncol = 1, heights = c(5, 2.5)
)

# Save all plots
ggsave(file.path(output_dir, "p_all.tiff"), p_all,
       width = 2200, height = 1800, dpi = 300, units = "px")
ggsave(file.path(output_dir, "p_qq.tiff"), p_qq,
       width = 1600, height = 1600, dpi = 300, units = "px")
ggsave(file.path(output_dir, "p_72h_bar_pvalue.tiff"), p_bar,
       width = 2400, height = 1800, dpi = 300, units = "px")
ggsave(file.path(output_dir, "p_72h_bar_asterisk_with_table.tiff"), combined_plot,
       width = 2400, height = 2600, dpi = 300, units = "px")

cat("\nAll plots saved successfully.\n")

cat("\n=== Analysis complete! ===\n")
cat(sprintf("Results saved to: %s\n", output_dir))