library(tidyverse)

# ----------------------------
# CONFIGURATION
# ----------------------------
ko_col <- "KO"   # <<< change this ONLY if your KO column has a different name

# ----------------------------
# Load GhostKOALA files
# ----------------------------
hs_df <- read_csv(
  "HS_GhostKOALA.csv",
  comment = "#",
  show_col_types = FALSE
)

hm_df <- read_csv(
  "HM_GhostKOALA.csv",
  comment = "#",
  show_col_types = FALSE
)

ha_df <- read_csv(
  "HA_GhostKOALA.csv",
  comment = "#",
  show_col_types = FALSE
)

# ----------------------------
# Sanity check: KO column
# ----------------------------
check_ko_col <- function(df, name) {
  if (!ko_col %in% colnames(df)) {
    stop(
      paste0(
        "ERROR: Column '", ko_col,
        "' not found in ", name,
        ". Available columns: ",
        paste(colnames(df), collapse = ", ")
      )
    )
  }
}

check_ko_col(hs_df, "HS_GhostKOALA.csv")
check_ko_col(hm_df, "HM_GhostKOALA.csv")
check_ko_col(ha_df, "HA_GhostKOALA.csv")

# ----------------------------
# Helper: count annotated proteins
# ----------------------------
count_annotated <- function(df) {
  sum(df[[ko_col]] != "-" & !is.na(df[[ko_col]]))
}

# ----------------------------
# Coverage summary
# ----------------------------
coverage_data <- tibble(
  Strain = factor(c("HS", "HM", "HA"), levels = c("HS", "HM", "HA")),
  Total = c(nrow(hs_df), nrow(hm_df), nrow(ha_df)),
  Annotated = c(
    count_annotated(hs_df),
    count_annotated(hm_df),
    count_annotated(ha_df)
  )
) %>%
  mutate(
    Unannotated = Total - Annotated,
    Pct = round(100 * Annotated / Total, 1)
  )

# ----------------------------
# Long format for plotting
# ----------------------------
coverage_long <- coverage_data %>%
  pivot_longer(
    cols = c(Annotated, Unannotated),
    names_to = "Status",
    values_to = "Count"
  ) %>%
  mutate(
    Status = factor(Status, levels = c("Unannotated", "Annotated"))
  )

# ----------------------------
# Plot
# ----------------------------
p_coverage <- ggplot(
  coverage_long,
  aes(x = Strain, y = Count, fill = Status)
) +
  geom_col(width = 0.7, color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold",
    size = 5
  ) +
  geom_text(
    data = coverage_data,
    aes(
      x = Strain,
      y = Total * 1.05,
      label = paste0(Pct, "%")
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Unannotated" = "#E57373",
      "Annotated" = "#2E7D32"
    ),
    name = ""
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "GhostKOALA Annotation Coverage",
    x = "Strain",
    y = "Number of Proteins"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

# ----------------------------
# Save figure
# ----------------------------
ggsave(
  filename = "ghostkoala_annotation_coverage.tiff",
  plot = p_coverage,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  device = "tiff"
)

cat("✓ Annotation coverage figure saved\n")
