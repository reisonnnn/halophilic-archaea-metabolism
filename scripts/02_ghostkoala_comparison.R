library(tidyverse)
library(httr)
library(jsonlite)

if (!require("eulerr", quietly = TRUE)) {
  install.packages("eulerr")
  library(eulerr)
}

project_root <- getwd()

dir.create(file.path(project_root, "results/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results/figures"), recursive = TRUE, showWarnings = FALSE)

load_ghostkoala <- function(filepath, strain_id) {
  data <- read.delim(filepath, header = FALSE, sep = "\t", 
                     col.names = c("protein_id", "KO"), 
                     stringsAsFactors = FALSE, fill = TRUE)
  data$KO[data$KO == ""] <- NA
  data$strain <- strain_id
  data <- data %>% filter(!is.na(KO))
  return(data)
}

fetch_kegg_info <- function(ko_list, batch_size = 10) {
  ko_list <- unique(ko_list)
  results <- data.frame(KO = ko_list, name = NA_character_, definition = NA_character_, 
                        stringsAsFactors = FALSE)
  
  cat("Fetching KEGG descriptions for", length(ko_list), "KOs...\n")
  
  for (i in seq(1, length(ko_list), by = batch_size)) {
    batch <- ko_list[i:min(i + batch_size - 1, length(ko_list))]
    
    for (ko in batch) {
      tryCatch({
        url <- paste0("https://rest.kegg.jp/get/", ko)
        response <- GET(url, timeout(10))
        
        if (status_code(response) == 200) {
          content <- content(response, "text", encoding = "UTF-8")
          lines <- strsplit(content, "\n")[[1]]
          
          name_line <- grep("^NAME", lines, value = TRUE)
          def_line <- grep("^DEFINITION", lines, value = TRUE)
          
          if (length(name_line) > 0) {
            results$name[results$KO == ko] <- gsub("^NAME\\s+", "", name_line[1])
          }
          if (length(def_line) > 0) {
            results$definition[results$KO == ko] <- gsub("^DEFINITION\\s+", "", def_line[1])
          }
        }
        Sys.sleep(0.1)
      }, error = function(e) {
        message("Error fetching ", ko, ": ", e$message)
      })
    }
    
    if (i %% 50 == 1) {
      cat("  Processed", min(i + batch_size - 1, length(ko_list)), "/", length(ko_list), "\n")
    }
  }
  
  return(results)
}

cat("Loading GhostKOALA annotations...\n")

hs <- load_ghostkoala(file.path(project_root, "annotation_outputs/GhostKOALA/HS_GhostKOALA.txt"), "HS")
hm <- load_ghostkoala(file.path(project_root, "annotation_outputs/GhostKOALA/HM_GhostKOALA.txt"), "HM")
ha <- load_ghostkoala(file.path(project_root, "annotation_outputs/GhostKOALA/HA_GhostKOALA.txt"), "HA")

cat("HS:", nrow(hs), "proteins with KO\n")
cat("HM:", nrow(hm), "proteins with KO\n")
cat("HA:", nrow(ha), "proteins with KO\n")

hs_kos <- unique(hs$KO)
hm_kos <- unique(hm$KO)
ha_kos <- unique(ha$KO)

core_kos <- Reduce(intersect, list(hs_kos, hm_kos, ha_kos))
hs_hm_only <- setdiff(intersect(hs_kos, hm_kos), ha_kos)
hs_ha_only <- setdiff(intersect(hs_kos, ha_kos), hm_kos)
hm_ha_only <- setdiff(intersect(hm_kos, ha_kos), hs_kos)
hs_unique <- setdiff(hs_kos, union(hm_kos, ha_kos))
hm_unique <- setdiff(hm_kos, union(hs_kos, ha_kos))
ha_unique <- setdiff(ha_kos, union(hs_kos, hm_kos))

cat("\nKO distribution:\n")
cat("  Core (all 3):", length(core_kos), "\n")
cat("  HS-HM shared:", length(hs_hm_only), "\n")
cat("  HS-HA shared:", length(hs_ha_only), "\n")
cat("  HM-HA shared:", length(hm_ha_only), "\n")
cat("  HS unique:", length(hs_unique), "\n")
cat("  HM unique:", length(hm_unique), "\n")
cat("  HA unique:", length(ha_unique), "\n")

all_kos <- unique(c(core_kos, hs_hm_only, hs_ha_only, hm_ha_only, hs_unique, hm_unique, ha_unique))
kegg_info <- fetch_kegg_info(all_kos)

make_table <- function(ko_list, kegg_info) {
  if (length(ko_list) == 0) {
    return(data.frame(KO = character(), name = character(), definition = character()))
  }
  df <- data.frame(KO = ko_list, stringsAsFactors = FALSE) %>%
    left_join(kegg_info, by = "KO") %>%
    arrange(KO)
  return(df)
}

cat("\nSaving tables...\n")

write.csv(make_table(core_kos, kegg_info),
          file.path(project_root, "results/tables/ghostkoala_core_all_three.csv"), row.names = FALSE)

write.csv(make_table(hs_hm_only, kegg_info),
          file.path(project_root, "results/tables/ghostkoala_shared_HS_HM.csv"), row.names = FALSE)

write.csv(make_table(hs_ha_only, kegg_info),
          file.path(project_root, "results/tables/ghostkoala_shared_HS_HA.csv"), row.names = FALSE)

write.csv(make_table(hm_ha_only, kegg_info),
          file.path(project_root, "results/tables/ghostkoala_shared_HM_HA.csv"), row.names = FALSE)

write.csv(make_table(hs_unique, kegg_info),
          file.path(project_root, "results/tables/ghostkoala_unique_HS.csv"), row.names = FALSE)

write.csv(make_table(hm_unique, kegg_info),
          file.path(project_root, "results/tables/ghostkoala_unique_HM.csv"), row.names = FALSE)

write.csv(make_table(ha_unique, kegg_info),
          file.path(project_root, "results/tables/ghostkoala_unique_HA.csv"), row.names = FALSE)

cat("Creating Venn diagram...\n")

euler_fit <- euler(list(
  "HS" = hs_kos,
  "HM" = hm_kos,
  "HA" = ha_kos
))

png(file.path(project_root, "results/figures/ghostkoala_venn_diagram.png"),
    width = 2400, height = 2400, res = 300)

plot(euler_fit,
     fills = list(fill = c("#E74C3C", "#3498DB", "#2ECC71"), alpha = 0.6),
     quantities = list(fontsize = 12, fontface = "bold"),
     labels = list(fontsize = 14, fontface = "bold"),
     edges = list(col = "white", lwd = 2),
     main = "GhostKOALA KO Overlap Between Strains")

dev.off()

summary_df <- data.frame(
  Category = c("Core (all 3)", "HS-HM shared", "HS-HA shared", "HM-HA shared",
               "HS unique", "HM unique", "HA unique"),
  KO_count = c(length(core_kos), length(hs_hm_only), length(hs_ha_only), 
               length(hm_ha_only), length(hs_unique), length(hm_unique), length(ha_unique))
)

write.csv(summary_df, file.path(project_root, "results/tables/ghostkoala_summary.csv"), row.names = FALSE)

cat("\nDone.\n")
cat("Tables saved in results/tables/\n")
cat("Venn diagram saved in results/figures/\n")