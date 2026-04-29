if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Function to check and install a library
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

install_if_missing("dplyr")
install_if_missing("readr")
install_if_missing("stringr")
install_if_missing("ggplot2")

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# Define file path
gff_file <- "genomic.gff"

# Read the GFF file (tab-delimited, skipping comment lines)
gff_data <- read_delim(
  gff_file,
  delim = "\t",
  comment = "#",
  col_names = c("seqid", "source", "type", "start", "end", 
                "score", "strand", "phase", "attributes")
)

# Extract only CDS features
cds_data <- gff_data %>%
  filter(type == "CDS") %>%
  select(attributes)

# Parse the 'attributes' column to extract fields
parsed_data <- cds_data %>%
  mutate(
    locus_tag = str_extract(attributes, "locus_tag=[^;]+") %>% str_remove("locus_tag="),
    protein_id = str_extract(attributes, "protein_id=[^;]+") %>% str_remove("protein_id="),
    product = str_extract(attributes, "product=[^;]+") %>% str_remove("product=")
  ) %>%
  select(locus_tag, protein_id, product)

# Print parsed locus_tags
print(parsed_data$locus_tag)

# Save the parsed data to a CSV file
write_csv(parsed_data, "gff_parsed.csv")

#Second part#
proteins <- read_csv("gff_parsed.csv")
proteins <- na.omit(proteins)

kegg <- read_csv("HM_GhostKOALA.csv", col_names = FALSE)
colnames(kegg) <- c("protein_id", "KO")
kegg[kegg == ""] <- NA
kegg_new <- na.omit(kegg)

combined_data <- merge(proteins, kegg_new, by = "protein_id", all = FALSE)
colnames(combined_data)[colnames(combined_data) == "locus_tag"] <- "gene"
write.csv(combined_data, "merged_with_ko.csv", row.names = FALSE)

# Load KO-to-pathway mapping
ko_to_pathway <- read.table("ko.txt", sep = "\t", header = FALSE, col.names = c("KO", "Pathway"))
head(ko_to_pathway)

# Load pathway descriptions
pathway_descriptions <- read.table("hmu.txt", sep = "\t", header = FALSE, col.names = c("Pathway", "Description"))
head(pathway_descriptions)

# Merge annotated genes with KO-to-pathway mapping
annotated_genes_new <- merge(combined_data, ko_to_pathway, by = "KO", all = FALSE)
annotated_genes_new <- merge(annotated_genes_new, pathway_descriptions, by = "Pathway", all = FALSE)
annotated_genes_new <- annotated_genes_new %>%
  mutate(across(everything(), ~str_replace_all(., "%", "-")))
write.csv(annotated_genes_new, file = "HM_annotated_pathway_genes.csv", row.names = FALSE)

# Count the number of genes per pathway description
pathway_counts <- annotated_genes_new %>%
  group_by(Description) %>%
  summarise(gene_count = n()) %>%
  arrange(desc(gene_count))

# Plot: Bar chart of gene counts per pathway
ggplot(pathway_counts, aes(x = reorder(Description, gene_count), y = gene_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Number of Genes per KEGG Pathway",
       x = "Pathway",
       y = "Gene Count") +
  theme_minimal()

ggsave("pathway_gene_counts.png", width = 14, height = 12)
