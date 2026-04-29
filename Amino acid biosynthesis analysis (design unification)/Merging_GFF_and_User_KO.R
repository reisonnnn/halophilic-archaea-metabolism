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

proteins <- read.csv("gff_parsed.csv")
proteins_new <-na.omit(proteins)

kegg <- read_csv("user_ko.csv", col_names = FALSE)
colnames(kegg) <- c("protein_id", "KO")

kegg[kegg == ""] <- NA
kegg_new <- na.omit(kegg)

combined_data <- merge(proteins, kegg_new, by = "protein_id", all = FALSE)
colnames(combined_data)[colnames(combined_data) == "locus_tag"] <- "gene"
write.csv(combined_data, "merged_with_ko.csv", row.names = FALSE)

# Load KO-to-pathway mapping
ko_to_pathway <- read.table("ko.txt", sep = "\t", header = FALSE, col.names = c("KO", "Pathway"))
# View the data to ensure it's in the correct format
head(ko_to_pathway)

# Load pathway descriptions
pathway_descriptions <- read.table("pathway.txt", sep = "\t", header = FALSE, col.names = c("Pathway", "Description"))
# View the data to ensure it's in the correct format
head(pathway_descriptions)

# Merge annotated genes with KO-to-pathway mapping
annotated_genes_new <- merge(combined_data, ko_to_pathway, by = "KO")

annotated_genes_new <- merge(annotated_genes_new, pathway_descriptions, by = "Pathway")

write.csv(annotated_genes_new, file="annotated_pathway_genes.csv", row.names=FALSE)