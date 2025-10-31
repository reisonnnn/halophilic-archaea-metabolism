# ==============================================================================
# Organize Your Downloaded Genome Files
# Works with YOUR folder structure
# ==============================================================================

library(fs)

project_root <- getwd()

cat("=========================================\n")
cat("ORGANIZING YOUR GENOME FILES\n")
cat("=========================================\n\n")

# Your strain information
strains <- data.frame(
  strain_id = c("KBTZ05", "KBTZ06", "KBTZ03"),
  genus = c("Halomicrobium", "Haloarcula", "Halobacterium"),
  assembly = c("GCA_041530035.1", "GCA_048595605.1", "GCA_037094595.1"),
  stringsAsFactors = FALSE
)

# Create organized directories
dir_create(file.path(project_root, "data/raw/genomes/KBTZ03"))
dir_create(file.path(project_root, "data/raw/genomes/KBTZ05"))
dir_create(file.path(project_root, "data/raw/genomes/KBTZ06"))
dir_create(file.path(project_root, "data/raw/eggnog"))

# Function to copy and rename files
organize_files <- function(strain_id, assembly_prefix) {
  
  cat("\n-----------------------------------\n")
  cat("Processing:", strain_id, "\n")
  cat("-----------------------------------\n")
  
  # Source: your download structure
  source_dir <- file.path(project_root, strain_id, "ncbi_dataset", "data", assembly_prefix)
  
  # Destination: organized structure
  dest_dir <- file.path(project_root, "data/raw/genomes", strain_id)
  
  if (!dir.exists(source_dir)) {
    cat("✗ Source not found:", source_dir, "\n")
    return(FALSE)
  }
  
  # Find and copy files
  files <- list.files(source_dir, full.names = TRUE)
  
  # Genome FASTA (ends with _genomic.fna or .fna)
  genome_file <- files[grepl("genomic\\.fna$", files)][1]
  if (!is.na(genome_file)) {
    dest <- file.path(dest_dir, paste0(strain_id, "_genome.fasta"))
    file.copy(genome_file, dest, overwrite = TRUE)
    size_mb <- file.info(dest)$size / (1024^2)
    cat(sprintf("✓ Genome: %.2f MB\n", size_mb))
  } else {
    cat("✗ Genome file not found\n")
  }
  
  # Protein FASTA
  protein_file <- files[grepl("protein\\.faa$", files)][1]
  if (!is.na(protein_file)) {
    # Copy to genomes folder
    dest <- file.path(dest_dir, paste0(strain_id, "_proteins.faa"))
    file.copy(protein_file, dest, overwrite = TRUE)
    
    # Copy to eggnog folder (for annotation)
    eggnog_dest <- file.path(project_root, "data/raw/eggnog", paste0(strain_id, "_proteins.faa"))
    file.copy(protein_file, eggnog_dest, overwrite = TRUE)
    
    # Count sequences
    lines <- readLines(dest, warn = FALSE)
    n_seqs <- sum(grepl("^>", lines))
    size_mb <- file.info(dest)$size / (1024^2)
    
    cat(sprintf("✓ Proteins: %d sequences (%.2f MB)\n", n_seqs, size_mb))
  } else {
    cat("✗ Protein file not found\n")
  }
  
  # GFF annotation
  gff_file <- files[grepl("genomic\\.gff$", files)][1]
  if (!is.na(gff_file)) {
    dest <- file.path(dest_dir, paste0(strain_id, "_annotation.gff"))
    file.copy(gff_file, dest, overwrite = TRUE)
    size_mb <- file.info(dest)$size / (1024^2)
    cat(sprintf("✓ GFF: %.2f MB\n", size_mb))
  } else {
    cat("✗ GFF file not found\n")
  }
  
  return(TRUE)
}

# Process all three strains
cat("Organizing files from your downloads...\n")

# KBTZ05 - Halomicrobium
organize_files("KBTZ05", "GCF_041530035.1")

# KBTZ06 - Haloarcula  
organize_files("KBTZ06", "GCA_048595605.1")

# KBTZ03 - Halobacterium
organize_files("KBTZ03", "GCA_037094595.1")

# Final summary
cat("\n=========================================\n")
cat("ORGANIZATION COMPLETE!\n")
cat("=========================================\n\n")

# Check organized structure
cat("Your organized files:\n")
cat("--------------------\n\n")

for (strain in c("KBTZ03", "KBTZ05", "KBTZ06")) {
  strain_dir <- file.path(project_root, "data/raw/genomes", strain)
  
  cat(strain, ":\n")
  
  if (dir.exists(strain_dir)) {
    files <- list.files(strain_dir)
    if (length(files) > 0) {
      for (f in files) {
        fpath <- file.path(strain_dir, f)
        size_mb <- file.info(fpath)$size / (1024^2)
        cat(sprintf("  ✓ %s (%.2f MB)\n", f, size_mb))
      }
    } else {
      cat("  ✗ No files\n")
    }
  }
  cat("\n")
}

# Check eggNOG folder
cat("Files ready for eggNOG-mapper:\n")
cat("------------------------------\n")
eggnog_dir <- file.path(project_root, "data/raw/eggnog")

if (dir.exists(eggnog_dir)) {
  eggnog_files <- list.files(eggnog_dir, pattern = "\\.faa$", full.names = TRUE)
  
  if (length(eggnog_files) > 0) {
    for (f in eggnog_files) {
      lines <- readLines(f, warn = FALSE)
      n_seqs <- sum(grepl("^>", lines))
      cat(sprintf("✓ %s: %d proteins\n", basename(f), n_seqs))
    }
    
    cat("\n✓ ALL READY FOR EGGNOG-MAPPER!\n")
  } else {
    cat("✗ No protein files found\n")
  }
}

cat("\n=========================================\n")
cat("Next: Run eggNOG-mapper annotation\n")
cat("=========================================\n")

