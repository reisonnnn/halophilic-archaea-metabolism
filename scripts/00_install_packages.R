

cat("============================================================================\n")
cat("Installing Required R Packages\n")
cat("============================================================================\n\n")

# List of required packages
required_packages <- c(
  "tidyverse",      # Data manipulation and visualization
  "readxl",         # Reading Excel files
  "VennDiagram",    # Creating Venn diagrams
  "ggplot2",        # Advanced plotting (part of tidyverse but listed for clarity)
  "gridExtra"       # Arranging multiple plots
)

# Function to check and install packages
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", package))
    install.packages(package, dependencies = TRUE, repos = "http://cran.rstudio.com/")
    if (require(package, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("✓ %s installed successfully\n\n", package))
      return(TRUE)
    } else {
      cat(sprintf("✗ Failed to install %s\n\n", package))
      return(FALSE)
    }
  } else {
    cat(sprintf("✓ %s already installed\n", package))
    return(TRUE)
  }
}

# Install packages
cat("Checking and installing packages...\n\n")
results <- sapply(required_packages, install_if_missing)

