# ___ _  _ ____ ___  ____ ____ _  _ _ ____ _  _ ____
#  |  |  | |__/ |__] |  | |___ |\ | | | __ |\/| |__|
#  |  |__| |  \ |__] |__| |___ | \| | |__] |  | |  |
#
# Input of single proteome and output of UNIPROT IDs, protein classification and
# characteristics for 2D-Gel support
#
library(tidyverse)
library(seqinr)

# Read in fasta file
df <- read.fasta('./data/test_proteome.fasta', seqtype = 'AA')

# Clean up UNIPROT Accession IDs and return a clean list
names <-
  attributes(df)$name %>% lapply(., function (x)
    str_extract(x, '\\|([^\\|]+)\\|')) %>%
  gsub('\\|', '', .)

# Create a df with the molecular weight and pI for each uniprot identifier
protein.df <- tibble(names)

# Construct molecular weight column
protein.df['mw'] <- as.numeric(0)
for (i in 1:length(df)) {
  protein.df$mw[i] <- pmw(df[[i]])
}

# Construct pI column
protein.df['pI'] <- as.numeric(0)
for (i in 1:length(df)) {
  protein.df$pI[i] <- computePI(df[[i]])
}

# Construct protein scatter plot of molecular weight vs pI
proteome.plot <- ggplot(protein.df, aes(x = mw, y = pI, fill = mw)) +
  geom_point(
    shape = 21,
    size = 4,
    alpha = 0.75,
    fill = 'darkorange'
  ) +
  theme_minimal() +
  scale_x_log10() +
  labs(title = 'Distribution of Individual Proteins Across Proteome',
       x = 'Molecular Weight (Da)',
       y = 'Isoelectric Point')
proteome.plot
