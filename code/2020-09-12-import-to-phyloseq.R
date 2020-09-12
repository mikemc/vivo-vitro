# Goal: Create a phyloseq object with a simplified feature set for analyzing
# bias
# - Import the abundance tables for A1, S1, and S2
# - Filter to a reduced set of features that correspond to just the inoculum
# and positive control strains 
# - Combine into a single phyloseq object, with sample data
# - Add taxonomy and tree

# Prep ------------------------------------------------------------------------

library(here)
library(tidyverse)
library(speedyseq)

# Strain data
strain_data <- here("output/strain-data", "strain-pheno-geno-tax.Rds") %>%
  readRDS %>%
  mutate(species_name = ncbi_species_display_name)

# Sample data for all 4 sequencing centers.
sam <- here("output", "sample-data", "sample-data.Rds") %>% readRDS

# GTDB tree + metadata
tree <- ape::read.tree(here("output", "strain-data",
    "gtdb-representatives.tree"))
tree$tip.label <- strain_data %>%
  {.$species_name[match(tree$tip.label, .$gtdb_genome_representative)]}

# GTDB taxonomy
rnks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
tax <- strain_data %>%
  select(species_name, gtdb_taxonomy) %>%
  separate(gtdb_taxonomy, into = rnks, sep = ";") %>%
  mutate(across(all_of(rnks), str_replace, "^[a-z]__", ""))

# Get initial phyloseq object -------------------------------------------------

# Will use the taxonomy table slot temporarily, to simplify merging with the
# shotgun data

## A1
# Get a map from ASV sequence to species name, for A1 inoculum and T. mobilis
# ASVs; note, excludes the fecal E. coli ASV
map.a1 <- here("output", "a1", "a1-asv-strain-assignments.csv") %>%
  read_csv %>%
  filter(str_detect(source, "inoculum|mobilis")) %>%
  transmute(
    feature_id = sequence, 
    species_name = str_glue("{genus} {species}") %>% as.character
  )
# Abundance table, subset to just the focal strains
abun.a1 <- here("output", "a1", "dada2", "seqtab-nochim.Rds") %>%
  readRDS %>%
  otu_table(taxa_are_rows = FALSE)
ps.a1 <- phyloseq(abun.a1, tax_table(map.a1))
# tax_table(ps.a1)

## S1 and S2
# Map from feature id to species name
map.s <- strain_data %>%
  select(feature_id = reference_genome, species_name)
# Shotgun data, subset to just the focal strains
abun.s <- here("output", "shotgun", "2020-07-09-count-matrix.Rds") %>%
  readRDS %>%
  otu_table(taxa_are_rows = FALSE)
ps.s <- phyloseq(abun.s, tax_table(map.s))
# Note: Includes host reads, but not Unmapped reads
# tax_table(ps.s)

# Combine tables and merge features by species name
ps <- merge_phyloseq(ps.a1, ps.s) %>% tax_glom("species_name")
taxa_names(ps) <- tax_table(ps)[,"species_name" %>% c]

# Create final phyloseq object ------------------------------------------------

# Add strain group to tax table facilitate downstream filtering
tax0 <- tax %>%
  left_join(strain_data %>% select(species_name, strain_group), 
    by = "species_name") %>%
  tax_table

# Merge with tree + replace tax
# Note: Adding tree removes the Host entry from the table
ps <- merge_phyloseq(ps, sample_data(sam), tree)
tax_table(ps) <- tax0

saveRDS(ps, here("output/phyloseq", "2020-09-12-phyloseq-focal-strains.Rds"))

# To filter to just the inoculum strains, run:
# ps %>% subset_taxa(strain_group == "Inoculum")
# To filter to just the inoculum strains that actually colonized the mice, run:
# ps %>% subset_taxa(strain_group == "Inoculum" & genus != "Faecalibacterium")
