library(tidyverse)
library(googlesheets4)
library(here)

# Data for the DNA samples ----------------------------------------------------

# One "dna sample" corresponds to one DNA extraction (before aliquoting onto
# plates for the four sequencing centers)

## Load from Google drive spreadsheet
sheet_url <- "https://docs.google.com/spreadsheets/d/1VYC4m4yk5K3AkXl2QOInaOuSla85ZkLHRDjAqFgabyA/"
# map from sample ids to specimen id and aliquot number
dna_sample_map <- read_sheet(sheet_url, sheet = 1, na = "NA",
  col_types = "cii") %>%
  rename(dna_sample_id = sample_id)
# information about individual specimens
specimen_data <- read_sheet(sheet_url, sheet = 2, na = "NA",
  col_types = "iciiDcciDii")
# aliquot_number to protocol used
aliquot_data <- read_sheet(sheet_url, sheet = 3, na = "NA",
  col_types = "ic")
# dna concentrations
qubit_results <- read_sheet(sheet_url, sheet = 4, na = "NA", 
  col_types = "cd") %>%
  rename(dna_sample_id = sample_id)
# which well on the plate + which plates the sample was aliquoted onto
plate_layout <- read_sheet(sheet_url, sheet = 6, na = "NA",
  col_types = "ccllll") %>%
  rename(dna_sample_id = sample_id)

## Save plain-text copies
write_csv(dna_sample_map, here("data/sample-data", "dna-sample-map.csv"))
write_csv(specimen_data, here("data/sample-data", "specimen-data.csv"))
write_csv(aliquot_data, here("data/sample-data", "aliquot-data.csv"))
write_csv(qubit_results, here("data/sample-data", "qubit-results.csv"))
write_csv(plate_layout, here("data/sample-data", "plate-layout.csv"))
