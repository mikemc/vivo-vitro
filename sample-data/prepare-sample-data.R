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
dir.create(here("sample-data"))
write_csv(dna_sample_map, here("sample-data", "dna-sample-map.csv"))
write_csv(specimen_data, here("sample-data", "specimen-data.csv"))
write_csv(aliquot_data, here("sample-data", "aliquot-data.csv"))
write_csv(qubit_results, here("sample-data", "qubit-results.csv"))
write_csv(plate_layout, here("sample-data", "plate-layout.csv"))

## Join into a single table + minor data cleaning
dna_sam <- dna_sample_map %>%
  left_join(specimen_data, by = "specimen_id") %>%
  left_join(aliquot_data, by = "aliquot_number") %>%
  left_join(qubit_results, by = "dna_sample_id") %>%
  left_join(plate_layout, by = "dna_sample_id") %>%
  # Change Qubit results from estimated DNA yield (total ng in 100 uL) to DNA
  # concentration (ng / uL) 
  mutate_at("dna_yield", ~ signif(. / 100, 5)) %>%
  rename(dna_conc = dna_yield) %>%
  # rename protocol var and shorten from "Protocol i" to "i"
  rename(extraction_protocol = protocol) %>%
  mutate_at("extraction_protocol", str_replace, "Protocol ", "") %>%
  # `number_of_pellets` is missing for fecal specimens, but can be set to 1
  mutate(
    number_of_pellets = case_when(
      specimen_type == "fecal" ~ 1L,
      TRUE ~ number_of_pellets
    )
  )

# Save plain-text and Rds versions (to preserve types)
write_csv(dna_sam, here("sample-data", "dna-sample-data.csv"))
saveRDS(dna_sam, here("sample-data", "dna-sample-data.Rds"))

# # Convert vars to factors (not run)
# dna_sam1 <- dna_sam %>%
#   mutate_at(
#     vars(specimen_id, aliquot_number, protocol, extraction_batch, 
#       specimen_type, host_species, host_sex), 
#     factor)
