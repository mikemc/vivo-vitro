
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

## Adjust variables
# Change Qubit results from estimated DNA yield (total ng in 100 uL) to DNA
# concentration (ng / uL) 
qubit_results <- qubit_results %>%
  mutate_at("dna_yield", ~ signif(. / 100, 5)) %>%
  rename(dna_conc = dna_yield)
# Shorten extraction protocol names from Protocol i to Pi
aliquot_data <- aliquot_data %>%
  mutate_at("protocol", str_replace, "Protocol ", "P")
# `number_of_pellets` is missing for mouse specimens, but should equal 1
specimen_data <- specimen_data %>%
  mutate(number_of_pellets = case_when(
      specimen_type == "fecal" ~ 1L,
      TRUE ~ number_of_pellets)
  )

## Save plain-text copies
dir.create(here("sample-data"))
write_csv(dna_sample_map, here("sample-data", "dna-sample-map.csv"))
write_csv(specimen_data, here("sample-data", "specimen-data.csv"))
write_csv(aliquot_data, here("sample-data", "aliquot-data.csv"))
write_csv(qubit_results, here("sample-data", "qubit-results.csv"))
write_csv(plate_layout, here("sample-data", "plate-layout.csv"))


## Join into a single table
dna_sam <- dna_sample_map %>%
  left_join(specimen_data, by = "specimen_id") %>%
  left_join(aliquot_data, by = "aliquot_number") %>%
  left_join(qubit_results, by = "dna_sample_id") %>%
  left_join(plate_layout, by = "dna_sample_id")

# Save plain-text and Rds versions (to preserve types)
write_csv(dna_sam, here("sample-data", "dna-sample-data.csv"))
saveRDS(dna_sam, here("sample-data", "dna-sample-data.Rds"))

# # Convert vars to factors (not run)
# dna_sam1 <- dna_sam %>%
#   mutate_at(
#     vars(specimen_id, aliquot_number, protocol, extraction_batch, 
#       specimen_type, host_species, host_sex), 
#     factor)

# Data for the sequenced samples ----------------------------------------------

# Here, "sample" corresponds to a sequenced sample, which came from one well on
# one of four plates. Sample names will take the form
# {center_id}_{dna_sample_id}

centers <- c("A1", "A2", "S1", "S2")

## Plate A1 (pilot data)

dna_sam <- readRDS(here("sample-data", "dna-sample-data.Rds"))

# Restrict to just the samples that were on the A1 plate; 
sam0 <- dna_sam %>%
  filter(A1) %>%
  select(-all_of(centers)) %>% 
  mutate(
    center_id = "A1",
    sample_id = str_glue("{center_id}_{dna_sample_id}") %>% as.character
  ) %>%
  select(-dna_sample_id) %>%
  select(sample_id, center_id, everything())

# Add the control samples used by center A1
sam.controls <- tibble(
  center_id = "A1",
  sample_id = c("A1_water_negc", "A1_zymo_posc"),
  specimen_type = c("water", "DNA mock"),
  # well info from "Final Plate Map" in spreadsheet from Center A1
  well = c("H12", "H11"),
  )
sam1 <- bind_rows(sam0, sam.controls)

# Add plate column and row
sam1 <- sam1 %>%
  mutate(
    row = str_sub(well, 1, 1) %>% as.ordered,
    column = str_sub(well, 2) %>% as.integer %>% as.ordered,
    )

write_csv(sam1, here("sample-data", "sample-data.csv"))
saveRDS(sam1, here("sample-data", "sample-data.Rds"))
