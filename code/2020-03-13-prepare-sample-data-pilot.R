library(tidyverse)
library(here)

# Data for the DNA samples ----------------------------------------------------

# One "dna sample" corresponds to one DNA extraction (before aliquoting onto
# plates for the four sequencing centers)

## Load component tables
# map from sample ids to specimen id and aliquot number
dna_sample_map <- read_csv(here("data/sample-data", "dna-sample-map.csv"),
  col_types = "cii")
# information about individual specimens
specimen_data <- read_csv(here("data/sample-data", "specimen-data.csv"),
  col_types = "iciiDcciDii")
# aliquot_number to protocol used
aliquot_data <- read_csv(here("data/sample-data", "aliquot-data.csv"),
  col_types = "ic")
# dna concentrations
qubit_results <- read_csv(here("data/sample-data", "qubit-results.csv"),
  col_types = "cd")
# which well on the plate + which plates the sample was aliquoted onto
plate_layout <- read_csv(here("data/sample-data", "plate-layout.csv"),
  col_types = "ccllll")

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
write_csv(dna_sam, here("output/sample-data", "pilot-dna-sample-data.csv"))
saveRDS(dna_sam, here("output/sample-data", "pilot-dna-sample-data.Rds"))
