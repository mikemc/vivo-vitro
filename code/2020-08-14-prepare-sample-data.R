library(tidyverse)
library(here)

# NOTE: Edited on 2021-06-20 to add the columns `volume_pbs` and
# `volume_per_aliquot`

# Data for the DNA samples ----------------------------------------------------

# One "dna sample" corresponds to one DNA extraction (before aliquoting onto
# plates for the four sequencing centers)

## Load component tables
# map from sample ids to specimen id and aliquot number
dna_sample_map <- read_csv(here("data/sample-data", "dna-sample-map.csv"),
  col_types = "cff"
) %>%
  mutate(across(where(is.factor), fct_inseq)) %>%
  glimpse
  
# information about individual specimens
specimen_data <- read_csv(here("data/sample-data", "specimen-data.csv"),
  col_types = cols(
    .default = col_factor(),
    collection_date = col_date(),
    number_of_pellets = col_integer(),
    number_of_aliquots = col_integer(),
    volume_pbs = col_double(),
    volume_per_aliquot = col_double()
  )
) %>%
  mutate(
    # `number_of_pellets` is missing for fecal specimens, but can be set to 1
    number_of_pellets = case_when(
      specimen_type == "fecal" ~ 1L,
      TRUE ~ number_of_pellets),
    # Order factors with numeric levels numerically
    across(c(specimen_id, collection_week, collection_day, extraction_batch),
      fct_inseq),
    # For other factors, capitalize levels and sort alphabetically
    across(c(specimen_type, host_species, host_sex), fct_relabel, 
      .fun = str_to_sentence),
    across(c(specimen_type, host_species, host_sex), fct_relevel, sort),
  ) %>%
  glimpse

# aliquot_number to protocol used
aliquot_data <- read_csv(here("data/sample-data", "aliquot-data.csv"),
  col_types = "fc"
) %>%
  # rename protocol var and shorten from "Protocol i" to "i"
  rename(extraction_protocol = protocol) %>%
  mutate(
    across(extraction_protocol, str_replace, "Protocol ", ""),
    across(extraction_protocol, as.factor),
    ) %>%
  mutate(across(where(is.factor), fct_inseq)) %>%
  glimpse

# which well on the plate + which plates the sample was aliquoted onto
plate_layout <- read_csv(here("data/sample-data", "plate-layout.csv")) %>%
  # Add plate row and column as ordered factors
  mutate(
    row = str_sub(well, 1, 1) %>% as.ordered,
    column = str_sub(well, 2) %>% as.integer %>% as.ordered,
  ) %>%
  # Adjust column order
  select(dna_sample_id, well, row, column, A1, A2, S1, S2) %>%
  glimpse

# DNA concentrations measured by us and centers S1 and S2
conc_us <- read_csv(here("data/sample-data", "qubit-results.csv")) %>%
  # Change Qubit results from estimated DNA yield (total ng in 100 uL) to DNA
  # concentration (ng / uL) 
  transmute(dna_sample_id, dna_conc = signif(dna_yield / 100, 5)) %>%
  glimpse
conc_s1 <- read_csv(here("data/sample-data", "s1-picogreen-results.csv")) %>%
  select(dna_sample_id = sample_id, dna_conc_s1) %>%
  mutate(
    across(dna_conc_s1, str_replace, "Undetected", "0"),
    across(dna_conc_s1, as.numeric)
  ) %>%
  glimpse
conc_s2 <- read_csv(here("data/sample-data", "s2-qubit-results.csv")) %>%
  mutate(
    across(dna_conc_s2, str_replace, "Too Low", "0"),
    across(dna_conc_s2, as.numeric),
    # Wells are in "A01" format; change to "A1" to join w/ plate layout
    across(well, str_replace, "(?<=[A-H])0", ""),
  ) %>%
  # Get sample_id from well
  left_join(plate_layout, by = "well") %>%
  select(dna_sample_id, dna_conc_s2) %>%
  glimpse

## Join into a single table + minor data cleaning
dna_sam <- dna_sample_map %>%
  left_join(specimen_data, by = "specimen_id") %>%
  left_join(aliquot_data, by = "aliquot_number") %>%
  left_join(conc_us, by = "dna_sample_id") %>%
  left_join(conc_s1, by = "dna_sample_id") %>%
  left_join(conc_s2, by = "dna_sample_id") %>%
  left_join(plate_layout, by = "dna_sample_id") %>%
  arrange(specimen_id, aliquot_number) %>%
  glimpse

# Save plain-text and Rds versions (to preserve types)
write_csv(dna_sam, here("output/sample-data", "dna-sample-data.csv"))
saveRDS(dna_sam, here("output/sample-data", "dna-sample-data.Rds"))

# Data for the sequenced samples ----------------------------------------------

# Get data frame with 1 row per sequenced sample
sam <- dna_sam %>%
  pivot_longer(c(A1, A2, S1, S2), 
    names_to = "center_id", values_to = "plated") %>%
  filter(plated) %>%
  select(-plated) %>%
  mutate(
    center_id = as.factor(center_id),
    # Set sample_id of form {center_id}_{dna_sample_id}
    sample_id = str_c(center_id, dna_sample_id, sep = "_"),
    ) %>%
  relocate(sample_id, center_id, .before = "dna_sample_id") %>%
  arrange(center_id, specimen_id, aliquot_number)

# Add the sample data for center A1's controls
a1_controls <- tibble(
  center_id = "A1" %>% as.factor,
  sample_id = c("A1_water_negc", "A1_zymo_posc"),
  specimen_type = c("Water", "DNA mock") %>% as.factor,
  # well info from "Final Plate Map" in spreadsheet from Center A1
  well = c("H12", "H11") %>% as.factor,
) %>%
  mutate(
    row = str_sub(well, 1, 1) %>% ordered(levels = levels(sam$row)),
    column = str_sub(well, 2) %>% ordered(levels = levels(sam$column))
  )
sam0 <- bind_rows(sam, a1_controls) %>%
  mutate(
    library_strategy = case_when(
      center_id %in% c("A1", "A2") ~ "Amplicon",
      center_id %in% c("S1", "S2") ~ "Shotgun",
      ) %>% as.factor,
    )

# Add DNA concentration measurements from the sequencing centers?

# Save plain-text and Rds versions (to preserve types)
write_csv(sam0, here("output/sample-data", "sample-data.csv"))
saveRDS(sam0, here("output/sample-data", "sample-data.Rds"))
