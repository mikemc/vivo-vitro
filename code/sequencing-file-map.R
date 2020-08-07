library(tidyverse)
library(here)
library(fs)

dna_sam <- readRDS(here("output/sample-data", "dna-sample-data.Rds"))

# Get file paths relative to project root for given seq center
get_center_files <- function(center) {
  base_path <- str_glue("data/{str_to_lower(center)}/reads/raw")
  path(
    base_path, 
    dir_ls(here(base_path), glob = "*.fastq.gz") %>% path_file
  )
}

# A1 --------------------------------------------------------------------------
center <- "A1"
fns <- get_center_files(center)
fns %>% path_file %>% head(2)
#> [1] "1-1_S28_L001_R1_001.fastq.gz" "1-1_S28_L001_R2_001.fastq.gz"
fns %>% path_file %>% tail(4)
#> [1] "Water-NEGC_S96_L001_R1_001.fastq.gz"
#> [2] "Water-NEGC_S96_L001_R2_001.fastq.gz"
#> [3] "Zymo-POSC_S95_L001_R1_001.fastq.gz" 
#> [4] "Zymo-POSC_S95_L001_R2_001.fastq.gz" 
a1_files <- tibble(center_id = center, path = fns) %>%
  mutate(
    # R1 or R2?
    read = str_extract(path, "(?<=_)R[1-2](?=_)"),
    # DNA sample id; undefined for A1 controls
    dna_sample_id = str_extract(path, "[0-9]{1,2}-[1-4]") %>%
      str_replace("-", "_"),
    # Sample id (prefix by center)
    sample_id = case_when(
      !is.na(dna_sample_id) ~ str_c(center, dna_sample_id, sep = "_"),
      str_detect(path, "Zymo") ~ "A1_zymo_posc",
      str_detect(path, "Water") ~ "A1_water_negc",
    )
  ) %>%
  # pivot to separate cols for R1 and R2
  mutate(across(read, str_to_lower)) %>%
  pivot_wider(names_from = read, values_from = path, names_prefix = "path_")

# A2 --------------------------------------------------------------------------

# TBD

# S1 --------------------------------------------------------------------------
# NOTE: S1 is the only center with single end sequencing
center <- "S1"
fns <- get_center_files(center)
fns %>% path_file %>% head(2)
#> [1] "10_1.fastq.gz" "10_2.fastq.gz"
s1_files <- tibble(center_id = center, path = fns) %>%
  mutate(
    # DNA sample id
    dna_sample_id = str_extract(path, "[0-9]{1,2}_[1-4]"),
    # Sample id (prefix by center)
    sample_id = str_c(center, dna_sample_id, sep = "_")
  ) %>%
  rename(path_r1 = path) %>%
  # Add dummy path_r2 col to facilitate binding with other centers' tibbles
  add_column(path_r2 = fs::path(NA))

# S2 --------------------------------------------------------------------------
# NOTE: S2 files are named by well instead of DNA sample id
center <- "S2"
fns <- get_center_files(center)
fns %>% path_file %>% head(2)
#> [1] "A01_S13_CIDNS200509-022_R1_001.fastq.gz"
#> [2] "A01_S13_CIDNS200509-022_R2_001.fastq.gz"
s2_files <- tibble(center_id = center, path = fns) %>%
  mutate(
    # R1 or R2?
    read = str_extract(path, "(?<=_)R[1-2](?=_)"),
    #  Well
    well = str_extract(path, "[A-H][0-1][0-9]") %>%
      # Convert well format from A01 to A1
      str_replace("(?<=^[A-H])0(?=[0-9])", "")
  ) %>%
  # Replace well with DNA sample id
  left_join(dna_sam %>% select(dna_sample_id, well), by = "well") %>%
  select(-well) %>%
  # Sample id (prefix by center)
  mutate(sample_id = str_c(center, dna_sample_id, sep = "_")) %>%
  # pivot to separate cols for R1 and R2
  mutate(across(read, str_to_lower)) %>%
  pivot_wider(names_from = read, values_from = path, names_prefix = "path_")

# Combine ---------------------------------------------------------------------

bind_rows(a1_files, s1_files %>% mutate(path_r2 = fs::path(NA)), s2_files)

ftb <- bind_rows(a1_files, s1_files, s2_files) %>%
  # Sort
  left_join(dna_sam %>% select(dna_sample_id, specimen_id, aliquot_number), 
    by = "dna_sample_id") %>%
  arrange(center_id, specimen_id, aliquot_number) %>%
  # final cols
  select(sample_id, center_id, dna_sample_id, path_r1, path_r2)
ftb %>% head(5)
ftb %>% tail(5)
# tail(ftb)
write_csv(ftb, here("output/sample-data", "sequencing-file-map.csv"))
