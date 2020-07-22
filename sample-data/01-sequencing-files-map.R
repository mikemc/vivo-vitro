library(tidyverse)
library(here)
library(fs)

dotenv::load_dot_env(here::here(".env"))
data_path <- file.path(Sys.getenv("DATA_PATH"))

dna_sam <- readRDS(here("sample-data", "dna-sample-data.Rds"))

# A1 files
x <- "data/a1/reads/raw"
fns <- dir_ls(here(x), glob = "*.fastq.gz")
fns %>% path_file %>% head(2)
#> [1] "1-1_S28_L001_R1_001.fastq.gz" "1-1_S28_L001_R2_001.fastq.gz"
a1_files <- tibble(center_id = "A1", full_path = fns) %>%
  mutate(
    file = path_file(full_path),
    # Path relative to project root folder
    path = path(x, file),
    # R1 or R2?
    read = str_extract(file, "(?<=_)R[1-2](?=_)"),
    dna_sample_id = str_extract(file, "[0-9]{1,2}-[1-4]") %>%
      str_replace("-", "_")
  )

# A2 TBD

# S1 files
x <- "data/s1/reads/raw"
fns <- dir_ls(here(x), glob = "*.fastq.gz")
fns %>% path_file %>% head(2)
#> [1] "10_1.fastq.gz" "10_2.fastq.gz"
s1_files <- tibble(center_id = "S1", full_path = fns) %>%
  mutate(
    file = path_file(full_path),
    # Path relative to project root folder
    path = path(x, file),
    dna_sample_id = str_extract(file, "[0-9]{1,2}_[1-4]")
  )

# S2 files --------------------------------------------------------------------
x <- "data/s2/reads/raw"
fns <- dir_ls(here(x), glob = "*.fastq.gz")
fns %>% path_file %>% head(2)
#> [1] "A01_S13_CIDNS200509-022_R1_001.fastq.gz"
#> [2] "A01_S13_CIDNS200509-022_R2_001.fastq.gz"
s2_files <- tibble(center_id = "S2", full_path = fns) %>%
  mutate(
    file = path_file(full_path),
    # Path relative to project root folder
    path = path(x, file),
    # R1 or R2?
    read = str_extract(file, "(?<=_)R[1-2](?=_)"),
    well = str_extract(file, "^[A-H][0-1][0-9]") %>%
      # Convert well format from A01 to A1
      str_replace("(?<=^[A-H])0(?=[0-9])", "")
  ) %>%
  left_join(dna_sam %>% select(dna_sample_id, well), by = "well") %>%
  select(-well)

# Combine ---------------------------------------------------------------------
# add well back in as well as a double check

ftb <- bind_rows(a1_files, s1_files, s2_files) %>%
  select(-full_path) %>%
  relocate(dna_sample_id, .after = center_id)
ftb

write_csv(ftb, here("sample-data", "sequencing-read-files.csv"))
