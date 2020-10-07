# Convert a blast6out-format file into a tibble saved in an Rds file
args <- commandArgs(trailingOnly = TRUE)
stopifnot(identical(length(args), 1L))

fn <- args[[1]]
rds_fn <- fs::path_ext_set(fn, ".Rds")

if (fs::file_exists(rds_fn))
  stop("Rds file already exists")

col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore")
col_types <- "ccdiiiiiiidd"

tb <- readr::read_tsv(fn, col_names = col_names, col_types = col_types)

saveRDS(tb, rds_fn)
fs::file_delete(fn)
