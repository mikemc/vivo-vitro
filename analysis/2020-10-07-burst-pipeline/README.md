## BURST-based shotgun analysis pipeline against v2 reference set

File descriptions:

- map-reads.org: Build a BURST database for the v2 reference set and map all S1
  and S2 samples with BURST using the ALLPATHs output option
- map-sample.sh: Sbatch submission script for mapping an individual sample
  (Fastq file)
- process-maps.Rmd: Process the BURST maps into count matrices
- burst-postprocessing-dev.Rmd: Notebook used to develop the functions for
  processing the burst maps
- output/: Location for BURST output (from map-reads.org) and the intermediate
  files `s1-counts.Rds` and `s2-counts.Rds` (from (process-maps.Rmd)
