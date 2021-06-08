### 2020-10-12-count-matrix.Rds

Count matrix from BURST output for bacterial reference genomes. 
Each row is a sample and each column is a reference genome, and counts denote the number of unambiguous insert maps to the reference.

[source Rmd file](./analysis/2020-10-07-burst-pipeline/process-maps.Rmd)

### 2021-06-08-count-matrix.Rds

Same as matrix in 2020-10-12-count-matrix.Rds except also includes the host read counts (first column), based on the host reads filtered in the v1 (BBMAP) pipeline.

[source Rmd file](./analysis/2020-10-07-burst-pipeline/add-host-reads.Rmd)
