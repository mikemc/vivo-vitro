## Data directory

The analysis files assume that sequencing reads are in an external directory
whose path is defined in the `DATA_PATH` variable of the `.env` file. This
folder should have the structure

```
DATA_PATH
├── a1
│   └── reads
├── a2
│   └── reads
├── references
├── s1
│   └── reads
├── s2
│   └── reads
└── sanger
```
<!-- edited from output of `tree -L 2` ran in `DATA_PATH` dir -->

where a1/, a2/, s1/, and s2/ house the data from the four sequencing centers,
sanger/ contains Sanger sequencing results for each strain, and references/
contains 16S and whole genome reference sequences obtained from public
databases.

## Analyses

### Shotgun pipeline

- [2020-06-15-bbmap-pipeline-test](./2020-06-15-bbmap-pipeline-test/) develops
  and tests the initial version of our pipeline on a subset of samples from S1.

- [2020-07-09-bbmap-pipeline](./2020-07-09-bbmap-pipeline/) runs the initial
  pipeline on all S1 and S2 samples, creating a table of mapped insert counts
  per reference genome in [output/](../output/)

- [2020-07-16-mash](./2020-07-16-mash/) uses the tool [Mash
  Screen](http:/dx.doi.org/10.1186/s13059-019-1841-x) to screen of a subset of
  samples against RefSeq and Genbank genomes, for use in identifying
  contaminant strains and validating our chosen references for our experimental
  strains.

- [2020-09-18-shotgun-testing](./2020-09-18-shotgun-testing/) 
  - Use Mash Screen results to choose new genomes; download from NCBI and save
    metadata
  - Extract 16S sequences
  - Test mapping pipelines: Simulate reads and map with BBSplit and BURST
  - Evaluate mapping performance on the simulated reads
