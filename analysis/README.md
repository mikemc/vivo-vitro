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
