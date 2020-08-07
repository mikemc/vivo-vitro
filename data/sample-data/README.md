
Sample data from extraction experiment, created from the [original source
(Google drive spreadsheet)](https://docs.google.com/spreadsheets/d/1VYC4m4yk5K3AkXl2QOInaOuSla85ZkLHRDjAqFgabyA/)
by the script [00-fetch-from-drive.R](./00-fetch-from-drive.R).

* dna-sample-map.csv: Map from DNA sample to specimen and aliquot number
* specimen-data.csv: Data for each specimen (specimen type, collection date,
  etc.)
* aliquot-data.csv: Map from aliquot number to extraction protocol (same for
  each specimen)
* plate-layout.csv: Map from DNA sample to plate well, with logicals
  indicating which plates the sample was added to
* qubit-results.csv: DNA yield in ng measured by Qubit assay by Angie Mordant

Additional data from the sequencing centers. These were created by manually
from files received from the centers by email.

* a1-index-qc.csv: QC stats and indexing information from the Center-A1 MiSeq
  run; exported from the "Index QC" sheet in the sequencing report from A1
* s1-picogreen-results.csv: DNA concentration (ng/ul) measured by Picogreen
  assay by Center S1
* s2-qubit-results.csv: DNA concentration (ng/ul) measured by Qubit assay by
  Center S2
