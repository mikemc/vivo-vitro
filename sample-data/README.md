# Sample data

Terminology: I use "sample" to correspond to an individual sequencing library;
these are the samples that correspond to rows in feature-abundance tables. I
use the term "DNA sample" to refer to the result of an individual DNA
extraction, before it was split into the four replicate plates for the four
sequencing centers. DNA sample IDs have the form `{specimen_id}_{aliquot_id}`;
Sample IDs have the form `{center_id}_{specimen_id}_{aliquot_id}`.

Where there are .Rds and .csv versions of a file, the data is identical except
that the Rds files include additional data formatting information (e.g. the
`row` and `column` variables are stored as ordered factors).

## DNA sample data for DNA extractions

* dna-sample-data.(csv/Rds): The DNA sample data table.

* dictionary.csv: "Data dictionary" with information about each variable.

### Component tables

The extraction sample data is derived from the following component tables.

* dna-sample-map.csv: Map from DNA sample ids to specimen id and aliquot
  number.

* specimen-data.csv: Information about each specimen.

* aliquot-data.csv: Map from aliquot number to extraction protocol.

* qubit-results.csv: DNA concentrations measured by Qubit assay.

* plate-layout.csv: Map from plated DNA samples to their corresponding well on
  the plate, plus columns for each sequencing center indicating whether the
  well was filled on that plate.

## QC data from sequencing centers

### a1-index-qc.csv

TODO

### s1-picogreen-results.csv: DNA concentrations measured by Center S1

S1 DNA concentrations were measured via PicoGreen assay.

Columns:

1. sample_id: The sample id's that we submitted to S1; these are actually the
   specimen ids.
1. well: Plate well in the format `A01`
1. dna_conc_s1: DNA concentrations measured by S1, in ng/uL

### s2-qubit-results.csv: DNA concentrations measured from Center S2

S2 DNA concentrations were measured via Qubit assay. 

Columns:

1. index: Index from 1:`n_samples`
1. well: Plate well in the format `A01`; also the sample names we submitted and
   that were used by S2
1. dna_conc_s2: DNA concentrations measured by S2, in ng/uL
1. dna_conc_submitted: DNA concentrations that we reported to S2, in ng/uL

## Sample data for sequencing results

* See the results/ folder for sample data to link with sequencing results.
