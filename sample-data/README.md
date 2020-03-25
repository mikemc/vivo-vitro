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

## Sample data for sequencing results

* See the results/ folder for sample data to link with sequencing results.
