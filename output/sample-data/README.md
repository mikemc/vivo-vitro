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
