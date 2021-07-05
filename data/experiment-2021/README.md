## Mock mixing volumes

Mixing volumes used for construction of the mock communities.
These values were originally generated in the [Experimental planning Rmd](/_posts/2021-06-22-planning-the-2021-experiment/planning-the-2021-experiment.Rmd).

Fields:

* `organism` - specification of the strain by "Genus species" except for /E. coli/ HS, for which "Genus species strain" is used to distinguish from /E. coli/ Nissle 1917. This format corresponds to the `display_name` field used elsewhere
* M1, M2, M3 - volume of culture added to the mock community of the given name, in mL, before being split into duplicates

## NanoDrop results

Measurement of DNA concentration of 6 samples to verify the extractions worked.

Fields:

* `dna_sample_id`, in the format `{base specimen}-{dilution step}-{aliquot number}`
* `dna_conc` = Measured concentration of DNA in ng/uL.

## OD600 results

Optical density (OD600) measurements of the cultures used to create the mock community mixtures.

Fields:

* `organism` - same as in the mixing-volumes table but followed by a three-letter abbreviation in parentheses
* `od600`; upper limit of detection of 2 is reached for several strains and denoted ">2.0"
