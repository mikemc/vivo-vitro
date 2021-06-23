
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Phyloseq objects

## 2020-09-12 Simplified phyloseq object for bias analysis

The phyloseq object stored in “2020-09-12-phyloseq-focal-strains.Rds”
was created by “/code/2020-09-12-import-to-phyloseq.R”. It contains a
phyloseq object with the current abundance data from A1, S1, and S2 for
a simplified feature set corresponding to just the main experimental
strains: the inoculum strains and the *T. mobilis* control.

``` r
library(here)
library(speedyseq)
library(tidyverse)

ps <- readRDS(here("output/phyloseq", "2020-09-12-phyloseq-focal-strains.Rds"))
ps
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:          [ 14 taxa and 278 samples ]:
#> sample_data() Sample Data:        [ 278 samples by 22 sample variables ]:
#> tax_table()   Taxonomy Table:     [ 14 taxa by 8 taxonomic ranks ]:
#> phy_tree()    Phylogenetic Tree:  [ 14 tips and 13 internal nodes ]:
#> taxa are columns
```

The taxonomy table contains the GTDB taxonomy and an additional column,
`strain_group`, that can be used to subset to just the inoculum strains,

``` r
ps0 <- ps %>% 
  subset_taxa(strain_group == "Inoculum")
```

Each feature (taxon) corresponds to a strain; however, the feature
“Escherichia coli” will also include reads from the second *Escherichia
coli* strain that colonized the mice. *Faecalibacterium prausnitzii* is
included, though it did not colonize the mice. All samples are included;
the `specimen_type` variable can be used to exclude the controls.
Therefore, prior to analyzing bias one may wish to subset as follows,

``` r
ps1 <- ps %>% 
  subset_samples(specimen_type %in% c("Fecal", "Inoculum")) %>%
  subset_taxa(
    strain_group == "Inoculum" & 
    genus != "Faecalibacterium" &
    genus != "Escherichia" 
  )
sample_data(ps1)$specimen_type %>% table
#> .
#>    Fecal Inoculum 
#>      174       90
taxa_names(ps1)
#>  [1] "Bacteroides uniformis"        "Bacteroides ovatus"          
#>  [3] "Bacteroides caccae"           "Bacteroides thetaiotaomicron"
#>  [5] "Barnesiella intestinihominis" "Akkermansia muciniphila"     
#>  [7] "Marvinbryantia formatexigens" "Eubacterium rectale"         
#>  [9] "Roseburia intestinalis"       "Clostridium symbiosum"       
#> [11] "Collinsella aerofaciens"
```

## 2020-11-28 Full results except host reads for the v2 pipeline

The phyloseq object in “2020-11-28-phyloseq.Rds” contains lightly
filtered results from bacterial profiling by the v2 amplicon and shotgun
pipelines for all sequencing centers. This file was created by the
script
“/analysis/2020-11-11-assign-strains/2020-11-28-phyloseq-import.Rmd”.
Counts of host reads are not included but will be in subsequent
versions.

``` r
ps <- readRDS(here("output/phyloseq", "2020-11-28-phyloseq.Rds"))
ps
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:          [ 126 taxa and 354 samples ]:
#> sample_data() Sample Data:        [ 354 samples by 22 sample variables ]:
#> tax_table()   Taxonomy Table:     [ 126 taxa by 10 taxonomic ranks ]:
#> refseq()      DNAStringSet:       [ 126 reference sequences ]
#> taxa are columns
```

Taxonomic features are ASVs (amplicon data) or reference accessions
(shotgun data). The taxonomy table includes Silva assignments for both
the ASVs and reference accessions (to genus level). It also includes
additional columns with information for the experimental, control, and
contaminant strains.

``` r
ps %>% tax_table %>% ps_tibble %>% glimpse
#> Rows: 126
#> Columns: 11
#> $ .otu         <chr> "v2_A1_ASV1", "v2_A1_ASV2", "v2_A1_ASV3", "v2_A1_ASV4", …
#> $ domain       <chr> "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteri…
#> $ phylum       <chr> "Bacteroidota", "Proteobacteria", "Verrucomicrobiota", "…
#> $ class        <chr> "Bacteroidia", "Gammaproteobacteria", "Verrucomicrobiae"…
#> $ order        <chr> "Bacteroidales", "Enterobacterales", "Verrucomicrobiales…
#> $ family       <chr> "Bacteroidaceae", "Enterobacteriaceae", "Akkermansiaceae…
#> $ genus        <chr> "Bacteroides", "Escherichia-Shigella", "Akkermansia", "B…
#> $ species      <chr> "Bacteroides ovatus", "Escherichia coli", "Akkermansia m…
#> $ organism     <chr> "Bacteroides ovatus", "Escherichia coli HS; Escherichia …
#> $ strain_group <chr> "Inoculum", "Inoculum; Mouse contaminant; Zymo", "Inocul…
#> $ source       <chr> "A1", "A1", "A1", "A1", "A1", "A1", "A1", "A1", "A1", "A…
ps %>% tax_table %>% ps_tibble %>% count(source)
#> # A tibble: 3 x 2
#>   source     n
#>   <chr>  <int>
#> 1 A1        37
#> 2 A2        69
#> 3 RefSeq    20
ps %>% tax_table %>% ps_tibble %>% count(strain_group)
#> # A tibble: 10 x 2
#>    strain_group                                         n
#>    <chr>                                            <int>
#>  1 Experimental control                                 3
#>  2 Inoculum                                            47
#>  3 Inoculum contaminant                                15
#>  4 Inoculum contaminant; Inoculum contaminant           4
#>  5 Inoculum contaminant; Zymo                           1
#>  6 Inoculum; Mouse contaminant; Zymo                    2
#>  7 Mouse contaminant                                    3
#>  8 Zymo                                                 6
#>  9 Zymo; Inoculum contaminant; Inoculum contaminant     1
#> 10 <NA>                                                44
```

The strain_group and source columns can be used to filter to specific
strains, while the species and organism columns are useful for filtering
as well as taxonomic merging. For example, the following code chunk
restricts to the inoculum strains and other *E. coli* strains, merges
taxonomic features from both amplicon and shotgun measurements to the
species level, and names the new taxa by species; it also restricts to
the non-control samples.

``` r
ps0 <- ps %>%
  filter_sample_data(
    specimen_type %in% c("Fecal", "Inoculum"),
  ) %>%
  filter_tax_table(
    strain_group == "Inoculum" | species == "Escherichia coli",
  ) %>%
  tax_glom("species") %>%
  mutate_tax_table(.otu = species)
ps0 %>% taxa_names
#>  [1] "Bacteroides ovatus"           "Escherichia coli"            
#>  [3] "Akkermansia muciniphila"      "Bacteroides uniformis"       
#>  [5] "Bacteroides thetaiotaomicron" "Clostridium symbiosum"       
#>  [7] "Roseburia intestinalis"       "Bacteroides caccae"          
#>  [9] "Faecalibacterium prausnitzii" "Collinsella aerofaciens"     
#> [11] "Marvinbryantia formatexigens" "Eubacterium rectale"         
#> [13] "Barnesiella intestinihominis"
```

This simplified phyloseq object is a natural starting point for
analyzing bias across all 4 sequencing centers, keeping in mind that

-   It lumps together reads from the two *E. coli* strains (HS and
    Nissle 1917)
-   It includes *Faecalibacterium prausnitzii*, which apparently did not
    colonize the mice and so we suspect any presence in the fecal
    samples is due to cross contamination

Thus depending on aims it can also make sense to remove one or both of
these species prior to analysis.

As noted, some light sample and taxonomic filtering of the has already
been done in the base phyloseq object; in particular,

-   Only samples with more than 1000 reads are included; this primarily
    removed samples from Center A2
-   ASVs whose phylum could not be assigned or were unexpectedly short
    were removed
-   Taxa (ASVs or reference genomes) with 100 or fewer reads were
    removed
