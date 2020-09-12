
<!-- README.md is generated from README.Rmd. Please edit that file -->

Phyloseq objects
================

2020-09-12 Simplified phyloseq object for bias analysis
-------------------------------------------------------

The phyloseq object stored in “2020-09-12-phyloseq-focal-strains.Rds”
was created by “code/2020-09-12-import-to-phyloseq.R”. It contains a
phyloseq object with the current abundance data from A1, S1, and S2 for
a simplified feature set corresponding to just the main experimental
strains: the inoculum strains and the *T. mobilis* control.

    library(here)
    library(speedyseq)

    ps <- readRDS(here("output/phyloseq", "2020-09-12-phyloseq-focal-strains.Rds"))
    ps
    #> phyloseq-class experiment-level object
    #> otu_table()   OTU Table:         [ 14 taxa and 278 samples ]
    #> sample_data() Sample Data:       [ 278 samples by 22 sample variables ]
    #> tax_table()   Taxonomy Table:    [ 14 taxa by 8 taxonomic ranks ]
    #> phy_tree()    Phylogenetic Tree: [ 14 tips and 13 internal nodes ]
    #> taxa are columns

The taxonomy table contains the GTDB taxonomy and an additional column,
`strain_group`, that can be used to subset to just the inoculum strains,

    ps0 <- ps %>% 
      subset_taxa(strain_group == "Inoculum")

Each feature (taxon) corresponds to a strain; however, the feature
“Escherichia coli” will also include reads from the second *Escherichia
coli* strain that colonized the mice. *Faecalibacterium prausnitzii* is
included, though it did not colonize the mice. All samples are included;
the `specimen_type` variable can be used to exclude the controls.
Therefore, prior to analyzing bias one may wish to subset as follows,

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
    #>  [1] "Bacteroides uniformis"        "Bacteroides ovatus"           "Bacteroides caccae"          
    #>  [4] "Bacteroides thetaiotaomicron" "Barnesiella intestinihominis" "Akkermansia muciniphila"     
    #>  [7] "Marvinbryantia formatexigens" "Eubacterium rectale"          "Roseburia intestinalis"      
    #> [10] "Clostridium symbiosum"        "Collinsella aerofaciens"
