# lefse output of oxygen

    Code
      mm_lefse
    Output
      microbiomeMarker-class inherited from phyloseq-class
      normalization method:              [ CPM ]
      microbiome marker identity method: [ lefse ]
      marker_table() Marker Table:       [ 12 microbiome markers with 5 variables ]
      otu_table()    OTU Table:          [ 276 taxa and  177 samples ]
      sample_data()  Sample Data:        [ 177 samples by  71 sample variables ]
      tax_table()    Taxonomy Table:     [ 276 taxa by 1 taxonomic ranks ]

---

    Code
      marker_table(mm_lefse)
    Output
                                                                                                                                     feature
      marker1                                                    k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae
      marker2                                k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium
      marker3        k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_s__
      marker4                                                                                        k__Bacteria|p__Firmicutes|c__Clostridia
      marker5                                                                       k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales
      marker6                             k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcaceae_g__
      marker7  k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcaceae_g__|s__Ruminococcaceae_g___s__
      marker8                                                                                           k__Bacteria|p__Firmicutes|c__Bacilli
      marker9                                                                        k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales
      marker10             k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__Streptococcus_s__
      marker11                                                   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae
      marker12                                  k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus
               enrich_group   ef_lda       pvalue         padj
      marker1       Healthy 4.993303 7.154793e-05 7.154793e-05
      marker2       Healthy 4.832765 5.914547e-04 5.914547e-04
      marker3       Healthy 4.830757 6.043983e-04 6.043983e-04
      marker4       Healthy 4.648541 7.176046e-04 7.176046e-04
      marker5       Healthy 4.648541 7.176046e-04 7.176046e-04
      marker6       Healthy 4.317335 6.990210e-03 6.990210e-03
      marker7       Healthy 4.317335 6.990210e-03 6.990210e-03
      marker8         Tumor 4.648541 7.176046e-04 7.176046e-04
      marker9         Tumor 4.514859 3.267911e-03 3.267911e-03
      marker10        Tumor 4.211031 4.371742e-03 4.371742e-03
      marker11        Tumor 4.122063 4.724670e-03 4.724670e-03
      marker12        Tumor 4.121855 4.812008e-03 4.812008e-03

