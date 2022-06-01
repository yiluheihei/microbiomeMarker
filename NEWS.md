# microbiomeMarker 1.3.2

+ fix error on subgroup in lefse, #62, #55

# microbiomeMarker 1.3.1 (2022-05-26)

+ Development version on Bioconductor.

# microbiomeMarker 1.2.1 (2022-05-26)

+ Confounder analysis.
+ Comparison of different methods.

# microbiomeMarker 1.2.0 (2022-04-27) 

+ Released on Bioconductor 3.15.

# microbiomeMarker 1.1.2 (2022-03-07)

+ Development version on Bioconductor
+ Use 3rd version of testthat to fix test error (use `expect_snapshot()` rather 
than `expect_known_ouput`).
+ Add two new arguments in `plot_heatmap()` `scale_by_row` and `annotation_col`
to improve heatmap viaualization, #52.
+ Set slot `marker_table` to `NULL` if no marker was identified.
+ Add new import function `import_picrust2()` to import prediction functional
table from PICRUSt2, and all DA functions support for PICRUSt2 output data.
+ Keep color consistent between legend and plot in cladogram, #42. 
+ Add a new argument `clade_label_font_size` in `plot_cladogram()` to specify 
font size of clade label, #49.

# microbiomeMarker  1.1.1 (2021-03-07)

+ Add a para `only_marker` in `plot_cladogram` to specify whether only show the 
markers or all features in the cladogram.
+ Fix a bug in `run_test_multiple_groups()`, error group names for enrich 
groups (2021-10-12, #48).
+ Fix a bug in `plot_abundance()`, error var name of effect size in
`marker_table` (2021-10-17, #47).

# microbiomeMarker 1.0.0 (2021-10-27)

+ Released on Bioconductor.

# microbiomeMarker 0.99.1 (2021-10-11)

+ Accepted by Bioconductor.

# microbiomeMarker 0.99.0 (2021-09-01)

+ Submitted to Bioconductor
