example_txt <- system.file("share", "sb", "trinotate_head.csv.xz", package = "hpgldata")

## Create a dataframe from the peculiarly formatted tsv/csv from trinotate.
trinotate_df <- load_trinotate_annotations(example_txt)
head(trinotate_df)

## Extract the gene ontology information from trinotate.
trinotate_go <- load_trinotate_go(example_txt,
                                  blast2go_column = "gene_ontology_blast")
head(trinotate_go[["go_table"]])
