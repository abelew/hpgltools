## Load some annotations from a gff file
example_gff <- system.file("share", "gas.gff", package = "hpgldata")

gas_iranges <- gff2irange(example_gff)
colnames(as.data.frame(gas_iranges))

gas_gr <- gff2gr(example_gff)
head(as.data.frame(gas_gr))

gas_df <- load_gff_annotations(example_gff, id_col = "locus_tag")
dim(gas_df)
