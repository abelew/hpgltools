## This function takes forever, let us only run it once in a while
## Perhaps I should have it create a savefile
if (interactive()) {
  pathfindr_kegg <- make_kegg_df("hsa")
}

## Load a pile of KEGGREST annotations.
coli_kegg <- load_kegg_annotations()
head(coli_kegg)
