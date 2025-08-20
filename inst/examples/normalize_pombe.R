
## Check if some input data is loaded; if not, load it!
pombe_norm <- get0("pombe_norm")
if (is.null(pombe_norm)) {
  pombe_norm <- normalize(pombe_se, transform = "log2", convert = "cpm",
                          filter = TRUE)
}
summary(pombe_norm)
