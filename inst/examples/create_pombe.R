
## Check if some input data is loaded; if not, load it!
pombe_se <- get0("pombe_se")
if (is.null(pombe_se)) {
  pombe_se <- make_pombe_se()
}
summary(pombe_se)
