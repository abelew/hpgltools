# #' Metadata sanitizers for an expt
# #' @export
#setMethod(
#  "sanitize_metadata", signature = signature(meta = "expt"),
#  definition = function(meta, columns = NULL, na_string = "notapplicable",
#                        lower = TRUE, punct = TRUE, factorize = "heuristic",
#                        max_levels = NULL, spaces = FALSE, numbers = NULL) {
#    old_meta <- pData(meta)
#    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
#                                  lower = lower, punct = punct,
#                                  factorize = factorize, max_levels = max_levels,
#                                  spaces = spaces, numbers = numbers)
#    pData(meta) <- new_meta
#    return(meta)
#  })
#
# #' Metadata sanitizers for an expressionset
# #' @export
#setMethod(
#  "sanitize_metadata", signature = signature(meta = "ExpressionSet"),
#  definition = function(meta, columns = NULL, na_string = "notapplicable",
#                        lower = TRUE, punct = TRUE, spaces = FALSE,
#                        numbers = NULL) {
#    old_meta <- pData(meta)
#    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
#                                  lower = lower, punct = punct, spaces = spaces,
#                                  numbers = numbers)
#    pData(meta) <- new_meta
#    return(meta)
#  })
#
# #' Metadata sanitizers for a Summarized Experiment.
# #' @export
#setMethod(
#  "sanitize_metadata", signature = signature(meta = "SummarizedExperiment"),
#  definition = function(meta, columns = NULL, na_string = "notapplicable",
#                        lower = TRUE, punct = TRUE, factorize = "heuristic",
#                        max_levels = NULL, spaces = FALSE, numbers = NULL) {
#    old_meta <- pData(meta)
#    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
#                                  lower = lower, punct = punct,
#                                  factorize = factorize, max_levels = max_levels,
#                                  spaces = spaces, numbers = numbers)
#    pData(meta) <- new_meta
#    return(meta)
#  })
