
#' A getter for the annotation databased used to create an expt/se.
#'
#' @param object One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom BiocGenerics annotation
#' @export
setMethod(
  "annotation", signature = signature(object = "expt"),
  definition = function(object) {
    BiocGenerics::annotation(object[["expressionset"]])
  })

#' A setter for the annotation database used to create an expt/se.
#'
#' @param object One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param value New annotation slot for the expt/se.
#' @importFrom Biobase annotation<-
#' @export
setMethod(
  "annotation<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    annotation(object[["expressionset"]]) <- value
    return(object)
  })

#' A getter to pull the assay data from an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
  "assay", signature = signature(x = "expt"),
  definition = function(x, withDimnames = TRUE, ...) {
    mtrx <- Biobase::exprs(x[["expressionset"]])
    return(mtrx)
  })

#' A setter to put the assay data into an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i specific samples to replace the data.
#' @param withDimnames I do not know.
#' @param ... Extra args, currently unused.
#' @param value New assay values to fill in the data structure.
#' @importFrom SummarizedExperiment assay<-
#' @export
setMethod(
  "assay<-", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::exprs(x[["expressionset"]]) <- value
    return(x)
  })

#' A getter to pull the assay data from an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param withDimnames I do not know.
#' @param ... Extra args!
#' @export
setMethod(
  "assay", signature = signature(x = "ExpressionSet"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::exprs(x)
  })

#' A setter to put the assay data into an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i Subset to replace.
#' @param withDimnames I do not know, I need to look this up.
#' @param ... Extra args.
#' @param value New values for the expressionset.
#' @export
setMethod(
  "assay<-", signature = signature(x = "ExpressionSet"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::exprs(x) <- value
    return(x)
  })

#' @export
batches <- function(se) {
  batches <- colData(se)[["batch"]]
  names(batches) <- sampleNames(se)
  return(batches)
}
setGeneric("batches")

#' @export
`batches<-` <- function(se, values) {
  colData(se)[["batch"]] <- values
  return(se)
}
setGeneric("batches<-")

#' Change the batches of an expt.
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify.
#' @param fact Batches to replace using this factor.
#' @param ids Specific samples to change.
#' @param ... Extra options are like spinach.
#' @return The original expt with some new metadata.
#' @seealso [create_expt()] [set_expt_conditions()] [Biobase]
#' @examples
#' \dontrun{
#'  expt = set_expt_batches(big_expt, factor = c(some,stuff,here))
#' }
#' @export
set_expt_batches <- function(expt, fact, ids = NULL, ...) {
  arglist <- list(...)
  original_batches <- pData(expt)[["batch"]]
  original_length <- length(original_batches)
  if (length(fact) == 1) {
    ## Assume it is a column in the design
    if (fact %in% colnames(pData(expt))) {
      fact <- pData(expt)[[fact]]
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  }

  if (length(fact) != original_length) {
    stop("The new factor of batches is not the same length as the original.")
  }
  pData(expt[["expressionset"]])[["batch"]] <- fact
  message("The number of samples by batch are: ")
  print(table(pData(expt)[["batch"]]))
  return(expt)
}
setGeneric("set_expt_batches")

#' @export
setMethod(
  "set_expt_batches", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, fact, ids = NULL, ...) {
    set_se_batches(expt, fact, ids, ...)
  })

#' @export
set_se_batches <- function(se, fact, ids = NULL, ...) {
  arglist <- list(...)
  original_batches <- colData(se)[["batch"]]
  original_length <- length(original_batches)
  if (length(fact) == 1) {
    ## Assume it is a column in the design
    if (fact %in% colnames(colData(se))) {
      fact <- colData(se)[[fact]]
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  }

  if (length(fact) != original_length) {
    stop("The new factor of batches is not the same length as the original.")
  }
  colData(se)[["batch"]] <- fact
  message("The number of samples by batch are: ")
  print(table(colData(se)[["batch"]]))
  return(se)
}

#' A getter to pull the sample data from an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param withDimnames Again, haven't looked it up yet.
#' @param ... Extra args.
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
  "colData", signature = signature(x = "expt"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::pData(x[["expressionset"]])
  })

#' A setter to put the sample data into an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i Subset to replace.
#' @param withDimnames indeed.
#' @param ... extra args.
#' @param value New Sample data for the expt.
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
  "colData<-", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::pData(x[["expressionset"]]) <- value
    return(x)
  })

#' A getter to pull the sample data from an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param withDimnames indeed.
#' @param ... extra args.
#' @export
setMethod(
  "colData", signature = signature(x = "ExpressionSet"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::pData(x)
  })

#' A setter to put the sample data into an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i Slice to replace.
#' @param withDimnames yes
#' @param ... args for the arglist
#' @param value New values for the expressionset.
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
  "colData<-", signature = signature(x = "ExpressionSet"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::pData(x) <- value
    return(x)
  })

setGeneric("colors", signature = signature(expt = "expt"),
           function(expt) standardGeneric("colors"))

#' A getter to pull the colors from an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "colors", signature = signature(expt = "expt"),
  definition = function(expt) {
    expt[["colors"]]
  })

#' A getter to pull the colors from a SummarizedExperiment.
#'
#' @param expt An expt.
#' @export
setMethod(
  "colors", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt) {
    S4Vectors::metadata(expt)[["colors"]]
  })

#' @export
`colors<-` <- function(expt, colors, ...) {
  set_se_colors(expt, colors, ...)
}
setGeneric("colors<-")

#' Get a named vector of colors by condition.
#'
#' Usually we give a vector of all samples by colors.  This just
#' simplifies that to one element each.  Currently only used in
#' combine_de_tables() but I think it will have use elsewhere.
#'
#' @param expt Expression from which to gather colors.
#' @return List of colors by condition.
#' @export
get_expt_colors <- function(expt, keep_underscore = TRUE) {
  all_colors <- colors(expt)
  condition_fact <- as.character(pData(expt)[["condition"]])
  if (isTRUE(keep_underscore)) {
    condition_fact <- gsub(pattern="[^_[:^punct:]]", replacement = "", x = condition_fact, perl = TRUE)
  } else {
    condition_fact <- gsub(x = condition_fact, pattern = "[[:punct:]]", replacement = "")
  }
  names(all_colors) <- condition_fact
  single_idx <- !duplicated(all_colors)
  all_colors <- all_colors[single_idx]
  return(all_colors)
}

#' @export
get_se_colors <- function(se, keep_underscore = TRUE) {
  all_colors <- colors(se)
  condition_fact <- as.character(colData(se)[["condition"]])
  if (isTRUE(keep_underscore)) {
    condition_fact <- gsub(pattern="[^_[:^punct:]]", replacement = "",
                           x = condition_fact, perl = TRUE)
  } else {
    condition_fact <- gsub(x = condition_fact, pattern = "[[:punct:]]", replacement = "")
  }
  names(all_colors) <- condition_fact
  single_idx <- !duplicated(all_colors)
  all_colors <- all_colors[single_idx]
  return(all_colors)
}

#' Change the colors of an expt
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param colors colors to replace
#' @param chosen_palette I usually use Dark2 as the RColorBrewer palette.
#' @param change_by Assuming a list is passed, cross reference by condition or sample?
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [set_expt_batches()] [RColorBrewer]
#' @examples
#' \dontrun{
#' unique(esmer_expt$design$conditions)
#' chosen_colors <- list(
#'    "cl14_epi" = "#FF8D59",
#'    "clbr_epi" = "#962F00",
#'    "cl14_tryp" = "#D06D7F",
#'    "clbr_tryp" = "#A4011F",
#'    "cl14_late" = "#6BD35E",
#'    "clbr_late" = "#1E7712",
#'    "cl14_mid" = "#7280FF",
#'    "clbr_mid" = "#000D7E")
#' esmer_expt <- set_expt_colors(expt = esmer_expt, colors = chosen_colors)
#' }
#' @export
set_expt_colors <- function(expt, colors = TRUE,
                            chosen_palette = "Dark2", change_by = "condition") {
  condition_factor <- as.factor(pData(expt)[["condition"]])

  ## Since I already have logic for named characters, just convert a list to one...
  if ("list" %in% class(colors)) {
    new_colors <- as.character(colors)
    names(new_colors) <- names(colors)
    colors <- new_colors
  }

  num_conditions <- length(levels(condition_factor))
  design <- pData(expt)
  num_samples <- nrow(design)
  sample_ids <- design[["sampleid"]]
  ## chosen_colors <- expt[["conditions"]]
  chosen_colors <- condition_factor
  chosen_names <- names(chosen_colors)
  sample_colors <- NULL
  if (is.null(colors) | isTRUE(colors)) {
    sample_colors <- sm(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else if (class(colors) == "factor") {
    if (change_by == "condition") {
      mesg("The new colors are a factor, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      colors_allocated <- names(colors) %in% levels(pData(expt)[["condition"]])
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        stop("Colors for the following categories are not being used: ",
             names(missing_colors), ".")
      }
      possible_conditions <- levels(pData(expt)[["condition"]])
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- pData(expt)[["condition"]] == cond
          missing_samples_by_cond <- rownames(pData(expt))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      mapping <- colors
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a factor, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "character") {
    if (is.null(names(colors))) {
      names(colors) <- levels(as.factor(expt[["conditions"]]))
    }
    if (change_by == "condition") {
      mesg("The new colors are a character, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- colors
      pd_factor <- as.factor(pData(expt)[["condition"]])
      possible_conditions <- levels(pd_factor)
      colors_allocated <- names(colors) %in% possible_conditions
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        warning("Colors for the following categories are not being used: ",
                names(missing_colors), ".")
      }
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- pData(expt)[["condition"]] == cond
          missing_samples_by_cond <- rownames(pData(expt))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a character, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "list") {
    if (change_by == "condition") {
      mesg("The new colors are a list, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- as.character(colors)
      names(mapping) <- names(colors)
      chosen_colors <- mapping[chosen_colors]
    } else if (change_by == "sample") {
      mesg("The new colors are a list, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- expt[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
        ## Set the condition for the changed samples to something unique.
        original_condition <- pData(expt)[sampleid, "condition"]
        changed_condition <- glue("{original_condition}{snum}")
        ## expt[["design"]][sampleid, "condition"] <- changed_condition
        tmp_pdata <- pData(expt)
        old_levels <- levels(tmp_pdata[["condition"]])
        new_levels <- c(old_levels, changed_condition)
        levels(tmp_pdata[["condition"]]) <- new_levels
        tmp_pdata[sampleid, "condition"] <- changed_condition
        pData(expt[["expressionset"]]) <- tmp_pdata
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (is.null(colors)) {
    mesg("Setting colors according to a color ramp.")
    colors <- sm(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    ## Check that all conditions are named in the color list:
    mapping <- setNames(colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else {
    warning("Number of colors provided does not match the number of conditions nor samples.")
    warning("Unsure of what to do, so choosing colors with RColorBrewer.")
    sample_colors <- suppressWarnings(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  }

  ## Catchall in case I forgot to set the names before now.
  names(chosen_colors) <- chosen_names
  expt[["colors"]] <- chosen_colors
  return(expt)
}

#' @export
set_se_colors <- function(se, colors = TRUE,
                            chosen_palette = "Dark2", change_by = "condition") {
  condition_factor <- as.factor(colData(se)[["condition"]])

  ## Since I already have logic for named characters, just convert a list to one...
  if ("list" %in% class(colors)) {
    new_colors <- as.character(colors)
    names(new_colors) <- names(colors)
    colors <- new_colors
  }

  num_conditions <- length(levels(condition_factor))
  design <- colData(se)
  num_samples <- nrow(design)
  sample_ids <- design[["sampleid"]]
  ## chosen_colors <- se[["conditions"]]
  chosen_colors <- condition_factor
  chosen_names <- names(chosen_colors)
  sample_colors <- NULL
  if (is.null(colors) | isTRUE(colors)) {
    sample_colors <- sm(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else if (class(colors) == "factor") {
    if (change_by == "condition") {
      mesg("The new colors are a factor, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      colors_allocated <- names(colors) %in% levels(colData(se)[["condition"]])
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        stop("Colors for the following categories are not being used: ",
             names(missing_colors), ".")
      }
      possible_conditions <- levels(colData(se)[["condition"]])
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- colData(se)[["condition"]] == cond
          missing_samples_by_cond <- rownames(colData(se))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      mapping <- colors
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a factor, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- se[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "character") {
    if (is.null(names(colors))) {
      names(colors) <- levels(as.factor(se[["conditions"]]))
    }
    if (change_by == "condition") {
      mesg("The new colors are a character, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- colors
      pd_factor <- as.factor(colData(se)[["condition"]])
      possible_conditions <- levels(pd_factor)
      colors_allocated <- names(colors) %in% possible_conditions
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        warning("Colors for the following categories are not being used: ",
                names(missing_colors), ".")
      }
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- colData(se)[["condition"]] == cond
          missing_samples_by_cond <- rownames(colData(se))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a character, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- se[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "list") {
    if (change_by == "condition") {
      mesg("The new colors are a list, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- as.character(colors)
      names(mapping) <- names(colors)
      chosen_colors <- mapping[chosen_colors]
    } else if (change_by == "sample") {
      mesg("The new colors are a list, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- se[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
        ## Set the condition for the changed samples to something unique.
        original_condition <- colData(se)[sampleid, "condition"]
        changed_condition <- glue("{original_condition}{snum}")
        ## se[["design"]][sampleid, "condition"] <- changed_condition
        tmp_pdata <- colData(se)
        old_levels <- levels(tmp_pdata[["condition"]])
        new_levels <- c(old_levels, changed_condition)
        levels(tmp_pdata[["condition"]]) <- new_levels
        tmp_pdata[sampleid, "condition"] <- changed_condition
        colData(se[["expressionset"]]) <- tmp_pdata
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (is.null(colors)) {
    mesg("Setting colors according to a color ramp.")
    colors <- sm(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    ## Check that all conditions are named in the color list:
    mapping <- setNames(colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else {
    warning("Number of colors provided does not match the number of conditions nor samples.")
    warning("Unsure of what to do, so choosing colors with RColorBrewer.")
    sample_colors <- suppressWarnings(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  }

  ## Catchall in case I forgot to set the names before now.
  names(chosen_colors) <- chosen_names
  metadata(se)[["colors"]] <- chosen_colors
  return(se)
}

#' A setter to put the colors into an expt.
#'
#' @param expt An expt.
#' @param value List of new colors.
#' @export
setMethod(
  "colors<-", signature = signature(expt = "expt"),
  definition = function(expt, value) {
    expt[["colors"]] <- value
    return(expt)
  })

#' A setter to put the colors into a SummarizedExperiment.
#'
#' @param expt A SummarizedExperiment.
#' @param value List of new colors.
#' @export
setMethod(
  "colors<-", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, value) {
    S4Vectors::metadata(expt)[["colors"]] <- value
    return(expt)
  })

#' @importFrom BiocGenerics conditions
#' @export
setMethod(
  "conditions", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, ...) {
    cond <- colData(object)[["condition"]]
    names(cond) <- sampleNames(object)
    return(cond)
  })

#' @export
`conditions<-` <- function(object, ..., value) {
  message("I think this should get pulled by the importFrom BiocGenerics.")
}
setGeneric("conditions<-")

#' @importFrom BiocGenerics `conditions<-`
#' @export
setMethod(
  "conditions<-", signature = signature(object = "expt"),
  definition = function(object, ..., value) {
    pData(object)[["condition"]] <- value
    new <- set_se_colors(object)
    return(new)
  })

#' @export
setMethod(
  "conditions<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, ..., value) {
    colData(object)[["condition"]] <- value
    se <- set_se_colors(object)
    return(se)
  })

#' Change the condition of an expt
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param fact Conditions to replace
#' @param ids Specific sample IDs to change.
#' @param prefix Add a prefix to the samples?
#' @param null_cell How to fill elements of the design which are null?
#' @param colors While we are here, set the colors.
#' @param ... Extra arguments are given to arglist.
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_batches()] [create_expt()]
#' @examples
#' \dontrun{
#'  expt = set_expt_conditions(big_expt, factor = c(some,stuff,here))
#' }
#' @export
set_expt_conditions <- function(expt, fact = NULL, ids = NULL,
                                prefix = NULL, null_cell = "null", colors = TRUE,
                                ...) {
  arglist <- list(...)
  if (!is.null(arglist[["factor"]])) {
    warning("I probably should change this argument to factor, but it is 'fact'.")
    fact <- arglist[["factor"]]
  }
  original_conditions <- pData(expt)[["condition"]]
  original_length <- length(original_conditions)
  original_num_conditions <- length(levels(as.factor(original_conditions)))
  new_expt <- expt  ## Explicitly copying expt to new_expt
  ## because when I run this as a function call() it seems to be not properly setting
  ## the conditions and I do not know why.
  fact_vector <- NULL
  fact_name <- "condition"
  if (!is.null(ids)) {
    ## Change specific id(s) to given condition(s).
    mesg("Setting condition for ids ", toString(ids), " to ", fact, ".")
    old_pdata <- pData(expt)
    old_cond <- as.character(old_pdata[["condition"]])
    names(old_cond) <- rownames(old_pdata)
    new_cond <- old_cond
    new_cond[ids] <- fact
    new_pdata <- old_pdata
    new_pdata[["condition"]] <- as.factor(new_cond)
    pData(new_expt[["expressionset"]]) <- new_pdata
    new_conditions <- as.character(new_expt[["conditions"]])
    names(new_conditions) <- names(new_expt[["conditions"]])
    new_conditions[ids] <- fact
    new_expt[["conditions"]] <- as.factor(new_conditions)
    ## new_expt[["design"]][["condition"]] <- new_cond
    fact_vector <- new_conditions
  } else if (length(fact) == 1) {
    fact_name <- fact
    ## Assume it is a column in the design
    if (fact %in% colnames(pData(expt))) {
      new_fact <- pData(expt)[[fact]]
      null_ids <- is.na(new_fact) | is.null(new_fact)
      ## Only do this if there are some null entries.
      if (sum(null_ids) > 0) {
        new_fact[null_ids] <- null_cell
      }
      if (!is.null(prefix)) {
        new_fact <- paste0(prefix, new_fact)
      }
      fact_vector <- new_fact
      new_expt[["conditions"]] <- new_fact
      pData(new_expt[["expressionset"]])[["condition"]] <- new_fact
      ## new_expt[["design"]][["condition"]] <- new_fact
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  } else if (length(fact) != original_length) {
    stop("The new factor of conditions is not the same length as the original.")
  } else {
    new_expt[["conditions"]] <- fact
    pData(new_expt[["expressionset"]])[["condition"]] <- fact
    ## new_expt[["design"]][["condition"]] <- fact
    fact_vector <- fact
  }

  message("The numbers of samples by condition are: ")
  print(table(pData(new_expt)[["condition"]]))
  condition_states <- levels(as.factor(pData(new_expt)[["condition"]]))
  if (class(colors)[1] == "list") {
    ## A list of colors may either be a color_choices list or
    ## a hash of states->color which could/should be a named vector.
    color_state_names <- names(colors)
    found_colors <- sum(color_state_names %in% condition_states)
    found_names <- sum(fact_name %in% color_state_names)
    ## In this first instance, the choices should be in this element.
    if (found_names > 0) {
      mesg("The colors appear to be a list delineated by state name.")
      colors <- colors[[fact]]
    } else if (found_colors > 0) {
      mesg("The colors appear to be a single list delineated by condition.")
    } else {
      message("A list of colors was provided, but element ", fact,
              " is not in it; using defaults")
      colors <- NULL
    }
  }
  new_expt <- set_expt_colors(new_expt, colors = colors)
  return(new_expt)
}

#' @export
set_se_conditions <- function(se, fact = NULL, ids = NULL, prefix = NULL,
                              null_cell = "null", colors = TRUE,
                              ...) {
  arglist <- list(...)
  if (!is.null(arglist[["factor"]])) {
    warning("I probably should change this argument to factor, but it is 'fact'.")
    fact <- arglist[["factor"]]
  }
  original_conditions <- colData(se)[["condition"]]
  original_length <- length(original_conditions)
  original_num_conditions <- length(levels(as.factor(original_conditions)))
  new_se <- se  ## Explicitly copying se to new_se
  ## because when I run this as a function call() it seems to be not properly setting
  ## the conditions and I do not know why.
  fact_vector <- NULL
  fact_name <- "condition"
  if (!is.null(ids)) {
    ## Change specific id(s) to given condition(s).
    mesg("Setting condition for ids ", toString(ids), " to ", fact, ".")
    old_pdata <- colData(se)
    old_cond <- as.character(old_pdata[["condition"]])
    names(old_cond) <- rownames(old_pdata)
    new_cond <- old_cond
    new_cond[ids] <- fact
    new_pdata <- old_pdata
    new_pdata[["condition"]] <- as.factor(new_cond)
    colData(new_se) <- new_pdata
  } else if (length(fact) == 1) {
    fact_name <- fact
    ## Assume it is a column in the design
    if (fact %in% colnames(colData(se))) {
      new_fact <- colData(se)[[fact]]
      null_ids <- is.na(new_fact) | is.null(new_fact)
      ## Only do this if there are some null entries.
      if (sum(null_ids) > 0) {
        new_fact[null_ids] <- null_cell
      }
      if (!is.null(prefix)) {
        new_fact <- paste0(prefix, new_fact)
      }
      fact_vector <- new_fact
      colData(new_se)[["condition"]] <- new_fact
      ## new_se[["design"]][["condition"]] <- new_fact
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  } else if (length(fact) != original_length) {
    stop("The new factor of conditions is not the same length as the original.")
  } else {
    colData(new_se)[["condition"]] <- fact
    ## new_se[["design"]][["condition"]] <- fact
    fact_vector <- fact
  }

  message("The numbers of samples by condition are: ")
  print(table(colData(new_se)[["condition"]]))
  condition_states <- levels(as.factor(colData(new_se)[["condition"]]))
  if (class(colors)[1] == "list") {
    ## A list of colors may either be a color_choices list or
    ## a hash of states->color which could/should be a named vector.
    color_state_names <- names(colors)
    found_colors <- sum(color_state_names %in% condition_states)
    found_names <- sum(fact_name %in% color_state_names)
    ## In this first instance, the choices should be in this element.
    if (found_names > 0) {
      mesg("The colors appear to be a list delineated by state name.")
      colors <- colors[[fact]]
    } else if (found_colors > 0) {
      mesg("The colors appear to be a single list delineated by condition.")
    } else {
      message("A list of colors was provided, but element ", fact,
              " is not in it; using defaults")
      colors <- NULL
    }
  }
  new_se <- set_se_colors(new_se, colors = colors)
  return(new_se)
}

#' A getter to pull the expression data from an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "exprs", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::exprs(object[["expressionset"]])
  })

#' A setter to put the expression data into an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "exprs<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(colnames(exprs(object)),
                           colnames(value))
    if (class(value)[1] == "data.frame") {
      value <- as.matrix(value)
    }
    exprs(object[["expressionset"]]) <- value
    return(object)
  })

#' A setter to put the expression data into an expt.
#'
#' @param object ExpressionSet to modify.
#' @param value New expression data.
#' @export
setMethod(
  "exprs<-", signature = signature(object = "ExpressionSet", value = "data.frame"),
  definition = function(object, value) {
    testthat::expect_equal(colnames(exprs(object)),
                           colnames(value))
    object <- as.matrix(exprs(value))
    exprs(object) <- value
    return(object)
  })

#' A getter to pull the expression data from a SummarizedExperiment.
#'
#' @param expt A SummarizedExperiment.
#' @export
setMethod(
  "exprs", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    SummarizedExperiment::assay(object)
  })

#' A setter to put the expression data to a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @export
setMethod(
  "exprs<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    testthat::expect_equal(colnames(exprs(object)),
                           colnames(value))
    if (class(value)[1] == "data.frame") {
      value <- as.matrix(value)
    }
    SummarizedExperiment::assay(object) <- value
    return(object)
  })

#' Change the factors (condition and batch) of an expt
#'
#' When exploring differential analyses, it might be useful to play with the
#' conditions/batches of the experiment.  Use this to make that easier.
#'
#' @param expt Expt to modify
#' @param condition New condition factor
#' @param batch New batch factor
#' @param ids Specific sample IDs to change.
#' @param table When set to 'metadata', use pData, otherwise fData.
#' @param class Set the data to this class by default.
#' @param columns Change these columns.
#' @param ... Arguments passed along (likely colors)
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [set_expt_batches()]
#' @examples
#' \dontrun{
#'  expt = set_expt_factors(big_expt, condition = "column", batch = "another_column")
#' }
#' @export
set_expt_factors <- function(expt, condition = NULL, batch = NULL, ids = NULL,
                             table = "metadata", class = "factor", columns = NULL, ...) {
  arglist <- list(...)
  if (!is.null(condition)) {
    expt <- set_expt_conditions(expt, fact = condition, ...)
  }
  if (!is.null(batch)) {
    expt <- set_expt_batches(expt, fact = batch, ...)
  }
  fd <- fData(expt)
  pd <- pData(expt)
  if (!is.null(columns)) {
    if (is.null(class)) {
      stop("If columns is set, then this assumes you want to set those columns to a given class.")
    }
    meta_columns <- colnames(pd)
    for (col in columns) {
      if (! col %in% meta_columns) {
        warning("The column: ", col, " is not in the metadata, skipping it.")
        next
      }
      mesg("Setting ", col, " to type ", class, ".")
      if (class == "factor") {
        if (table == "metadata") {
          pd[[col]] <- as.factor(pd[[col]])
        } else {
          fd[[col]] <- as.factor(fd[[col]])
        }
      } else if (class == "character") {
        if (table == "metadata") {
          pd[[col]] <- as.character(pd[[col]])
        } else {
          fd[[col]] <- as.character(fd[[col]])
        }
      } else if (class == "numeric") {
        if (table == "metadata") {
          pd[[col]] <- as.numeric(pd[[col]])
        } else {
          fd[[col]] <- as.numeric(fd[[col]])
        }
      } else {
        stop("I do not know this class.")
      }
    }
  }
  pData(expt[["expressionset"]]) <- pd
  fData(expt[["expressionset"]]) <- fd
  return(expt)
}

#' A getter to pull the gene annotation data from an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "fData", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::fData(object[["expressionset"]])
  })

#' A setter to put the gene annotation data into an expt.
#'
#' @param expt An expt.
#' @export
setMethod(
  "fData<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(fData(object)),
                           rownames(value))
    fData(object[["expressionset"]]) <- value
    return(object)
  })

#' A getter to pull the gene annotation data from a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @export
setMethod(
  "fData", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    SummarizedExperiment::rowData(object)
  })

#' A setter to put the gene annotation data into a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @export
setMethod(
  "fData<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(fData(object)),
                           rownames(value))
    SummarizedExperiment::rowData(object) <- value
    return(object)
  })

#' Switch the gene names of an expressionset using a column from fData.
#'
#' I am not sure if set_expt_genenames() is smart enough to check for
#' missing values.  It definitely handles duplicates.
#'
#' @param expt Current expressionSet.
#' @param new_column Column from the gene annotations containing the
#'  new gene IDs.
#' @return The expressionset with swapped out IDs.
#' @export
set_expt_genename_column <- function(expt, new_column) {
  start_df <- fData(expt)
  start_df[["start"]] <- rownames(start_df)
  start_df[["end"]] <- start_df[[new_column]]
  start_df <- start_df[, c("start", "end")]
  new <- set_expt_genenames(expt, start_df)
  return(new)
}

#' Change the gene names of an expt.
#'
#' I want to change all the gene names of a big expressionset to the
#' ortholog groups.  But I want to also continue using my expts.
#' Ergo this little function.
#'
#' @param expt Expt to modify
#' @param ids Specific sample IDs to change.
#' @param ... Extra arguments are given to arglist.
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [create_expt()]
#' @examples
#' \dontrun{
#'  expt = set_expt_conditions(big_expt, factor = c(some,stuff,here))
#' }
#' @export
set_expt_genenames <- function(expt, ids = NULL, ...) {
  arglist <- list(...)
  expr <- expt[["expressionset"]]

  ## Make sure the order of the IDs stays consistent.
  current_ids <- rownames(exprs(expr))
  if (class(ids) == "data.frame") {    ## our_column contains the IDs from my species.
    our_column <- NULL
    ## their_column contains the IDs from the species we want to map against.
    their_column <- NULL
    ## Grab the first ID in the first column.
    ## We will explicitly assume there are 2 columns in the data frame.
    test_first <- sum(ids[[1]] %in% current_ids)
    test_second <- sum(ids[[2]] %in% current_ids)
    if (test_first > 0) {
      mesg("Found: ", test_first, " ids in common using the first column of the IDs.")
      our_column <- colnames(ids)[1]
      their_column <- colnames(ids)[2]
    } else if (test_second > 0) {
      mesg("Found: ", test_second, " ids in common using the second column of the IDs.")
      our_column <- colnames(ids)[2]
      their_column <- colnames(ids)[1]
    } else {
      stop("Unable to match the IDs.")
    }
    mesg("Our column is: ", our_column, ", their column is: ", their_column, ".")
    ## Now the job is to ensure that the ordering is maintained.
    ## We need therefore to merge the rownames of the current IDs into the
    ## data frame of the ID mapping between our species.
    ## In addition, we must keep all of the IDs from our species (all.x = TRUE).
    exprs_id_df <- as.data.frame(current_ids)
    reordered <- merge(exprs_id_df, ids, by.x = "current_ids", by.y = our_column, all.x = TRUE)
    ## This merge should give a NA in the other column if there is no gene in the
    ## other species for a gene in our species.  In addition, it should give us
    ## duplicate IDs where a single gene from our species maps against multiple
    ## genes from the other species.
    ## Take care of the first scenario:
    na_ids <- is.na(reordered[[their_column]])
    reordered[na_ids, ] <- reordered[na_ids, "current_ids"]
    ## Then get rid of the duplicates.
    dup_ids <- duplicated(reordered[["current_ids"]])
    reordered <- reordered[!dup_ids, ]
    ## Now we should have a data frame with the same number of rows as our expressionset.
    ## The second column of this should contain the IDs from the other species when possible
    ## and a copy of the ID from our species when it is not.

    ## One final caveat: some of our new IDs may be duplicated in this (multigene families present
    ## in the other species),
    ## I will make.names() them and report how many there are, this might not be the correct way
    ## to handle this scenario!
    dup_ids <- sum(duplicated(reordered[[their_column]]))
    mesg("There are ", dup_ids, " duplicated IDs in the ", colnames(reordered)[2], " column.")
    ids <- make.names(reordered[[their_column]], unique = TRUE)
  }

  rownames(expr) <- ids
  expt[["expressionset"]] <- expr

  if (!is.null(expt[["tximport"]])) {
    rownames(expt[["tximport"]][["raw"]][["abundance"]]) <- ids
    rownames(expt[["tximport"]][["raw"]][["counts"]]) <- ids
    rownames(expt[["tximport"]][["raw"]][["length"]]) <- ids
    rownames(expt[["tximport"]][["scaled"]][["abundance"]]) <- ids
    rownames(expt[["tximport"]][["scaled"]][["counts"]]) <- ids
    rownames(expt[["tximport"]][["scaled"]][["length"]]) <- ids
  }
  return(expt)
}

#' Getter for library sizes in an expt.
#'
#' @param x input
#' @export
libsize <- function(x) {
  x[["libsize"]]
}
setGeneric("libsize")

#' Get the library sizes of a summarized experiment.
#' @export
setMethod(
  "libsize", signature = signature(x = "SummarizedExperiment"),
  definition = function(x) {
    meta <- S4Vectors::metadata(x)
    meta[["libsize"]]
  })

#' Setter for library sizes in an expt.
#' @export
`libsize<-` <- function(x, vector) {
  message("Using the generic?")
  x[["libsize"]] <- vector
  return(x)
}
setGeneric("libsize<-")

#' @export
setMethod(
  "libsize<-", signature = signature(x = "SummarizedExperiment"),
  definition = function(x, vector) {
    message("Using the se method.")
    meta <- S4Vectors::metadata(x)
    meta[["libsize"]] <- vector
    S4Vectors::metadata(x) <- meta
    return(x)
  })

#' A getter to pull the notes an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "notes", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::notes(object[["expressionset"]])
  })

#' A getter to pull the experimental metadata from an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::pData(object[["expressionset"]])
  })

#' A setter to put the experimental metadata into an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(pData(object)),
                           rownames(value))
    pData(object[["expressionset"]]) <- value
    return(object)
  })

#' A getter to pull the experimental metadata from a SummarizedExperiment.
#'
#' This is essentially synonymous with colData, except I cannot seem
#' to remember that function when I am working; so I just added
#' another signature to pData.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    SummarizedExperiment::colData(object)
  })

#' A setter to put the experimental metadata into a SummarizedExperiment.
#'
#' This is essentially synonymous with colData, except I cannot seem
#' to remember that function when I am working; so I just added
#' another signature to pData.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    testthat::expect_equal(
      rownames(SummarizedExperiment::colData(object)),
      rownames(value))
    SummarizedExperiment::colData(object) <- value
    return(object)
  })

#' A getter of the gene information from an expt, synonymous with fData().
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod(
  "rowData", signature = signature(x = "expt"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::fData(x[["expressionset"]])
  })

#' A setter to put the gene information into an expt.
#' @export
setMethod(
  "rowData<-", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    x <- Biobase::fData(x[["expressionset"]]) <- value
    return(x)
  })

#' A getter of the gene information from an ExpressionSet, synonymous with fData().
#' @export
setMethod(
  "rowData", signature = signature(x = "ExpressionSet"),
  definition = function(x, withDimnames = TRUE, ...) {
    Biobase::fData(x)
  })

#' A getter to get the samples names from an expt.
#' @export
setMethod(
  "sampleNames", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::sampleNames(object[["expressionset"]])
  })

#' A setter to put the samples names into an expt.
#' @export
setMethod(
  "sampleNames<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    set_expt_samplenames(object, value)
  })

#' A getter to get the samples names from a SummarizedExperiment.
#' @export
setMethod(
  "sampleNames", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    BiocGenerics::colnames(object)
  })

#' A setter to put the samples names into a SummarizedExperiment.
#' @export
setMethod(
  "sampleNames<-", signature = signature(object = "SummarizedExperiment", value = "character"),
  definition = function(object, value) {
    BiocGenerics::colnames(object) <- value
    return(object)
  })

#' Change the sample names of an expt.
#'
#' Sometimes one does not like the hpgl identifiers, so provide a way to change
#' them on-the-fly.
#'
#' @param expt Expt to modify
#' @param newnames New names, currently only a character vector.
#' @return expt Send back the expt with some new metadata
#' @seealso [set_expt_conditions()] [set_expt_batches()]
#' @examples
#' \dontrun{
#'  expt = set_expt_samplenames(expt, c("a","b","c","d","e","f"))
#' }
#' @export
set_expt_samplenames <- function(expt, newnames) {
  if (length(newnames) == 1) {
    ## assume it is a factor in the metadata and act accordingly.
    mesg("Using the column: ", newnames, " to rename the samples.")
    newer_names <- make.names(pData(expt)[[newnames]], unique=TRUE)
    result <- set_expt_samplenames(expt, newer_names)
    return(result)
  }

  new_expt <- expt
  ## oldnames <- rownames(new_expt[["design"]])
  oldnames <- sampleNames(new_expt)
  newnames <- make.unique(newnames)
  newnote <- glue("Sample names changed from: {toString(oldnames)} \\
                   to: {toString(newnames)} at: {date()}
")
  ## Things to modify include: batches, conditions
  names(batches(new_expt)) <- newnames
  new_colors <- colors(new_expt)
  names(new_colors) <- newnames
  colors(new_expt) <- new_colors
  ##newdesign <- new_expt[["design"]]
  ##newdesign[["oldnames"]] <- rownames(newdesign)
  ##rownames(newdesign) <- newnames
  ##newdesign[["sampleid"]] <- newnames
  ##new_expt[["design"]] <- newdesign
  new_expressionset <- new_expt[["expressionset"]]
  Biobase::sampleNames(new_expressionset) <- newnames
  pData(new_expressionset)[["sampleid"]] <- newnames
  new_expt[["expressionset"]] <- new_expressionset
  names(new_expt[["libsize"]]) <- newnames
  new_expt[["samplenames"]] <- newnames
  return(new_expt)
}

#' Get the state from an expt.
#'
#' @param expt One of my slightly modified ExpressionSets.
#' @return List with the methods used to modify the data (if any).
setGeneric("state",
  function(expt) standardGeneric("state"),
  signature = signature(expt = "expt"))

#' Get the state of the data in an expt.
#'
#' @param expt Experiment containing the state.
state <- function(expt) {
  return(expt[["state"]])
}

#' Extract the state of an expt vis a vis normalization.
#' @export
setMethod(
  "state", signature = signature(expt = "expt"),
  definition = function(expt) {
    expt[["state"]]
  })

#' Set the state of the data in an expt.
#'
#' @param expt Experiment requiring a state update.
#' @param value New state!
`state<-` <- function(expt, value) {
  expt[["state"]] <- value
  return(expt)
}
setGeneric("state<-")

#' Put the current state into an expt.
#' @export
setMethod(
  "state<-", signature = signature(expt = "expt"),
  definition = function(expt, value) {
    expt[["state"]] <- value
    return(expt)
  })

#' Get the state from a SummarizedExperiment.
#' @export
setMethod(
  "state", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt) {
    S4Vectors::metadata(expt)[["state"]]
  })

#' Put the state into a SummarizedExperiment.
#' @export
setMethod(
  "state<-", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, value) {
    S4Vectors::metadata(expt)[["state"]] <- value
    return(expt)
  })


#' Print a string describing what happened to this data.
#'
#' Sometimes it is nice to have a string like: log2(cpm(data)) describing what
#' happened to the data.
#'
#' @param expt The expressionset.
#' @param transform How was it transformed?
#' @param convert How was it converted?
#' @param norm How was it normalized?
#' @param filter How was it filtered?
#' @param batch How was it batch-corrected?
#' @return An expression describing what has been done to this data.
#' @seealso [create_expt()] [normalize_expt()]
#' @export
what_happened <- function(expt = NULL, transform = "raw", convert = "raw",
                          norm = "raw", filter = "raw", batch = "raw",
                          impute = "raw") {
  if (is.null(transform)) {
    transform <- "raw"
  }
  if (is.null(batch)) {
    batch <- "raw"
  }
  if (is.null(convert)) {
    convert <- "raw"
  }
  if (is.null(norm)) {
    norm <- "raw"
  }
  if (is.null(filter)) {
    filter <- "raw"
  }
  if (is.null(impute)) {
    impute <- "raw"
  }

  if (!is.null(expt)) {
    current <- state(expt)
    if (!is.null(current[["transform"]])) {
      transform <- current[["transform"]]
    }
    if (!is.null(current[["batch"]])) {
      batch <- current[["batch"]]
    }
    if (!is.null(current[["conversion"]])) {
      convert <- current[["conversion"]]
    }
    if (!is.null(current[["normalization"]])) {
      norm <- current[["normalization"]]
    }
    if (!is.null(current[["impute"]])) {
      norm <- current[["impute"]]
    }
    if (!is.null(current[["filter"]])) {
      filter <- current[["filter"]]
    }
  }
  ## Short circuit if nothing was done.
  if (transform == "raw" && batch == "raw" &&
        convert == "raw" && norm == "raw" &&
        filter == "raw") {
    what <- "raw(data)"
    return(what)
  }

  what <- ""
  if (impute != "raw") {
    what <- glue("{what}{impute}(")
  }
  if (transform != "raw") {
    what <- glue("{what}{transform}(")
  }
  if (batch != "raw") {
    if (isTRUE(batch)) {
      what <- glue("{what}batch-correct(")
    } else {
      what <- glue("{what}{batch}(")
    }
  }
  if (convert != "raw") {
    what <- glue("{what}{convert}(")
  }
  if (norm != "raw") {
    what <- glue("{what}{norm}(")
  }
  if (filter != "raw") {
    what <- glue("{what}{filter}(")
  }
  what <- glue("{what}data")
  if (transform != "raw") {
    what <- glue("{what})")
  }
  if (batch != "raw") {
    what <- glue("{what})")
  }
  if (convert != "raw") {
    what <- glue("{what})")
  }
  if (norm != "raw") {
    what <- glue("{what})")
  }
  if (filter != "raw") {
    what <- glue("{what})")
  }
  return(what)
}
