## I thought I was starting to get a decent handle on S4
## and then I started getting the following error apropros of nothing:
## ! Failed to load R/zzz_attributes.R
## Caused by error in `setMethod()`:
## ! no existing definition for function ‘conditions<-’
## I my mind, the original definition of conditions<- is in BiocGenerics
## and I have an explicit import of it in 01_hpgltools.R
## Let us see what happens if I made an explicit import here
## and use the roxygen include tag.
## Trying suggestions from: https://stackoverflow.com/questions/69846121/no-definition-for-classes-from-another-package

## Also, I have a mess of color functions which need to be sorted out.
## I think I get it; when I started getting a bunch of new weirdo errors about
## missing function declarations in BiocGenerics/SummarizedExperiment/Biobase
## I had deleted the NAMESPACE file in order to regenerate it
## I think therefore that those generics were not getting loaded.

#' @include 01_hpgltools.R
#' @importFrom BiocGenerics conditions conditions<-
NULL

#' Get the batch column from a se.
#'
#' @param se Input summarized experiment.
#' @example inst/examples/attributes_se.R
#' @export
batches <- function(se) {
  batches <- colData(se)[["batch"]]
  names(batches) <- sampleNames(se)
  return(batches)
}
setGeneric("batches")

#' Add a batch column to a se.
#'
#' @param se Summarized Experiment to modify.
#' @param value vector of batches.
#' @export
`batches<-` <- function(se, value) {
  colData(se)[["batch"]] <- value
  return(se)
}
setGeneric("batches<-")

#' While I am struggling with S4 dispatch, here is a garbage generic.
#'
#' @param x Input
#' @export
get_colors <- function(x) {
  message("hpgltools generic to extract colors.")
}
setGeneric("get_colors")

#' Add colors to a dataset
#'
#' @param x Object to modify
#' @param ... extra arguments.
#' @param value vector of colors
#' @export
`colors<-` <- function(x, ..., value) {
  set_se_colors(x, value, ...)
}
setGeneric("colors<-")

#' Get a named vector of colors by condition.
#'
#' Usually we give a vector of all samples by colors.  This just
#' simplifies that to one element each.  Currently only used in
#' combine_de_tables() but I think it will have use elsewhere.
#'
#' @param expt Expression from which to gather colors.
#' @param fact Use this metadata column to set the colors.
#' @param levels When not null, colors may be set to arbitrary samples.
#' @return List of colors by condition.
#' @export
define_expt_colors <- function(expt, fact = "condition", levels = NULL) {
  all_colors <- get_colors(expt)
  names(all_colors) <- rownames(pData(expt))
  conditions_by_sample <- conditions(expt)
  if (is.null(levels)) {
    condition_fact <- levels(droplevels(as.factor(pData(expt)[[fact]])))
    colors_by_condition <- as.character(condition_fact)
    names(colors_by_condition) <- as.character(condition_fact)
  } else {
    colors_by_condition <- as.character(levels)
    names(colors_by_condition) <- as.character(levels)
  }
  for (f in seq_along(colors_by_condition)) {
    element <- colors_by_condition[f]
    ## See if there are more than 1 color per sample type
    choice_idx <- conditions_by_sample == element
    potential_colors <- unique(all_colors[choice_idx])
    if (length(potential_colors) == 1) {
      colors_by_condition[element] <- potential_colors
    } else if (length(potential_colors) > 1) {
      warning("Condition: ", element, " has multiple colors, chosing the first.")
      colors_by_condition[element] <- potential_colors[1]
    } else {
      warning("Condition: ", element, " has no color, setting it to red.")
      colors_by_condition[element] <- "#dd0000"
    }
  }
  return(colors_by_condition)
}

#' Get the colors from a summarized experiment.
#'
#' @param se Input se
#' @param keep_underscore Sanitize the columns for underscores?
#' @example inst/examples/attributes_se.R
#' @export
get_se_colors <- function(se, keep_underscore = TRUE) {
  all_colors <- get_colors(se)
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

#' Extract library sizes.
#'
#' @param x input
#' @example inst/examples/attributes_se.R
#' @export
libsize <- function(x) {
  message("libsize() extracts the relative sizes of every sample in an experiment.")
}
setGeneric("libsize")

#' Setter for library sizes in an expt.
#'
#' @param x Starting data
#' @param ... extra args.
#' @param value New library sizes
#' @export
`libsize<-` <- function(x, ..., value) {
  message("On occasion one might wish to set the relative sizes of samples in an experiment.")
  return(x)
}
setGeneric("libsize<-")

#' Remove/keep specifically named genes from a data structure.
#'
#' I find subsetting weirdly confusing.  Hopefully this function will allow one
#' to include/exclude specific genes/families based on string comparisons.
#'
#' @param input Data structure to filter.
#' @param invert The default is to remove the genes with the semantic strings.
#'  Keep them when inverted.
#' @param topn Take the topn most abundant genes rather than a text based heuristic.
#' @param semantic Character list of strings to search for in the annotation
#'  data.
#' @param semantic_column Column in the annotations to search.
#' @return A presumably smaller expt.
#' @seealso [Biobase]
#' @export
semantic_filter <- function(input, invert = FALSE, topn = NULL,
                            semantic = c("mucin", "sialidase", "RHS", "MASP", "DGF", "GP63"),
                            semantic_column = "description") {
  mtrx <- assay(input)
  annots <- rowData(input)
  if (isTRUE(invert)) {
    new_annots <- data.frame()
    new_mtrx <- data.frame()
  } else {
    new_annots <- annots
    new_mtrx <- mtrx
  }
  start_rows <- nrow(mtrx)
  numbers_removed <- 0
  if (is.null(topn)) {
    for (string in semantic) {
      idx <- NULL
      if (isTRUE(invert)) {
        ## Keep the rows which match the ~7 strings above.
        ## For these, we will re-grep the full table each time and just add the matches.
        type <- "Kept"
        if (semantic_column == "rownames") {
          idx <- grepl(pattern = string, x = rownames(annots))
        } else {
          idx <- grepl(pattern = string, x = annots[, semantic_column])
        }
        message("Hit ", sum(idx), " genes for term ", string, ".")
        ## Then, after grepping, just append the matched rows to the new annotations and matrix.
        tmp_annots <- annots[idx, ]
        tmp_mtrx <- mtrx[idx, ]
        new_annots <- rbind(new_annots, tmp_annots)
        new_mtrx <- rbind(new_mtrx, tmp_mtrx)
      } else {
        type <- "Removed"
        ## In the case of removals, I need to only grep what is left after each iteration.
        if (semantic_column == "rownames") {
          idx <- grepl(pattern = string, x = rownames(new_annots))
        } else {
          idx <- grepl(pattern = string, x = new_annots[, semantic_column])
        }
        mesg("Hit ", sum(idx), " genes for term ", string, ".")
        idx <- ! idx
        ## So, we take the index of stuff to keep, and just subset on that index.
        new_annots <- new_annots[idx, ]
        new_mtrx <- new_mtrx[idx, ]
      }
    } ## End for loop
    end_rows <- nrow(new_mtrx)
    lost_rows <- start_rows - end_rows
    message("semantic_expt_filter(): Removed ", lost_rows, " genes.")
  } else {
    ## Instead of a string based sematic filter, take the topn most abundant
    medians <- rowMedians(mtrx)
    new_order <- order(medians, decreasing = TRUE)
    reordered <- mtrx[new_order, ]
    subset <- rownames(head(reordered, n = topn))
    new_annots <- annots[subset, ]
  }

  keepers <- rownames(new_annots)
  new_data <- input[keepers, ]
  new_libsizes <- colSums(assay(new_data))
  return(new_data)
}
setGeneric("semantic_filter")

setMethod(
  "semantic_filter", signature = signature(input = "expt"),
  definition = function(input, invert = FALSE, topn = NULL,
                        semantic = c("mucin", "sialidase", "RHS", "MASP",
                                     "DGF", "GP63"),
                        semantic_column = "description") {
    semantic_expt_filter(input, invert = invert, topn = topn,
                         semantic = semantic,
                         semantic_column = semantic_column)
  })

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

#' Set the batches for a summarized experiment.
#'
#' @param se Input se
#' @param fact factor of new batches
#' @param ids specific IDs to change
#' @export
set_se_batches <- function(se, fact, ids = NULL) {
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

#' Set the colors of a summarized experiment.
#'
#' @param se Input se
#' @param colors Set of colors to add.
#' @param chosen_palette If colors is TRUE, use this palette to set colors.
#' @param change_by Use this factor to set the colors.
#' @importFrom SummarizedExperiment metadata
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

#' Set conditions to a se.
#'
#' @param se Input se
#' @param fact Factor of conditions
#' @param ids Set of ids to change.
#' @param prefix Prefix of each sample name
#' @param null_cell If a cell is null, what to change it to?
#' @param colors Set the colors as well?
#' @param ... Arbitrary arguments.
#' @export
set_conditions <- function(se, fact = NULL, ids = NULL, prefix = NULL,
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
setGeneric("set_conditions")

set_se_conditions <- function(...) {
  set_conditions(...)
}

setMethod(
  "set_conditions", signature = signature(se = "expt"),
  definition = function(se, fact = NULL, ids = NULL, prefix = NULL,
                        null_cell = "null", colors = TRUE,
                        ...) {
    set_expt_conditions(se, fact = fact, ids = ids,
                        prefix = prefix, null_cell = null_cell,
                        colors = colors, ...)
  })



set_factors <- function(se, condition = NULL, batch = NULL, ids = NULL,
                        table = "metadata", class = "factor",
                        columns = NULL, ...) {
  arglist <- list(...)
  if (!is.null(condition)) {
    se <- set_conditions(se, fact = condition, ...)
  }
  if (!is.null(batch)) {
    se <- set_batches(se, fact = batch, ...)
  }
  fd <- rowData(se)
  pd <- colData(se)
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
  colData(se) <- pd
  rowData(se) <- fd
  return(se)
}
setGeneric("set_factors")

setMethod(
  "set_factors", signature = signature(se = "expt"),
  definition = function(se, condition = NULL, batch = NULL, ids = NULL,
                        table = "metadata", class = "factor", columns = NULL, ...) {
    set_expt_factors(se, condition = condition, batch = batch,
                     ids = ids, table = table, class = class,
                     columns = columns, ...)
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

#' Set the genenames of a SE
#' @export
set_se_genenames <- function(se, ids = NULL, column = NULL, ...) {
  arglist <- list(...)
  current_ids <- rownames(assay(se))
  if (is.null(column) && is.null(ids)) {
    stop("Nothing was provided to change the IDs.")
  } else if (is.null(column)) {
    rownames(se) <- ids
  } else if (is.null(ids)) {
    new_ids <- make.names(rowData(se)[[column]], unique = TRUE)
    rownames(se) <- new_ids
  } else {
    message("Both a set of IDs and a column was provided, I am going to use the IDs.")
    rownames(se) <- ids
  }
  return(se)
}

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
  new_colors <- get_colors(new_expt)
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

#' Get the state of the data in an expt.
#'
#' @param input Experiment containing the state.
state <- function(input) {
  message("I am responsible for getting state from input data.")
}
setGeneric("state")

#' Set the state of the data in an expt.
#'
#' @param input Experiment requiring a state update.
#' @param value New state!
`state<-` <- function(input, value) {
  message("I am responsible for setting the state of input.")
  return(input)
}
setGeneric("state<-")

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
#' @param impute Was the data imputed?
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

########################################
## Generics Live here
########################################

#' Get the colors from an expt.
#'
#' @param expt One of my slightly modified ExpressionSets.
#' @return List with the methods used to modify the data (if any).
setGeneric("colors", signature = signature(expt = "expt"),
           function(expt) standardGeneric("colors"))

########################################
## Methods live here
########################################

# #' A getter for the annotation databased used to create an expt/se.
# # '
# #' @param object One of my various expressionset analogs, expt,
# #'  expressionSet, or summarizedExperiment.
# #' @importFrom BiocGenerics annotation
# #' @export
#setMethod(
#  "annotation", signature = signature(object = "expt"),
#  definition = function(object) {
#    BiocGenerics::annotation(object[["expressionset"]])
#  })

# #' A setter for the annotation database used to create an expt/se.
# #'
# #' @param object One of my various expressionset analogs, expt,
# #'  expressionSet, or summarizedExperiment.
# #' @param value New annotation slot for the expt/se.
# #' @importFrom Biobase annotation<-
# #' @export
#setMethod(
#  "annotation<-", signature = signature(object = "expt"),
#  definition = function(object, value) {
#    annotation(object[["expressionset"]]) <- value
#    return(object)
#  })

#'
NULL

# #' If you mess up the NAMESPACE file, the following becomes necessary
# #'
# #' message("I am from SummarizedExperiment and am explicitly imported, wtf.")
# #' @param x The SummarizedExperiment input
# #' @param i undef
# #' @param withDimnames undef
# #' @param ... extra args.
# #' @importFrom SummarizedExperiment assay
# #' @export
# assay <- function(x, i, withDimnames = TRUE, ...) {
#    SummarizedExperiment::assay(x, i, withDimnames = withDimnames, ...)
# }

#' A getter to pull the assay data from an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i I am guessing a subsetter
#' @param withDimnames I do not know.
#' @param ... Extra args!
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
  "assay", signature = signature(x = "ExpressionSet"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    Biobase::exprs(x)
  })

#' A getter to pull the assay data from an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i I am guessing a subsetter
#' @param withDimnames I do not know.
#' @param ... Extra args!
#' @export
setMethod(
  "assay", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    mtrx <- Biobase::exprs(x[["expressionset"]])
    return(mtrx)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from SummarizedExperiment and am explicitly imported, wtf.")
#' @param x The SummarizedExperiment input
#' @param i undef
#' @param withDimnames undef
#' @param ... extra args.
#' @param value New value.
#' @importFrom SummarizedExperiment assay<-
#' @export
`assay<-` <- function(x, i, withDimnames = TRUE, ..., value) {
  SummarizedExperiment::assay(x, i, withDimnames = withDimnames, ...) <- value
  return(x)
}

#' A setter to put the assay data into an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i Subset to replace.
#' @param withDimnames I do not know, I need to look this up.
#' @param ... Extra args.
#' @param value New values for the expressionset.
#' @importFrom SummarizedExperiment assay<-
#' @export
setMethod(
  "assay<-", signature = signature(x = "ExpressionSet"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::exprs(x) <- value
    return(x)
  })

#' A setter to put the assay data into an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i specific samples to replace the data.
#' @param withDimnames I do not know.
#' @param ... Extra args, currently unused.
#' @param value New assay values to fill in the data structure.
#' @export
setMethod(
  "assay<-", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ..., value) {
    Biobase::exprs(x[["expressionset"]]) <- value
    return(x)
  })

set_batches <- function(se, ...) {
  se <- set_se_batches(se, ...)
  return(se)
}
setGeneric("set_batches")

setMethod(
  "set_batches", signature = signature(se = "expt"),
  definition = function(se, fact, ids = NULL) {
    set_expt_batches(se, fact = fact, ids = ids)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from SummarizedExperiment and am explicitly imported, wtf.")
#' @param x The SummarizedExperiment input
#' @param i undef
#' @param withDimnames undef
#' @param ... extra args.
#' @importFrom SummarizedExperiment colData
#' @export
colData <- function(x, i, withDimnames = TRUE, ...) {
  SummarizedExperiment::colData(x, i, withDimnames = withDimnames, ...)
}

#' A getter to pull the sample data from an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i not j
#' @param withDimnames indeed.
#' @param ... extra args.
#' @export
setMethod(
  "colData", signature = signature(x = "ExpressionSet"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    Biobase::pData(x)
  })

#' A getter to pull the sample data from an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param i not j
#' @param withDimnames Again, haven't looked it up yet.
#' @param ... Extra args.
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
  "colData", signature = signature(x = "expt"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    Biobase::pData(x[["expressionset"]])
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from SummarizedExperiment and am explicitly imported, wtf.")
#' @param x The SummarizedExperiment input
#' @param ... extra args.
#' @param value New value.
#' @importFrom SummarizedExperiment colData<-
#' @export
`colData<-` <- function(x, ..., value) {
  SummarizedExperiment::colData(x, ...) <- value
  return(x)
}

#' A setter to put the sample data into an ExpressionSet.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param ... args for the arglist
#' @param value New values for the expressionset.
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
  "colData<-", signature = signature(x = "ExpressionSet"),
  definition = function(x, ..., value) {
    Biobase::pData(x) <- value
    return(x)
  })

setMethod(
  "colData<-", signature = signature(x = "SummarizedExperiment", value = "data.frame"),
  definition = function(x, ..., value) {
    value <- DataFrame(value)
    message("Recasting the data.frame to DataFrame.")
    SummarizedExperiment::colData(x, ...) <- value
    return(x)
  })

#' A setter to put the sample data into an expt.
#'
#' @param x One of my various expressionset analogs, expt,
#'  expressionSet, or summarizedExperiment.
#' @param ... extra args.a
#' @param value New Sample data for the expt.
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
  "colData<-", signature = signature(x = "expt"),
  definition = function(x, ..., value) {
    Biobase::pData(x[["expressionset"]]) <- value
    return(x)
  })

#' A setter to put the colors into an x.
#'
#' @param x An x.
#' @param ... Extra args
#' @param value List of new colors.
#' @export
setMethod(
  "colors<-", signature = signature(x = "expt"),
  definition = function(x, ..., value) {
    x[["colors"]] <- value
    return(x)
  })

#' A setter to put the colors into a SummarizedExperiment.
#'
#' @param x A SummarizedExperiment.
#' @param ... extra args.
#' @param value List of new colors.
#' @example inst/examples/attributes_se.R
#' @export
setMethod(
  "colors<-", signature = signature(x = "SummarizedExperiment"),
  definition = function(x, ..., value) {
    S4Vectors::metadata(x)[["colors"]] <- value
    return(x)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from BiocGenerics and am explicitly imported, wtf.")
#' @param object Input object
#' @param ... extra args
#' @importFrom BiocGenerics conditions
#' @export
conditions <- function(object, ...) {
  BiocGenerics::conditions(object, ...)
}

#' A getter to pull the conditions from an expt.
#'
#' @param object Input expt
#' @param ... extra args
#' @export
setMethod(
  "conditions", signature = signature(object = "expt"),
  definition = function(object, ...) {
    cond <- pData(object)[["condition"]]
    names(cond) <- rownames(pData(object))
    return(cond)
  })

#' Get the experimental conditions from a SE
#'
#' @param object Input SE
#' @param ... extra args.
#' @export
setMethod(
  "conditions", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, ...) {
    cond <- colData(object)[["condition"]]
    names(cond) <- sampleNames(object)
    return(cond)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from BiocGenerics and am explicitly imported, wtf.")
#' @param object Input object
#' @param ... extra args
#' @param value New value.
#' @importFrom BiocGenerics conditions<-
#' @export
`conditions<-` <- function(object, ..., value) {
  BiocGenerics::conditions(object, ...) <- value
  return(object)
}

#' Add experimental conditions to an expt.
#'
#' @param object Output expt
#' @param ... extra args
#' @param value vector of new conditions
#' @export
setMethod(
  "conditions<-", signature = signature(object = "expt"),
  definition = function(object, ..., value) {
    pData(object)[["condition"]] <- value
    new <- set_se_colors(object)
    return(new)
  })

#' A setter to put the conditions into a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @param ... arbitrary arguments
#' @param value List of new conditions.
#' @export
setMethod(
  "conditions<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, ..., value) {
    colData(object)[["condition"]] <- value
    se <- set_se_colors(object)
    return(se)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @importFrom Biobase exprs
#' @export
exprs <- function(object) {
  Biobase::exprs(object)
}

#' A getter to pull the expression data from an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "exprs", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::exprs(object[["expressionset"]])
  })

#' A getter to pull the expression data from a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @export
setMethod(
  "exprs", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    SummarizedExperiment::assay(object)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @param value new value
#' @importFrom Biobase exprs<-
#' @export
`exprs<-` <- function(object, value) {
  message("I am from Biobase and am explicitly imported, wtf.")
  object <- Biobase::exprs(object, value)
  return(object)
}

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

#' A setter to put the expression data into an expt.
#'
#' @param object An expt.
#' @param value New expression data.
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

#' A setter to put the expression data to a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @param value New expression data.
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


#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @importFrom Biobase fData
#' @export
fData <- function(object) {
  Biobase::fData(object)
}

#' A getter to pull the gene annotation data from an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "fData", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::fData(object[["expressionset"]])
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

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @param value new value
#' @importFrom Biobase fData<-
#' @export
`fData<-` <- function(object, value) {
  Biobase::fData(object) <- value
  return(object)
}

#' A setter to put the gene annotation data into an expt.
#'
#' @param object An expt.
#' @param value New annotations for the expressionset.
#' @export
setMethod(
  "fData<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(fData(object)),
                           rownames(value))
    fData(object[["expressionset"]]) <- value
    return(object)
  })

#' A setter to put the gene annotation data into a SummarizedExperiment.
#'
#' @param object A SummarizedExperiment.
#' @param value New annotations for the se.
#' @export
setMethod(
  "fData<-", signature = signature(object = "SummarizedExperiment"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(fData(object)),
                           rownames(value))
    SummarizedExperiment::rowData(object) <- value
    return(object)
  })

#' Get the colors from any arbitrary object, I am still trying to figure out
#' where my S4 dispatch is going wrong; this resides in BiocGenerics.
#'
#' @param x Object from which to get the colors.
setMethod(
  "get_colors", signature = signature(x = "ANY"),
  definition = function(x) {
    message("I assumed the Generic would pick this up.")
  })

#' A getter to pull the colors from an expt.
#'
#' @param x An x.
#' @export
setMethod(
  "get_colors", signature = signature(x = "expt"),
  definition = function(x) {
    x[["colors"]]
  })

#' A getter to pull the colors from a SummarizedExperiment.
#'
#' @param x An x.
#' @export
setMethod(
  "get_colors", signature = signature(x = "SummarizedExperiment"),
  definition = function(x) {
    S4Vectors::metadata(x)[["colors"]]
  })

#' I keep messing with my S4 dispatch of object attributes.
#'
#' @param input dimensional object with color information.
get_input_colors <- function(input) {
  message("I want to get away from expt-specific stuff.")
}
setGeneric("get_input_colors")

#' Extract colors from an expt.
#'
#' @param input input expt.
setMethod(
  "get_input_colors", signature = signature(input = "expt"),
  definition = function(input) {
    define_expt_colors(input)
  })

#' Extract colors from a SE.
#'
#' @param input input SE.
setMethod(
  "get_input_colors", signature = signature(input = "SummarizedExperiment"),
  definition = function(input) {
    metadata(input)[["colors"]]
  })

#' A getter to pull the library sizes from an expt.
#'
#' @param x An expt.
#' @export
setMethod(
  "libsize", signature = signature(x = "expt"),
  definition = function(x) {
    x[["libsize"]]
  })

#' Get the library sizes of a summarized experiment.
#'
#' @param x a summarized experiment.
#' @export
setMethod(
  "libsize", signature = signature(x = "SummarizedExperiment"),
  definition = function(x) {
    meta <- S4Vectors::metadata(x)
    meta[["libsize"]]
  })

#' Setter for library sizes in an expt.
#'
#' @param x expt to add library sizes
#' @param ... extra args
#' @param value new library sizes
#' @export
setMethod(
  "libsize<-", signature = signature(x = "expt", value = "ANY"),
  definition = function(x, ..., value) {
    x[["libsize"]] <- value
    return(x)
  })

#' Setter for library sizes in a se.
#'
#' @param x se to add library sizes
#' @param ... extra args
#' @param value new library sizes
#' @export
setMethod(
  "libsize<-", signature = signature(x = "SummarizedExperiment", value = "ANY"),
  definition = function(x, ..., value) {
    meta <- S4Vectors::metadata(x)
    meta[["libsize"]] <- vector
    S4Vectors::metadata(x) <- meta
    return(x)
  })

libsize_factor <- function(object) {
  libsizes <- libsize(object)
  max_libsize <- max(libsizes)
  norm_factor <- libsizes / max_libsize
  names(norm_factor) <- sampleNames(object)
  return(norm_factor)
}

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @importFrom Biobase notes
#' @export
notes <- function(object) {
  Biobase::notes(object)
}

#' A getter to pull the notes an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "notes", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::notes(object[["expressionset"]])
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @importFrom Biobase pData
#' @export
pData <- function(object) {
  Biobase::pData(object)
}

#' A getter to pull the experimental metadata from an expt.
#'
#' @param object An expt.
#' @export
setMethod(
  "pData", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::pData(object[["expressionset"]])
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

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @param value new value
#' @importFrom Biobase pData<-
#' @export
`pData<-` <- function(object, value) {
  message("I am from Biobase and am explicitly imported, wtf.")
  Biobase::pData(object) <- value
  return(object)
}

#' A setter to put the experimental metadata into an expt.
#'
#' @param object An expt.
#' @param value New metadata.
#' @export
setMethod(
  "pData<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    testthat::expect_equal(rownames(pData(object)),
                           rownames(value))
    pData(object[["expressionset"]]) <- value
    return(object)
  })

#' A setter to put the experimental metadata into a SummarizedExperiment.
#'
#' This is essentially synonymous with colData, except I cannot seem
#' to remember that function when I am working; so I just added
#' another signature to pData.
#'
#' @param object An expt.
#' @param value New metadata.
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

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from SummarzedExperiment and am explicitly imported, wtf.")
#' @param x Input object
#' @param use.names new value
#' @param ... extra args
#' @importFrom SummarizedExperiment rowData
#' @export
rowData <- function(x, use.names = TRUE, ...) {
  SummarizedExperiment::rowData(x, use.names = use.names, ...)
}

#' A getter of the gene information from an ExpressionSet, synonymous with fData().
#'
#' @param x Input
#' @param use.names yes
#' @param ... them too!
#' @export
setMethod(
  "rowData", signature = signature(x = "ExpressionSet"),
  definition = function(x, use.names = TRUE, ...) {
    Biobase::fData(x)
  })

#' A getter of the gene information from an expt, synonymous with fData().
#'
#' @param x Input
#' @param use.names Use those names...
#' @param ... them too!
#' @export
setMethod(
  "rowData", signature = signature(x = "expt"),
  definition = function(x, use.names = TRUE, ...) {
    Biobase::fData(x[["expressionset"]])
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from SummarzedExperiment and am explicitly imported, wtf.")
#' @param x Input object
#' @param ... extra args
#' @param value new value.
#' @importFrom SummarizedExperiment rowData<-
#' @export
`rowData<-` <- function(x, ..., value) {
  SummarizedExperiment::rowData(x, use.names = use.names, ...) <- value
  return(x)
}

#' A setter to put the gene information into an expt.
#'
#' @param x Input
#' @param ... them too!
#' @param value New annotations to put into the input
#' @export
setMethod(
  "rowData<-", signature = signature(x = "expt"),
  definition = function(x, ..., value) {
    x <- Biobase::fData(x[["expressionset"]]) <- value
    return(x)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @importFrom Biobase sampleNames
#' @export
sampleNames <- function(object) {
  Biobase::sampleNames(object)
}

#' A getter to get the samples names from an expt.
#'
#' @param object Input
#' @export
setMethod(
  "sampleNames", signature = signature(object = "expt"),
  definition = function(object) {
    Biobase::sampleNames(object[["expressionset"]])
  })

#' A getter to get the samples names from a SummarizedExperiment.
#'
#' @param object Input
#' @export
setMethod(
  "sampleNames", signature = signature(object = "SummarizedExperiment"),
  definition = function(object) {
    BiocGenerics::colnames(object)
  })

#' If you mess up the NAMESPACE file, the following becomes necessary
#'
#' message("I am from Biobase and am explicitly imported, wtf.")
#' @param object Input object
#' @param value new value
#' @importFrom Biobase sampleNames<-
#' @export
`sampleNames<-` <- function(object, value) {
  Biobase::sampleNames(object) <- value
  return(object)
}


#' A setter to put the samples names into an expt.
#'
#' @param object Input
#' @param value New names.
#' @export
setMethod(
  "sampleNames<-", signature = signature(object = "expt"),
  definition = function(object, value) {
    set_expt_samplenames(object, value)
  })

#' A setter to put the samples names into a SummarizedExperiment.
#'
#' @param object Input
#' @param value new names.
#' @export
setMethod(
  "sampleNames<-", signature = signature(object = "SummarizedExperiment", value = "character"),
  definition = function(object, value) {
    BiocGenerics::colnames(object) <- value
    return(object)
  })

#' set the batches of a SE
#'
#' @param expt Input summarized experiment.
#' @param fact Factor to use.
#' @param ids Or specific IDs
#' @param ... Other arguments.
#' @export
setMethod(
  "set_expt_batches", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, fact, ids = NULL, ...) {
    set_se_batches(expt, fact, ids, ...)
  })

#' Extract the state of an expt vis a vis normalization.
#'
#' @param input Input expt.
#' @export
setMethod(
  "state", signature = signature(input = "expt"),
  definition = function(input) {
    input[["state"]]
  })

#' Get the state from a SummarizedExperiment.
#'
#' @param input Input summarized experiment.
#' @export
setMethod(
  "state", signature = signature(input = "SummarizedExperiment"),
  definition = function(input) {
    S4Vectors::metadata(input)[["state"]]
  })

#' Put the current state into an expt.
#'
#' @param input Input expt
#' @param value New state.
#' @export
setMethod(
  "state<-", signature = signature(input = "expt"),
  definition = function(input, value) {
    input[["state"]] <- value
    return(input)
  })

#' Put the state into a SummarizedExperiment.
#'
#' @param input Input summarized experiment.
#' @param value new state.
#' @export
setMethod(
  "state<-", signature = signature(input = "SummarizedExperiment"),
  definition = function(input, value) {
    S4Vectors::metadata(input)[["state"]] <- value
    return(input)
  })

#' A predicate to check for the various tabular formats.
#'
#' I might want an argument to recast the input to a specific datatype.
#' @param x Datum to check.
tabularp <- function(x) {
  tabular <- FALSE
  if ("data.frame" %in% class(x) || "DFrame" %in% class(x) ||
        "tbl_df" %in% class(x) || "tibble" %in% class(x) ||
          "data.table" %in% class(x) || "matrix" %in% class(x)) {
    tabular <- TRUE
  }
  return(tabular)
}

txinfo <- function(se) {
  tximport_info <- S4Vectors::metadata(se)[["tximport"]]
  return(tximport_info)
}
setGeneric("txinfo")

setMethod(
  "txinfo", signature = signature(se = "expt"),
  definition = function(se) {
    tximport_info <- se[["tximport"]]
    return(tximport_info)
  })

setMethod(
  "txinfo", signature = signature(se = "expt"),
  definition = function(se) {
    tximport_info <- se[["tximport"]]
    return(tximport_info)
  })

`txinfo<-` <- function(se, value) {
  S4Vectors::metadata(se)[["tximport"]] <- value
  return(se)
}
setGeneric("txinfo<-")

setMethod(
  "txinfo<-", signature = signature(se = "expt"),
  definition = function(se, value) {
    se[["tximport"]] <- value
    return(se)
  })
