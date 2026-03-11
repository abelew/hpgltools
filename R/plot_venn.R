## plot_venn: Some functions to assist with venn diagrams.
## I think much of this may be removed in lieu of upset.

#' Give the intersection list produced by vennerable easier-to-read names.
#'
#' @param venn Vennerable data structure.
#' @param lst List of names which are legible.
rename_vennerable_intersections <- function(venn, lst) {
  intersects <- venn@IntersectionSets
  list_names <- names(lst)
  for (i in seq(from = 2, to = length(intersects) - 1)) {
    characters <- names(intersects)[i]
    characters <- strsplit(x = characters, split = "")[[1]]
    new_name <- ""
    for (c in seq_along(characters)) {
      char <- characters[c]
      if (char == "1") {
        list_name <- list_names[c]
        new_name <- glue("{new_name}{list_name}, ")
      }
    }
    new_name <- gsub(pattern = ", $", replacement = "", x = new_name)
    names(intersects)[i] <- new_name
  } ## Iterating through every intersection
  names(intersects)[1] <- "none"
  names(intersects)[length(intersects)] <- "all"
  return(intersects)
}

#' Extract rows of interesting information from vennerable.
#'
#' @param tables The set of tables associated with the intersections.
#' @param intersections The intersections themselves.
get_vennerable_rows <- function(tables, intersections) {
  ## Skip the 'none' table.
  int_tables <- intersections
  int_names <- names(intersections)
  ## Skip 'none'
  for (t in seq(from = 2, to = length(intersections))) {
    int_name <- int_names[t]
    chosen_table_name <- strsplit(x = int_name, split = ", ")[[1]][1]
    chosen_table <- tables[[chosen_table_name]]
    chosen_rows <- intersections[[t]]
    rows <- chosen_table[chosen_rows, ]
    int_tables[[t]] <- rows
  }
  return(int_tables)
}
