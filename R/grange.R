gr_from_to <- function(gr, from, to, column = "gene_id",
                       type = "gene", padding = 100) {
  meta <- as.data.frame(mcols(gr))
  meta_na <- is.na(meta[[column]])
  if (sum(meta_na) > 0) {
    message("Dropping ", sum(meta_na), " rows with NA from column ", column, ".")
    gr <- gr[!meta_na, ]
    meta <- mcols(gr)
  }

  beginning_id <- meta[[column]] == from
  if (sum(beginning_id, na.rm = TRUE) != 1) {
    stop("This only works if we get one gene for the beginning.")
  }
  beginning_gr <- gr[beginning_id, ]
  end_id <- meta[[column]] == to
  if (sum(end_id) != 1) {
    stop("This only works if we get one gene for the end.")
  }
  end_gr <- gr[end_id, ]
  wanted_begin <- start(beginning_gr)
  wanted_end <- end(end_gr)
  wanted_idx <- start(gr) >= wanted_begin &
    end(gr) <= wanted_end
  mesg("This range results in ", sum(wanted_idx), " genes in the subset.")
  region <- gr[wanted_idx, ]
  mcols(region)[, "orientation"] <- TRUE
  minus_idx <- strand(region) == "-" | strand(region) == -1
  mcols(region)[minus_idx, "orientation"] <- FALSE
  mcols(region)[, "gene_name"] <- mcols(region)[, "gene_id"]
  mcols(region)[, "gene_biotype"] <- "gene"
  return(region)
}
