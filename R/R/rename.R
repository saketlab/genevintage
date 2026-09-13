# Shared by the matrix/data-frame branches of geneid2name() and correct_genenames().

#' Validate row identifiers against a matrix's row count
#' @noRd
.check_row_ids <- function(x, ids) {
  if (is.null(ids)) {
    stop("`x` has no row names; pass `row_ids` instead.", call. = FALSE)
  }
  if (length(ids) != nrow(x)) {
    stop(sprintf("`row_ids` has %d entries but `x` has %d rows.", length(ids), nrow(x)),
      call. = FALSE
    )
  }
}

#' Replace a matrix's row identifiers with the names `convert` computes
#'
#' Resolve and validate `row_ids`, then assign the converted names as plain
#' row labels so conversion metadata stays off the matrix.
#' @noRd
.rename_rows <- function(x, row_ids, convert) {
  rid <- if (is.null(row_ids)) rownames(x) else row_ids
  .check_row_ids(x, rid)
  rownames(x) <- as.vector(convert(rid))
  x
}
