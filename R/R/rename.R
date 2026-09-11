# Matrix in, same matrix out, gene names on the rows.

#' Replace a matrix's row identifiers with gene names
#'
#' Detects and converts the row identifiers, then puts the names back on the
#' same object. Dimensions never change and rows never move: a row whose
#' identifier has no gene name keeps the identifier.
#'
#' Works on anything with `rownames()`: a base matrix, a sparse `Matrix` or a
#' data frame, because only the dimnames are touched.
#'
#' @param x A matrix or data frame with gene identifiers as row names.
#' @param ids Row identifiers. Defaults to `rownames(x)`.
#' @param mapping Skip detection entirely and use this `id`/`name` table, as
#'   accepted by [gene_names()]. Useful when the mapping is already in hand.
#' @param unique Disambiguate names that would otherwise collide, by appending
#'   the identifier. On by default: duplicate row names are legal in a matrix
#'   but break almost everything downstream.
#' @inheritParams geneid2name
#'
#' @return `x` with new row names. Nothing else about the object changes: no
#'   extra attributes, so a sparse `Matrix` stays a valid S4 object through
#'   subsetting, binding and arithmetic. Keep `rownames(x)` yourself if you
#'   need to undo the rename.
#'
#' @seealso [geneid2name()] for the vector form.
#' @export
#' @examples
#' \dontrun{
#' m <- matrix(1:6,
#'   nrow = 3,
#'   dimnames = list(c(
#'     "ENSG00000141510", "ENSG00000284733",
#'     "ENSG00000012048"
#'   ), c("s1", "s2"))
#' )
#' rename_rows(m)
#' #>                 s1 s2
#' #> TP53             1  4
#' #> ENSG00000284733  2  5
#' #> BRCA1            3  6
#' }
rename_rows <- function(x, ids = rownames(x), species = NULL, release = NULL,
                        assembly = NULL, source = NULL, mapping = NULL,
                        unique = TRUE, quiet = FALSE) {
  if (is.null(ids)) {
    stop("`x` has no row names; pass `ids` instead.", call. = FALSE)
  }
  if (length(ids) != nrow(x)) {
    stop(sprintf("`ids` has %d entries but `x` has %d rows.", length(ids), nrow(x)),
      call. = FALSE
    )
  }
  nm <- if (is.null(mapping)) {
    geneid2name(ids,
      species = species, release = release, assembly = assembly,
      source = source, unique = unique, quiet = quiet
    )
  } else {
    gene_names(ids, mapping, unique = unique, warn = !quiet)
  }
  rownames(x) <- as.vector(nm)
  x
}
