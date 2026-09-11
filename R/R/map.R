# ID -> gene name, keeping the ID when no name exists.

#' A name that is really a name
#'
#' References signal "no symbol" by repeating the identifier, not by leaving the
#' field blank, as Ensembl and CellRanger both do. Blanks and `NA` count too.
#' @noRd
.is_name <- function(name, id) !is.na(name) & nzchar(name) & name != id

#' Convert gene identifiers to gene names
#'
#' Returns one value per input, in the input's order. An identifier with no
#' name in `mapping` is returned unchanged, so the result is always usable as
#' row labels.
#'
#' Version suffixes are ignored when matching, so `ENSG00000141510.16` resolves
#' against an unversioned mapping and the other way round.
#'
#' @param ids Character vector of gene identifiers.
#' @param mapping A data frame with `id` and `name` columns, or a named
#'   character vector of names indexed by identifier.
#' @param unique Make the result unique by appending the identifier to names
#'   that would otherwise collide. Off by default, since it changes names.
#' @param warn Warn when names collide or when little of the input mapped.
#'
#' @return A character vector the same length as `ids`, carrying a `mapped`
#'   attribute: the logical vector of which inputs found a real name.
#'
#' @seealso [fetch_mapping()] to obtain `mapping` for a species and release.
#' @export
#' @examples
#' m <- data.frame(
#'   id = c("ENSG00000141510", "ENSG00000284733"),
#'   name = c("TP53", "ENSG00000284733")
#' ) # second has no symbol
#' gene_names(c("ENSG00000141510.16", "ENSG00000284733", "ENSG00000999999"), m)
#' #> "TP53" "ENSG00000284733" "ENSG00000999999"
gene_names <- function(ids, mapping, unique = FALSE, warn = TRUE) {
  ids <- .as_ids(ids)

  if (is.character(mapping) && !is.null(names(mapping))) {
    mapping <- data.frame(id = names(mapping), name = unname(mapping))
  }
  if (!all(c("id", "name") %in% names(mapping))) {
    stop("`mapping` needs `id` and `name` columns.", call. = FALSE)
  }

  # factor columns would otherwise flow through as level codes
  mapping$id <- as.character(mapping$id)
  mapping$name <- as.character(mapping$name)

  keep <- .is_name(mapping$name, mapping$id)
  lut <- stats::setNames(mapping$name[keep], .bare(mapping$id[keep]))
  # duplicates agreeing on the name are normal; disagreeing ones are ambiguous
  dup <- duplicated(names(lut))
  if (any(dup) && warn) {
    # surviving de-duplication on (id, name) means more than one distinct name
    key <- paste(names(lut), lut, sep = "\r")
    conflict <- unique(names(lut)[dup & !duplicated(key)])
    if (length(conflict)) {
      warning(sprintf(
        "%d id(s) have conflicting names in `mapping`, e.g. %s; the first is used.",
        length(conflict), conflict[1]
      ), call. = FALSE)
    }
  }
  lut <- lut[!dup]

  hit <- lut[.bare(ids)]
  mapped <- unname(!is.na(hit))
  out <- .coalesce(unname(hit), ids) # the retention rule

  if (unique) {
    # suffix collisions with the id; repeated input ids still collide after that
    coll <- duplicated(out) | duplicated(out, fromLast = TRUE)
    out[coll & mapped] <- paste0(out[coll & mapped], "_", ids[coll & mapped])
    out <- make.unique(out, sep = "_")
  }

  if (warn) {
    n_dup <- sum(duplicated(out))
    if (n_dup > 0 && !unique) {
      warning(sprintf(
        "%d duplicate names; pass `unique = TRUE` to disambiguate.",
        n_dup
      ), call. = FALSE)
    }
    if (length(mapped) && mean(mapped) < 0.5) {
      warning(sprintf(
        "only %.0f%% of ids mapped -- wrong species or release?",
        100 * mean(mapped)
      ), call. = FALSE)
    }
  }

  attr(out, "mapped") <- mapped
  out
}
