# ID -> gene name, keeping the ID when no name exists.

#' Whether a name is present and differs from the identifier
#'
#' Repeated identifiers, blanks and `NA` are treated as missing symbols.
#' @noRd
.is_name <- function(name, id) !is.na(name) & nzchar(name) & name != id

#' The collision warning, classed so a caller can catch just this one
#'
#' The `genevintage_duplicate_names` class lets callers muffle collision
#' warnings when they disambiguate names downstream.
#' @noRd
.duplicate_names_condition <- function(n) {
  warningCondition(
    sprintf("%d duplicate names; pass `unique = TRUE` to disambiguate.", n),
    class = "genevintage_duplicate_names"
  )
}

#' Add a suffix to duplicate names, then make all names unique
#'
#' Only resolved entries (`ok`) in a collision receive the `tag` suffix.
#' [make.unique()] handles the remaining collisions and the unmapped IDs.
#' @noRd
.suffix_collisions <- function(out, ok, tag) {
  coll <- (duplicated(out) | duplicated(out, fromLast = TRUE)) & ok
  out[coll] <- paste0(out[coll], "_", tag[coll])
  make.unique(out, sep = "_")
}

#' Convert gene identifiers to gene names
#'
#' This function returns one value for each input, in the order of the input. It
#' returns an identifier with no name in `mapping` unchanged, so you can
#' use the result as row labels.
#'
#' The function ignores version suffixes when it matches.
#' `ENSG00000141510.16` resolves against an unversioned mapping. A versioned
#' mapping also resolves an unversioned identifier.
#'
#' @param ids Character vector of gene identifiers.
#' @param mapping A data frame with `id` and `name` columns, or a named
#'   character vector of names indexed by identifier.
#' @param unique Make the result unique. The function appends the identifier to
#'   names that otherwise collide. This is off by default, because it changes
#'   names.
#' @param warn Give a warning when names collide or when little of the input
#'   mapped. The collision warning has class `genevintage_duplicate_names`.
#'   A caller that disambiguates downstream can catch and muffle only that
#'   warning.
#'
#' @return A character vector with the same length as `ids`, with a `mapped`
#'   attribute. This attribute is the logical vector of the inputs that found a
#'   real name.
#'
#' @seealso [fetch_mapping()] to get `mapping` for a species and release.
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
    out <- .suffix_collisions(out, mapped, ids)
  }

  if (warn) {
    n_dup <- sum(duplicated(out))
    if (n_dup > 0 && !unique) {
      warning(.duplicate_names_condition(n_dup))
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
