# ID parsing and scheme detection. Everything downstream takes the parsed frame.

ENS_RX <- "^(ENS[A-Z]{0,4}?)([EGTPR])(\\d{11})(?:\\.(\\d+))?(_[A-Z0-9_]+)?$"

#' Identifier schemes, most specific first
#' @noRd
SCHEMES <- c(
  ensembl  = ENS_RX,
  flybase  = "^FBgn\\d{7}$",
  wormbase = "^WBGene\\d{8}$",
  tair     = "^AT[1-5CM]G\\d{5}(\\.\\d+)?$",
  sgd      = "^(Y[A-P][LR]\\d{3}[WC](-[A-Z])?|Q0\\d{3,4}|R\\d{4}[WC])$",
  xenbase  = "^Xetrov\\d+.*$",
  refseq   = "^[NX][MRP]_\\d+(\\.\\d+)?$",
  locus    = "^LOC\\d+$",
  entrez   = "^\\d+$"
)

#' Identifiers as a character vector, or a clear refusal
#'
#' Three entry points took the same argument three ways and reported the same
#' problem with three messages.
#' @noRd
.as_ids <- function(ids) {
  if (is.character(ids)) {
    return(ids)
  }
  if (is.factor(ids)) {
    return(as.character(ids))
  }
  stop("`ids` must be a character vector.", call. = FALSE)
}

#' Whether each element looks like a gene identifier of any known scheme
#'
#' Per element, not per vector: a caller may hand over identifiers and gene
#' names mixed together, and deciding for the whole vector loses one or the
#' other silently.
#' @noRd
.looks_like_id <- function(x) grepl(.ANY_SCHEME, .bare(x), perl = TRUE)

# The schemes as one alternation. Nine separate grepl passes over a matrix's
# row names cost three times as much as a single anchored one.
.ANY_SCHEME <- paste0("^(?:", paste(sub("\\$$", "", sub("^\\^", "", SCHEMES)),
  collapse = "|"
), ")$")

#' The numeric part of an identifier
#'
#' One definition, used both when the index is built and when input is scored
#' against it. The `qmax` feasibility test compares the two, so they must not
#' drift apart.
#' @noRd
.id_number <- function(x, bare = .bare(x)) {
  b <- bare
  out <- rep(NA_real_, length(b))
  # digits only; si.dkey.219e24 would otherwise read as 2.19e26
  ok <- grepl("^[^0-9]+[0-9]+$", b)
  out[ok] <- as.numeric(sub("^[^0-9]+", "", b[ok]))
  out
}

#' The version suffix of an identifier, or NA
#'
#' The scalar counterpart of [parse_ens()]'s `id_version` column, for callers
#' that want only the versions and not a full parse frame.
#' @noRd
.id_version <- function(x) {
  x <- .normalize_id(x)
  v <- rep(NA_integer_, length(x))
  hit <- grepl("\\.\\d+(_[A-Z0-9_]+)?$", x)
  v[hit] <- suppressWarnings(as.integer(sub("^.*\\.(\\d+)(_[A-Z0-9_]+)?$", "\\1", x[hit])))
  v
}

#' Whether a scheme's identifiers are assigned in increasing numeric order
#'
#' Exactly "`.id_number()` yields a usable number", so the flag and the
#' extraction cannot disagree. TAIR (`AT1G01010`) and SGD (`YDR387C`) embed the
#' chromosome, so neither yields one and `qmax` tells you nothing about them.
#' @noRd
.numeric_ordered <- function(ids) mean(!is.na(.id_number(ids))) > 0.9

#' Strip pipeline-added decoration from an identifier
#'
#' Some pipelines paste the symbol onto the ID (`ENSG00000141510!TP53`,
#' `ENSMMUG00000000001__CCNF`). `_PAR_Y` is a real Ensembl tag, not decoration.
#' @noRd
.normalize_id <- function(x) {
  x <- trimws(x)
  # perl: these run on every id of every call, and TRE is about twice the cost
  x <- sub("[!|].*$", "", x, perl = TRUE)
  sub("__.*$", "", x, perl = TRUE)
}

#' First non-missing of two vectors, element by element
#'
#' Not `ifelse()`: that takes its type and length from the test, so it returns
#' `logical(0)` for empty input where the type of `x` is meant. Seeding with the
#' fallback and overwriting the hits is type-stable by construction.
#' @noRd
.coalesce <- function(x, y) {
  na <- is.na(x)
  x[na] <- y[na]
  x
}

# whole-object fallback; .coalesce() is the element-wise one. base R has this
# from 4.4, and Depends is 4.1.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Bare identifier: decoration and any version suffix removed
#' @noRd
.bare <- function(x) sub("\\.\\d+(_[A-Z0-9_]+)?$", "", .normalize_id(x), perl = TRUE)

#' Which identifier scheme a vector of IDs uses
#'
#' @param ids Character vector of identifiers.
#' @param min_purity Smallest matching fraction to accept a scheme.
#' @return A list with `scheme`, `purity` and the per-scheme fractions.
#' @examples
#' guess_scheme(c("ENSG00000141510", "ENSG00000012048"))$scheme
#' guess_scheme(c("TP53", "BRCA1"))$scheme
#'
#' # purity is how a mixed vector announces itself
#' guess_scheme(c("ENSG00000141510", "TP53"))$purity
#' @export
guess_scheme <- function(ids, min_purity = 0.5) {
  ids <- .normalize_id(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  # no ids means no scheme; mean() over nothing is NaN and which.max() empty
  if (!length(ids)) {
    return(list(
      scheme = "symbol_or_mixed", purity = 0,
      all = stats::setNames(numeric(length(SCHEMES)), names(SCHEMES))
    ))
  }
  frac <- vapply(SCHEMES, function(rx) mean(grepl(rx, ids)), numeric(1))
  best <- which.max(frac)
  list(
    scheme = if (frac[best] >= min_purity) names(SCHEMES)[best] else "symbol_or_mixed",
    purity = unname(frac[best]), all = frac
  )
}

#' Parse Ensembl identifiers into their parts
#'
#' Handles unversioned and versioned IDs in one vector, `_PAR_Y`-style tags,
#' and non-Ensembl input (returned as `NA` rows).
#'
#' @param x Character vector of identifiers.
#' @return A data frame with one row per input: `id`, `prefix`, `feature`,
#'   `num`, `id_version`, `tag`.
#' @examples
#' parse_ens(c(
#'   "ENSG00000141510.16", "ENSMUSG00000059552",
#'   "ENSG00000182378.14_PAR_Y", "TP53"
#' ))
#' @export
parse_ens <- function(x) {
  x <- .normalize_id(x)
  # regexpr(perl = TRUE) exposes capture groups directly, with no per-element list
  m <- regexpr(ENS_RX, x, perl = TRUE)
  ok <- !is.na(m) & m != -1L
  n <- length(x)
  out <- data.frame(
    id = x,
    prefix = rep(NA_character_, n), feature = rep(NA_character_, n),
    num = rep(NA_real_, n), id_version = rep(NA_integer_, n),
    tag = rep(NA_character_, n), stringsAsFactors = FALSE
  )
  if (!any(ok)) {
    return(out)
  }

  st <- attr(m, "capture.start")[ok, , drop = FALSE]
  len <- attr(m, "capture.length")[ok, , drop = FALSE]
  grp <- function(i) {
    g <- substring(x[ok], st[, i], st[, i] + len[, i] - 1L)
    ifelse(len[, i] == 0L, NA_character_, g)
  }
  out$prefix[ok] <- grp(1)
  out$feature[ok] <- grp(2)
  out$num[ok] <- as.numeric(grp(3))
  out$id_version[ok] <- suppressWarnings(as.integer(grp(4)))
  out$tag[ok] <- grp(5)
  out
}
