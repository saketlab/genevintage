# Gene symbols Excel silently reinterprets as dates or numbers, and the way back.

# old and new HGNC spellings listed (SEPT1/SEPTIN1, MARCH1/MARCHF1, DEC1/DELEC1); annotation decides which applies
EXCEL_MONTHS <- list(
  jan = c("JAN"), feb = c("FEB"), mar = c("MARCH", "MARCHF", "MAR"),
  apr = c("APR"), may = c("MAY"), jun = c("JUN"), jul = c("JUL"),
  aug = c("AUG"), sep = c("SEPT", "SEPTIN", "SEP"), oct = c("OCT"),
  nov = c("NOV"), dec = c("DEC", "DELEC")
)

# no inner groups (each shifts every capture index); longest alternative first, or "mar" wins over "march"
.MONTH_RX <- paste("jan", "feb", "march|mar", "apr", "may", "june|jun",
  "july|jul", "aug", "september|sept|sep", "oct", "nov", "dec",
  sep = "|"
)

# Excel truncates sci notation to ~3 significant digits; unrecoverable, no candidates generated
.SCI_RX <- "^([0-9])\\.([0-9]+)[Ee]\\+?([0-9]+)$"

# leading apostrophe (forced-text) and CSV quoting don't belong to the symbol
.excel_clean <- function(x) {
  x <- gsub("^['\"]+|['\"]+$", "", trimws(x), perl = TRUE)
  # a datetime export carries a clock the symbol never had
  sub("[ T][0-9]{1,2}:[0-9]{2}(:[0-9]{2})?$", "", x, perl = TRUE)
}

# Separators differ by locale: 1-Mar, 1/Mar, 01.Sep, "Sep 15".
.SEP <- "[-/. ]"

# recognition + parsing grammar in one table, stated once rather than twice
.mon <- function(s) match(substr(tolower(s), 1, 3), tolower(month.abb))
.dmy <- function(p) list(c(.mon(p[3]), as.integer(p[2])))

DAMAGE <- list(
  list(
    kind = "date", rx = paste0("^([0-9]{1,2})", .SEP, "(", .MONTH_RX, ")$"),
    md = .dmy
  ),
  list(
    kind = "date", rx = paste0("^(", .MONTH_RX, ")", .SEP, "([0-9]{1,2})$"),
    md = function(p) list(c(.mon(p[2]), as.integer(p[3])))
  ),
  list(kind = "date", rx = paste0(
    "^([0-9]{1,2})", .SEP, "(", .MONTH_RX, ")",
    .SEP, "[0-9]{2,4}$"
  ), md = .dmy),
  list(
    kind = "date", rx = "^[0-9]{4}-([0-9]{2})-([0-9]{2})$",
    md = function(p) list(as.integer(p[2:3]))
  ),
  list(
    kind = "date", rx = "^([0-9]{1,2})[-/.]([0-9]{1,2})[-/.][0-9]{2,4}$",
    md = function(p) {
      a <- as.integer(p[2])
      b <- as.integer(p[3])
      # day/month and month/day are both plausible; offer each valid reading
      Filter(Negate(is.null), list(if (b <= 12) c(b, a), if (a <= 12) c(a, b)))
    }
  ),
  list(
    kind = "date_serial", rx = "^([0-9]{5})$",
    md = function(p) {
      n <- as.integer(p[2])
      # Excel's phantom Feb 29 1900 makes every serial <=60 off by one; not a plausible symbol
      if (n <= 60L) {
        return(list())
      }
      d <- as.Date(n, origin = "1899-12-30")
      list(c(as.integer(format(d, "%m")), as.integer(format(d, "%d"))))
    }
  ),
  list(kind = "sci", rx = .SCI_RX, md = function(p) list())
)

# cheap gate: no separator/digits/exponent rules out every shape (almost every row name)
.MAYBE_RX <- "[-/. ]|^[0-9]+$|[0-9][Ee]"

.excel_kind <- function(x) {
  x <- .excel_clean(x)
  out <- rep(NA_character_, length(x))
  maybe <- which(grepl(.MAYBE_RX, x, perl = TRUE))
  if (!length(maybe)) {
    return(out)
  }
  todo <- x[maybe]
  hit <- rep(NA_character_, length(todo))
  for (r in DAMAGE) {
    open <- is.na(hit)
    if (!any(open)) break
    hit[open][grepl(r$rx, todo[open], ignore.case = TRUE, perl = TRUE)] <- r$kind
  }
  out[maybe] <- hit
  out
}

# slash dates read both ways (locale unknown); annotation settles which
.excel_dates <- function(x) {
  lapply(.excel_clean(x), function(s) {
    for (r in DAMAGE) {
      p <- regmatches(s, regexec(r$rx, s, ignore.case = TRUE, perl = TRUE))[[1]]
      if (length(p)) {
        return(r$md(p))
      }
    }
    list()
  })
}

.excel_candidates <- function(x) {
  x <- .excel_clean(x) # match and extract from the same string
  dates <- .excel_dates(x)
  lapply(seq_along(x), function(i) {
    p <- regmatches(x[i], regexec(.SCI_RX, x[i]))[[1]]
    if (length(p)) {
      digits <- paste0(p[2], p[3])
      e <- suppressWarnings(as.integer(p[4])) - (nchar(digits) - 1L)
      if (is.na(e) || e < 0) {
        return(character())
      }
      # RIKEN clone names are the common casualty and Ensembl suffixes them Rik
      stem <- paste0(digits, "E", formatC(e, width = 2, flag = "0"))
      return(c(stem, paste0(stem, "Rik")))
    }
    unlist(lapply(dates[[i]], function(md) {
      # a day of 0 or 32 is not a date, so it is not this kind of damage either
      if (anyNA(md) || md[1] < 1 || md[1] > 12 || md[2] < 1 || md[2] > 31) {
        return(character())
      }
      paste0(EXCEL_MONTHS[[md[1]]], md[2])
    }))
  })
}

#' Repair gene symbols a spreadsheet has converted to dates or numbers
#'
#' Excel silently reinterprets some gene symbols: `SEPT9` becomes `9-Sep`,
#' `MARCH1` becomes `1-Mar`, `2310009E13` becomes `2.31E+13`. The damage is
#' widespread in published supplementary files (Ziemann et al. 2016;
#' Abeysooriya et al. 2021) and is not recoverable from the file alone, because
#' several symbols can produce the same date.
#'
#' The annotation resolves it. Every symbol that could have produced the token
#' is generated, and only those the caller's own species and release actually
#' contain are offered, so the answer is the symbol their matrix should hold.
#'
#' That matters more than it sounds: HGNC renamed the worst offenders in 2020,
#' so `9-Sep` is `SEPT9` against Ensembl release 95 and `SEPTIN9` against
#' release 116. A fixed lookup table cannot be right for both.
#'
#' Tokens that are not damaged are returned unchanged, so the result drops
#' straight back into `rownames()`. A token with no surviving candidate, or with
#' more than one, is also returned unchanged and reported.
#'
#' A damaged token carries no species of its own, so pass `mapping`, or both
#' `species` and `release`, whenever `x` contains one. Choosing the release
#' explicitly also decides historical names such as `SEPT9` versus `SEPTIN9`;
#' left to detect, a symbol-only vector would silently default to the newest
#' release instead of refusing to guess.
#'
#' @inheritParams mito_genes
#' @param x Gene names, some of which a spreadsheet may have converted.
#' @param unique Make the result unique by appending the damaged token to a
#'   repaired name that collides with another (as [gene_names()] does), then
#'   falling back to [make.unique()] for any collision that remains. Off by
#'   default, since it changes names.
#' @return A character vector the same length and order as `x`, carrying a
#'   `corrected` attribute (logical), an `ambiguous` attribute naming tokens
#'   that matched more than one symbol, and an `unresolved` attribute naming
#'   damaged tokens no symbol in this annotation produces.
#' @references
#' Ziemann M, Eren Y, El-Osta A (2016). Gene name errors are widespread in the
#' scientific literature. *Genome Biology* 17:177.
#'
#' Abeysooriya M, Soria M, Kasu MS, Ziemann M (2021). Gene name errors: Lessons
#' not learned. *PLOS Computational Biology* 17(7):e1008984.
#' @seealso [gene_names()], [gene_ids()].
#' @export
#' @examples
#' \dontrun{
#' # against a 2019 vintage these come back as SEPT9 and MARCH1
#' correct_names(c("9-Sep", "1-Mar", "TP53"), species = "human", release = 95)
#'
#' # against a current one, as SEPTIN9 and MARCHF1
#' correct_names(c("9-Sep", "1-Mar", "TP53"), species = "human", release = 116)
#' }
correct_names <- function(x, species = NULL, release = NULL, assembly = NULL,
                          source = NULL, mapping = NULL, unique = FALSE,
                          quiet = FALSE) {
  x <- .as_ids(x)
  out <- x
  corrected <- rep(FALSE, length(x))
  ambiguous <- character()
  unresolved <- character()
  hit <- which(!is.na(.excel_kind(x)))

  if (length(hit)) {
    if (is.null(mapping) && (is.null(species) || is.null(release))) {
      stop("gene symbols alone cannot fix Excel damage: pass `mapping`, or both `species` and `release`.",
        call. = FALSE
      )
    }

    loc <- .locate_annotation(x, species, release, assembly, source, quiet,
      mapping = mapping
    )
    known <- loc$map$name[!is.na(loc$map$name)]
    lower <- tolower(known)

    cand <- .excel_candidates(x[hit])
    for (j in seq_along(hit)) {
      real <- base::unique(known[match(tolower(cand[[j]]), lower, nomatch = 0L)])
      if (length(real) == 1L) {
        out[hit[j]] <- real
        corrected[hit[j]] <- TRUE
      } else if (length(real) > 1L) {
        ambiguous <- c(ambiguous, x[hit[j]])
      } else {
        unresolved <- c(unresolved, x[hit[j]])
      }
    }

    if (!quiet) {
      if (any(corrected)) {
        message(sprintf(
          "repaired %d of %d damaged name%s.", sum(corrected),
          length(hit), if (length(hit) == 1L) "" else "s"
        ))
      }
      if (length(unresolved)) {
        message(sprintf(
          "%d left alone: no symbol in this annotation produces %s.",
          length(unresolved), unresolved[1]
        ))
      }
      if (length(ambiguous)) {
        warning(sprintf(
          "%d token(s) match more than one symbol and were left alone, e.g. %s.",
          length(ambiguous), ambiguous[1]
        ), call. = FALSE)
      }
    }
  }

  if (unique) {
    out <- .suffix_collisions(out, corrected, x)
  }

  attr(out, "corrected") <- corrected
  attr(out, "ambiguous") <- ambiguous
  attr(out, "unresolved") <- unresolved
  out
}

#' Repair Excel-damaged row identifiers in place
#'
#' [correct_names()] applied to a matrix's row names, and put back on the same
#' object. Dimensions never change and rows never move: a token with no
#' surviving candidate, or with more than one, keeps its original value.
#'
#' Works on anything with `rownames()`: a base matrix, a sparse `Matrix` or a
#' data frame, because only the dimnames are touched.
#'
#' @param x A matrix or data frame with (possibly Excel-damaged) gene symbols
#'   as row names.
#' @param ids Row identifiers. Defaults to `rownames(x)`.
#' @param unique Disambiguate names that would otherwise collide, by appending
#'   the damaged token. On by default: duplicate row names are legal in a
#'   matrix but break almost everything downstream.
#' @inheritParams correct_names
#'
#' @return `x` with new row names. Nothing else about the object changes: no
#'   extra attributes, so a sparse `Matrix` stays a valid S4 object through
#'   subsetting, binding and arithmetic. Keep `rownames(x)` yourself if you
#'   need to see which tokens were ambiguous or unresolved; call
#'   [correct_names()] directly for that.
#'
#' @seealso [correct_names()] for the vector form and its attributes,
#'   [rename_rows()] for the ID-to-name case.
#' @export
#' @examples
#' \dontrun{
#' m <- matrix(1:4, nrow = 2, dimnames = list(c("9-Sep", "TP53"), c("s1", "s2")))
#' correct_rows(m, species = "human", release = 116)
#' }
correct_rows <- function(x, ids = rownames(x), species = NULL, release = NULL,
                         assembly = NULL, source = NULL, mapping = NULL,
                         unique = TRUE, quiet = FALSE) {
  .check_row_ids(x, ids)
  nm <- correct_names(ids,
    species = species, release = release, assembly = assembly,
    source = source, mapping = mapping, unique = unique, quiet = quiet
  )
  rownames(x) <- as.vector(nm)
  x
}
