# spellings Excel reads as a month
EXCEL_MONTHS <- list(
  jan = c("JAN"), feb = c("FEB"), mar = c("MARCH", "MARC", "MAR"),
  apr = c("APR"), may = c("MAY"), jun = c("JUN"), jul = c("JUL"),
  aug = c("AUG"), sep = c("SEPT", "SEP"), oct = c("OCT"),
  nov = c("NOV"), dec = c("DEC")
)
.EXCEL_SYMBOL_RX <- paste0("^(", paste(unlist(EXCEL_MONTHS), collapse = "|"), ")[0-9]{1,2}$")

# HGNC, MGI and ZFIN renames, per species; data-raw/build_excel_renames.R
.excel_renames <- .once(function() {
  utils::read.delim(.extdata("excel_renames.tsv"), comment.char = "#", colClasses = "character")
})

# a locale month can be a typed symbol (ago-01 from AGO1)
.LOCALE_MONTHS <- c(
  maart = 3L, mars = 3L, mrt = 3L, mai = 5L, mei = 5L, agosto = 8L, ago = 8L,
  oktober = 10L, okt = 10L, dez = 12L, dic = 12L
)
.LOCALE_LANGUAGES <- c(
  maart = "nl", mrt = "nl", mei = "nl", mars = "fr,sv,nb", mai = "fr,de,nb,da",
  agosto = "es,it,pt", ago = "es,it,pt", oktober = "de,nl,sv,nb,da",
  okt = "de,nl,sv,nb,da", dez = "de,pt", dic = "es,it"
)
.ENGLISH_MONTHS <- c(
  jan = 1L, feb = 2L, march = 3L, mar = 3L, apr = 4L, may = 5L, june = 6L, jun = 6L,
  july = 7L, jul = 7L, aug = 8L, september = 9L, sept = 9L, sep = 9L, oct = 10L,
  nov = 11L, dec = 12L
)
.MONTHS <- c(.LOCALE_MONTHS, .ENGLISH_MONTHS)
.MONTH_RX <- paste(names(.MONTHS), collapse = "|")
.ENGLISH_MONTH_RX <- paste(names(.ENGLISH_MONTHS), collapse = "|")
.LOCALE_MONTH_RX <- paste(names(.LOCALE_MONTHS), collapse = "|")

.SCI_RX <- "^([0-9])\\.([0-9]+)[Ee]\\+?([0-9]+)$"

# forced-text apostrophes and CSV quotes are not part of the symbol
.excel_clean <- function(x) {
  x <- gsub("^['\"]+|['\"]+$", "", trimws(x), perl = TRUE)
  # a datetime export adds midnight; any other time marks a real timestamp
  sub("(?i)[ T](0{1,2}:00(:00)?|12:00(:00)? ?AM)$", "", x, perl = TRUE)
}

# separators differ by locale: 1-Mar, 1/Mar, 01.Sep, Sep 15
.SEP <- "[-/. ]"

.mon <- function(s) unname(.MONTHS[tolower(s)])
.mdy <- function(p) list(c(.mon(p[2]), as.integer(p[3])))
.dmy <- function(p) list(c(.mon(p[3]), as.integer(p[2])))

DAMAGE <- list(
  list(
    form = "day_month", rx = paste0("^([0-9]{1,2})", .SEP, "(", .MONTH_RX, ")$"),
    md = .dmy
  ),
  list(
    form = "month_day", rx = paste0("^(", .ENGLISH_MONTH_RX, ")", .SEP, "([0-9]{1,2})$"),
    md = .mdy
  ),
  # Excel pads this day; unpadded mars.1 and AGO-1 are make.unique() labels
  list(
    form = "month_day", rx = paste0("^(", .LOCALE_MONTH_RX, ")", .SEP, "([0-9]{2})$"),
    md = .mdy
  ),
  list(form = "day_month_year", rx = paste0(
    "^([0-9]{1,2})", .SEP, "(", .MONTH_RX, ")",
    .SEP, "[0-9]{2,4}$"
  ), md = .dmy),
  list(
    form = "date_day", rx = "^(?:19|20)[0-9]{2}[-/.]([0-9]{1,2})[-/.]([0-9]{1,2})$",
    md = function(p) list(as.integer(p[2:3]))
  ),
  list(
    form = "numeric_date", rx = "^([0-9]{1,2})[-/.]([0-9]{1,2})[-/.][0-9]{2,4}$",
    md = function(p) {
      a <- as.integer(p[2])
      b <- as.integer(p[3])
      # locale unknown; offer both readings
      Filter(Negate(is.null), list(if (b <= 12) c(b, a), if (a <= 12) c(a, b)))
    }
  ),
  list(
    form = "serial", rx = "^([0-9]{5})$",
    md = function(p) {
      n <- as.integer(p[2])
      # serials <=60 predate or are Excel's phantom 29 Feb 1900
      if (n <= 60L) {
        return(list())
      }
      d <- as.Date(n, origin = "1899-12-30")
      list(c(as.integer(format(d, "%m")), as.integer(format(d, "%d"))))
    }
  ),
  list(form = "scientific", rx = .SCI_RX, md = function(p) list())
)

# cheap prefilter before the damage regexes
.MAYBE_RX <- "[-/. ]|^[0-9]+$|[0-9][Ee]"

# first matching DAMAGE rule wins; expects cleaned tokens
.excel_kind <- function(x) {
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
    hit[open][grepl(r$rx, todo[open], ignore.case = TRUE, perl = TRUE)] <- r$form
  }
  out[maybe] <- hit
  out
}

.local_words <- function(s) {
  intersect(strsplit(tolower(s), paste0(.SEP, "|[0-9]+"))[[1]], names(.LOCALE_MONTHS))
}

.languages <- function(x) {
  word <- vapply(x, function(s) .local_words(s)[1], "", USE.NAMES = FALSE)
  unname(.LOCALE_LANGUAGES[word])
}

# Entrez IDs and array coordinates hold no symbols; one token alone is an explicit request
.NUMERIC_RX <- "^(?=.*[0-9])[0-9./-]+$"
.numeric_list <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  length(x) > 1L && mean(grepl(.NUMERIC_RX, x, perl = TRUE)) > 0.5
}

# slash dates read both ways (locale unknown); annotation decides
.excel_dates <- function(x) {
  lapply(x, function(s) {
    for (r in DAMAGE) {
      p <- regmatches(s, regexec(r$rx, s, ignore.case = TRUE, perl = TRUE))[[1]]
      if (length(p)) {
        return(r$md(p))
      }
    }
    list()
  })
}

# Some Excel files store SEPT3 as the date 2003-09-01 and SEPT5 as 2005-09-01.
# In these dates, the gene number is in the year, not in the day.
# If a list has both dates, we read the gene number from the year.
# A list with only 2003-09-01 does not show this, so we read that date as SEPT1.
.month_year <- function(x) {
  rx <- "^20([0-9]{2})[-/.](0?[1-9]|1[0-2])[-/.]([0-9]{1,2})$"
  at <- grep(rx, x)
  year <- as.integer(sub(rx, "\\1", x[at]))
  month <- as.integer(sub(rx, "\\2", x[at]))
  day <- as.integer(sub(rx, "\\3", x[at]))
  # other days in the same year and month mean day-style dates, as in a merge
  span <- paste(year, month)
  one <- day == 1L & !span %in% span[day != 1L]
  if (any(tapply(year[one], month[one], function(v) length(unique(v)) > 1L))) {
    x[at[one]] <- sprintf("%d-%s", year[one], month.abb[month[one]])
  }
  x
}

# the year reading loses the month word; offer every locale's name
.local_names <- function(x) {
  lapply(.excel_dates(x), function(d) {
    sprintf("%s%d", names(.LOCALE_MONTHS)[.LOCALE_MONTHS == d[[1]][1]], d[[1]][2])
  })
}

# current symbols the species' authority gave these old ones
.renamed <- function(x, species) {
  ren <- .excel_renames()
  ren$symbol[ren$species == species & ren$previous %in% toupper(x)]
}

.excel_candidates <- function(x, species) {
  dates <- .excel_dates(x)
  lapply(seq_along(x), function(i) {
    p <- regmatches(x[i], regexec(.SCI_RX, x[i]))[[1]]
    if (length(p)) {
      digits <- paste0(p[2], p[3])
      e <- suppressWarnings(as.integer(p[4])) - (nchar(digits) - 1L)
      if (is.na(e) || e < 0) {
        return(character())
      }
      # RIKEN clone names; Ensembl suffixes them Rik
      stem <- paste0(digits, "E", formatC(e, width = 2, flag = "0"))
      return(c(stem, paste0(stem, "Rik")))
    }
    local <- .local_words(x[i])
    unlist(lapply(dates[[i]], function(md) {
      if (anyNA(md) || md[1] < 1 || md[1] > 12 || md[2] < 1 || md[2] > 31) {
        return(character())
      }
      cand <- paste0(EXCEL_MONTHS[[md[1]]], md[2])
      c(cand, sprintf("%s%d", local, md[2]), .renamed(cand, species))
    }))
  })
}

#' Repair gene symbols corrupted by a spreadsheet
#'
#' Excel changes some gene symbols: `SEPT9` becomes `9-Sep`, `MARCH1` becomes
#' `1-Mar`, and `2310009E13` becomes `2.31E+13`. Different symbols can give the
#' same token, so the function compares candidates with the species and
#' release that you supply. For example, `9-Sep` becomes `SEPT9` in Ensembl
#' release 95 and `SEPTIN9` in release 116.
#'
#' The function does not change a gene name that has no damage, so you can use the
#' result in `rownames()`. `MARCH1` and `MARC1` both become `1-Mar`. If a
#' candidate is present without damage in `x`, the function removes that
#' candidate.
#'
#' The function compares candidates with `mapping` if you supply it. If not, it
#' uses `species` at `release`. Without `release`, it uses the newest release,
#' which has the current HGNC names (`SEPTIN9`, not `SEPT9`). But if `x` has
#' stable IDs for at least half of its names, it finds the release from these
#' IDs. Without `species`, it finds the species from stable IDs in `x`, and it
#' refuses a vector with only symbols.
#'
#' The nomenclature authority of the species gives more renames: HGNC for human,
#' MGI for mouse, and ZFIN for zebrafish (`SEP15` to `SELENOF`, `sept9` to
#' `septin9a`). If you supply `mapping` without `species`, the function uses
#' human nomenclature.
#'
#' The function reads month names that Excel writes in other locales (`mars-02`,
#' `ago-01`). Such a token can also be a symbol that a person typed. For example,
#' `mars-02` can be `MARCHF2` or `MARS2`. The function reports it as ambiguous
#' unless one of the two is already present.
#'
#' For a matrix, sparse `Matrix`, or data frame, the function changes only the
#' row names. It keeps the dimensions, the row order, and other attributes. To
#' get the `corrected`, `ambiguous`, `unresolved`, and `skipped` attributes,
#' call `correct_genenames()` on `rownames(x)`.
#'
#' @inheritParams mito_genes
#' @param x Gene names, or a matrix, sparse `Matrix`, or data frame with gene
#'   names as row names. Some of these names can have spreadsheet damage.
#' @param unique Make the result unique. If a repaired name is the same as a
#'   different name, the function adds the damaged token to it, as
#'   [gene_names()] does. Then [make.unique()] changes any remaining duplicates.
#'   Defaults to `TRUE` when `x` is a matrix, sparse `Matrix`, or data frame, and
#'   to `FALSE` otherwise.
#' @param row_ids When `x` is a matrix, sparse `Matrix`, or data frame, the row
#'   identifiers to correct. Defaults to `rownames(x)`. Set this when `x` has no
#'   row names.
#' @return If `x` is a matrix, sparse `Matrix` or data frame: the same object
#'   with corrected row names and none of the attributes below.
#'
#'   If `x` is a plain vector: a character vector with the same length and order
#'   as `x`. The vector has these attributes:
#'   * `corrected`: logical.
#'   * `ambiguous`: tokens that match more than one symbol.
#'   * `unresolved`: damaged tokens that no symbol in this annotation can give.
#'   * `skipped`: tokens that look damaged, in a list that holds numbers.
#'   * `category`: the form of each damaged token, `NA` for other tokens. The
#'     values are `day_month` (`2-Sep`), `month_day` (`Sep-02`), `day_month_year`
#'     (`2-Sep-2019`), `numeric_date` (`9/2/2019`), `date_day` (`2020-09-02`),
#'     `date_year` (`2003-09-01`, where the year gives the gene number), `serial`
#'     (`42248`) and `scientific` (`2.31E+13`).
#'   * `locale`: the languages that use the month word of a token, for example
#'     `"nl"` for `mei-01`, `NA` if there is no such word. In the year reading,
#'     only a repaired gene shows the language.
#' @references
#' Ziemann M, Eren Y, El-Osta A (2016). Gene name errors are widespread in the
#' scientific literature. *Genome Biology* 17:177.
#'
#' Abeysooriya M, Soria M, Kasu MS, Ziemann M (2021). Gene name errors: Lessons
#' not learned. *PLOS Computational Biology* 17(7):e1008984.
#' @seealso [gene_names()], [gene_ids()], [geneid2name()].
#' @export
#' @examples
#' hs <- read.delim(system.file("extdata", "example_mapping.tsv", package = "genevintage"))
#' correct_genenames(c("9-Sep", "1-Dec", "TP53"), mapping = hs)
#' m <- matrix(1:4, nrow = 2, dimnames = list(c("9-Sep", "TP53"), c("s1", "s2")))
#' correct_genenames(m, mapping = hs)
#' @examplesIf interactive()
#' # release 95 resolves these as SEPT9 and DEC1
#' correct_genenames(c("9-Sep", "1-Dec", "TP53"), species = "human", release = 95)
#'
#' # release 116 resolves these as SEPTIN9 and DELEC1
#' correct_genenames(c("9-Sep", "1-Dec", "TP53"), species = "human", release = 116)
#'
#' # no release: the newest, so current spellings
#' correct_genenames(c("9-Sep", "1-Dec"), species = "human")
#'
#' correct_genenames(c("1-Mar", "MTARC1"), species = "human", release = 116)
correct_genenames <- function(x, species = NULL, release = NULL, assembly = NULL,
                              source = NULL, mapping = NULL, unique, quiet = FALSE,
                              row_ids = NULL) {
  if (!is.null(dim(x))) {
    if (missing(unique)) unique <- TRUE
    return(.rename_rows(x, row_ids, function(rid) {
      correct_genenames(rid,
        species = species, release = release, assembly = assembly,
        source = source, mapping = mapping, unique = unique, quiet = quiet
      )
    }))
  }
  if (missing(unique)) unique <- FALSE
  x <- .as_ids(x)
  out <- x
  corrected <- rep(FALSE, length(x))
  ambiguous <- character()
  unresolved <- character()
  clean <- .excel_clean(x)
  kind <- .excel_kind(clean)
  damaged <- which(!is.na(kind))
  category <- kind
  locale <- rep(NA_character_, length(x))
  hit <- damaged
  numeric_list <- .numeric_list(clean)
  if (numeric_list) hit <- damaged[!grepl(.NUMERIC_RX, clean[damaged], perl = TRUE)]
  skipped <- x[setdiff(damaged, hit)]

  if (length(hit)) {
    loc <- .locate_annotation(x, species, release, assembly, source, quiet,
      mapping = mapping
    )
    known <- loc$map$name[!is.na(loc$map$name)]
    lower <- tolower(known)
    typed <- damaged[tolower(clean[damaged]) %in% lower] # worm ago-1 and sep-1 are real symbols
    hit <- setdiff(hit, typed)
    category[typed] <- NA_character_
    # a plain number that is no gene means serials are numbers too
    serial <- kind[hit] == "serial"
    plain <- clean[is.na(kind) & grepl("^[0-9]+$", clean)]
    if (any(serial) && any(!plain %in% lower)) {
      skipped <- c(skipped, x[hit[serial]])
      hit <- hit[!serial]
    }

    nomen <- loc$species %||% "homo_sapiens"
    read <- .month_year(clean)[hit]
    shifted <- read != clean[hit]
    cand <- .excel_candidates(read, nomen)
    cand[shifted] <- Map(c, cand[shifted], .local_names(read[shifted]))
    category[hit[shifted]] <- "date_year"
    locale[hit] <- .languages(clean[hit])
    # a year reading that names no gene falls back to the day
    if (any(shifted) && !any(tolower(unlist(cand[shifted])) %in% lower)) {
      cand[shifted] <- .excel_candidates(clean[hit][shifted], nomen)
    }
    named <- known[match(tolower(unlist(cand)), lower)]
    named <- split(named, factor(rep(seq_along(cand), lengths(cand)), seq_along(cand)))
    # a gene present undamaged, under either spelling, is not the one a token lost
    rest <- x[setdiff(seq_along(x), hit)]
    present <- intersect(tolower(unlist(named)), tolower(c(rest, .renamed(rest, nomen))))
    for (j in seq_along(hit)) {
      real <- base::unique(named[[j]][!is.na(named[[j]])])
      absent <- real[!tolower(real) %in% present]
      if (length(absent)) real <- absent
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

  if (!quiet && length(skipped)) {
    message(sprintf(
      "%d left alone: %s, so not read as gene names, e.g. %s.",
      length(skipped), if (numeric_list) "mostly numbers" else "the list has plain numbers",
      skipped[1]
    ))
  }

  if (unique) {
    out <- .suffix_collisions(out, corrected, x)
  }

  attr(out, "corrected") <- corrected
  attr(out, "ambiguous") <- ambiguous
  attr(out, "unresolved") <- unresolved
  attr(out, "skipped") <- skipped
  # the year reading lost the month word; take the language from the repaired gene
  year <- which(corrected & category %in% "date_year")
  locale[year] <- .languages(out[year])
  attr(out, "category") <- category
  attr(out, "locale") <- locale
  out
}
