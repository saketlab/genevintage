# Species detection from gene symbols using a release-aware symbol index.

#' Build a compact species/release gene-symbol index from GTF mappings
#'
#' @param mappings A named list of mapping data frames returned by
#'   [fetch_mapping()], or one such data frame with `species` and `release`
#'   columns.
#' @param min_release Optional minimum release number.
#' @return A data frame with `species`, `release`, `assembly`, `source`, and
#'   `symbol`.
#' @export
build_symbol_index <- function(mappings, min_release = NULL) {
  if (is.data.frame(mappings)) mappings <- list(mappings)
  if (!is.list(mappings) || !length(mappings)) stop("mappings must be a non-empty list or data frame.", call. = FALSE)
  out <- lapply(mappings, function(m) {
    if (!is.data.frame(m) || !all(c("name") %in% names(m))) stop("Each mapping needs a `name` column.", call. = FALSE)
    n <- nrow(m)
    col <- function(k, default = NA_character_) {
      if (!k %in% names(m)) {
        return(rep(default, n))
      }
      as.character(m[[k]])
    }
    data.frame(
      species = col("species"), release = col("release"),
      assembly = col("assembly"), source = col("source", "ensembl"),
      symbol = toupper(trimws(as.character(m$name))), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[!is.na(out$symbol) & nzchar(out$symbol), , drop = FALSE]
  if (!is.null(min_release)) {
    rn <- .rel_num(out$release)
    out <- out[is.na(rn) | rn >= min_release, , drop = FALSE]
  }
  unique(out)
}


#' Fetch GTF mappings and build a symbol index
#'
#' @param species Species names accepted by [resolve_species()].
#' @param releases Ensembl/GENCODE release numbers to fetch.
#' @param source Annotation family.
#' @return A symbol index suitable for [detect_species_names()].
#' @export
build_symbol_index_from_gtf <- function(species, releases, source = "ensembl") {
  species <- vapply(species, resolve_species, character(1), quiet = TRUE)
  maps <- list()
  for (sp in species) {
    for (rel in releases) {
      m <- tryCatch(fetch_mapping(sp, rel, source = source), error = function(e) NULL)
      if (!is.null(m)) {
        m$species <- sp
        m$release <- as.character(rel)
        m$source <- source
        maps[[paste(sp, rel, source, sep = "|")]] <- m
      }
    }
  }
  if (!length(maps)) stop("No requested GTF mappings could be fetched.", call. = FALSE)
  build_symbol_index(maps)
}

#' Detect species from gene symbols
#'
#' Scores each indexed species/release by informative symbol overlap. Common
#' symbols are downweighted and Excel-prone date-like symbols are excluded.
#' Results below the overlap or score threshold have `status = "insufficient"`;
#' remaining results are `"ambiguous"` unless the top candidate meets the margin.
#' @param symbols Character vector of gene names.
#' @param index A symbol index from [build_symbol_index()].
#' @param min_overlap Minimum informative symbols required.
#' @param min_score Minimum weighted overlap fraction.
#' @param min_margin Required score gap between the top two candidates.
#' @return Ranked data frame with `species`, `release`, `score`, `overlap`,
#'   `margin`, and `status`.
#' @export
detect_species_names <- function(symbols, index, min_overlap = 20L,
                                 min_score = 0.25, min_margin = 0.05) {
  symbols <- toupper(trimws(as.character(symbols)))
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!is.data.frame(index) || !all(c("species", "symbol") %in% names(index))) stop("index must contain species and symbol columns.", call. = FALSE)
  # Date-like tokens are valid symbols, but are too error-prone to provide
  # positive evidence after spreadsheet conversion (e.g. SEPT1/MARCH1).
  excel <- grepl(paste0("^(", paste(unlist(EXCEL_MONTHS), collapse = "|"), ")[0-9]+$"), symbols)
  informative <- symbols[!excel]
  if (!length(informative)) {
    return(data.frame())
  }
  freq <- tabulate(match(index$symbol, informative), nbins = length(informative))
  w <- 1 / log1p(freq)
  w[!is.finite(w)] <- 1
  names(w) <- informative
  sum_w <- sum(w)
  for (col in c("release", "assembly", "source")) {
    if (!col %in% names(index)) index[[col]] <- NA_character_
  }
  key <- paste(index$species, index$release, index$assembly, index$source, sep = "\r")
  symbols_by_group <- split(index$symbol, key)
  meta <- index[match(names(symbols_by_group), key), c("species", "release", "assembly", "source")]
  rows <- Map(function(syms, sp, rel, asm, src) {
    hit <- intersect(informative, syms)
    data.frame(
      species = sp, release = rel, assembly = asm, source = src,
      overlap = length(hit), score = sum(w[hit]) / sum_w, stringsAsFactors = FALSE
    )
  }, symbols_by_group, meta$species, meta$release, meta$assembly, meta$source)
  out <- do.call(rbind, rows)
  out <- out[order(-out$score, -out$overlap, out$species), , drop = FALSE]
  if (!nrow(out)) {
    return(out)
  }
  out$margin <- c(if (nrow(out) > 1) out$score[[1]] - out$score[[2]] else out$score[[1]], rep(NA_real_, nrow(out) - 1))
  out$status <- ifelse(out$overlap < min_overlap | out$score < min_score, "insufficient",
    ifelse(seq_len(nrow(out)) == 1 & out$margin >= min_margin, "selected", "ambiguous")
  )
  rownames(out) <- NULL
  out
}
