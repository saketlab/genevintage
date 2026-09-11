# Cached mappings and homologies live under tools::R_user_dir(). CRAN's policy
# for that location asks that the contents be actively managed, which needs a
# way to see and remove them.

# One place knows the cache layout, and both writers go through it.
#
# "-" separates fields: no source, species, assembly or release contains one,
# while "_" appears in two of them (homo_sapiens, gencode_all) and species run
# to three tokens. Separating on a character the fields can contain forced the
# reader to carry the source vocabulary and a field count.
#
#   mapping-<source>-<species>-<assembly>-<release>.rds
#   orthologs-<species>-<to>-<release>-<homology>.rds
#   xrefs-<db>-<species>-<assembly>-<release>.rds
.CACHE_FIELDS <- list(
  mapping = c("source", "species", "assembly", "release"),
  orthologs = c("species", "to", "release", "homology"),
  xrefs = c("db", "species", "assembly", "release")
)

.cache_path <- function(kind, ...) {
  file.path(.cache_dir(), paste0(paste(c(kind, ...), collapse = "-"), ".rds"))
}

# .tmp marks a write that did not finish. anything else is not ours and is
# dropped, so a file somebody left here is never reported nor deleted.
.cache_parts <- function(f) {
  p <- strsplit(sub("\\.rds(\\.tmp)?$", "", basename(f)), "-", fixed = TRUE)
  kind <- vapply(p, function(x) x[1L], "")
  ours <- kind %in% names(.CACHE_FIELDS) &
    lengths(p) == lengths(.CACHE_FIELDS)[match(kind, names(.CACHE_FIELDS))] + 1L
  ours[is.na(ours)] <- FALSE
  f <- f[ours]
  p <- p[ours]
  kind <- kind[ours]
  # a field by name, for the rows whose kind has one
  field <- function(nm) {
    vapply(seq_along(p), function(j) {
      i <- match(nm, .CACHE_FIELDS[[kind[j]]])
      if (is.na(i)) NA_character_ else p[[j]][i + 1L]
    }, "")
  }
  data.frame(
    file = basename(f), kind = kind, species = field("species"),
    to = field("to"), release = field("release"),
    path = f, stringsAsFactors = FALSE, row.names = NULL
  )
}

#' Inspect or clear the download cache
#'
#' Mapping tables and ortholog pairs are streamed once and cached under
#' `tools::R_user_dir("genevintage", "cache")`, so a species you have used once
#' works offline afterwards. This reports what is cached and removes it.
#'
#' The location follows R's own convention and can be moved by setting the
#' `R_USER_CACHE_DIR` environment variable before the package is used.
#'
#' An interrupted download leaves a `.rds.tmp` file behind. Those are listed as
#' `partial` and removed like anything else, so nothing the package writes is
#' invisible to the user or beyond their reach.
#'
#' @param clear Remove cached files. `TRUE` removes all of them; a number
#'   removes those not used in that many days.
#' @return A data frame with one row per cached file (`file`, `kind`, `species`,
#'   `to`, `release`, `bytes`, `modified`), invisibly when `clear` is not
#'   `FALSE`. `to` is the second species of an ortholog pair and
#'   is `NA` for a mapping. Zero rows when nothing is cached.
#'
#' @examples
#' gene_cache() # what is cached, and how large
#' \dontrun{
#' gene_cache(clear = 90) # drop anything untouched for 90 days
#' gene_cache(clear = TRUE) # drop all of it
#' }
#' @export
gene_cache <- function(clear = FALSE) {
  d <- .cache_dir()
  # list.files() returns character(0) for a directory that does not exist
  parts <- .cache_parts(list.files(d, pattern = "\\.rds(\\.tmp)?$", full.names = TRUE))
  info <- file.info(parts$path)
  out <- data.frame(parts[setdiff(names(parts), "path")],
    bytes = info$size, modified = info$mtime,
    stringsAsFactors = FALSE, row.names = NULL
  )
  out$kind[endsWith(out$file, ".tmp")] <- "partial"
  if (identical(clear, FALSE)) {
    return(out)
  }

  # "everything" is just a threshold no file can be younger than
  days <- if (isTRUE(clear)) {
    -Inf
  } else if (is.numeric(clear)) {
    clear
  } else {
    stop("`clear` must be TRUE, FALSE or a number of days.", call. = FALSE)
  }
  stale <- difftime(Sys.time(), info$mtime, units = "days") > days

  if (any(stale)) {
    gone <- parts$path[stale]
    unlink(gone)
    # the in-session memos hold the same tables by full path; leaving them would
    # keep serving what was just deleted
    for (memo in list(.mapping_memo, .ortholog_memo, .xref_memo)) {
      rm(list = intersect(gone, ls(memo)), envir = memo)
    }
    # leave no empty directory behind when the last file goes
    if (!length(list.files(d, all.files = TRUE, no.. = TRUE))) unlink(d, recursive = TRUE)
  }
  message(sprintf(
    "removed %d of %d cached file%s (%s).", sum(stale), nrow(out),
    if (nrow(out) == 1L) "" else "s",
    format(structure(sum(info$size[stale]), class = "object_size"),
      units = "auto"
    )
  ))
  invisible(out[stale, , drop = FALSE])
}
