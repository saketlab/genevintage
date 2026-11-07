# Cache filenames use "-" between fields so names containing "_" stay intact.
#
#   mapping-<source>-<species>-<assembly>-<release>.rds
#   lengths-<source>-<species>-<assembly>-<release>.rds
#   transcripts-<source>-<species>-<assembly>-<release>.rds
#   orthologs-<species>-<to>-<release>-<homology>.rds
#   xrefs-<db>-<species>-<assembly>-<release>.rds
.CACHE_FIELDS <- list(
  mapping = c("source", "species", "assembly", "release"),
  lengths = c("source", "species", "assembly", "release"),
  transcripts = c("source", "species", "assembly", "release"),
  orthologs = c("species", "to", "release", "homology"),
  xrefs = c("db", "species", "assembly", "release")
)

.cache_path <- function(kind, ...) {
  file.path(.cache_dir(), paste0(paste(c(kind, ...), collapse = "-"), ".rds"))
}

# .tmp marks an interrupted write. Only recognised kinds with the expected
# field count are listed and eligible for removal.
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

#' Show or clear the download cache
#'
#' The package streams mapping tables and ortholog pairs once. It caches them
#' under `tools::R_user_dir("genevintage", "cache")`, so a species that you
#' used once works offline after that. This function reports the contents of the
#' cache and can remove them.
#'
#' To move the location, set the `R_USER_CACHE_DIR` environment variable
#' before you use the package.
#'
#' An interrupted write leaves a `.rds.tmp` file. The function lists these files
#' as `partial`, and you can remove them.
#'
#' @param clear Remove cached files. `TRUE` removes all of them. A number
#'   removes the files with a modification time older than that number of days.
#' @return A data frame with one row for each cached file (`file`, `kind`,
#'   `species`, `to`, `release`, `bytes`, `modified`). The function returns it
#'   invisibly when `clear` is not `FALSE`. `to` is the second species of an
#'   ortholog pair and is `NA` for a mapping. The data frame has zero rows when
#'   the cache is empty. When `clear` is not `FALSE`, it has only the removed
#'   files.
#'
#' @examples
#' # a scratch cache, so the example leaves yours alone
#' old <- Sys.getenv("R_USER_CACHE_DIR")
#' Sys.setenv(R_USER_CACHE_DIR = tempfile())
#' gene_cache() # what is cached, and how large
#' gene_cache(clear = 90) # drop files older than 90 days
#' gene_cache(clear = TRUE) # drop all of it
#' Sys.setenv(R_USER_CACHE_DIR = old)
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

  # -Inf makes every finite file age eligible for removal
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
    for (memo in list(
      .mapping_memo, .ortholog_memo, .xref_memo, .lengths_memo,
      .transcript_memo
    )) {
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
