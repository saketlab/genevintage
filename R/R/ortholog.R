# Orthologs, from Ensembl Compara. Same shape as the mapping path: stream, cache
# the small result, convert offline.

.ortholog_memo <- new.env(parent = emptyenv())

# select the primary assembly when a release carries multiple assemblies
.mapping_for <- function(species, release) {
  row <- .primary_row(species, release)
  fetch_mapping(species, release,
    assembly = if (nrow(row)) row$assembly else NULL
  )
}

.named <- function(ids, species, release) {
  as.vector(gene_names(ids, .mapping_for(species, release), warn = FALSE))
}

# Gene names back to stable identifiers
#
# a symbol can name several genes, so keep every match and warn the caller
.ids_from_names <- function(x, species, release) {
  m <- .mapping_for(species, release)
  lut <- m[.is_name(m$name, m$id), c("id", "name"), drop = FALSE]
  lut$id <- .bare(lut$id)

  # split only matching rows to limit lookup work to the requested names
  keep <- lut$name %in% x
  nm <- lut$name[keep]
  hits <- split(lut$id[keep], factor(nm, levels = unique(nm)))[x]
  # try case-insensitive matching for symbols with no exact match
  miss <- vapply(hits, is.null, logical(1))
  if (any(miss)) {
    lx <- tolower(x[miss])
    lname <- tolower(lut$name)
    lk <- lname %in% lx
    lower <- lname[lk]
    lo <- split(lut$id[lk], factor(lower, levels = unique(lower)))
    hits[miss] <- lo[lx]
  }
  n <- lengths(hits)
  if (any(n == 0)) {
    warning(sprintf(
      "%d name(s) matched no gene in %s release %s, e.g. %s.",
      sum(n == 0), species, release, x[n == 0][1]
    ), call. = FALSE)
  }
  if (any(n > 1)) {
    k <- which.max(n)
    warning(sprintf(
      "%d name(s) match more than one gene (%s matches %d); all are included.",
      sum(n > 1), x[k], n[k]
    ), call. = FALSE)
  }
  unlist(hits, use.names = FALSE)
}

# retain the homology columns in empty results so callers can access them
.empty_homology <- function(names = FALSE) {
  out <- data.frame(
    id = character(), name = character(),
    ortholog_id = character(), ortholog_name = character(),
    homology_type = character(), identity = numeric(),
    ortholog_identity = numeric(), high_confidence = logical(),
    stringsAsFactors = FALSE
  )
  if (!names) out[c("name", "ortholog_name")] <- NULL
  out
}

# derive column order from the empty result schema
.HOM_COLS <- names(.empty_homology(TRUE))

# read from the other species' file, a one2many becomes a many2one
.swap_homology <- function(x) {
  ifelse(x == "ortholog_one2many", "ortholog_many2one",
    ifelse(x == "ortholog_many2one", "ortholog_one2many", x)
  )
}

#' Stream one species' orthologs against one other species
#'
#' Filters the Compara stream in chunks to retain orthologs for the requested
#' target species, limiting memory to the selected pair.
#'
#' @param url A Compara `homologies.tsv.gz` URL.
#' @param to Ensembl species name to keep, e.g. `"mus_musculus"`.
#' @param swap Reverse the pair when reading it from the target species' file,
#'   including the direction of one-to-many relationships.
#' @param chunk Lines to decompress per iteration.
#' @return A data frame of `id`, `ortholog_id`, `homology_type`, `identity`,
#'   `ortholog_identity`, `high_confidence`.
#' @examples
#' \dontrun{
#' url <- paste0(
#'   "https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/",
#'   "homologies/mus_musculus/",
#'   "Compara.116.protein_default.homologies.tsv.gz"
#' )
#' head(stream_homologies(url, to = "homo_sapiens"))
#' }
#' @export
stream_homologies <- function(url, to, swap = FALSE, chunk = 200000) {
  need <- c(
    "gene_stable_id", "homology_type", "homology_gene_stable_id",
    "homology_species", "identity", "homology_identity",
    "is_high_confidence"
  )
  # cheap fixed-string prefilter before the far more expensive split
  tag <- paste0("\t", to, "\t")
  j <- NULL
  nc <- NULL
  # read and remove the header from the first chunk
  k <- .stream_filter(url, function(x) {
    if (is.null(j)) {
      cols <- strsplit(x[1], "\t", fixed = TRUE)[[1]]
      miss <- setdiff(need, cols)
      if (length(miss)) {
        stop("homology file is missing column(s): ", paste(miss, collapse = ", "),
          call. = FALSE
        )
      }
      j <<- match(need, cols)
      nc <<- length(cols)
      x <- x[-1]
    }
    x[grepl(tag, x, fixed = TRUE)]
  }, chunk)
  if (is.null(j)) stop("empty homology file: ", url, call. = FALSE)
  if (!length(k)) {
    return(.empty_homology())
  }

  # reshape tab-separated fields into rows
  f <- matrix(unlist(strsplit(k, "\t", fixed = TRUE), use.names = FALSE),
    ncol = nc, byrow = TRUE
  )
  # the prefilter matches the species anywhere on the line
  f <- f[f[, j[4]] == to, , drop = FALSE]
  num <- function(v) suppressWarnings(as.numeric(replace(v, v == "NULL", NA)))
  # id, ortholog_id, identity, ortholog_identity; straight or reversed
  p <- if (swap) c(3L, 1L, 6L, 5L) else c(1L, 3L, 5L, 6L)
  out <- data.frame(
    id                = f[, j[p[1]]],
    ortholog_id       = f[, j[p[2]]],
    # a one2many read from the other side is a many2one
    homology_type     = if (swap) .swap_homology(f[, j[2]]) else f[, j[2]],
    identity          = num(f[, j[p[3]]]),
    ortholog_identity = num(f[, j[p[4]]]),
    high_confidence   = f[, j[7]] == "1",
    stringsAsFactors  = FALSE
  )
  # paralogs share the file; only orthology crosses species
  out[startsWith(out$homology_type, "ortholog"), , drop = FALSE]
}

# prefer an Ensembl row with a missing frozen_release value for Compara
.primary_row <- function(species, release) {
  fp <- .fp()
  row <- fp[fp$species == species & fp$release == as.character(release) &
    fp$source == "ensembl", ]
  if (!nrow(row)) {
    return(row)
  }
  row[order(!is.na(row$frozen_release)), ][1, ]
}

# Which file holds one pair
#
# try source collections, then target collections with the pair reversed,
# because Compara can record a pair in either species' files
.collections <- function(root, release, kind, species) {
  dir_url <- sprintf(
    "%s/release-%s/tsv/ensembl-compara/homologies/%s/",
    root, release, species
  )
  idx <- .listing(dir_url, sprintf(
    "Compara has no homologies for %s at release %s (%s). Older releases predate the homology dumps; try a newer one, and see ortholog_species().",
    species, release, dir_url
  ))
  # no release in the pattern; plants r63 ships Compara.116.protein_default...
  rx <- sprintf("^Compara\\.[0-9]+\\.%s_[a-z]+\\.homologies\\.tsv\\.gz$", kind)
  f <- grep(rx, idx, value = TRUE)
  if (!length(f)) {
    stop(sprintf(
      "Compara publishes no %s homologies for %s at release %s.",
      kind, species, release
    ), call. = FALSE)
  }
  # the default collection holds most pairs
  paste0(dir_url, f[order(!grepl("_default\\.", f))])
}

.homologies_for_pair <- function(root, release, kind, species, to) {
  # a species may appear in its partner's file even when its own directory is absent
  err <- character()
  side <- function(sp, swap) {
    cs <- tryCatch(.collections(root, release, kind, sp),
      error = function(e) {
        err <<- c(err, conditionMessage(e))
        NULL
      }
    )
    for (u in cs) {
      m <- stream_homologies(u, if (swap) species else to, swap = swap)
      if (nrow(m)) {
        return(m)
      }
    }
    NULL
  }
  m <- side(species, FALSE)
  if (!is.null(m)) {
    return(m)
  }
  m <- side(to, TRUE)
  if (!is.null(m)) {
    return(m)
  }
  # propagate the source error if both collection lookups fail
  if (length(err) == 2L) stop(err[1], call. = FALSE)
  # successful lookups with no matching pair return an empty homology table
  .empty_homology()
}

#' Species Compara publishes homologies for
#'
#' Lists the species with Compara files for a release and division so callers
#' can check availability before streaming a pair.
#'
#' @param release Ensembl release number. Under Ensembl Genomes, that division's
#'   own release number.
#' @param division Which division's Compara to list. Defaults to vertebrates.
#'
#' @return A character vector of species names.
#' @export
#' @examples
#' \dontrun{
#' "mus_musculus" %in% ortholog_species(116)
#' }
ortholog_species <- function(release, division = "vertebrates") {
  fp <- .fp()
  root <- unique(fp$root[fp$division == division & fp$source == "ensembl"])
  if (!length(root)) stop("unknown division: ", division, call. = FALSE)
  dir_url <- sprintf("%s/release-%s/tsv/ensembl-compara/homologies/", root[1], release)
  idx <- .listing(dir_url, sprintf(
    "Compara has no homologies at %s release %s.",
    division, release
  ))
  sort(grep("^[a-z]+_[a-z_]+$", idx, value = TRUE))
}

#' Ortholog table for one species pair and Ensembl release
#'
#' Streamed from Ensembl Compara on first use, then cached, so later calls are
#' offline. Compara is Ensembl's, so `release` is an Ensembl release even when
#' the identifiers came from GENCODE.
#'
#' @param species Ensembl species the identifiers come from.
#' @param to Species to find orthologs in, in any form [resolve_species()] accepts.
#' @param release Ensembl release number.
#' @param kind `"protein"` (the default) or `"ncrna"`: Compara publishes the two
#'   trees separately, and non-coding orthologs live only in the second.
#' @param refresh Re-download even when a cached copy exists.
#'
#' @return A data frame of `id`, `ortholog_id`, `homology_type`, `identity`,
#'   `ortholog_identity`, `high_confidence`.
#' @examples
#' \dontrun{
#' # first call per pair streams a Compara file and caches it; then offline
#' head(fetch_orthologs("human", to = "mouse", release = 116))
#' }
#' @export
fetch_orthologs <- function(species, to, release, kind = c("protein", "ncrna"),
                            refresh = FALSE) {
  kind <- match.arg(kind)
  species <- resolve_species(species)
  to <- resolve_species(to)
  if (identical(species, to)) {
    stop("`to` must differ from the source species.", call. = FALSE)
  }
  fp <- .fp()
  row <- .primary_row(species, release)
  if (!nrow(row)) {
    stop(sprintf("no Ensembl index row for %s release %s.", species, release),
      call. = FALSE
    )
  }
  trow <- fp[fp$species == to & fp$source == "ensembl", ]
  if (!nrow(trow)) {
    stop(sprintf("'%s' is not a species in the index.", to), call. = FALSE)
  }
  # require both species to belong to the same division's Compara
  if (!identical(row$division, trow$division[1])) {
    stop(sprintf(
      "%s (%s) and %s (%s) are in different Ensembl divisions; Compara does not relate them.",
      species, row$division, to, trow$division[1]
    ), call. = FALSE)
  }

  f <- .cache_path("orthologs", species, to, release, kind)
  if (!refresh) {
    hit <- .cache_get(f, .ortholog_memo, c("id", "ortholog_id"))
    if (!is.null(hit)) {
      return(hit)
    }
  }

  # use the division's root to locate Compara homologies
  .cache_put(
    f, .ortholog_memo,
    .homologies_for_pair(row$root, release, kind, species, to)
  )
}

#' Find orthologs of gene identifiers in another species
#'
#' Detects the species from identifiers and defaults to its newest indexed
#' Ensembl release. Fetches that release's Compara orthologs against `to`
#' (streamed once, then cached) and returns the orthologous identifiers with
#' their gene names.
#'
#' A gene can have several orthologs or none. With `one2one = TRUE` (the
#' default) only one-to-one orthologs are returned, so
#' there is at most one row per input identifier and the result can be used
#' directly as a lookup. Set it to `FALSE` to see one-to-many and many-to-many
#' relationships as well.
#'
#' @param ids Character vector of gene identifiers or gene names. A symbol can
#'   label multiple genes, so all matches are included. Supply `species` for
#'   input containing only names.
#' @param to Species to find orthologs in, in any form [resolve_species()]
#'   accepts: `"mus_musculus"`, `"Mouse"` or `"mouse"`.
#' @param species,release Skip detection by naming them.
#' @param one2one Keep only one-to-one orthologs. See details.
#' @param names Also look up gene names for both sides of the pair, which costs
#'   a streamed mapping for each species. Both are named from the same Ensembl
#'   release, so the two halves of a row describe one vintage.
#' @param kind Passed to [fetch_orthologs()].
#' @param quiet Suppress the message reporting what was detected.
#'
#' @return A data frame with one row per ortholog relationship: `id` and `name`
#'   for the query gene, `ortholog_id` and `ortholog_name` for its counterpart
#'   in `to` (the two name columns when `names`), then `homology_type`,
#'   `identity`, `ortholog_identity` and `high_confidence`. Identifiers with no
#'   ortholog in `to` are absent, so the result is usually shorter than `ids`.
#'
#' @seealso [geneid2name()] for names within one species.
#' @export
#' @examples
#' \dontrun{
#' orthologs(c("ENSG00000141510", "ENSG00000012048"), to = "mus_musculus")
#' }
orthologs <- function(ids, to, species = NULL, release = NULL, one2one = TRUE,
                      names = TRUE, kind = c("protein", "ncrna"), quiet = FALSE) {
  kind <- match.arg(kind)
  to <- resolve_species(to)
  if (!is.null(species)) species <- resolve_species(species)
  ids <- .as_ids(ids)
  # classify each element so mixed inputs retain identifiers and names;
  # input containing only names requires an explicit species
  is_id <- .looks_like_id(ids)
  by_name <- !all(is_id)
  if (by_name && !any(is_id) && is.null(species)) {
    stop("gene names do not say which species they belong to; pass `species`.",
      call. = FALSE
    )
  }
  if (is.null(species)) {
    cand <- detect_species(ids[is_id])
    if (!nrow(cand)) stop("cannot identify the species from these ids.", call. = FALSE)
    species <- cand$species[1]
  }
  if (is.null(release)) {
    # default to the newest indexed Ensembl release for Compara lookups
    release <- .newest_release(species)
    if (!quiet) message(sprintf("%s -> %s, Ensembl release %s", species, to, release))
  }

  if (by_name) {
    ids <- c(ids[is_id], .ids_from_names(ids[!is_id], species, release))
    if (!length(ids)) {
      return(.empty_homology(names))
    }
  }

  h <- fetch_orthologs(species, to, release, kind = kind)
  if (one2one) h <- h[h$homology_type == "ortholog_one2one", , drop = FALSE]

  # match bare identifiers against Compara's unversioned stable IDs
  bare <- .bare(ids)
  keep <- h[h$id %in% bare, , drop = FALSE]
  # report the caller's own identifiers back, versions and all. when the
  # input was names, the stable id is the useful thing to return, and the name
  # comes back in its own column anyway
  if (!by_name) keep$id <- ids[match(keep$id, bare)]

  if (names) {
    if (!nrow(keep)) {
      return(.empty_homology(TRUE))
    }
    # name both sides from the same Ensembl release; preserve IDs lacking symbols
    keep$name <- .named(keep$id, species, release)
    keep$ortholog_name <- .named(keep$ortholog_id, to, release)
    keep <- keep[, .HOM_COLS]
  }
  rownames(keep) <- NULL
  keep
}
