# Cross-references from Ensembl TSV files under pub/release-N/tsv/.

# match the database suffix because filenames use the division's release number
.XREF_FILE <- "\\.%s\\.tsv\\.gz$"

#' Stream one cross-reference TSV and keep the gene-level pairs
#'
#' This function removes duplicates of gene-level pairs from transcript rows
#' during the stream. This limits memory to distinct pairs.
#'
#' @param url A `.tsv.gz` cross-reference URL.
#' @param db The `db_name` to keep, for example `"EntrezGene"`.
#' @param chunk Lines to decompress in each iteration.
#' @return A data frame of `id` and `xref`.
#' @export
#' @examplesIf interactive()
#' stream_xrefs(paste0(
#'   "https://ftp.ensembl.org/pub/release-116/tsv/saccharomyces_cerevisiae/",
#'   "Saccharomyces_cerevisiae.R64-1-1.63.entrez.tsv.gz"
#' ))
stream_xrefs <- function(url, db = "EntrezGene", chunk = 200000) {
  seen <- character()
  k <- .stream_filter(url, function(x) {
    x <- x[!startsWith(x, "gene_stable_id")]
    if (!length(x)) {
      return(character())
    }
    f <- strsplit(x, "\t", fixed = TRUE)
    col <- function(i) vapply(f, function(z) if (length(z) >= i) z[i] else NA_character_, "")
    xref <- col(4L)
    keep <- col(5L) == db & !is.na(xref) & nzchar(xref)
    if (!any(keep)) {
      return(character())
    }
    pair <- paste(col(1L)[keep], xref[keep], sep = "\t")
    pair <- pair[!duplicated(pair) & !pair %in% seen]
    seen <<- c(seen, pair)
    pair
  }, chunk)
  if (!length(k)) {
    return(data.frame(
      id = character(), xref = character(),
      stringsAsFactors = FALSE
    ))
  }
  f <- strsplit(k, "\t", fixed = TRUE)
  data.frame(
    id = vapply(f, `[`, "", 1L), xref = vapply(f, `[`, "", 2L),
    stringsAsFactors = FALSE
  )
}

#' Ensembl gene identifiers and their cross-references, cached
#'
#' This function streams the cross-reference TSV that Ensembl publishes for a
#' species and release. It keeps the gene-level pairs and caches them. Later
#' calls work offline.
#'
#' @param species Species in any form that [resolve_species()] accepts.
#' @param release Ensembl release number for the cross-reference files.
#' @param db Cross-reference database. `"entrez"` is the only indexed database.
#' @param refresh Fetch the table again and overwrite the cached copy.
#' @return A data frame of `id` and `xref`, one row for each distinct pair.
#' @seealso [entrez_ids()] for the lookup that uses this table.
#' @export
#' @examplesIf interactive()
#' head(fetch_xrefs("yeast", release = 116))
fetch_xrefs <- function(species, release, db = "entrez", refresh = FALSE) {
  species <- resolve_species(species)
  release <- as.character(release)
  fp <- .fp()
  row <- fp[fp$species == species & fp$source == "ensembl" &
    fp$release == release, , drop = FALSE]
  if (!nrow(row)) {
    stop(sprintf("no Ensembl index row for %s release %s.", species, release),
      call. = FALSE
    )
  }
  row <- .primary_row(species, release)

  f <- .cache_path("xrefs", db, species, row$assembly, release)
  if (!refresh) {
    hit <- .cache_get(f, .xref_memo, c("id", "xref"))
    if (!is.null(hit)) {
      return(hit)
    }
  }

  dir_url <- sprintf("%s/release-%s/tsv/%s/", row$annotation_root, release, species)
  idx <- .listing(dir_url, paste("cannot reach", dir_url))
  cand <- grep(sprintf(.XREF_FILE, db), idx, value = TRUE)
  if (!length(cand)) {
    stop(sprintf(
      "%s publishes no %s cross-references at release %s.",
      species, db, release
    ), call. = FALSE)
  }
  .cache_put(f, .xref_memo, stream_xrefs(paste0(dir_url, cand[1]), .XREF_DB[[db]]))
}

# what the db_name column calls each of them
.XREF_DB <- list(entrez = "EntrezGene")

.xref_memo <- new.env(parent = emptyenv())

#' NCBI Entrez identifiers for a set of genes
#'
#' This function takes gene identifiers or gene names and returns their Entrez
#' identifiers. It detects the species in the same way as the rest of the
#' package. It streams cross-references once for each species and release and
#' then caches them.
#'
#' A gene can have several Entrez identifiers or none. Each row is a pair.
#' `entrez` is `NA` where Ensembl records no cross-reference.
#'
#' Ensembl publishes cross-references against Ensembl releases. If identifiers
#' belong to a GENCODE release, the function looks them up in the newest indexed
#' Ensembl release for that species.
#'
#' @inheritParams mito_genes
#' @param ids Gene identifiers, or gene names. For names, you must give
#'   `species`, because a name does not tell which species it came from.
#' @return A data frame of `input`, `id`, `name` and `entrez`, one row for each
#'   pair.
#' @seealso [fetch_xrefs()], [gene_names()].
#' @export
#' @examplesIf interactive()
#' entrez_ids(c("ENSG00000141510", "ENSG00000012048"))
#' entrez_ids(c("TP53", "BRCA1"), species = "human")
entrez_ids <- function(ids, species = NULL, release = NULL, assembly = NULL,
                       source = NULL, mapping = NULL, quiet = FALSE) {
  loc <- .locate_annotation(ids, species, release, assembly, source, quiet,
    mapping = mapping
  )
  hit <- .restrict(loc$map, ids)
  if (!nrow(hit)) {
    return(data.frame(
      input = character(), id = character(), name = character(),
      entrez = character(), stringsAsFactors = FALSE
    ))
  }
  # a GENCODE vintage has no cross-reference file of its own, and Ensembl
  # publishes tsv/ only from release 85; stable ids persist into the newest
  rel <- if (identical(loc$source, "ensembl") && .rel_num(loc$release) >= 85) {
    loc$release
  } else {
    .newest_release(loc$species)
  }
  x <- fetch_xrefs(loc$species, rel)

  out <- merge(hit[, c("input", "id", "name")], x, by = "id", all.x = TRUE)
  names(out)[names(out) == "xref"] <- "entrez"
  out <- out[order(match(out$input, .as_ids(ids)), out$entrez),
    c("input", "id", "name", "entrez"),
    drop = FALSE
  ]
  row.names(out) <- NULL
  out
}
