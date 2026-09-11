# Cross-references to other identifier systems. Ensembl publishes one TSV per
# species and release under pub/release-N/tsv/, so this needs no MySQL dump and
# no biomaRt.

# Do not put the release in the pattern. Under Ensembl Genomes the file carries
# that division's release number, so yeast at vertebrates 116 ships
# Saccharomyces_cerevisiae.R64-1-1.63.entrez.tsv.gz.
.XREF_FILE <- "\\.%s\\.tsv\\.gz$"

#' Stream one cross-reference TSV, keeping the gene-level pairs
#'
#' The file lists one row per transcript, so a gene with twelve transcripts
#' appears twelve times. Deduplicating as the stream goes keeps memory
#' proportional to the gene count.
#'
#' @param url A `.tsv.gz` cross-reference URL.
#' @param db The `db_name` to keep, e.g. `"EntrezGene"`.
#' @param chunk Lines to decompress per iteration.
#' @return A data frame of `id` and `xref`.
#' @export
stream_xrefs <- function(url, db = "EntrezGene", chunk = 200000) {
  seen <- character()
  k <- .stream_filter(url, function(x) {
    x <- x[!startsWith(x, "gene_stable_id")]
    if (!length(x)) {
      return(character())
    }
    f <- strsplit(x, "\t", fixed = TRUE)
    col <- function(i) vapply(f, function(z) if (length(z) >= i) z[i] else NA_character_, "")
    keep <- col(5L) == db & !is.na(col(4L)) & nzchar(col(4L))
    if (!any(keep)) {
      return(character())
    }
    pair <- paste(col(1L)[keep], col(4L)[keep], sep = "\t")
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
#' Streams the cross-reference TSV Ensembl publishes for a species and release,
#' keeps the gene-level pairs and caches them, so later calls are offline.
#'
#' @param species Species in any form [resolve_species()] accepts.
#' @param release Ensembl release. Cross-references are published against
#'   Ensembl releases, so a GENCODE release number is not one.
#' @param db Cross-reference database. `"entrez"` is the only one indexed.
#' @param refresh Refetch and overwrite the cached copy.
#' @return A data frame of `id` and `xref`, one row per distinct pair.
#' @seealso [entrez_ids()] for the lookup this feeds.
#' @export
#' @examples
#' \dontrun{
#' head(fetch_xrefs("yeast", release = 116))
#' }
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
#' Takes gene identifiers or gene names and returns their Entrez identifiers,
#' detecting the species the way the rest of the package does. Cross-references
#' are streamed once per species and release, then cached.
#'
#' Entrez is not a relabelling of Ensembl. A gene may carry several Entrez
#' identifiers or none, so a row is a pair rather than a translation, and
#' `entrez` is `NA` where Ensembl records no cross-reference. Roughly a third of
#' human genes are in that position, mostly novel long non-coding RNAs.
#'
#' Cross-references are published against Ensembl releases. Identifiers dated to
#' a GENCODE release are looked up in the newest indexed Ensembl release for
#' that species, since a cross-reference is a relationship between stable
#' identifiers rather than a vintage.
#'
#' @inheritParams mito_genes
#' @param ids Gene identifiers, or gene names. Names need `species`, since a
#'   name does not say which species it came from.
#' @return A data frame of `input`, `id`, `name` and `entrez`, one row per pair.
#' @seealso [fetch_xrefs()], [gene_names()].
#' @export
#' @examples
#' \dontrun{
#' entrez_ids(c("ENSG00000141510", "ENSG00000012048"))
#' entrez_ids(c("TP53", "BRCA1"), species = "human")
#' }
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
  # a GENCODE vintage has no cross-reference file of its own
  rel <- if (identical(loc$source, "ensembl")) {
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
