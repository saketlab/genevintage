# The other direction: names and Entrez ids back to stable identifiers, and the
# whole annotation record for a set of genes.

#' Gene identifiers for gene names or Entrez identifiers
#'
#' The reverse of [gene_names()]. Published marker lists and enrichment results
#' carry symbols or Entrez identifiers; count matrices carry Ensembl ones.
#'
#' Symbols and Entrez identifiers can each match multiple genes, so the table
#' returns every match. Unmatched inputs are retained with `id = NA`.
#'
#' @inheritParams mito_genes
#' @param x Gene names, or Entrez identifiers.
#' @param from What `x` holds: `"name"`, `"entrez"`, or `"auto"` to decide from
#'   the input. All-digit input is Entrez; anything else is a name.
#' @return A data frame of `input`, `id` and `name`, one row per match.
#' @seealso [gene_names()] for the forward direction, [entrez_ids()].
#' @export
#' @examples
#' \dontrun{
#' gene_ids(c("TP53", "BRCA1"), species = "human")
#' gene_ids(c("7157", "672"), from = "entrez", species = "human")
#' }
gene_ids <- function(x, from = c("auto", "name", "entrez"), species = NULL,
                     release = NULL, assembly = NULL, source = NULL,
                     mapping = NULL, quiet = FALSE) {
  from <- match.arg(from)
  x <- .as_ids(x)
  if (from == "auto") from <- if (length(x) && all(grepl("^[0-9]+$", x))) "entrez" else "name"

  # Entrez lookup through a supplied mapping requires its own entrez column
  if (from == "entrez" && !is.null(mapping) && !"entrez" %in% names(mapping)) {
    stop("`mapping` carries no cross-references; give it an `entrez` column, ",
      "or name `species` and `release` instead of `mapping`.",
      call. = FALSE
    )
  }
  loc <- .locate_annotation(if (from == "name") x else NULL, species, release,
    assembly, source, quiet,
    mapping = mapping
  )
  map <- loc$map

  if (from == "name") {
    hit <- .restrict(map, x)[, c("input", "id", "name"), drop = FALSE]
  } else if (!is.null(mapping)) {
    xr <- data.frame(
      id = mapping$id, xref = as.character(mapping$entrez),
      stringsAsFactors = FALSE
    )
    xr <- xr[!is.na(xr$xref) & xr$xref %in% x, , drop = FALSE]
  } else {
    rel <- if (identical(loc$source, "ensembl")) {
      loc$release
    } else {
      .newest_release(loc$species)
    }
    xr <- fetch_xrefs(loc$species, rel)
    xr <- xr[xr$xref %in% x, , drop = FALSE]
  }
  if (from == "entrez") {
    hit <- data.frame(
      input = xr$xref, id = xr$id,
      name = map$name[match(xr$id, map$id)],
      stringsAsFactors = FALSE
    )
  }
  .keep_unmatched(hit, x)
}

# pad unmatched tokens with NA row indices to preserve column types;
# grouped positions expand multiple matches in input order
.keep_unmatched <- function(hit, x) {
  u <- unique(x)
  pos <- match(hit$input, u)
  o <- order(pos)
  cnt <- tabulate(pos, nbins = length(u))
  cstart <- cumsum(c(0L, cnt))[seq_along(cnt)] # 0-based block start in `o`

  xi <- match(x, u)
  len <- pmax(cnt[xi], 1L) # an unmatched token still gets a row
  k <- rep(seq_along(x), len)
  i <- xi[k]
  idx <- rep(NA_integer_, length(k))
  got <- cnt[i] > 0L
  idx[got] <- o[cstart[i[got]] + sequence(len)[got]]

  out <- hit[idx, , drop = FALSE]
  out$input <- x[k]
  row.names(out) <- NULL
  out
}

#' Everything the annotation knows about a set of genes
#'
#' Symbol, biotype and coordinates in one table, for genes given as identifiers
#' or as names. The caller's own token comes back in `input`, so the result
#' joins straight onto a results table or a feature frame.
#'
#' Unmatched inputs are retained with `NA` columns. Names matching several
#' genes yield a row for each.
#'
#' @inheritParams mito_genes
#' @return A data frame of `input`, `id`, `name`, `chr`, `biotype` and, where
#'   the annotation carries them, `start`, `end`, `strand` and `span`. `span` is
#'   the genomic extent from first to last base. See [stream_gtf()].
#' @seealso [gene_names()] for one value per input.
#' @export
#' @examples
#' \dontrun{
#' ann <- annotate_genes(rownames(counts))
#' res$biotype <- ann$biotype[match(rownames(res), ann$input)]
#' }
annotate_genes <- function(ids, species = NULL, release = NULL, assembly = NULL,
                           source = NULL, mapping = NULL, quiet = FALSE) {
  ids <- .as_ids(ids)
  loc <- .locate_annotation(ids, species, release, assembly, source, quiet,
    mapping = mapping
  )
  hit <- .restrict(loc$map, ids)
  cols <- c(
    "input", "id", "name", "chr", "biotype",
    intersect(c("start", "end", "strand", "span"), names(hit))
  )
  .keep_unmatched(hit[, cols, drop = FALSE], ids)
}
