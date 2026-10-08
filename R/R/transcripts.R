# Transcript identifiers to their genes, and transcript-level counts to gene level.

.transcript_memo <- new.env(parent = emptyenv())

# coverage below which older assemblies are tried
TRANSCRIPT_COVERAGE <- 0.9
# dated releases tried, one per annotation family
TRANSCRIPT_TRIES <- 3L
# past this distance ids are too few to date; they'd pick the smallest release
TRANSCRIPT_MAX_DIST <- 1

#' Transcript-to-gene records from a GTF, streamed without writing it to disk
#'
#' Keeps `transcript` lines. GTFs before Ensembl release 75 have none, so a
#' chunk without them is read from its exon and CDS lines, which carry the same
#' `transcript_id`, `gene_id` and `gene_name`; duplicates are dropped.
#'
#' @param url A `.gtf.gz` URL.
#' @param chunk Lines to decompress per iteration.
#' @return A data frame of `transcript_id`, `id` (the gene) and `name` (the gene
#'   name, `NA` where the annotation gives none). Identifiers keep whatever
#'   version suffix the GTF writes: GENCODE versions them, Ensembl does not.
#' @seealso [fetch_transcripts()] for the cached table of one release.
#' @examples
#' \dontrun{
#' stream_transcripts(paste0(
#'   "https://ftp.ensembl.org/pub/release-116/gtf/",
#'   "saccharomyces_cerevisiae/",
#'   "Saccharomyces_cerevisiae.R64-1-1.116.gtf.gz"
#' ))
#' }
#' @export
stream_transcripts <- function(url, chunk = 100000) {
  # shrink lines as they stream; the attribute column is most of each line
  k <- .stream_filter(url, function(x) {
    t <- x[grepl("\ttranscript\t", x, fixed = TRUE)]
    if (!length(t)) t <- x[grepl('transcript_id "', x, fixed = TRUE)]
    t <- t[!startsWith(t, "#")]
    tid <- .gtf_attr(t, "transcript_id")
    keep <- !is.na(tid) & !duplicated(tid)
    t <- t[keep]
    paste(tid[keep], .gtf_attr(t, "gene_id"), .gtf_attr(t, "gene_name"), sep = "\t")
  }, chunk)
  if (!length(k)) stop("no transcript records in ", url, call. = FALSE)
  # pre-r75 exon lines repeat a transcript across chunks
  k <- k[!duplicated(sub("\t.*", "", k))]
  f <- matrix(unlist(strsplit(k, "\t", fixed = TRUE)), ncol = 3L, byrow = TRUE)
  name <- f[, 3L]
  name[name == "NA"] <- NA_character_ # paste() writes a missing gene_name as "NA"
  data.frame(
    transcript_id = f[, 1L], id = f[, 2L], name = name,
    stringsAsFactors = FALSE
  )
}

#' Transcript-to-gene table for one species and release
#'
#' Streamed from the annotation's GTF on first use, then cached, so later calls
#' are offline.
#'
#' @inheritParams fetch_mapping
#' @return A data frame of `transcript_id`, `id` and `name`, as
#'   [stream_transcripts()] returns it.
#' @seealso [transcript2gene()] to map a set of transcript identifiers.
#' @examples
#' \dontrun{
#' head(fetch_transcripts("chicken", release = 116))
#' }
#' @export
fetch_transcripts <- function(species, release, assembly = NULL, source = NULL,
                              refresh = FALSE) {
  .transcript_table(.annotation_source(species, release, assembly, source), refresh)
}

#' fetch_transcripts() for an already-resolved annotation source
#' @noRd
.transcript_table <- function(a, refresh = FALSE) {
  f <- .cache_path("transcripts", a$row$source, a$species, a$row$assembly, a$release)
  tab <- if (refresh) NULL else .cache_get(f, .transcript_memo, c("transcript_id", "id", "name"))
  if (is.null(tab)) tab <- .cache_put(f, .transcript_memo, stream_transcripts(.annotation_url(a)))
  tab
}

#' The releases worth trying for transcripts of one species
#'
#' For species outside the transcript index. Stable identifiers persist within
#' an assembly, so the newest Ensembl release of each assembly stands for it,
#' newest assembly first. Chicken renumbered its identifiers with GRCg7b, so an
#' older assembly is worth a second fetch.
#' @noRd
.transcript_candidates <- function(species) {
  fp <- .fp()
  fp <- fp[fp$species == species & fp$source == "ensembl", , drop = FALSE]
  if (!nrow(fp)) {
    stop(sprintf("no Ensembl releases indexed for %s.", species), call. = FALSE)
  }
  # frozen rows (GRCh37) and relocated annotations are tried last
  fp <- fp[order(!is.na(fp$frozen_release), -.rel_num(fp$release)), , drop = FALSE]
  fp <- fp[!duplicated(fp$assembly), , drop = FALSE]
  fp[, c("release", "assembly"), drop = FALSE]
}

#' The best-scoring release for a set of transcript identifiers
#'
#' The best few, best first. `NULL` when the species is outside the
#' transcript index, no release is consistent with the ids, or none fits them
#' (too few ids to date), so the caller falls back to coverage.
#' @noRd
.detected_transcript_release <- function(ids, species, assembly = NULL, source = NULL) {
  tx <- .tx_fp()
  if (!any(tx$species == species & (is.null(source) | tx$source == source)) || !.is_transcript(ids)) {
    return(NULL)
  }
  rel <- suppressWarnings(detect_release(ids, species = species))
  if (!is.null(source)) rel <- rel[rel$source == source, , drop = FALSE]
  if (!is.null(assembly)) rel <- rel[rel$assembly == as.character(assembly), , drop = FALSE]
  if (isTRUE(attr(rel, "numeric_ordered"))) rel <- rel[rel$feasible, , drop = FALSE]
  if (!nrow(rel) || rel$dist[1] > TRANSCRIPT_MAX_DIST) {
    return(NULL)
  }
  # qmax pins a filtered set's vintage but not its family; GENCODE CHR drops scaffolds
  rel <- rel[!duplicated(rel$source), c("release", "assembly", "source"), drop = FALSE]
  utils::head(rel, TRANSCRIPT_TRIES)
}

#' Bring a caller's transcript table to `transcript_id`, `gene_id`, `gene_name`
#'
#' Accepts [fetch_transcripts()] output (`id`, `name`) or [transcript2gene()]
#' output (`gene_id`, `gene_name`).
#' @noRd
.as_transcript_table <- function(mapping) {
  if (!is.data.frame(mapping)) stop("`mapping` must be a data frame.", call. = FALSE)
  pick <- function(a, b) if (a %in% names(mapping)) a else if (b %in% names(mapping)) b else NA
  gid <- pick("gene_id", "id")
  gnm <- pick("gene_name", "name")
  if (!"transcript_id" %in% names(mapping) || is.na(gid)) {
    stop("`mapping` needs `transcript_id` and `id` (or `gene_id`) columns. ",
      "See fetch_transcripts() for the shape.",
      call. = FALSE
    )
  }
  tab <- data.frame(
    transcript_id = as.character(mapping$transcript_id),
    gene_id = as.character(mapping[[gid]]),
    gene_name = if (is.na(gnm)) NA_character_ else as.character(mapping[[gnm]]),
    stringsAsFactors = FALSE
  )
  tab[!is.na(tab$transcript_id) & !is.na(tab$gene_id), , drop = FALSE]
}

#' Match ids against a transcript table, ignoring versions
#' @noRd
.match_transcripts <- function(ids, tab) {
  # GENCODE files PAR_Y copies under their own gene; match the tag first
  hit <- match(.unversioned(ids), .unversioned(tab$transcript_id))
  miss <- is.na(hit)
  if (any(miss)) hit[miss] <- match(.bare(ids[miss]), .bare(tab$transcript_id))
  name <- tab$gene_name[hit]
  gid <- tab$gene_id[hit]
  # unnamed genes keep their id, as in gene_names()
  fill <- !is.na(gid) & !.is_name(name, gid)
  name[fill] <- gid[fill]
  data.frame(
    input = ids, transcript_id = tab$transcript_id[hit],
    gene_id = gid, gene_name = name, stringsAsFactors = FALSE
  )
}

#' Map transcript identifiers to their genes
#'
#' Detects the species from the identifier prefix (`ENSMUST`, `ENSGALT`, ...),
#' fetches the annotation's transcript table (streamed once, then cached) and
#' returns each transcript's gene identifier and gene name. Version suffixes
#' are ignored when matching, so `ENSMUST00000178537.2` resolves against an
#' unversioned annotation and the other way round.
#'
#' Without `release`, [detect_release()] dates the transcript identifiers
#' against a transcript fingerprint index of Ensembl and GENCODE for the
#' commonly quantified species. `source` and `assembly` narrow the candidates.
#'
#' When the best-scoring release leaves any id unmapped, the best release of
#' each other annotation family (Ensembl, GENCODE, GENCODE ALL) is fetched too
#' and the best-covering one kept: a filtered transcriptome (cDNA only) dates by
#' its newest id but leaves its family open.
#'
#' A species outside that index, a set too small to date (a handful of ids),
#' or dated releases explaining less than 90% of the input fall back to
#' coverage: the newest Ensembl release of each assembly is tried,
#' newest first, and the best-covering one kept. Chicken
#' needs this: GRCg7b (release 107 on) renumbered every identifier, so GRCg6a
#' transcripts only map at release 106 or earlier.
#'
#' @inheritParams geneid2name
#' @param ids Character vector of transcript identifiers.
#' @param release Annotation release. Without it, chosen as above.
#' @param mapping Skip fetching and use this table: [fetch_transcripts()]
#'   output, or a data frame with `transcript_id`, `gene_id` and `gene_name`.
#'
#' @return A data frame with one row per input, in order: `input`,
#'   `transcript_id` (as the annotation writes it), `gene_id` and `gene_name`.
#'   Unmapped inputs have `NA` in the last three. A gene without a name has its
#'   identifier as `gene_name`. Attributes `species`, `source`, `release` and
#'   `assembly` record the annotation used.
#' @seealso [aggregate_transcripts()] to sum a transcript-level matrix to genes.
#' @export
#' @examples
#' \dontrun{
#' transcript2gene(c("ENSMUST00000178537.2", "ENSMUST00000178862.2"))
#' transcript2gene(c("ENSGALT00000000003", "ENSGALT00000000004"))
#' }
transcript2gene <- function(ids, species = NULL, release = NULL, assembly = NULL,
                            source = NULL, mapping = NULL, quiet = FALSE) {
  ids <- .as_ids(ids)
  if (!is.null(mapping)) {
    return(.match_transcripts(ids, .as_transcript_table(mapping)))
  }
  if (!is.null(species)) species <- resolve_species(species)
  if (is.null(species)) {
    cand <- detect_species(ids)
    if (!nrow(cand)) stop("cannot identify the species from these ids.", call. = FALSE)
    species <- cand$species[1]
  }
  cand <- if (!is.null(release)) {
    data.frame(release = as.character(release), assembly = assembly %||% NA_character_, source = source %||% NA_character_)
  } else {
    .detected_transcript_release(ids, species, assembly, source)
  }
  n_dated <- if (is.null(release)) NROW(cand) else 0L
  # coverage fallback knows only Ensembl rows
  if (is.null(release) && (is.null(source) || identical(source, "ensembl"))) {
    fb <- .transcript_candidates(species)
    if (!is.null(assembly)) fb <- fb[fb$assembly == as.character(assembly), , drop = FALSE]
    fb$source <- "ensembl"
    cand <- rbind(cand, fb[!paste(fb$release, fb$assembly) %in% paste(cand$release, cand$assembly), , drop = FALSE])
  }
  if (!NROW(cand)) {
    stop(sprintf("no %s release of %s to try; pass `release`.", source %||% "Ensembl", species), call. = FALSE)
  }

  best <- NULL
  for (i in seq_len(nrow(cand))) {
    asm <- if (is.na(cand$assembly[i])) NULL else cand$assembly[i]
    src <- if (is.na(cand$source[i])) NULL else cand$source[i]
    a <- .annotation_source(species, cand$release[i], asm, src)
    out <- .match_transcripts(ids, .as_transcript_table(.transcript_table(a)))
    cov <- if (length(ids)) mean(!is.na(out$gene_id)) else 0
    if (is.null(best) || cov > best$cov) {
      best <- list(out = out, cov = cov, a = a, release = cand$release[i])
    }
    # a dated release must map every id; the fallback needs 90%
    if (cov >= if (i <= n_dated) 1 else TRANSCRIPT_COVERAGE) break
    if (i == n_dated && best$cov >= TRANSCRIPT_COVERAGE) break
  }
  out <- best$out
  if (!quiet) {
    message(sprintf(
      "%s, %s release %s (assembly %s): %d of %d transcripts mapped",
      species, best$a$row$source, best$release, best$a$row$assembly,
      sum(!is.na(out$gene_id)), length(ids)
    ))
  }
  structure(out,
    species = species, source = best$a$row$source,
    release = best$release, assembly = best$a$row$assembly
  )
}

#' The unmapped-transcripts condition, classed so a caller can catch just this one
#' @noRd
.unmapped_transcripts_condition <- function(miss, n) {
  warningCondition(
    sprintf(
      "only %d of %d transcripts mapped to a gene, e.g. %s -- wrong species or release?",
      n - length(miss), n, miss[1]
    ),
    class = "genevintage_unmapped_transcripts"
  )
}

#' Sum the rows of a matrix within groups, keeping its class
#' @noRd
.sum_rows <- function(x, group) {
  lev <- unique(group)
  if (inherits(x, "Matrix")) {
    g <- Matrix::sparseMatrix(
      i = match(group, lev), j = seq_along(group), x = 1,
      dims = c(length(lev), length(group))
    )
    out <- g %*% x
    dimnames(out) <- list(lev, colnames(x))
    return(out)
  }
  if (is.data.frame(x) && !all(vapply(x, is.numeric, NA))) {
    stop("every column of `x` must be numeric to sum transcripts.", call. = FALSE)
  }
  rowsum(x, group, reorder = FALSE)
}

#' Sum transcript-level counts to gene level
#'
#' Maps each row's transcript identifier with [transcript2gene()] and sums the
#' rows of each gene. Rows come back in the order their gene was first seen. Columns,
#' column names and the matrix class (dense, sparse `Matrix`, data frame) are
#' kept.
#'
#' Summing is right for counts and for TPM; it is not for effective lengths.
#'
#' @inheritParams transcript2gene
#' @param x A numeric matrix, sparse `Matrix` or data frame whose rows are
#'   transcripts.
#' @param row_ids Overrides which identifiers to map. Defaults to `rownames(x)`;
#'   set this when `x` has none.
#' @param label Label gene rows by `"id"` (the gene identifier) or `"name"`
#'   (the gene name, with the identifier appended where two genes share one).
#' @param unmapped `"keep"` leaves unmapped rows under their own identifier, as
#'   [geneid2name()] keeps an identifier it cannot name; `"drop"` removes them.
#' @param quiet Suppress the messages reporting the annotation and how many rows
#'   mapped.
#'
#' @return `x` with one row per gene. A warning of class
#'   `genevintage_unmapped_transcripts` is raised when fewer than half the rows
#'   mapped; call [transcript2gene()] and pass its result as `mapping` to see
#'   which ones.
#' @seealso [transcript2gene()], [geneid2name()].
#' @export
#' @examples
#' tx <- data.frame(
#'   transcript_id = c("ENSMUST00000000001", "ENSMUST00000000002"),
#'   id = "ENSMUSG00000000001", name = "Gnai3"
#' )
#' m <- matrix(1:6, nrow = 3, dimnames = list(
#'   c("ENSMUST00000000001.4", "ENSMUST00000000002.1", "ENSMUST00000999999.1"),
#'   c("s1", "s2")
#' ))
#' aggregate_transcripts(m, mapping = tx)
#'
#' \dontrun{
#' aggregate_transcripts(salmon_counts) # species and release from the ids
#' }
aggregate_transcripts <- function(x, species = NULL, release = NULL, assembly = NULL,
                                  source = NULL, mapping = NULL, row_ids = NULL,
                                  label = c("id", "name"),
                                  unmapped = c("keep", "drop"), quiet = FALSE) {
  label <- match.arg(label)
  unmapped <- match.arg(unmapped)
  if (is.null(dim(x))) {
    stop("`x` must be a matrix, sparse `Matrix` or data frame; ",
      "use transcript2gene() for a vector of ids.",
      call. = FALSE
    )
  }
  rid <- row_ids %||% rownames(x)
  .check_row_ids(x, rid)
  tab <- transcript2gene(rid,
    species = species, release = release, assembly = assembly,
    source = source, mapping = mapping, quiet = quiet
  )
  ok <- !is.na(tab$gene_id)
  if (!quiet && any(!ok)) {
    message(sprintf(
      "%d of %d transcripts have no gene; %s.", sum(!ok), length(ok),
      if (unmapped == "keep") "kept under their own ids" else "dropped"
    ))
  }
  if (length(ok) && mean(ok) < 0.5) {
    warning(.unmapped_transcripts_condition(rid[!ok], length(ok)))
  }
  gid <- .coalesce(tab$gene_id, rid)
  lab <- .coalesce(tab$gene_name, rid)
  if (unmapped == "drop") {
    x <- x[ok, , drop = FALSE]
    gid <- gid[ok]
    lab <- lab[ok]
  }
  # group on id so genes sharing a name stay apart
  out <- .sum_rows(x, gid)
  if (label == "name") {
    named <- rownames(out) %in% tab$gene_id[ok]
    rownames(out) <- .suffix_collisions(lab[match(rownames(out), gid)], named, rownames(out))
  }
  out
}
