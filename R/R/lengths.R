# exon-union gene lengths, the FPKM and TPM denominator

.lengths_memo <- new.env(parent = emptyenv())

#' Total length of the union of each gene's exons
#'
#' Exons that overlap count once. GTF coordinates are 1-based and closed.
#'
#' @param gene,start,end Parallel vectors, one element per exon.
#' @return A named integer vector of lengths, one per gene, in first-seen order.
#' @noRd
.exon_union_length <- function(gene, start, end) {
  if (!length(gene)) {
    return(stats::setNames(integer(), character()))
  }
  o <- order(gene, start, end, method = "radix")
  g <- gene[o]
  s <- start[o]
  # a new block starts at each gene, or at an exon past every end seen so far
  reach <- stats::ave(end[o], g, FUN = cummax)
  first <- c(TRUE, g[-1L] != g[-length(g)])
  new_block <- first | s > c(-Inf, reach[-length(reach)])
  # reach at a block's last exon is that block's end
  last <- c(which(new_block)[-1L] - 1L, length(s))
  out <- rowsum(reach[last] - s[new_block] + 1, g[new_block], reorder = FALSE)[, 1]
  ids <- unique(gene)
  stats::setNames(as.integer(out[ids]), ids)
}

#' Stream exon-union lengths from a GTF and write nothing to disk
#'
#' This function keeps only `exon` lines. It reads gene identifiers from the
#' `gene_id` attribute, so a GENCODE GTF gives versioned identifiers, as its
#' `gene` lines do in [stream_gtf()].
#'
#' @param url A `.gtf.gz` URL.
#' @param chunk Lines to decompress in each iteration.
#' @return A data frame of `id` and `length` (bases).
#' @examplesIf interactive()
#' stream_exon_lengths(paste0(
#'   "https://ftp.ensembl.org/pub/release-116/gtf/",
#'   "saccharomyces_cerevisiae/",
#'   "Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz"
#' ))
#' @export
stream_exon_lengths <- function(url, chunk = 100000) {
  # shrink lines as they stream; the attribute column is most of each line
  k <- .stream_filter(url, function(x) {
    x <- x[grepl("\texon\t", x, fixed = TRUE)]
    if (!length(x)) {
      return(x)
    }
    f <- strsplit(x, "\t", fixed = TRUE)
    id <- .gtf_attr(x, "gene_id")
    keep <- !is.na(id)
    paste(id, vapply(f, `[`, "", 4L), vapply(f, `[`, "", 5L), sep = "\t")[keep]
  }, chunk)
  if (!length(k)) stop("no exon records in ", url, call. = FALSE)
  f <- matrix(unlist(strsplit(k, "\t", fixed = TRUE)), ncol = 3L, byrow = TRUE)
  start <- as.numeric(f[, 2L])
  end <- as.numeric(f[, 3L])
  if (anyNA(start) || anyNA(end) || any(end < start)) {
    stop("exon coordinates could not be read.", call. = FALSE)
  }
  len <- .exon_union_length(f[, 1L], start, end)
  data.frame(id = names(len), length = unname(len), stringsAsFactors = FALSE)
}

#' Exon-union gene lengths for one species and release
#'
#' @inheritParams fetch_mapping
#' @param ids Optional gene identifiers or gene names to return lengths for, in
#'   order. Versioned identifiers match their unversioned record. The function
#'   matches names through [fetch_mapping()]. A name that several genes share
#'   gets `NA`.
#'
#' @return Without `ids`, a data frame of `id` and `length` for each gene in the
#'   annotation. With `ids`, a named integer vector with the length of `ids`,
#'   and `NA` where there is no match.
#' @seealso [fetch_mapping()] for genomic span, [stream_exon_lengths()].
#' @examplesIf interactive()
#' gene_lengths("human", release = 116, ids = c("ENSG00000141510", "TP53"))
#' @export
gene_lengths <- function(species, release, ids = NULL, assembly = NULL,
                         source = NULL, refresh = FALSE) {
  a <- .annotation_source(species, release, assembly, source)
  f <- .cache_path("lengths", a$row$source, a$species, a$row$assembly, a$release)
  tab <- if (refresh) NULL else .cache_get(f, .lengths_memo, c("id", "length"))
  if (is.null(tab)) tab <- .cache_put(f, .lengths_memo, stream_exon_lengths(.annotation_url(a)))
  if (is.null(ids)) {
    return(tab)
  }
  ids <- .as_ids(ids)
  bare <- .bare(tab$id)
  hit <- match(.bare(ids), bare)
  # unmatched entries may be gene names
  miss <- is.na(hit)
  if (any(miss)) {
    m <- fetch_mapping(a$species, release, assembly = a$row$assembly, source = a$row$source)
    m <- m[.is_name(m$name, m$id), ]
    shared <- m$name[duplicated(m$name)]
    m <- m[!m$name %in% shared, ]
    hit[miss] <- match(.bare(m$id[match(ids[miss], m$name)]), bare)
  }
  stats::setNames(tab$length[hit], ids)
}
