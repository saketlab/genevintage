# Mapping tables, streamed from the Ensembl GTF. No disk, no biomaRt.

# read once per session; rename_rows() over many matrices would re-inflate the
# same 60k-row table each call
.mapping_memo <- new.env(parent = emptyenv())

# memory, then disk, then the network.
# a half-written file survives a crash, so trust a cached value only when it
# carries the columns it should; discarding is safe, nothing here is unrebuildable
.cache_get <- function(f, memo, cols) {
  hit <- memo[[f]]
  if (!is.null(hit)) {
    return(hit)
  }
  if (!file.exists(f)) {
    return(NULL)
  }
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.data.frame(x) && all(cols %in% names(x))) {
    memo[[f]] <- x
    return(x)
  }
  unlink(f)
  NULL
}

# write then rename, so an interrupted save cannot leave a half-written cache
.cache_put <- function(f, memo, x) {
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(f, ".tmp")
  saveRDS(x, tmp) # gzip: 10x faster to write, 200KB larger
  file.rename(tmp, f)
  memo[[f]] <- x
  x
}

# not created here; .cache_put() makes it, so a query never leaves a directory
# behind in the user's filespace
.cache_dir <- function() tools::R_user_dir("genevintage", "cache")

# Ensembl's FTP host refuses the occasional connection under load; the failure
# is not sticky
.attempt <- function(f, tries = 3L, wait = 2) {
  for (i in seq_len(tries)) {
    v <- try(f(), silent = TRUE)
    if (!inherits(v, "try-error")) {
      return(v)
    }
    if (i < tries) Sys.sleep(wait * i)
  }
  v
}

# ---- protocol -------------------------------------------------------------
#
# FTP preferred, probed once per session so a blocked port 21 doesn't cost a
# timeout on every fetch. options(genevintage.protocol=) overrides;
# "ftp" disables the fallback.

.as_ftp <- function(url) sub("^https?://", "ftp://", url)
.as_https <- function(url) sub("^ftp://", "https://", url)

.ftp_works <- local({
  cache <- NULL
  function(url) {
    if (!is.null(cache)) {
      return(cache)
    }
    # short leash; a blocked port should cost seconds, not a timeout
    old <- options(timeout = 8)
    on.exit(options(old), add = TRUE)
    root <- sub("(^ftp://[^/]+/[^/]+/).*$", "\\1", .as_ftp(url))
    v <- try(suppressWarnings(readLines(root, n = 1L, warn = FALSE)), silent = TRUE)
    cache <<- !inherits(v, "try-error")
    cache
  }
})

#' The URLs to try, in order
#'
#' One entry when the protocol is pinned or FTP is unreachable, two when FTP is
#' worth trying and HTTPS is there to catch it.
#' @noRd
.url_candidates <- function(url) {
  pref <- getOption("genevintage.protocol", "auto")
  if (identical(pref, "https")) {
    return(.as_https(url))
  }
  if (identical(pref, "ftp")) {
    return(.as_ftp(url))
  }
  if (grepl("^ftp://", url) || .ftp_works(url)) {
    return(unique(c(.as_ftp(url), .as_https(url))))
  }
  .as_https(url)
}

#' Do something with a URL, over whichever protocol answers
#'
#' Only the last candidate is retried: a protocol that is blocked fails on the
#' first attempt and there is nothing to wait for, while a transient refusal on
#' the one that does work is worth a second and third go.
#' @noRd
.by_protocol <- function(url, f) {
  cand <- .url_candidates(url)
  for (i in seq_along(cand)) {
    v <- .attempt(function() f(cand[i]), tries = if (i == length(cand)) 3L else 1L)
    if (!inherits(v, "try-error")) {
      return(v)
    }
  }
  v
}

# ---- directory listings ---------------------------------------------------

#' Entry names from a directory index, whichever protocol served it
#'
#' HTTPS returns an HTML index and FTP returns `ls -l` lines, so the two are
#' normalised to bare names here rather than every caller learning both shapes.
#' @noRd
.listing_names <- function(x) {
  x <- x[nzchar(trimws(x))]
  html <- grepl("<[a-zA-Z/!]", x)
  if (any(grepl("href=", x, fixed = TRUE))) {
    n <- unlist(regmatches(x, gregexpr('(?<=href=")[^"?][^"]*', x, perl = TRUE)))
  } else if (any(html)) {
    # markup but no links: an error page. raw HTML is not a filename.
    return(character())
  } else {
    # ls -l: nine fields, the ninth the name, which may contain spaces
    perms <- "^[-dlbcps][-rwxstST]{9}[.+]?"
    lsl <- grepl(paste0(perms, "\\s+(\\S+\\s+){7}\\S"), x)
    # in an NLST listing a permissions string means malformed, not a filename
    n <- if (any(lsl)) sub("^(\\S+\\s+){8}", "", x[lsl]) else x[!grepl(perms, x)]
    # symlink lines end "name -> target"
    n <- sub(" -> .*$", "", n)
  }
  n <- sub("/$", "", trimws(n))
  unique(n[nzchar(n) & !grepl("^[./]", n)])
}

#' A directory's entry names, or a message naming what could not be reached
#' @noRd
.listing <- function(dir_url, what = dir_url) {
  idx <- .by_protocol(dir_url, function(u) readLines(u, warn = FALSE))
  if (inherits(idx, "try-error")) stop(what, call. = FALSE)
  .listing_names(idx)
}

# keep only the lines the filter wants; peak memory tracks what is kept, not
# the size of the file
.stream_filter <- function(url, keep, chunk) {
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  con <- .by_protocol(url, function(u) gzcon(url(u, open = "rb")))
  if (inherits(con, "try-error")) {
    stop("cannot open ", url, call. = FALSE)
  }
  on.exit(close(con), add = TRUE)

  parts <- list()
  i <- 0L
  repeat {
    x <- readLines(con, n = chunk, warn = FALSE)
    if (!length(x)) break
    x <- keep(x)
    if (!length(x)) next
    i <- i + 1L
    parts[[i]] <- x
  }
  unlist(parts, use.names = FALSE)
}

#' Pull gene records out of a GTF without ever writing it to disk
#'
#' Decompresses the HTTP stream in chunks and keeps only `gene` lines, so peak
#' memory tracks the gene count rather than the size of the file.
#'
#' @param url A `.gtf.gz` URL.
#' @param chunk Lines to decompress per iteration.
#' @return A data frame of `id`, `name`, `chr`, `start`, `end`, `strand`,
#'   `span`, `id_version` and `biotype`. `span` is the genomic extent, not the
#'   exonic length TPM and FPKM require.
#' @examples
#' \dontrun{
#' # any Ensembl-shaped GTF; nothing is written to disk
#' stream_gtf(paste0(
#'   "https://ftp.ensembl.org/pub/release-116/gtf/",
#'   "saccharomyces_cerevisiae/",
#'   "Saccharomyces_cerevisiae.R64-1-1.116.gtf.gz"
#' ))
#' }
#' @export
stream_gtf <- function(url, chunk = 100000) {
  # fixed = TRUE; the regex form rescans every line for the same result
  k <- .stream_filter(
    url, function(x) grep("\tgene\t", x, value = TRUE, fixed = TRUE),
    chunk
  )
  if (!length(k)) {
    # releases 50-74 have no gene rows. exon and CDS carry the same attributes,
    # so take the first line per gene, deduplicating as the stream goes.
    seen <- character()
    k <- .stream_filter(url, function(x) {
      x <- x[!startsWith(x, "#")]
      if (!length(x)) {
        return(x)
      }
      # lazy and PCRE; the greedy TRE form backtracks from the end of every line
      id <- sub('.*?gene_id "([^"]*)".*', "\\1", x, perl = TRUE)
      new <- !duplicated(id) & !(id %in% seen)
      seen <<- c(seen, id[new])
      x[new]
    }, chunk)
  }
  if (!length(k)) stop("no gene records in ", url, call. = FALSE)

  # anchor on the separator; unanchored ".*gene_id" also matches havana_gene_id
  attr_of <- function(key) {
    # lookbehind, not a leading anchor; the first attribute follows a tab
    rx <- sprintf('(?<![A-Za-z_])%s "([^"]*)"', key)
    m <- regexpr(rx, k, perl = TRUE)
    out <- rep(NA_character_, length(k))
    hit <- m != -1L
    if (any(hit)) {
      st <- attr(m, "capture.start")[hit, 1]
      len <- attr(m, "capture.length")[hit, 1]
      out[hit] <- substring(k[hit], st, st + len - 1L)
    }
    out
  }
  # GTF columns 1, 4, 5, 7
  f <- strsplit(k, "\t", fixed = TRUE)
  col <- function(i) vapply(f, `[`, "", i)
  chr <- col(1L)
  start <- suppressWarnings(as.integer(col(4L)))
  end <- suppressWarnings(as.integer(col(5L)))
  strand <- col(7L)
  ids <- attr_of("gene_id")
  if (mean(is.na(ids)) > 0.01) {
    stop("gene_id could not be read from ", url,
      " -- the attribute format is not what was expected.",
      call. = FALSE
    )
  }
  # GENCODE ships no gene_version; its version is in the id (ENSG00000223972.5)
  ver <- as.integer(attr_of("gene_version"))
  if (all(is.na(ver))) ver <- .id_version(ids)
  # GENCODE writes gene_type where Ensembl writes gene_biotype; without the
  # fallback every biotype-driven function returns nothing on GENCODE
  bio <- attr_of("gene_biotype")
  if (all(is.na(bio))) bio <- attr_of("gene_type")
  data.frame(
    id = ids,
    name = attr_of("gene_name"), # NA for ~45% of human genes
    chr = chr,
    start = start,
    end = end,
    strand = strand,
    # genomic extent, first to last base. TPM and FPKM divide by the union of a
    # gene's exons, which streaming only gene lines never reads.
    span = end - start + 1L,
    id_version = ver,
    biotype = bio,
    stringsAsFactors = FALSE
  )
}

#' Where GENCODE keeps one release's GTF
#'
#' GENCODE has no `release-NN/gtf/<species>/` tree to list: the path is
#' predictable, the species are only human and mouse, and each release ships two
#' gene sets: the primary assembly, and the scaffold- and patch-inclusive one.
#' Shared with `data-raw/build_gencode.R` so the layout is stated once.
#' @noRd
.gencode_url <- function(root, species, release, source = "gencode") {
  dir <- c(homo_sapiens = "Gencode_human", mus_musculus = "Gencode_mouse")[species]
  if (is.na(dir)) {
    stop("GENCODE covers human and mouse only, not ", species, ".", call. = FALSE)
  }
  scope <- if (identical(source, "gencode_all")) ".chr_patch_hapl_scaff" else ""
  sprintf(
    "%s/%s/release_%s/gencode.v%s%s.annotation.gtf.gz",
    root, dir, release, release, scope
  )
}

#' Mapping table for one species and Ensembl release
#'
#' Streamed from the Ensembl GTF on first use, then cached, so later calls are
#' offline. The GTF is the only source that carries `gene_version`, and unlike
#' the biomaRt archives it exists for every release.
#'
#' @param species Species, in any form [resolve_species()] accepts:
#'   `"homo_sapiens"`, `"Human"` or `"human"`.
#' @param release Ensembl release number.
#' @param assembly Assembly suffix as the index records it, e.g. `"38"` or
#'   `"37"`. Together with `species`, `source` and `release` it selects the index
#'   row that carries the annotation's FTP root and, where an assembly is
#'   frozen, the release its annotation actually stopped at.
#' @param source Annotation family: `"ensembl"`, or `"gencode"` /
#'   `"gencode_all"` for GENCODE's primary and scaffold-inclusive gene sets.
#'   Needed only where one release number exists in more than one of them.
#' @param refresh Re-download even when a cached copy exists.
#'
#' @return A data frame of `id`, `name`, `chr`, `start`, `end`, `strand`,
#'   `span`, `id_version` and `biotype`, suitable as the `mapping` argument of
#'   [gene_names()]. `span` is the genomic extent, not the exonic length TPM and
#'   FPKM require. See [stream_gtf()].
#' @examples
#' \dontrun{
#' m <- fetch_mapping("yeast", release = 116)
#' head(m)
#'
#' # human release 112 ships two assemblies, so one must be named
#' fetch_mapping("human", release = 112, assembly = "38")
#' }
#' @export
fetch_mapping <- function(species, release, assembly = NULL, source = NULL,
                          refresh = FALSE) {
  # annotation_root and frozen_release carry the relocation as data; human
  # GRCh37 is both relocated and frozen at r87
  species <- resolve_species(species)
  fp <- .fp()
  row <- fp[fp$species == species & fp$release == as.character(release), ]
  if (!is.null(source)) row <- row[row$source == source, ]
  if (!is.null(assembly)) row <- row[row$assembly == as.character(assembly), ]
  if (!nrow(row)) {
    stop(
      sprintf(
        "no index row for %s release %s%s%s.", species, release,
        if (is.null(source)) "" else paste0(" source ", source),
        if (is.null(assembly)) "" else paste0(" assembly ", assembly)
      ),
      call. = FALSE
    )
  }
  if (nrow(row) > 1L) {
    # GENCODE's two gene sets share a release and assembly; naming the
    # assemblies alone would read "2 assemblies (38, 38)"
    if (length(unique(row$source)) > 1L) {
      stop(sprintf(
        "%s release %s exists in %d annotation sources (%s); pass `source`. detect_release(ids, species = \"%s\") reports which one the ids came from.",
        species, release, length(unique(row$source)),
        paste(unique(row$source), collapse = ", "), species
      ), call. = FALSE)
    }
    stop(
      sprintf(
        "%s release %s has %d assemblies (%s); pass `assembly`. detect_release(ids, species = \"%s\") reports which one the ids came from.",
        species, release, nrow(row), paste(row$assembly, collapse = ", "), species
      ),
      call. = FALSE
    )
  }
  rel <- if (is.na(row$frozen_release)) row$release else row$frozen_release

  f <- .cache_path("mapping", row$source, species, row$assembly, rel)
  if (!refresh) {
    hit <- .cache_get(f, .mapping_memo, c("id", "name", "chr", "biotype", "span"))
    if (!is.null(hit) && all(is.na(hit$biotype))) hit <- NULL # pre-gene_type cache
    if (!is.null(hit)) {
      return(hit)
    }
  }

  if ("annotation_url" %in% names(row) && !is.na(row$annotation_url) &&
      nzchar(row$annotation_url)) {
    # Old archives use filenames such as gencode_v4.annotation.GRCh37.gtf.gz.
    url <- row$annotation_url
  } else if (startsWith(row$source, "gencode")) {
    url <- .gencode_url(row$annotation_root, species, rel, row$source)
  } else {
    dir_url <- sprintf("%s/release-%s/gtf/%s/", row$annotation_root, rel, species)
    idx <- .listing(dir_url, paste("cannot reach", dir_url))
    # the plain <Species>.<assembly>.<release>.gtf.gz, not abinitio/chr/patch
    cand <- grep("\\.gtf\\.gz$", idx, value = TRUE)
    cand <- grep("abinitio|\\.chr\\.|patch", cand, value = TRUE, invert = TRUE)
    if (!length(cand)) stop("no GTF found for ", species, " release ", rel, call. = FALSE)
    url <- paste0(dir_url, cand[1])
  }

  .cache_put(f, .mapping_memo, stream_gtf(url))
}
