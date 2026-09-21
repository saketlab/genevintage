# Species and release detection from a bare vector of identifiers.

# how far qmax may fall below an observed id number before a release is ruled
# out; Ensembl retires top-numbered genes occasionally
QMAX_TOLERANCE <- 0.999

#' Read once, keep for the session
#'
#' Cache the first non-NULL result for subsequent calls.
#' @noRd
.once <- function(f) {
  v <- NULL
  function() {
    if (is.null(v)) v <<- f()
    v
  }
}

#' A data file that ships in inst/extdata
#' @noRd
.extdata <- function(name) {
  f <- system.file("extdata", name, package = "genevintage")
  if (!nzchar(f)) f <- file.path("inst/extdata", name) # fallback for sourced code
  f
}

#' The shipped fingerprint index
#' @noRd
.fp <- .once(function() {
  f <- .extdata("fingerprints.rds")
  if (!file.exists(f)) {
    stop("fingerprint index not found; the install is incomplete.", call. = FALSE)
  }
  readRDS(f)
})

# compare release numbers numerically after removing GENCODE letter tags
.rel_num <- function(x) {
  # Historical GENCODE releases include 2b, 3b and 3c as well as mouse M tags.
  suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", x)))
}

# Use the newest indexed release for stable-identifier lookups.
.newest_release <- function(species, source = "ensembl") {
  fp <- .fp()
  r <- fp$release[fp$species == species & fp$source == source]
  if (!length(r)) {
    stop(sprintf("no %s releases indexed for %s.", source, species), call. = FALSE)
  }
  r[which.max(.rel_num(r))]
}

#' Where the index violates its own ordering invariant
#'
#' Identifier numbers are assigned upward, so `qmax` should not fall as releases
#' advance within one species and assembly. Ensembl retires genes at the top of
#' the range occasionally, so only a fall of more than 1% is reported.
#' Only numerically ordered schemes are checked.
#'
#' @param fp A fingerprint index. Defaults to the shipped one.
#' @return A character vector of violation messages, empty when the index is sound.
#' @examples
#' # report ordering violations in the shipped index
#' qmax_violations()
#' @export
qmax_violations <- function(fp = .fp()) {
  fp <- fp[fp$ordered, ]
  if (!nrow(fp)) {
    return(character())
  }
  key <- paste(fp$source, fp$species, fp$assembly)
  # sort by group and release so neighbouring qmax values can be compared
  o <- order(key, .rel_num(fp$release))
  k <- key[o]
  q <- fp$qmax[o]
  bad <- c(FALSE, k[-1] == k[-length(k)] & diff(q) < -0.01 * q[-length(q)])
  if (!any(bad)) {
    return(character())
  }
  # one message per (species, assembly), naming every release that fell
  lab <- paste(fp$species[o], fp$assembly[o])[bad]
  rel <- fp$release[o][bad]
  ag <- tapply(rel, factor(lab, levels = unique(lab)), paste, collapse = ", ")
  unname(sprintf(
    "%s: qmax falls sharply at release %s -- check the assembly labels",
    names(ag), ag
  ))
}

#' The distinct species rows, derived once from the index
#' @noRd
.known_species <- .once(function() {
  # merge annotation sources sharing a species, prefix and scheme
  x <- .fp()[, c("species", "division", "prefix", "scheme", "source")]
  x <- x[order(x$source != "ensembl"), ] # keep the Ensembl division
  x <- x[!duplicated(paste(x$species, x$prefix, x$scheme)), ]
  x$source <- NULL
  x
})

#' The pattern that identifies one species' identifiers
#'
#' Every scheme but Ensembl belongs to a single species, so its own regex from
#' [guess_scheme()] is the whole answer. Ensembl shares one scheme across many
#' species, distinguished by the stem before the feature letter, so only that
#' case needs the prefix.
#' @noRd
.species_pattern <- function(scheme, prefix) {
  if (identical(scheme, "ensembl")) {
    sprintf("^%s[EGTPR]\\d+$", sub("[EGTPR]$", "", prefix))
  } else if (scheme %in% names(SCHEMES)) {
    SCHEMES[[scheme]]
  } else {
    # a scheme SCHEMES does not know; its prefix identifies it
    sprintf("^\\Q%s\\E[0-9]+$", prefix)
  }
}

#' Which species a set of identifiers belongs to
#'
#' Matches identifier prefixes to indexed species.
#'
#' @param ids Character vector of gene identifiers.
#' @return A data frame of candidate species ordered by the share of `ids` they
#'   explain, with columns `species`, `common`, `division`, `prefix`, `scheme`
#'   and `frac`.
#'   Rows can tie at the same `frac`: dog breeds and mouse strains share an
#'   identifier space. Check for ties when that distinction matters.
#' @examples
#' detect_species(c("ENSMUSG00000051951", "ENSMUSG00000089699"))
#'
#' # the prefix identifies bonobo
#' detect_species("ENSPPAG00000021109")[1, c("species", "common", "frac")]
#' @export
detect_species <- function(ids) {
  ids <- .bare(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  known <- .known_species()
  # empty input returns no candidates, avoiding undefined fractions
  if (!length(ids)) known <- known[0, ]
  # Ensembl patterns differ only in the stem, so one extraction and a table
  # answers all of them at once
  frac <- numeric(nrow(known))
  ens <- known$scheme == "ensembl"
  if (any(ens)) {
    m <- regexpr("^ENS[A-Z]{0,4}?(?=[EGTPR][0-9]+$)", ids, perl = TRUE)
    hit <- m != -1L
    tb <- table(substring(ids[hit], 1L, attr(m, "match.length")[hit]))
    f <- as.numeric(tb[sub("[EGTPR]$", "", known$prefix[ens])]) / length(ids)
    f[is.na(f)] <- 0
    frac[ens] <- f
  }
  # unknown schemes are all "^<prefix><digits>$", so one pass answers them
  gen <- !ens & !(known$scheme %in% names(SCHEMES))
  if (any(gen)) {
    stem <- sub("[0-9]+$", "", ids)
    stem[!grepl("^[^0-9]+[0-9]+$", ids)] <- NA_character_
    f <- as.numeric(table(stem)[known$prefix[gen]]) / length(ids)
    f[is.na(f)] <- 0
    frac[gen] <- f
  }
  for (i in which(!ens & !gen)) {
    frac[i] <- mean(grepl(.species_pattern(known$scheme[i], known$prefix[i]), ids))
  }
  known$frac <- frac
  # rank only species with a matching identifier
  known <- known[frac > 0, ]

  # prefixes are shared between species, so frac ties and a name tie-break would
  # decide by spelling; rank on reach, then gene-count fit
  num <- .id_number(ids, bare = ids) # ids are already bare here
  obs_max <- suppressWarnings(max(num, na.rm = TRUE))
  obs_n <- length(unique(ids))
  fp <- .fp()
  fp <- fp[fp$species %in% known$species, ]
  qmax_by <- vapply(
    split(fp$qmax, fp$species),
    function(q) suppressWarnings(max(q, na.rm = TRUE)), numeric(1)
  )
  fit_by <- vapply(
    split(abs(log1p(fp$n) - log1p(obs_n)), fp$species),
    function(d) suppressWarnings(min(d, na.rm = TRUE)), numeric(1)
  )
  q <- qmax_by[known$species]
  reach <- !is.finite(obs_max) | !is.finite(q) | q >= obs_max * QMAX_TOLERANCE
  fit <- unname(fit_by[known$species])

  # rank by fraction, prefix length, reach, gene-count fit and species name
  known <- known[order(
    -known$frac, -nchar(known$prefix), !reach, fit,
    known$species
  ), ]
  known$common <- .common_name(known$species)
  known <- known[, c("species", "common", setdiff(names(known), c("species", "common")))]
  rownames(known) <- NULL
  known
}

#' Which Ensembl release a set of identifiers came from
#'
#' Scores every release the index knows for this species, marks infeasible
#' those that could not have produced the input, and ranks the rest by fit. See
#' `vignette("how-it-works")` for the scoring.
#'
#' @param ids Character vector of gene identifiers.
#' @param species Species name, e.g. `"homo_sapiens"`. Detected when `NULL`.
#'
#' A best match at the oldest or newest indexed release raises a warning: the
#' index covers a finite span, and a set from outside it can only be scored
#' against the nearest edge.
#'
#' @return A data frame of candidate releases, best first, with `release`,
#'   `source` (`ensembl`, `gencode` for GENCODE's primary-assembly gene set, or
#'   `gencode_all` for its scaffold- and patch-inclusive one), `assembly`,
#'   `feasible`, `dist` and `n_index`. Each assembly has its own row within a
#'   release. All candidates are returned because adjacent releases can have
#'   indistinguishable fingerprints.
#' @examples
#' ids <- readLines(system.file("extdata", "example_ids.txt.gz",
#'   package = "genevintage"
#' ))
#' head(detect_release(ids), 3)
#' @export
detect_release <- function(ids, species = NULL) {
  ids <- unique(ids[nzchar(ids)])
  if (!is.null(species)) species <- resolve_species(species)
  if (is.null(species)) {
    cand <- detect_species(ids)
    if (!nrow(cand)) stop("no known species matches these identifiers.", call. = FALSE)
    species <- cand$species[1]
  }
  fp <- .fp()
  ref <- fp[fp$species == species, ]
  if (!nrow(ref)) stop("no fingerprints for species '", species, "'.", call. = FALSE)

  # score only this species' ids; a stray symbol like TP53 contributes 53
  pat <- .species_pattern(ref$scheme[1], ref$prefix[1])
  bare <- .bare(ids)
  own <- grepl(pat, bare)
  if (!any(own)) {
    stop(sprintf("none of the identifiers match %s.", species), call. = FALSE)
  }
  ids <- ids[own]
  bare <- bare[own]

  ordered <- all(ref$ordered)

  # reuse the bare identifiers when extracting their numeric parts
  num <- .id_number(ids, bare = bare)
  num <- num[is.finite(num)]
  if (ordered && !length(num)) {
    stop("no numeric identifiers to fingerprint.", call. = FALSE)
  }
  ver <- .id_version(ids)

  obs_n <- length(ids)
  obs_q90 <- if (length(num)) stats::quantile(num, 0.9, na.rm = TRUE) else NA_real_
  obs_qmax <- if (length(num)) max(num) else NA_real_
  obs_vmax <- suppressWarnings(max(ver, na.rm = TRUE))

  # a release can't predate its own ids/versions; allow QMAX_TOLERANCE slack
  # since Ensembl sometimes renumbers downward, and skip the check where that
  # breaks a species' own qmax monotonicity.
  feasible <- rep(TRUE, nrow(ref))
  if (ordered && !length(qmax_violations(ref))) {
    feasible <- !is.na(ref$qmax) & ref$qmax >= obs_qmax * QMAX_TOLERANCE
  }
  # Ensembl Genomes records no gene versions; an NA must not veto a release
  if (is.finite(obs_vmax)) {
    feasible <- feasible & (is.na(ref$vmax) | ref$vmax >= obs_vmax)
  }

  dist <- abs(log1p(ref$n) - log1p(obs_n))
  if (ordered) dist <- dist + abs(log1p(ref$q90) - log1p(obs_q90))

  out <- data.frame(
    species = species, source = ref$source,
    release = ref$release, assembly = ref$assembly,
    feasible = feasible, dist = round(dist, 4), n_index = ref$n
  )
  out <- out[order(!out$feasible, out$dist), ]
  rownames(out) <- NULL
  attr(out, "numeric_ordered") <- ordered
  if (!ordered) {
    out$feasible <- NA # feasibility is unknown
    message(sprintf(
      "%s identifiers are not numerically ordered; ranking by gene count alone.",
      ref$scheme[1]
    ))
  }
  # a best match at either end may mean the true release lies outside the index.
  # each source numbers its own releases, so the edge is that source's edge.
  best <- out$release[1]
  span <- ref$release[ref$source == out$source[1]]
  span <- span[order(.rel_num(span))]
  span <- span[c(1L, length(span))]
  if (best %in% span) {
    warning(
      sprintf(
        "best match is %s release %s, the %s of the indexed range (%s-%s); the true release may lie outside it.",
        out$source[1], best, if (best == span[1]) "oldest" else "newest", span[1], span[2]
      ),
      call. = FALSE
    )
  }
  if (ordered && !any(feasible)) {
    warning("no release is consistent with these ids -- likely a merged set.",
      call. = FALSE
    )
  }
  out
}
