# The end-to-end path: identifiers in, gene names out.

#' Convert gene identifiers to gene names, detecting species and release
#'
#' Detects the species from the identifier prefixes, ranks the Ensembl releases
#' consistent with the identifier set, fetches that release's mapping (streamed
#' once, then cached) and converts. Identifiers with no gene name are returned
#' unchanged.
#'
#' @param ids Character vector of gene identifiers.
#' @param species Skip species detection by naming it. Any form
#'   [resolve_species()] accepts: `"homo_sapiens"`, `"Homo sapiens"`,
#'   `"Human"` or `"human"`.
#' @param release Skip release detection by naming it.
#' @param assembly Assembly suffix, e.g. `"37"` for human GRCh37. Detected
#'   alongside the release when both are left `NULL`.
#' @param source Annotation family: `"ensembl"`, or `"gencode"` /
#'   `"gencode_all"` for GENCODE's primary and scaffold-inclusive gene sets.
#'   Detected alongside the release.
#' @param mapping Skip detection and use this `id`/`name` table, as accepted
#'   by [gene_names()].
#' @param unique Passed to [gene_names()]. Defaults to `TRUE` when `ids` is a
#'   matrix, sparse `Matrix` or data frame, and to `FALSE` otherwise.
#' @param quiet Suppress the message reporting what was detected.
#' @param row_ids When `ids` is a matrix, sparse `Matrix` or data frame,
#'   overrides which identifiers to detect and map. Defaults to `rownames(ids)`;
#'   set this when `ids` has none.
#'
#' @return If `ids` is a plain vector: a character vector the same length as
#'   `ids`, carrying `mapped`, `species`, `source`, `release` and `assembly`
#'   attributes. If `ids` is a matrix, sparse `Matrix` or data frame: the same
#'   object with only its row names replaced by gene names; dimensions, row
#'   order and other attributes are preserved.
#'
#' @seealso [detect_species()], [detect_release()], [gene_names()].
#' @export
#' @examples
#' \dontrun{
#' ids <- c("ENSG00000141510", "ENSG00000284733", "ENSG00000012048")
#' geneid2name(ids)
#'
#' m <- matrix(1:6, nrow = 3, dimnames = list(ids, c("s1", "s2")))
#' geneid2name(m)
#' #>                 s1 s2
#' #> TP53             1  4
#' #> ENSG00000284733  2  5
#' #> BRCA1            3  6
#' }
geneid2name <- function(ids, species = NULL, release = NULL, assembly = NULL,
                        source = NULL, mapping = NULL, unique, quiet = FALSE,
                        row_ids = NULL) {
  if (!is.null(dim(ids))) {
    if (missing(unique)) unique <- TRUE
    return(.rename_rows(ids, row_ids, function(rid) {
      geneid2name(rid,
        species = species, release = release, assembly = assembly,
        source = source, mapping = mapping, unique = unique, quiet = quiet
      )
    }))
  }
  if (missing(unique)) unique <- FALSE
  ids <- .as_ids(ids)
  if (!is.null(mapping)) {
    return(gene_names(ids, mapping, unique = unique, warn = !quiet))
  }
  if (!is.null(species)) species <- resolve_species(species)
  if (is.null(species)) {
    cand <- detect_species(ids)
    if (!nrow(cand)) stop("cannot identify the species from these ids.", call. = FALSE)
    species <- cand$species[1]
    if (cand$frac[1] < 0.5) {
      warning(sprintf("only %.0f%% of ids look like %s.", 100 * cand$frac[1], species),
        call. = FALSE
      )
    }
    # report tied species sharing an identifier prefix, such as dog breeds
    tied <- sum(cand$frac >= cand$frac[1] - 1e-9)
    if (tied > 1) {
      message(sprintf(
        "%d species share these identifiers (%s); using %s.",
        tied, paste(utils::head(cand$species, 3), collapse = ", "), species
      ))
    }
    if (nrow(cand) > 1 && cand$frac[2] >= 0.1 && cand$prefix[2] != cand$prefix[1]) {
      warning(sprintf(
        "ids look like a mixture: %s (%.0f%%) and %s (%.0f%%).",
        cand$species[1], 100 * cand$frac[1],
        cand$species[2], 100 * cand$frac[2]
      ), call. = FALSE)
    }
  }
  if (is.null(release)) {
    rel <- detect_release(ids, species = species)
    # a named source or assembly wins over the best-scoring row
    if (!is.null(source)) {
      keep <- rel$source == source
      if (!any(keep)) {
        stop(sprintf("no %s release is indexed for %s.", source, species), call. = FALSE)
      }
      rel <- rel[keep, ]
    }
    if (!is.null(assembly)) {
      keep <- rel$assembly == as.character(assembly)
      if (!any(keep)) {
        stop(sprintf("no release of %s on assembly %s is indexed.", species, assembly),
          call. = FALSE
        )
      }
      rel <- rel[keep, ]
    }
    # nothing feasible means no release explains the set
    if (isTRUE(attr(rel, "numeric_ordered")) && !any(rel$feasible, na.rm = TRUE)) {
      stop("no release is consistent with these ids; pass `release` explicitly ",
        "if you want to map them anyway.",
        call. = FALSE
      )
    }
    best <- rel[1, ] # detect_release returns them in preference order
    release <- best$release
    assembly <- best$assembly
    source <- best$source
    if (!quiet) {
      feas <- rel[!is.na(rel$feasible) & rel$feasible, ]
      n_close <- if (nrow(feas)) sum(feas$dist <= feas$dist[1] * 1.1) else 0L
      message(sprintf(
        "Best fingerprint match: %s, %s release %s (assembly %s)%s", species, source, release, assembly,
        if (n_close > 1) sprintf(" -- %d candidates score within 10%%", n_close) else ""
      ))
    }
  }
  m <- fetch_mapping(species, release, assembly = assembly, source = source)
  out <- gene_names(ids, m, unique = unique)
  attr(out, "species") <- species
  attr(out, "source") <- source
  attr(out, "release") <- release
  attr(out, "assembly") <- assembly
  out
}
