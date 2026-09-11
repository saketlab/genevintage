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
#'   Detected alongside the release. Human and mouse counts are usually produced
#'   against GENCODE, whose gene set is the smaller of the two.
#' @param unique Passed to [gene_names()].
#' @param quiet Suppress the message reporting what was detected.
#'
#' @return A character vector the same length as `ids`, carrying `mapped`,
#'   `species`, `source`, `release` and `assembly` attributes.
#'
#' @seealso [detect_species()], [detect_release()], [gene_names()].
#' @export
#' @examples
#' \dontrun{
#' ids <- c("ENSG00000141510", "ENSG00000284733", "ENSG00000012048")
#' geneid2name(ids)
#' }
geneid2name <- function(ids, species = NULL, release = NULL, assembly = NULL,
                        source = NULL, unique = FALSE, quiet = FALSE) {
  ids <- .as_ids(ids)
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
    # a second species with a real share means a merged set
    # dog breeds all use ENSCAFG; report the tie rather than a determination
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
        "%s, %s release %s (assembly %s)%s", species, source, release, assembly,
        if (n_close > 1) sprintf(" -- %d fit equally well", n_close) else ""
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
