#' Convert gene identifiers to gene names with species and release detection
#'
#' This function detects the species from the identifier prefixes. It ranks the
#' Ensembl releases that are consistent with the identifier set. It fetches the
#' mapping of the best release and converts the identifiers.
#'
#' @param ids Character vector of gene identifiers.
#' @param species Give this value to skip species detection. Any form that
#'   [resolve_species()] accepts: `"homo_sapiens"`, `"Homo sapiens"`, `"Human"`
#'   or `"human"`.
#' @param release Give this value to skip release detection.
#' @param assembly Assembly suffix, for example `"37"` for human GRCh37. The
#'   function detects it together with the release when both are `NULL`.
#' @param source Annotation family: `"ensembl"`, `"gencode"` for the primary
#'   gene set of GENCODE, or `"gencode_all"` for its gene set with scaffolds.
#'   The function detects it together with the release.
#' @param mapping Skip detection and use this `id`/`name` table, in the form
#'   that [gene_names()] accepts.
#' @param unique The function passes this to [gene_names()]. The default is
#'   `TRUE` when `ids` is a matrix, sparse `Matrix` or data frame, and `FALSE`
#'   otherwise.
#' @param quiet Do not show the message that reports the detected values. With
#'   `mapping`, also do not give the warnings from [gene_names()]. Without
#'   `mapping`, these warnings stay, because they can show a wrong detection.
#' @param row_ids When `ids` is a matrix, sparse `Matrix` or data frame, this
#'   argument overrides the identifiers to detect and map. The default is
#'   `rownames(ids)`. Set this argument when `ids` has no row names.
#'
#' @return If `ids` is a plain vector, a character vector of the same length as
#'   `ids`, with `mapped`, `species`, `source`, `release` and `assembly`
#'   attributes. If `ids` is a matrix, sparse `Matrix` or data frame, the same
#'   object with only its row names replaced by gene names. The function keeps
#'   the dimensions, row order and other attributes.
#'
#' @seealso [detect_species()], [detect_release()], [gene_names()].
#' @export
#' @examples
#' # an excerpt of human Ensembl 116; drop `mapping` to detect and fetch the real one
#' hs <- read.delim(system.file("extdata", "example_mapping.tsv", package = "genevintage"))
#' ids <- c("ENSG00000141510", "ENSG00000284616", "ENSG00000012048")
#' geneid2name(ids, mapping = hs) # the unnamed lncRNA keeps its id
#'
#' m <- matrix(1:6, nrow = 3, dimnames = list(ids, c("s1", "s2")))
#' geneid2name(m, mapping = hs)
#' @examplesIf interactive()
#' geneid2name(ids) # species and release detected, then fetched
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
