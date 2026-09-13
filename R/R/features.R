# Gene sets selected by sequence, biotype and symbol patterns.

# Mitochondrial sequence names vary by species: MT for vertebrates, Mito for
# yeast, MtDNA for C. elegans, Mt for Arabidopsis,
# mitochondrion_genome for Drosophila, chrM for UCSC-style references.
MITO_SEQNAMES <- "^(chr)?(MT|MtDNA|Mito|mitochondrion(_genome)?|M)$"

# RPL/RPS in human and yeast, Rpl/Rps in mouse, RpL/RpS in fly, rpl-/rps- in
# worm, so match case-insensitively. the mitoribosome takes an M prefix.
RIBO_PROTEIN_RX <- "^RP[SL]"
MITORIBO_PROTEIN_RX <- "^MRP[SL]"

# the S6 kinases (RPS6KA1..RPS6KL1) and RPS19BP1 match the prefix and are
# protein-coding, so a name rule removes them.
RIBO_PROTEIN_NOT_RX <- "^RP[SL][0-9]*K|BP[0-9]*$"

# Ensembl calls the mitochondrial set Mt_tRNA and the nuclear one tRNA.
TRNA_BIOTYPES <- "^(Mt_)?tRNA$"

# The RNA components of the ribosome. rRNA_pseudogene belongs with the
# pseudogenes.
RRNA_BIOTYPES <- "^(Mt_)?rRNA$"

# letters restricted to ABDEGMQZ to exclude HBEGF/HBS1L/HBP1; the three tail
# shapes cover human (HBA1), mouse (Hba-a1) and zebrafish (hbaa1) spellings.
HEMOGLOBIN_RX <- "^HB[ABDEGMQZ]([0-9]*|[A-Za-z][0-9]*|-[A-Za-z0-9]+)$"

# H2- (mouse) keeps its hyphen to avoid matching the H2A histones; HLA-
# (human) and mhc1/mhc2 (zebrafish) round out the per-clade spellings.
MHC_RX <- "^(HLA-|H2-[A-Za-z]|MHC[12])"

# match lncRNA biotype spellings across species and releases
LNCRNA_BIOTYPES <- paste0(
  "^(lncRNA|lincRNA|antisense(_RNA)?|macro_lncRNA|bidirectional_promoter_lncRNA",
  "|sense_intronic|sense_overlapping|processed_transcript",
  "|3prime_overlapping_ncrna)$"
)

#' The sex chromosomes present in one annotation
#'
#' Treat X, Y, Z and W as sex chromosomes. If IX or XI is present, exclude X
#' because Roman-numbered annotations use it as an autosome name.
#' @noRd
.sex_seqnames <- function(seqnames) {
  present <- unique(seqnames)
  roman <- any(c("IX", "XI") %in% present)
  want <- if (roman) c("Y", "Z", "W") else c("X", "Y", "Z", "W")
  bare <- sub("^chr", "", present)
  present[bare %in% want]
}

#' Which annotation to read a gene set out of
#'
#' With identifiers in hand, date them the way `geneid2name()` does. Without,
#' the caller must name the species and gets its newest indexed Ensembl release.
#' @noRd
.locate_annotation <- function(ids, species, release, assembly, source, quiet,
                               mapping = NULL) {
  # validate the supplied annotation's required columns and return it directly
  if (!is.null(mapping)) {
    if (!is.data.frame(mapping)) {
      stop("`mapping` must be a data frame.", call. = FALSE)
    }
    miss <- setdiff(c("id", "name", "chr", "biotype"), names(mapping))
    if (length(miss)) {
      stop("`mapping` needs column(s): ", paste(miss, collapse = ", "),
        ". See fetch_mapping() for the shape.",
        call. = FALSE
      )
    }
    return(list(
      species = if (is.null(species)) "the supplied mapping" else species,
      source = source, release = release, assembly = assembly,
      map = mapping
    ))
  }
  if (!is.null(species)) species <- resolve_species(species)
  if (is.null(ids) && is.null(species)) {
    stop("name a `species`, or pass `ids` to detect one from.", call. = FALSE)
  }
  # identify stable IDs once; gene symbols alone require a release fallback
  x <- if (is.null(ids)) character() else .as_ids(ids)
  own <- x[.looks_like_id(x)]
  datable <- length(own) > 0L
  if (is.null(species)) {
    if (!datable) {
      stop("gene names do not say which species they belong to; pass `species`.",
        call. = FALSE
      )
    }
    species <- detect_species(own)$species[1]
  }

  if (is.null(release)) {
    if (!datable) {
      release <- .newest_release(species)
      source <- "ensembl"
    } else {
      best <- suppressWarnings(detect_release(own, species = species))[1, ]
      release <- best$release
      assembly <- assembly %||% best$assembly
      source <- source %||% best$source
    }
  }
  if (!quiet) {
    message(sprintf("%s, %s release %s", species, source %||% "ensembl", release))
  }
  # return the selected annotation metadata with its mapping
  list(
    species = species, source = source %||% "ensembl", release = release,
    assembly = assembly,
    map = fetch_mapping(species, release, assembly = assembly, source = source)
  )
}

#' Pick the rows of an annotation a caller's tokens refer to
#'
#' Match each token as an identifier and as a gene name, preferring the ID
#' match. Return the caller's token in `input` for row subsetting and retain
#' the annotation's stable identifier in `id`.
#' @noRd
.restrict <- function(map, ids) {
  if (is.null(ids)) {
    return(map)
  }
  ids <- .as_ids(ids)
  by_id <- match(.bare(map$id), .bare(ids))
  by_nm <- match(map$name, ids)
  hit <- .coalesce(by_id, by_nm)
  keep <- map[!is.na(hit), , drop = FALSE]
  keep$input <- ids[hit[!is.na(hit)]]
  keep
}

#' Mitochondrial genes
#'
#' Selects genes by mitochondrial sequence names in the annotation for use in
#' mitochondrial count fractions. Sequence-name matching supports the naming
#' conventions of different species.
#'
#' @param ids Optional character vector of gene identifiers or gene names. When
#'   given, the species and release are detected from them and only these genes
#'   are considered. When omitted, `species` must be named and the whole
#'   mitochondrial gene set for its newest indexed release is returned.
#' @param species,release,assembly,source Skip detection by naming them.
#'   `species` accepts anything [resolve_species()] does.
#' @param mapping Supply an annotation directly to bypass detection and
#'   downloading. A data frame with at least `id`, `name`, `chr` and `biotype`,
#'   as [fetch_mapping()] returns; supports local GTFs and custom builds.
#' @param quiet Suppress the message reporting what was detected.
#'
#' @return A data frame of `id`, `name`, `chr` and `biotype`, one row per
#'   mitochondrial gene. When `ids` were given, an `input` column comes first
#'   holding the caller's own token, identifier or gene name as passed, so
#'   `rownames(m) %in% mito_genes(rownames(m))$input` subsets a
#'   matrix directly whichever kind its row names are.
#'
#' @seealso [sex_genes()] for the sex chromosomes.
#' @export
#' @examples
#' \dontrun{
#' mito_genes(species = "human") # the mitochondrial genes
#' mt <- mito_genes(rownames(counts)) # just those in your matrix
#' colSums(counts[mt$input, ]) / colSums(counts) # the mitochondrial fraction
#' }
mito_genes <- function(ids = NULL, species = NULL, release = NULL,
                       assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  .gene_subset(ids, species, release, assembly, source, quiet,
    mapping = mapping,
    function(sub, all, sp) grepl(MITO_SEQNAMES, sub$chr, ignore.case = TRUE)
  )
}

#' Sex-chromosome genes
#'
#' Sex chromosome names vary by clade: mammals use X and Y, birds Z and W,
#' and C. elegans X. Annotations without recognised sex chromosomes return
#' no rows.
#'
#' @inheritParams mito_genes
#' @param which Restrict to one chromosome: `"X"`, `"Y"`, `"Z"` or `"W"`.
#'   Defaults to every sex chromosome the annotation has.
#'
#' @return A data frame of `id`, `name`, `chr` and `biotype`, one row per gene,
#'   shaped as [mito_genes()] returns.
#'
#' @seealso [mito_genes()].
#' @export
#' @examples
#' \dontrun{
#' sex_genes(species = "human", which = "Y") # Y-linked genes
#' table(sex_genes(rownames(counts))$chr) # X vs Y in your matrix
#' }
sex_genes <- function(ids = NULL, which = NULL, species = NULL, release = NULL,
                      assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  .gene_subset(ids, species, release, assembly, source, quiet, mapping = mapping, function(sub, all, sp) {
    # use the whole annotation so Roman-numbered chromosomes can identify X
    keep <- .sex_seqnames(all$chr)
    if (!is.null(which)) {
      which <- toupper(which)
      bad <- setdiff(which, sub("^chr", "", keep))
      if (length(bad)) {
        stop(
          sprintf(
            "%s has no chromosome %s; it has %s.", sp,
            paste(bad, collapse = ", "),
            if (length(keep)) paste(keep, collapse = ", ") else "none"
          ),
          call. = FALSE
        )
      }
      keep <- keep[sub("^chr", "", keep) %in% which]
    }
    sub$chr %in% keep
  })
}


#' One gene set, however it is picked
#'
#' `pick` receives the caller's subset, the whole annotation (for questions that
#' need the full picture, like which chromosomes exist) and the species.
#' @noRd
.gene_subset <- function(ids, species, release, assembly, source, quiet,
                         mapping = NULL, pick) {
  a <- .locate_annotation(ids, species, release, assembly, source, quiet, mapping)
  sub <- .restrict(a$map, ids)
  out <- sub[pick(sub, a$map, a$species), , drop = FALSE]
  rownames(out) <- NULL
  cols <- c(if (!is.null(ids)) "input", "id", "name", "chr", "biotype")
  # carry the coordinates through when the annotation has them
  cols <- c(cols, intersect(c("start", "end", "strand", "span"), names(out)))
  out[, cols]
}

#' Group Ensembl biotypes by coding, pseudogene and immune status
#'
#' Test `pseudogene` before the immunoglobulin prefix so `IG_V_pseudogene`
#' lands with pseudogenes. Missing biotypes, TEC and artifact entries are
#' `other`; remaining biotypes are `noncoding`.
#' @noRd
.biotype_class <- function(x) {
  # classify distinct biotypes once when the input exceeds 100 entries
  if (length(x) > 100L) {
    u <- unique(x)
    return(.biotype_class(u)[match(x, u)])
  }
  ifelse(is.na(x) | !nzchar(x), "other",
    ifelse(x == "protein_coding", "coding",
      ifelse(grepl("pseudogene", x), "pseudogene",
        ifelse(grepl("^(IG|TR)_", x), "immune",
          ifelse(x %in% c("TEC", "artifact"), "other", "noncoding")
        )
      )
    )
  )
}

BIOTYPE_CLASSES <- c("coding", "noncoding", "pseudogene", "immune", "other")

#' What biotypes an annotation contains
#'
#' Use this to see what is there before filtering: the vocabulary differs by
#' species and drifts between releases.
#'
#' @inheritParams mito_genes
#' @return A data frame of `biotype`, `class` and `n`, commonest first.
#' @seealso [genes_by_biotype()], [coding_genes()], [noncoding_genes()].
#' @export
#' @examples
#' \dontrun{
#' gene_biotypes(species = "human")
#' gene_biotypes(rownames(counts)) # only what your matrix holds
#' }
gene_biotypes <- function(ids = NULL, species = NULL, release = NULL,
                          assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  a <- .locate_annotation(ids, species, release, assembly, source, quiet, mapping)
  b <- .restrict(a$map, ids)$biotype
  tb <- sort(table(b, useNA = "ifany"), decreasing = TRUE)
  out <- data.frame(
    biotype = names(tb), class = .biotype_class(names(tb)),
    n = as.integer(tb), stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

#' Genes of a given biotype, or class of biotypes
#'
#' @inheritParams mito_genes
#' @param biotype One or more Ensembl biotypes (`"lncRNA"`, `"miRNA"`,
#'   `"snoRNA"`, ...) or class names: `"coding"`, `"noncoding"`, `"pseudogene"`,
#'   `"immune"` (immunoglobulin and T-cell receptor segments) or `"other"`.
#'   [gene_biotypes()] lists what an annotation actually has.
#'
#' @return A data frame of `id`, `name`, `chr` and `biotype`, shaped as
#'   [mito_genes()] returns.
#' @seealso [gene_biotypes()] to see the vocabulary first.
#' @export
#' @examples
#' \dontrun{
#' genes_by_biotype("lncRNA", species = "human")
#' genes_by_biotype(c("miRNA", "snoRNA"), rownames(counts))
#' genes_by_biotype("pseudogene", rownames(counts)) # a whole class
#' }
genes_by_biotype <- function(biotype, ids = NULL, species = NULL, release = NULL,
                             assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  if (!is.character(biotype) || !length(biotype)) {
    stop("`biotype` must name at least one biotype or class.", call. = FALSE)
  }
  .gene_subset(ids, species, release, assembly, source, quiet, mapping = mapping, function(sub, all, sp) {
    cls <- intersect(biotype, BIOTYPE_CLASSES)
    exact <- setdiff(biotype, BIOTYPE_CLASSES)
    unknown <- setdiff(exact, all$biotype)
    if (length(unknown)) {
      stop(sprintf(
        "%s has no biotype %s. See gene_biotypes().", sp,
        paste(unknown, collapse = ", ")
      ), call. = FALSE)
    }
    sub$biotype %in% exact | .biotype_class(sub$biotype) %in% cls
  })
}

#' Protein-coding genes
#'
#' Selects the `protein_coding` biotype.
#' Immunoglobulin and T-cell receptor segments encode protein too but carry
#' their own biotypes; ask for `"immune"` through [genes_by_biotype()] to add
#' them.
#'
#' @inheritParams mito_genes
#' @return A data frame of `id`, `name`, `chr` and `biotype`.
#' @seealso [noncoding_genes()], [genes_by_biotype()].
#' @export
#' @examples
#' \dontrun{
#' nrow(coding_genes(species = "human"))
#' counts <- counts[coding_genes(rownames(counts))$input, ]
#' }
coding_genes <- function(ids = NULL, species = NULL, release = NULL,
                         assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  genes_by_biotype("coding", ids, species, release, assembly, source, mapping, quiet)
}

#' Non-coding RNA genes
#'
#' Genes transcribing functional RNA: lncRNA, miRNA, snRNA, snoRNA, rRNA and
#' the rest. Pseudogenes are a separate class; ask for `"pseudogene"` through
#' [genes_by_biotype()] when you want them.
#'
#' @inheritParams mito_genes
#' @return A data frame of `id`, `name`, `chr` and `biotype`.
#' @seealso [coding_genes()], [genes_by_biotype()].
#' @export
#' @examples
#' \dontrun{
#' table(noncoding_genes(rownames(counts))$biotype)
#' }
noncoding_genes <- function(ids = NULL, species = NULL, release = NULL,
                            assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  genes_by_biotype("noncoding", ids, species, release, assembly, source, mapping, quiet)
}


#' Ribosomal genes
#'
#' `which` selects the ribosomal set:
#'
#' * `"protein"` (the default): the genes encoding the ribosome's proteins,
#'   the set behind a single-cell ribosomal-content metric.
#' * `"rrna"`: the RNA components, by biotype (`rRNA` and `Mt_rRNA`).
#' * `"mito"`: the mitoribosome's proteins, named with an `M` prefix. These are
#'   nuclear genes; [mito_genes()] does not return them.
#' * `"all"`: every one of the above.
#'
#' The `"^RP[SL]"` symbol pattern is restricted to `protein_coding` genes to
#' exclude pseudogenes and lncRNAs. A name rule also excludes S6 kinases
#' (`RPS6KA1` to `RPS6KL1`) and `RPS19BP1`.
#'
#' @inheritParams mito_genes
#' @param which Which ribosomal set to return. See details.
#'
#' @return A data frame of `id`, `name`, `chr` and `biotype`, shaped as
#'   [mito_genes()] returns.
#'
#' @seealso [mito_genes()], [genes_by_biotype()].
#' @export
#' @examples
#' \dontrun{
#' rb <- ribosomal_genes(rownames(counts))
#' colSums(counts[rb$input, ]) / colSums(counts) # ribosomal content
#' ribosomal_genes(species = "human", which = "rrna")
#' }
ribosomal_genes <- function(ids = NULL, which = c("protein", "rrna", "mito", "all"),
                            species = NULL, release = NULL, assembly = NULL,
                            source = NULL, mapping = NULL, quiet = FALSE) {
  which <- match.arg(which)
  .gene_subset(ids, species, release, assembly, source, quiet, mapping = mapping, function(sub, all, sp) {
    nm <- ifelse(is.na(sub$name), "", sub$name)
    coding <- !is.na(sub$biotype) & sub$biotype == "protein_coding"
    hit <- rep(FALSE, nrow(sub))
    if (which %in% c("protein", "all")) {
      # select coding RP-prefixed genes, excluding kinase and binding-protein names
      hit <- hit | (coding & grepl(RIBO_PROTEIN_RX, nm, ignore.case = TRUE) &
        !grepl(RIBO_PROTEIN_NOT_RX, nm, ignore.case = TRUE))
    }
    if (which %in% c("mito", "all")) {
      hit <- hit | (coding & grepl(MITORIBO_PROTEIN_RX, nm, ignore.case = TRUE))
    }
    if (which %in% c("rrna", "all")) {
      hit <- hit | (!is.na(sub$biotype) & grepl(RRNA_BIOTYPES, sub$biotype))
    }
    hit
  })
}


#' Transfer RNA genes
#'
#' Ensembl names them `Mt_tRNA` on the mitochondrial chromosome and `tRNA` in
#' the nucleus. Which of those an annotation actually contains varies: the
#' vertebrate GTFs carry only the mitochondrial tRNAs, since nuclear tRNA genes
#' are predicted separately and lie outside the gene set. Invertebrate, fungal
#' and plant annotations include the nuclear ones. The `chr` and `biotype`
#' columns say which kind each row is.
#'
#' @inheritParams mito_genes
#' @return A data frame of `id`, `name`, `chr` and `biotype`, shaped as
#'   [mito_genes()] returns.
#' @seealso [ribosomal_genes()], [genes_by_biotype()].
#' @export
#' @examples
#' \dontrun{
#' trna_genes(species = "human") # the mitochondrial tRNAs
#' table(trna_genes(species = "yeast")$chr) # nuclear, by chromosome
#' }
trna_genes <- function(ids = NULL, species = NULL, release = NULL,
                       assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  .gene_subset(ids, species, release, assembly, source, quiet,
    mapping = mapping,
    function(sub, all, sp) {
      !is.na(sub$biotype) & grepl(TRNA_BIOTYPES, sub$biotype)
    }
  )
}


#' Ribosomal RNA genes
#'
#' The RNA components of the ribosome, by biotype: `rRNA` in the nucleus and
#' `Mt_rRNA` in the mitochondrion.
#'
#' `rRNA_pseudogene` belongs with the pseudogenes. Ask for it with
#' `genes_by_biotype("rRNA_pseudogene")`.
#'
#' @inheritParams mito_genes
#' @return A data frame of `id`, `name`, `chr` and `biotype`, shaped as
#'   [mito_genes()] returns.
#' @seealso [ribosomal_genes()] for the ribosome's *proteins*, [trna_genes()].
#' @export
#' @examples
#' \dontrun{
#' rrna_genes(species = "human")
#' table(rrna_genes(rownames(counts))$biotype)
#' }
rrna_genes <- function(ids = NULL, species = NULL, release = NULL,
                       assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  .gene_subset(ids, species, release, assembly, source, quiet,
    mapping = mapping,
    function(sub, all, sp) {
      !is.na(sub$biotype) & grepl(RRNA_BIOTYPES, sub$biotype)
    }
  )
}

#' Hemoglobin subunit genes
#'
#' The set behind a haemoglobin-contamination metric in single-cell work, beside
#' [mito_genes()] and [ribosomal_genes()]. Named per clade: `HBA1`/`HBB` in
#' human, `Hba-a1`/`Hbb-bs` in mouse, `hbaa1`/`hbba1` in zebrafish.
#'
#' Subunit-letter matching excludes `HBEGF`, `HBS1L` and `HBP1`.
#' Requiring `protein_coding` also excludes haemoglobin pseudogenes such as
#' `HBAP1`, `HBBP1` and `HBZP1`.
#'
#' @inheritParams mito_genes
#' @return A data frame of `id`, `name`, `chr` and `biotype`, shaped as
#'   [mito_genes()] returns.
#' @seealso [mito_genes()], [ribosomal_genes()].
#' @export
#' @examples
#' \dontrun{
#' hb <- hemoglobin_genes(rownames(counts))
#' colSums(counts[hb$input, ]) / colSums(counts) # haemoglobin fraction
#' }
hemoglobin_genes <- function(ids = NULL, species = NULL, release = NULL,
                             assembly = NULL, source = NULL, mapping = NULL,
                             quiet = FALSE) {
  .gene_subset(ids, species, release, assembly, source, quiet,
    mapping = mapping,
    function(sub, all, sp) {
      nm <- ifelse(is.na(sub$name), "", sub$name)
      !is.na(sub$biotype) & sub$biotype == "protein_coding" &
        grepl(HEMOGLOBIN_RX, nm, ignore.case = TRUE)
    }
  )
}

#' Major histocompatibility complex genes
#'
#' The MHC, which human calls HLA. Every clade names it differently: `HLA-A` in
#' human, `H2-K1` in mouse, `mhc1uba` in zebrafish, so matching `"^HLA-"` finds
#' nothing outside primates.
#'
#' Requires `protein_coding` and an MHC family name pattern. The patterns
#' exclude regional genes such as `B2M`, `TAP1`, `TAP2` and `MICA`.
#' Annotations without matching genes return an empty result.
#'
#' @inheritParams mito_genes
#' @return A data frame of `id`, `name`, `chr` and `biotype`, shaped as
#'   [mito_genes()] returns.
#' @seealso [genes_by_biotype()] for the immunoglobulin and T-cell receptor
#'   segments, which are a different question.
#' @export
#' @examples
#' \dontrun{
#' mhc_genes(species = "human") # coding HLA genes
#' mhc_genes(species = "mouse") # H2-K1, H2-D1, H2-Aa, ...
#' }
mhc_genes <- function(ids = NULL, species = NULL, release = NULL,
                      assembly = NULL, source = NULL, mapping = NULL,
                      quiet = FALSE) {
  .gene_subset(ids, species, release, assembly, source, quiet,
    mapping = mapping,
    function(sub, all, sp) {
      nm <- ifelse(is.na(sub$name), "", sub$name)
      hit <- !is.na(sub$biotype) & sub$biotype == "protein_coding" &
        grepl(MHC_RX, nm, ignore.case = TRUE)
      if (!any(hit) && !quiet) {
        message(
          "no MHC genes in this annotation for ", sp,
          "; only jawed vertebrates have one."
        )
      }
      hit
    }
  )
}

#' Long non-coding RNA genes
#'
#' Matches lncRNA biotype spellings across species and releases, including
#' `lncRNA`, `lincRNA`, `antisense`, `processed_transcript` and
#' `3prime_overlapping_ncrna`.
#'
#' Generic `ncRNA` entries are excluded. If the annotation contains them and
#' the subset has no lncRNA matches, a message points to
#' `genes_by_biotype("ncRNA")`.
#'
#' @inheritParams mito_genes
#' @return A data frame of `id`, `name`, `chr` and `biotype`, shaped as
#'   [mito_genes()] returns.
#' @seealso [noncoding_genes()] for every non-coding class, [gene_biotypes()].
#' @export
#' @examples
#' \dontrun{
#' nrow(lncrna_genes(species = "human"))
#' table(lncrna_genes(species = "zebrafish")$biotype) # antisense and lincRNA
#' }
lncrna_genes <- function(ids = NULL, species = NULL, release = NULL,
                         assembly = NULL, source = NULL, mapping = NULL, quiet = FALSE) {
  .gene_subset(ids, species, release, assembly, source, quiet,
    mapping = mapping,
    function(sub, all, sp) {
      hit <- !is.na(sub$biotype) & grepl(LNCRNA_BIOTYPES, sub$biotype)
      if (!any(hit) && !quiet && any(all$biotype == "ncRNA", na.rm = TRUE)) {
        message(sprintf(
          "%s files non-coding RNA under a generic 'ncRNA' biotype and does not mark lncRNA separately; genes_by_biotype(\"ncRNA\") returns those %d genes.",
          sp, sum(all$biotype == "ncRNA", na.rm = TRUE)
        ))
      }
      hit
    }
  )
}
