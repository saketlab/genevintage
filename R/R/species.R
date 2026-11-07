.species_tbl <- .once(function() {
  f <- .extdata("species_names.rds")
  if (!file.exists(f)) {
    # return the species-table schema when the bundled name table is absent
    return(data.frame(
      species = character(), common = character(),
      taxonomy_id = integer(), division = character(),
      stringsAsFactors = FALSE
    ))
  }
  readRDS(f)
})

#' A name from the user, or a message that gives the reason it is not valid
#'
#' Validate scalar names at public entry points. `.norm_name()` also handles
#' vectors from the species table.
#' @noRd
.as_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(sprintf("`%s` must be a single, non-missing, non-empty name.", arg),
      call. = FALSE
    )
  }
  trimws(x)
}

#' Which rows of the species table a normalised name matches
#' @noRd
.name_match <- function(tbl, n) {
  grepl(n, .norm_name(tbl$common), fixed = TRUE) |
    grepl(n, .norm_name(tbl$species), fixed = TRUE)
}

# lowercase, replace punctuation and collapse whitespace for name matching
.norm_name <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub(" +", " ", x))
}

# map common shorthand names to their Ensembl species identifiers
SPECIES_ALIASES <- c(
  human = "homo_sapiens", mouse = "mus_musculus", rat = "rattus_norvegicus",
  fly = "drosophila_melanogaster", fruitfly = "drosophila_melanogaster",
  `fruit fly` = "drosophila_melanogaster",
  worm = "caenorhabditis_elegans", nematode = "caenorhabditis_elegans",
  `c elegans` = "caenorhabditis_elegans",
  yeast = "saccharomyces_cerevisiae", `budding yeast` = "saccharomyces_cerevisiae",
  arabidopsis = "arabidopsis_thaliana", `thale cress` = "arabidopsis_thaliana",
  zebrafish = "danio_rerio", chicken = "gallus_gallus", pig = "sus_scrofa",
  dog = "canis_lupus_familiaris", cow = "bos_taurus", cattle = "bos_taurus",
  frog = "xenopus_tropicalis", `african clawed frog` = "xenopus_tropicalis",
  chimp = "pan_troglodytes", chimpanzee = "pan_troglodytes",
  macaque = "macaca_mulatta", `rhesus macaque` = "macaca_mulatta",
  rhesus = "macaca_mulatta", horse = "equus_caballus", sheep = "ovis_aries",
  rabbit = "oryctolagus_cuniculus", cat = "felis_catus",
  `guinea pig` = "cavia_porcellus", medaka = "oryzias_latipes",
  `zebra finch` = "taeniopygia_guttata", opossum = "monodelphis_domestica",
  platypus = "ornithorhynchus_anatinus", `naked mole rat` = "heterocephalus_glaber_female"
)

#' Resolve a species name
#'
#' This function accepts the Ensembl name (`"homo_sapiens"`), the scientific
#' name (`"Homo sapiens"`) or the Ensembl display name (`"Zebrafish"`,
#' `"Giant panda"`). It also accepts a common short name (`"human"`, `"mouse"`,
#' `"rat"`, `"fly"`, `"yeast"`). The match ignores case, punctuation and
#' underscores.
#'
#' The function accepts and reports an unambiguous partial match. It refuses an
#' ambiguous partial match and shows the candidates. The `"mouse"` alias selects
#' `mus_musculus` even when strain genomes also match the word.
#'
#' @param x A species name in one of the accepted forms.
#' @param quiet Do not show the message that reports a partial match.
#'
#' @return The Ensembl species name, for example `"homo_sapiens"`.
#' @export
#' @examples
#' resolve_species("human")
#' resolve_species("Zebrafish")
resolve_species <- function(x, quiet = FALSE) {
  x <- .as_name(x, "species")
  tbl <- .species_tbl()
  if (x %in% tbl$species) {
    return(x)
  } # already an Ensembl name

  n <- .norm_name(x)
  # An Ensembl name with the underscores knocked out ("homo sapiens") is also
  # the scientific name, so this covers both.
  hit <- tbl$species[.norm_name(tbl$species) == n]
  if (length(hit) == 1L) {
    return(hit)
  }

  hit <- tbl$species[.norm_name(tbl$common) == n]
  if (length(hit) == 1L) {
    return(hit)
  }

  # [[ throws when the name is absent; [ gives NA
  a <- unname(SPECIES_ALIASES[n])
  if (!is.na(a) && a %in% tbl$species) {
    return(a)
  }

  # last resort: an unambiguous substring
  part <- .name_match(tbl, n)
  hit <- tbl$species[part]
  if (length(hit) == 1L) {
    if (!quiet) message(sprintf("'%s' -> %s (%s)", x, hit, tbl$common[part]))
    return(hit)
  }
  if (length(hit) > 1L) {
    stop(sprintf(
      "'%s' matches %d species (%s%s); name one exactly.", x, length(hit),
      paste(utils::head(hit, 4), collapse = ", "),
      if (length(hit) > 4) ", ..." else ""
    ), call. = FALSE)
  }
  stop(sprintf("'%s' is not a species in the index. See species_table().", x),
    call. = FALSE
  )
}

#' Species the index covers, with their common names
#'
#' @param pattern Optional. Keep only species whose Ensembl name or common name
#'   matches. The match ignores case and punctuation.
#'
#' @return A data frame of `species`, `common`, `taxonomy_id` and `division`.
#' @export
#' @examples
#' head(species_table("mouse"))
species_table <- function(pattern = NULL) {
  tbl <- .species_tbl()
  if (!is.null(pattern)) {
    tbl <- tbl[.name_match(tbl, .norm_name(.as_name(pattern, "pattern"))), ]
    rownames(tbl) <- NULL
  }
  tbl
}

#' The common name of a species, or the Ensembl name when there is none
#' @noRd
.common_name <- function(species) {
  tbl <- .species_tbl()
  .coalesce(tbl$common[match(species, tbl$species)], species)
}
