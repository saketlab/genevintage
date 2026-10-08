# Naming a species: Ensembl name, scientific name, display name, or shorthand.

test_that("every accepted form of a name resolves to the Ensembl one", {
  skip_if(nrow(species_table()) == 0, "no species table")
  expect_equal(resolve_species("homo_sapiens"), "homo_sapiens")
  expect_equal(resolve_species("Homo sapiens"), "homo_sapiens")
  expect_equal(resolve_species("Human"), "homo_sapiens")
  expect_equal(resolve_species("human"), "homo_sapiens")
  expect_equal(resolve_species("HUMAN"), "homo_sapiens")
  expect_equal(resolve_species("Zebrafish"), "danio_rerio")
})

test_that("shorthands Ensembl's own display names do not cover still work", {
  skip_if(nrow(species_table()) == 0, "no species table")
  # "Norway rat - BN/NHsdMcwi", "Drosophila melanogaster - (Fruit fly)" and
  # "Caenorhabditis elegans (Nematode, N2)" are the published names; nobody
  # types those
  expect_equal(resolve_species("rat"), "rattus_norvegicus")
  expect_equal(resolve_species("fly"), "drosophila_melanogaster")
  expect_equal(resolve_species("worm"), "caenorhabditis_elegans")
  expect_equal(resolve_species("yeast"), "saccharomyces_cerevisiae")
})

test_that("'mouse' is the species, not one of its strains", {
  skip_if(nrow(species_table()) == 0, "no species table")
  # twenty-odd mus_musculus_* strain genomes have "Mouse" in their display name
  expect_equal(resolve_species("mouse"), "mus_musculus")
  expect_gt(nrow(species_table("mouse")), 5)
})

test_that("an ambiguous name is refused with its candidates, not guessed", {
  skip_if(nrow(species_table()) == 0, "no species table")
  expect_error(resolve_species("macaca"), "matches 3 species")
  expect_error(resolve_species("mus"), "name one exactly")
})

test_that("an unknown name says so", {
  expect_error(resolve_species("not_a_species"), "not a species in the index")
  expect_error(resolve_species(c("a", "b")), "single, non-missing, non-empty")
  expect_error(resolve_species(""), "single, non-missing, non-empty")
})

test_that("an unambiguous partial match is accepted and reported", {
  skip_if(nrow(species_table()) == 0, "no species table")
  expect_message(r <- resolve_species("Giant panda"), NA) # exact, no message
  expect_equal(r, "ailuropoda_melanoleuca")
  expect_message(resolve_species("panda"), "ailuropoda_melanoleuca")
})

test_that("detect_species reports the common name alongside the Ensembl one", {
  skip_no_index()
  d <- detect_species(c("ENSG00000141510", "ENSG00000012048"))
  expect_true(all(c("species", "common") %in% names(d)))
  expect_equal(d$species[1], "homo_sapiens")
  expect_equal(d$common[1], "Human")
})

test_that("common names are accepted wherever a species is", {
  skip_no_index()
  ids <- c("ENSG00000141510", "ENSG00000012048")
  # two ids sit at the edge of the indexed range, which warns; not the point here
  expect_identical(
    suppressWarnings(detect_release(ids, species = "human")),
    suppressWarnings(detect_release(ids, species = "homo_sapiens"))
  )
  expect_error(fetch_orthologs("human", "human", 116), "must differ")
})

test_that("a pattern that cannot filter is refused, not silently obeyed", {
  skip_if(nrow(species_table()) == 0, "no species table")
  # grepl() takes the first element of a longer pattern, errors on an empty one,
  # and returns NA for a missing one -- which indexed the table with NA and
  # manufactured rows of NA that read as real species
  for (p in list(character(0), NA_character_, c("mouse", "rat"), "", 42)) {
    expect_error(species_table(p), "single, non-missing, non-empty")
  }
})
