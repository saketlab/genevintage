# The whole path: ids in, names out, nothing reordered.
test_that("geneid2name detects, fetches and converts", {
  skip_on_cran()
  skip_if_offline()
  skip_no_index()
  skip_if(!"saccharomyces_cerevisiae" %in% .fp()$species)

  ids <- c("YDL246C", "YDR387C", "NOT_A_GENE")
  out <- geneid2name(ids, species = "saccharomyces_cerevisiae", release = 63, quiet = TRUE)
  expect_length(out, 3)
  expect_equal(as.vector(out)[3], "NOT_A_GENE")
  expect_equal(attr(out, "species"), "saccharomyces_cerevisiae")
})

test_that("an explicit release skips detection entirely", {
  skip_on_cran()
  skip_if_offline()
  out <- geneid2name("YDL246C",
    species = "saccharomyces_cerevisiae",
    release = 63, quiet = TRUE
  )
  expect_equal(attr(out, "release"), 63)
})

test_that("a non-Ensembl scheme resolves end to end without an explicit species", {
  skip_on_cran()
  skip_if_offline()
  skip_no_index()
  # this is what used to fail: TAIR and SGD were undetectable
  expect_equal(detect_species("YDL246C")$species[1], "saccharomyces_cerevisiae")
  expect_equal(detect_species("AT1G01010")$species[1], "arabidopsis_thaliana")
})

test_that("an explicitly named source is honoured, not overwritten by detection", {
  skip_no_index()
  skip_if_not(
    file.exists(test_path("pipelines", "gdc_gencode_v22.txt.gz")),
    "no pipeline fixture"
  )
  ids <- readLines(gzfile(test_path("pipelines", "gdc_gencode_v22.txt.gz")), warn = FALSE)
  # these ids score best against GENCODE; a caller asking for Ensembl means it
  r <- suppressWarnings(suppressMessages(detect_release(ids, species = "human")))
  expect_equal(r$source[1], "gencode")
  expect_error(
    geneid2name(head(ids, 100),
      species = "human", source = "no_such_source",
      quiet = TRUE
    ),
    "no no_such_source release is indexed"
  )
})
