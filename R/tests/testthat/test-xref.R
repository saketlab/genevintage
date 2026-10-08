# Entrez cross-references. The parsing is exercised offline against a local
# file; the network path is covered by test-fetch.R's guards.

xref_file <- function(..., env = parent.frame()) {
  f <- tempfile(fileext = ".tsv.gz")
  withr::defer(unlink(f), envir = env)
  hdr <- c(
    "gene_stable_id", "transcript_stable_id", "protein_stable_id",
    "xref", "db_name", "info_type", "source_identity", "xref_identity",
    "linkage_type"
  )
  rows <- vapply(list(...), paste, character(1), collapse = "\t")
  writeLines(c(paste(hdr, collapse = "\t"), rows), gzfile(f))
  paste0("file://", f)
}

test_that("one row per transcript collapses to one row per pair", {
  u <- xref_file(
    c("ENSG1", "ENST1", "ENSP1", "111", "EntrezGene", "DEPENDENT", "-", "-", "-"),
    c("ENSG1", "ENST2", "ENSP2", "111", "EntrezGene", "DEPENDENT", "-", "-", "-"),
    c("ENSG2", "ENST3", "ENSP3", "222", "EntrezGene", "DEPENDENT", "-", "-", "-")
  )
  x <- stream_xrefs(u)
  expect_equal(nrow(x), 2L)
  expect_equal(x$id, c("ENSG1", "ENSG2"))
  expect_equal(x$xref, c("111", "222"))
})

test_that("a gene with several Entrez ids keeps all of them", {
  # Entrez is not a relabelling of Ensembl; one gene can carry several
  u <- xref_file(
    c("ENSG1", "ENST1", "ENSP1", "111", "EntrezGene", "DEPENDENT", "-", "-", "-"),
    c("ENSG1", "ENST1", "ENSP1", "999", "EntrezGene", "DEPENDENT", "-", "-", "-")
  )
  x <- stream_xrefs(u)
  expect_equal(nrow(x), 2L)
  expect_setequal(x$xref, c("111", "999"))
})

test_that("other cross-reference databases in the same file are ignored", {
  u <- xref_file(
    c("ENSG1", "ENST1", "ENSP1", "111", "EntrezGene", "DEPENDENT", "-", "-", "-"),
    c("ENSG2", "ENST2", "ENSP2", "P12345", "Uniprot/SWISSPROT", "DIRECT", "-", "-", "-")
  )
  x <- stream_xrefs(u)
  expect_equal(x$id, "ENSG1")
})

test_that("a file with no rows of the wanted database returns the right shape", {
  u <- xref_file(
    c("ENSG1", "ENST1", "ENSP1", "P1", "Uniprot/SWISSPROT", "DIRECT", "-", "-", "-")
  )
  x <- stream_xrefs(u)
  expect_s3_class(x, "data.frame")
  expect_equal(nrow(x), 0L)
  expect_named(x, c("id", "xref"))
})

test_that("dedup spans chunk boundaries", {
  # the same pair split across two reads must not survive twice
  u <- xref_file(
    c("ENSG1", "ENST1", "ENSP1", "111", "EntrezGene", "DEPENDENT", "-", "-", "-"),
    c("ENSG1", "ENST2", "ENSP2", "111", "EntrezGene", "DEPENDENT", "-", "-", "-"),
    c("ENSG1", "ENST3", "ENSP3", "111", "EntrezGene", "DEPENDENT", "-", "-", "-")
  )
  expect_equal(nrow(stream_xrefs(u, chunk = 1L)), 1L)
  expect_equal(nrow(stream_xrefs(u, chunk = 2L)), 1L)
})

test_that("the cache names an xref file by its own fields", {
  d <- withr::local_tempdir()
  withr::local_envvar(c(R_USER_CACHE_DIR = d))
  cd <- tools::R_user_dir("genevintage", "cache")
  dir.create(cd, recursive = TRUE)
  saveRDS(data.frame(id = "x"), file.path(cd, "xrefs-entrez-homo_sapiens-38-116.rds"))
  x <- gene_cache()
  expect_equal(x$kind, "xrefs")
  expect_equal(x$species, "homo_sapiens")
  expect_equal(x$release, "116")
})

test_that("entrez_ids works from a supplied mapping, offline", {
  m <- data.frame(
    id = c("ENSG1", "ENSG2"), name = c("TP53", "BRCA1"),
    chr = "17", biotype = "protein_coding", stringsAsFactors = FALSE
  )
  fake <- data.frame(
    id = c("ENSG1", "ENSG1"), xref = c("7157", "999"),
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    fetch_xrefs = function(...) fake,
    .newest_release = function(...) "116"
  )

  out <- entrez_ids(c("TP53", "BRCA1"), species = "human", mapping = m, quiet = TRUE)
  # every match is kept, and a gene with no cross-reference is NA rather than absent
  expect_equal(nrow(out), 3L)
  expect_setequal(out$entrez[out$input == "TP53"], c("7157", "999"))
  expect_true(is.na(out$entrez[out$input == "BRCA1"]))
  expect_named(out, c("input", "id", "name", "entrez"))
})

test_that("an input matching no gene is absent rather than a row of NA", {
  m <- data.frame(
    id = "ENSG1", name = "TP53", chr = "17",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    fetch_xrefs = function(...) {
      data.frame(
        id = "ENSG1", xref = "7157",
        stringsAsFactors = FALSE
      )
    },
    .newest_release = function(...) "116"
  )
  out <- entrez_ids(c("TP53", "NOT_A_GENE"),
    species = "human", mapping = m,
    quiet = TRUE
  )
  expect_equal(out$input, "TP53")
})

test_that("an input matching nothing returns the empty shape without fetching", {
  m <- data.frame(
    id = "ENSG1", name = "TP53", chr = "17",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  local_mocked_bindings(fetch_xrefs = function(...) stop("must not fetch"))
  out <- entrez_ids("NOT_A_GENE", species = "human", mapping = m, quiet = TRUE)
  expect_equal(nrow(out), 0L)
  expect_named(out, c("input", "id", "name", "entrez"))
  expect_type(out$entrez, "character")
})

test_that("an Ensembl vintage keeps its release; a GENCODE one falls back to newest", {
  m <- data.frame(
    id = "ENSG1", name = "TP53", chr = "17",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  x <- data.frame(id = "ENSG1", xref = "7157", stringsAsFactors = FALSE)
  asked <- NULL

  # cross-references are published against Ensembl releases, so a GENCODE
  # vintage has no file of its own
  local_mocked_bindings(
    .locate_annotation = function(...) {
      list(
        species = "homo_sapiens",
        source = "ensembl", release = "112",
        assembly = "38", map = m
      )
    },
    .newest_release = function(...) stop("must not need the newest release"),
    fetch_xrefs = function(species, release, ...) {
      asked <<- release
      x
    }
  )
  entrez_ids("TP53", species = "human", quiet = TRUE)
  expect_equal(asked, "112")

  local_mocked_bindings(
    .locate_annotation = function(...) {
      list(
        species = "homo_sapiens",
        source = "gencode", release = "44",
        assembly = "38", map = m
      )
    },
    .newest_release = function(...) "116",
    fetch_xrefs = function(species, release, ...) {
      asked <<- release
      x
    }
  )
  entrez_ids("TP53", species = "human", quiet = TRUE)
  expect_equal(asked, "116")
})

test_that("rows come back in the caller's order", {
  m <- data.frame(
    id = c("ENSG1", "ENSG2"), name = c("TP53", "BRCA1"),
    chr = "17", biotype = "protein_coding", stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    fetch_xrefs = function(...) {
      data.frame(
        id = c("ENSG2", "ENSG1"),
        xref = c("672", "7157"),
        stringsAsFactors = FALSE
      )
    },
    .newest_release = function(...) "116"
  )
  out <- entrez_ids(c("BRCA1", "TP53"), species = "human", mapping = m, quiet = TRUE)
  expect_equal(out$input, c("BRCA1", "TP53"))
  expect_equal(out$entrez, c("672", "7157"))
  expect_equal(row.names(out), c("1", "2"))
})
