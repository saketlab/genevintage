# The end-to-end path, exercised offline. Its network tests skip on CRAN, which
# left the headline function at 18% coverage; mocking the fetch covers the
# branch logic without a server.

HS <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000284733")

map_of <- function(id, name) {
  data.frame(
    id = id, name = name, chr = "1", biotype = "protein_coding",
    span = 1L, stringsAsFactors = FALSE
  )
}

# a ranked table shaped like detect_release()'s, with the ordering attribute
rel_table <- function(..., ordered = TRUE) {
  x <- do.call(rbind, lapply(list(...), function(r) {
    data.frame(
      source = r[[1]], release = r[[2]], assembly = r[[3]],
      feasible = r[[4]], dist = r[[5]], stringsAsFactors = FALSE
    )
  }))
  attr(x, "numeric_ordered") <- ordered
  x
}

test_that("a named species and release skip detection entirely", {
  seen <- NULL
  local_mocked_bindings(
    detect_species = function(...) stop("must not detect the species"),
    detect_release = function(...) stop("must not detect the release"),
    fetch_mapping = function(species, release, assembly = NULL, source = NULL, ...) {
      seen <<- list(species = species, release = release)
      map_of(HS[1:2], c("TP53", "BRCA1"))
    }
  )
  out <- geneid2name(HS, species = "Human", release = 116, quiet = TRUE)

  expect_equal(seen$species, "homo_sapiens") # resolved, not passed through
  expect_equal(seen$release, 116)
  # one value per input, in order, with the id retained where there is no name
  expect_equal(as.vector(out), c("TP53", "BRCA1", "ENSG00000284733"))
  expect_equal(attr(out, "mapped"), c(TRUE, TRUE, FALSE))
  expect_equal(attr(out, "species"), "homo_sapiens")
})

test_that("identifiers no species explains are refused before anything is fetched", {
  local_mocked_bindings(
    fetch_mapping = function(...) stop("must not fetch")
  )
  expect_error(
    geneid2name(c("not_an_id", "also_not"), quiet = TRUE),
    "cannot identify the species"
  )
})

test_that("a minority of recognisable identifiers warns twice, at both steps", {
  local_mocked_bindings(
    fetch_mapping = function(...) map_of(HS[1], "TP53")
  )
  w <- character()
  withCallingHandlers(
    geneid2name(c(HS[1], "junk1", "junk2", "junk3"), release = 116, quiet = TRUE),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  # detection says the ids are not mostly this species; naming says they did not map
  expect_match(w[1], "of ids look like homo_sapiens")
  expect_match(w[2], "of ids mapped")
})

test_that("a named source picks its own release, not the best-scoring row", {
  seen <- NULL
  local_mocked_bindings(
    detect_release = function(...) {
      rel_table(
        list("gencode", "44", "38", TRUE, 0.01), # ranks first
        list("ensembl", "112", "38", TRUE, 0.02)
      )
    },
    fetch_mapping = function(species, release, assembly = NULL, source = NULL, ...) {
      seen <<- list(release = release, source = source, assembly = assembly)
      map_of(HS[1], "TP53")
    }
  )
  out <- geneid2name(HS[1], species = "human", source = "ensembl", quiet = TRUE)
  expect_equal(seen$source, "ensembl")
  expect_equal(seen$release, "112")
  expect_equal(attr(out, "source"), "ensembl")
})

test_that("a source or assembly the index does not hold is an error, not a silent fallback", {
  local_mocked_bindings(
    detect_release = function(...) rel_table(list("ensembl", "116", "38", TRUE, 0.01)),
    fetch_mapping = function(...) stop("must not fetch")
  )
  expect_error(
    geneid2name(HS[1], species = "human", source = "refseq", quiet = TRUE),
    "no refseq release is indexed"
  )
  expect_error(
    geneid2name(HS[1], species = "human", assembly = "37", quiet = TRUE),
    "no release of homo_sapiens on assembly 37"
  )
})

test_that("an assembly is matched however it is typed", {
  seen <- NULL
  local_mocked_bindings(
    detect_release = function(...) {
      rel_table(
        list("ensembl", "112", "38", TRUE, 0.01),
        list("ensembl", "112", "37", TRUE, 0.05)
      )
    },
    fetch_mapping = function(species, release, assembly = NULL, ...) {
      seen <<- assembly
      map_of(HS[1], "TP53")
    }
  )
  geneid2name(HS[1], species = "human", assembly = 37, quiet = TRUE) # numeric
  expect_equal(seen, "37")
})

test_that("no feasible release refuses rather than mapping against the least bad", {
  local_mocked_bindings(
    detect_release = function(...) {
      rel_table(
        list("ensembl", "116", "38", FALSE, 0.01),
        list("ensembl", "112", "38", NA, 0.02)
      )
    },
    fetch_mapping = function(...) stop("must not fetch")
  )
  expect_error(
    geneid2name(HS[1], species = "human", quiet = TRUE),
    "no release is consistent with these ids"
  )
})

test_that("an unordered scheme maps anyway, since feasibility cannot be judged", {
  local_mocked_bindings(
    detect_release = function(...) {
      rel_table(
        list("ensembl", "116", "4", NA, 0.01),
        ordered = FALSE
      )
    },
    fetch_mapping = function(...) map_of("YDL246C", "SOR2")
  )
  expect_silent(out <- geneid2name("YDL246C", species = "yeast", quiet = TRUE))
  expect_equal(as.vector(out), "SOR2")
})

test_that("the report names the vintage and says when releases candidates score within 10%", {
  local_mocked_bindings(
    detect_release = function(...) {
      rel_table(
        list("ensembl", "116", "38", TRUE, 0.0100),
        list("ensembl", "115", "38", TRUE, 0.0105)
      )
    }, # within 10%
    fetch_mapping = function(...) map_of(HS[1], "TP53")
  )
  expect_message(
    geneid2name(HS[1], species = "human"),
    "Best fingerprint match: homo_sapiens, ensembl release 116 \\(assembly 38\\).*2 candidates score within 10%"
  )
})

test_that("quiet suppresses the report but not the attributes", {
  local_mocked_bindings(
    detect_release = function(...) rel_table(list("ensembl", "116", "38", TRUE, 0.01)),
    fetch_mapping = function(...) map_of(HS[1], "TP53")
  )
  expect_silent(out <- geneid2name(HS[1], species = "human", quiet = TRUE))
  expect_equal(attr(out, "release"), "116")
  expect_equal(attr(out, "assembly"), "38")
})

test_that("unique= is passed through to the naming step", {
  local_mocked_bindings(
    fetch_mapping = function(...) map_of(HS[1:2], c("SAME", "SAME"))
  )
  out <- geneid2name(HS[1:2],
    species = "human", release = 116, unique = TRUE,
    quiet = TRUE
  )
  expect_equal(length(unique(as.vector(out))), 2L)
})

test_that("a supplied mapping bypasses detection for a plain vector too", {
  local_mocked_bindings(
    detect_species = function(...) stop("must not detect the species"),
    detect_release = function(...) stop("must not detect the release"),
    fetch_mapping = function(...) stop("must not fetch")
  )
  local_map <- data.frame(id = HS[1], name = "TP53")
  out <- geneid2name(HS[1:2], mapping = local_map, quiet = TRUE)
  expect_equal(as.vector(out), c("TP53", HS[2]))
})
