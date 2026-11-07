# Which annotation a call uses, and which inputs stop it.

# A fake fetch_mapping() that records its arguments.
recorder <- function(map = hs_now()) {
  seen <- new.env()
  list(seen = seen, fetch = function(species, release, assembly = NULL, source = NULL) {
    seen$species <- species
    seen$release <- release
    seen$source <- source
    map
  })
}

test_that("undatable symbols with a species use the newest release", {
  skip_no_index()
  rec <- recorder()
  local_mocked_bindings(fetch_mapping = rec$fetch)
  for (tok in c("5", "7157", "LOC100506145", "NM_000546", "NR_046018.2", "XM_011524000", "XR_1234.1")) {
    x <- c(tok, "TP53", "9-Sep")
    expect_no_error(out <- correct_genenames(x, species = "human", quiet = TRUE))
    expect_equal(as.vector(out), c(tok, "TP53", "SEPTIN9"), label = tok)
    expect_equal(rec$seen$release, .newest_release("homo_sapiens"), label = tok)
    expect_equal(rec$seen$source, "ensembl", label = tok)
  }
})

test_that("a mix of undatable symbols does not stop the repair", {
  skip_no_index()
  local_mocked_bindings(fetch_mapping = function(...) hs_now())
  x <- c("5", "LOC100506145", "NM_000546", "TP53", "9-Sep", "7157", "MYC", "ACTB", "EGFR", "KRAS")
  out <- correct_genenames(x, species = "human", quiet = TRUE)
  expect_equal(as.vector(out)[5], "SEPTIN9")
  expect_equal(as.vector(out)[-5], x[-5])
})

test_that("the fallback is the same for mouse", {
  skip_no_index()
  rec <- recorder(mm_now())
  local_mocked_bindings(fetch_mapping = rec$fetch)
  out <- correct_genenames(c("LOC100506145", "9-Sep"), species = "mouse", quiet = TRUE)
  expect_equal(as.vector(out)[2], "Septin9")
  expect_equal(rec$seen$species, "mus_musculus")
  expect_equal(rec$seen$release, .newest_release("mus_musculus"))
})

test_that("the release is reported unless quiet", {
  skip_no_index()
  local_mocked_bindings(fetch_mapping = function(...) hs_now())
  msg <- capture_messages(correct_genenames(c("5", "9-Sep"), species = "human"))
  want <- paste0("homo_sapiens, ensembl release ", .newest_release("homo_sapiens"))
  expect_match(msg, want, all = FALSE)
})

test_that("an explicit release is kept", {
  skip_no_index()
  rec <- recorder()
  local_mocked_bindings(fetch_mapping = rec$fetch)
  correct_genenames(c("5", "9-Sep"), species = "human", release = 95, quiet = TRUE)
  expect_equal(rec$seen$release, 95)
})

test_that("gene sets accept undatable symbols with a species", {
  skip_no_index()
  m <- data.frame(
    id = c("g1", "g2", "g3", "g4"), name = c("MT-ND1", "5", "LOC100506145", "XIST"),
    chr = c("MT", "1", "2", "X"), biotype = c("protein_coding", "protein_coding", "lncRNA", "lncRNA")
  )
  local_mocked_bindings(fetch_mapping = function(...) m)
  ids <- c("5", "LOC100506145", "MT-ND1", "XIST")
  expect_equal(mito_genes(ids, species = "human", quiet = TRUE)$input, "MT-ND1")
  expect_equal(sex_genes(ids, species = "human", quiet = TRUE)$input, "XIST")
  expect_equal(mito_genes(c("5", "LOC100506145"), species = "human", quiet = TRUE)$name, character())
  expect_setequal(coding_genes(ids, species = "human", quiet = TRUE)$input, c("5", "MT-ND1"))
})

test_that("undatable symbols without a species are refused", {
  skip_no_index()
  expect_error(correct_genenames(c("5", "LOC1", "9-Sep"), quiet = TRUE), "pass `species`")
  expect_error(mito_genes(c("5", "LOC1"), quiet = TRUE), "pass `species`")
})

test_that("stable IDs still date the annotation", {
  skip_no_index()
  rec <- recorder()
  local_mocked_bindings(fetch_mapping = rec$fetch)
  ids <- c("ENSG00000141510", "ENSG00000012048")
  best <- suppressWarnings(detect_release(ids, species = "homo_sapiens"))[1, ]
  correct_genenames(c(ids, "9-Sep"), species = "human", quiet = TRUE)
  expect_equal(rec$seen$release, best$release)
  expect_equal(rec$seen$source, best$source)
})

test_that("stable IDs give the species when none is named", {
  skip_no_index()
  rec <- recorder()
  local_mocked_bindings(fetch_mapping = rec$fetch)
  out <- correct_genenames(c("ENSG00000141510", "9-Sep"), quiet = TRUE)
  expect_equal(rec$seen$species, "homo_sapiens")
  expect_equal(as.vector(out)[2], "SEPTIN9")
})

test_that("digits, LOC and RefSeq cannot date a release", {
  skip_no_index()
  expect_false(any(.datable(c("5", "7157", "LOC100506145", "NM_000546", "XR_1234.1", "TP53", "9-Sep"))))
  expect_true(all(.datable(c(
    "ENSG00000141510", "ENSG00000141510.16", "ENSMUSG00000051951", "FBgn0000008", "WBGene00000001"
  ))))
})

test_that("the newest release is the highest number", {
  skip_no_index()
  fp <- .fp()
  r <- fp$release[fp$species == "homo_sapiens" & fp$source == "ensembl"]
  expect_equal(as.numeric(.newest_release("homo_sapiens")), max(as.numeric(r)))
  expect_error(.newest_release("not_a_species"), "no ensembl releases")
})

test_that("identifiers of no use to a species are refused", {
  skip_no_index()
  expect_error(detect_release("TP53", species = "human"), "none of the identifiers match homo_sapiens")
  expect_error(detect_release("ENSMUSG00000051951", species = "human"), "none of the identifiers match")
  expect_error(detect_release("TP53"), "no known species")
})

test_that("a stray symbol does not change the score of a species' own IDs", {
  skip_no_index()
  ids <- c("ENSG00000141510", "ENSG00000012048")
  a <- suppressWarnings(detect_release(ids, species = "human"))
  b <- suppressWarnings(detect_release(c(ids, "TP53", "5"), species = "human"))
  expect_equal(a, b)
})

test_that("IDs of another species are refused in a repair", {
  skip_no_index()
  local_mocked_bindings(fetch_mapping = function(...) stop("must not fetch"))
  expect_error(
    correct_genenames(c("9-Sep", "ENSMUSG00000051951"), species = "human", quiet = TRUE),
    "none of the identifiers match"
  )
})

test_that("a mapping avoids the location step", {
  local_mocked_bindings(
    detect_release = function(...) stop("must not detect"),
    fetch_mapping = function(...) stop("must not fetch")
  )
  x <- c("ENSMUSG00000051951", "5", "9-Sep")
  expect_equal(as.vector(correct_genenames(x, mapping = hs_now(), quiet = TRUE))[3], "SEPTIN9")
})

test_that("a mapping is checked before it is used", {
  expect_error(correct_genenames("9-Sep", mapping = "x"), "must be a data frame")
  expect_error(correct_genenames("9-Sep", mapping = data.frame(id = "g1")), "needs column")
})
