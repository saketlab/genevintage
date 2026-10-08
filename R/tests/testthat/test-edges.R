# Empty, blank and degenerate inputs. These are the cases that reach a user as a
# cryptic error rather than an answer, and none of them need the network.

test_that("an empty identifier vector is an empty answer, not an error", {
  expect_length(gene_names(character(), data.frame(id = "ENSG1", name = "TP53")), 0L)
  expect_equal(attr(
    gene_names(character(), data.frame(id = "ENSG1", name = "TP53")),
    "mapped"
  ), logical(0))
  expect_equal(nrow(parse_ens(character())), 0L)
})

test_that("nothing to go on yields no scheme and no species", {
  # mean() over nothing is NaN, which made which.max() empty and the comparison
  # after it a zero-length condition
  expect_equal(guess_scheme(character())$purity, 0)
  expect_equal(guess_scheme(character())$scheme, "symbol_or_mixed")
  expect_equal(guess_scheme(c("", ""))$purity, 0)
  expect_equal(guess_scheme(c(NA_character_, ""))$purity, 0)
})

test_that("an empty input proposes no species rather than all of them", {
  skip_no_index()
  # NaN > 0 is NA, which let 21 species through as if they explained nothing
  for (v in list(character(), "", c("", ""), NA_character_)) {
    d <- detect_species(v)
    expect_equal(nrow(d), 0L)
    expect_true(all(c("species", "common", "frac") %in% names(d)))
  }
})

test_that("gene_names keeps its promises on a mapping that maps nothing", {
  m <- data.frame(id = c("ENSG1", "ENSG2"), name = c("ENSG1", "ENSG2"))
  out <- gene_names(c("ENSG1", "ENSG2"), m, warn = FALSE)
  # both repeat their id, which is how "no name" is written; both are retained
  expect_equal(as.vector(out), c("ENSG1", "ENSG2"))
  expect_equal(attr(out, "mapped"), c(FALSE, FALSE))
})

test_that("releases are ordered numerically, not as strings", {
  skip_no_index()
  # "100" sorts before "50" as text; every species would look like a violation
  expect_equal(genevintage:::.rel_num(c("100", "50", "M37")), c(100, 50, 37))
  expect_lt(length(qmax_violations()), 10L)
})
