# An explicit mapping keeps these offline; detection is exercised elsewhere.
local_map <- data.frame(
  id   = c("ENSG00000141510", "ENSG00000012048", "ENSG00000284733", "ENSG00000999998"),
  name = c("TP53", "BRCA1", "ENSG00000284733", "TP53") # 3rd unnamed, 4th collides
)
sparse_ids <- c("ENSG00000141510", "ENSG00000284733", "ENSG00000012048")
mk_sparse <- function() {
  Matrix::Matrix(c(1, 0, 0, 0, 2, 0),
    nrow = 3, sparse = TRUE,
    dimnames = list(sparse_ids, c("s1", "s2"))
  )
}

test_that("row names become gene names, dimensions and order untouched", {
  m <- mk(c("ENSG00000141510", "ENSG00000284733", "ENSG00000012048"))
  out <- geneid2name(m, mapping = local_map, quiet = TRUE)
  expect_equal(dim(out), dim(m))
  expect_equal(rownames(out), c("TP53", "ENSG00000284733", "BRCA1"))
  expect_equal(unname(out[, "s1"]), unname(m[, "s1"])) # values never move
})

test_that("colliding names are disambiguated by default", {
  out <- geneid2name(mk(c("ENSG00000141510", "ENSG00000999998")),
    mapping = local_map, quiet = TRUE
  )
  expect_equal(anyDuplicated(rownames(out)), 0L)
  expect_true(all(grepl("^TP53", rownames(out))))
})

test_that("an unmapped id keeps its identifier", {
  out <- geneid2name(mk(c("ENSG00000284733", "ENSG00000000000")),
    mapping = local_map, quiet = TRUE
  )
  expect_equal(rownames(out), c("ENSG00000284733", "ENSG00000000000"))
})

test_that("nothing is smuggled onto the object", {
  out <- geneid2name(mk(c("ENSG00000141510", "ENSG00000284733")),
    mapping = local_map, quiet = TRUE
  )
  expect_setequal(names(attributes(out)), c("dim", "dimnames"))
})

test_that("a mismatched or missing id vector is refused, not guessed", {
  m <- mk(c("ENSG00000141510", "ENSG00000012048"))
  expect_error(geneid2name(m, row_ids = "ENSG00000141510", mapping = local_map), "2 rows")
  dimnames(m) <- NULL
  expect_error(geneid2name(m, mapping = local_map), "no row names")
})

test_that("a sparse matrix comes back sparse, valid, and unmoved", {
  skip_if_not_installed("Matrix")
  sm <- mk_sparse()
  out <- geneid2name(sm, mapping = local_map, quiet = TRUE)

  expect_s4_class(out, "dgCMatrix") # not densified
  expect_true(methods::validObject(out)) # still a legal S4 object
  expect_equal(rownames(out), c("TP53", "ENSG00000284733", "BRCA1"))
  expect_equal(as.numeric(out[, "s2"]), c(0, 2, 0))
  expect_identical(out@x, sm@x) # storage untouched
  expect_identical(out@i, sm@i)
})

test_that("a renamed sparse matrix still survives normal Matrix operations", {
  skip_if_not_installed("Matrix")
  sm <- mk_sparse()
  out <- geneid2name(sm, mapping = local_map, quiet = TRUE)
  # these are where a stray attribute on an S4 object shows up as breakage
  expect_true(methods::validObject(out[1:2, , drop = FALSE]))
  expect_true(methods::validObject(rbind(out, out)))
  expect_true(methods::validObject(Matrix::t(out)))
  expect_equal(dim(out %*% matrix(1, ncol(out), 1)), c(3L, 1L))
})

test_that("data frames work too", {
  df <- data.frame(
    s1 = 1:2, s2 = 3:4,
    row.names = c("ENSG00000141510", "ENSG00000012048")
  )
  expect_equal(
    rownames(geneid2name(df, mapping = local_map, quiet = TRUE)),
    c("TP53", "BRCA1")
  )
})
