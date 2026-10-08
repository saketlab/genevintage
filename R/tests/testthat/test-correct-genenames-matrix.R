new_nom <- data.frame(
  id = paste0("g", 1:3),
  name = c("SEPTIN9", "MARCHF1", "TP53"), chr = "1",
  biotype = "protein_coding", stringsAsFactors = FALSE
)

test_that("row names are repaired, dimensions and values untouched", {
  m <- mk(c("9-Sep", "1-Mar", "TP53"))
  out <- correct_genenames(m, mapping = new_nom, quiet = TRUE)
  expect_equal(dim(out), dim(m))
  expect_equal(rownames(out), c("SEPTIN9", "MARCHF1", "TP53"))
  expect_equal(unname(out[, "s1"]), unname(m[, "s1"]))
})

test_that("an unresolved token keeps its damaged value", {
  out <- correct_genenames(mk(c("9-Sep", "15-Sep")), mapping = new_nom, quiet = TRUE)
  expect_equal(rownames(out), c("SEPTIN9", "15-Sep"))
})

test_that("colliding repairs are disambiguated with the causing token", {
  m <- mk(c("9-Sep", "9-September"))
  out <- correct_genenames(m, mapping = new_nom, quiet = TRUE)
  expect_equal(anyDuplicated(rownames(out)), 0L)
  expect_true(all(grepl("^SEPTIN9", rownames(out))))
})

test_that("unique = FALSE leaves collisions in place", {
  out <- correct_genenames(mk(c("9-Sep", "9-September")),
    mapping = new_nom, unique = FALSE, quiet = TRUE
  )
  expect_equal(rownames(out), c("SEPTIN9", "SEPTIN9"))
})

test_that("two identical unresolved tokens are disambiguated, not self-suffixed", {
  # unrepaired, so no self-suffix like 15-Sep_15-Sep
  out <- correct_genenames(mk(c("15-Sep", "15-Sep")), mapping = new_nom, quiet = TRUE)
  expect_equal(rownames(out), c("15-Sep", "15-Sep_1"))
})

test_that("nothing is smuggled onto the object", {
  out <- correct_genenames(mk(c("9-Sep", "TP53")), mapping = new_nom, quiet = TRUE)
  expect_setequal(names(attributes(out)), c("dim", "dimnames"))
})

test_that("a mismatched or missing id vector is refused, not guessed", {
  m <- mk(c("9-Sep", "TP53"))
  expect_error(correct_genenames(m, row_ids = "9-Sep", mapping = new_nom), "2 rows")
  dimnames(m) <- NULL
  expect_error(correct_genenames(m, mapping = new_nom), "no row names")
})

test_that("a damaged row without mapping or species is refused", {
  expect_error(correct_genenames(mk(c("9-Sep", "TP53"))), "species")
})

test_that("a sparse matrix comes back sparse, valid, and unmoved", {
  skip_if_not_installed("Matrix")
  sm <- Matrix::Matrix(c(1, 0, 0, 0, 2, 0),
    nrow = 3, sparse = TRUE,
    dimnames = list(c("9-Sep", "1-Mar", "TP53"), c("s1", "s2"))
  )
  out <- correct_genenames(sm, mapping = new_nom, quiet = TRUE)
  expect_s4_class(out, "dgCMatrix")
  expect_true(methods::validObject(out))
  expect_equal(rownames(out), c("SEPTIN9", "MARCHF1", "TP53"))
  expect_identical(out@x, sm@x)
})

test_that("data frames work too", {
  df <- data.frame(s1 = 1:2, s2 = 3:4, row.names = c("9-Sep", "TP53"))
  expect_equal(
    rownames(correct_genenames(df, mapping = new_nom, quiet = TRUE)),
    c("SEPTIN9", "TP53")
  )
})
