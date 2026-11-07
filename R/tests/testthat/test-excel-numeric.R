# A list of mostly numbers holds identifiers or coordinates, not gene names.

test_that("Entrez-like lists stay unchanged", {
  x <- c("42248", "7157", "672", "1956")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_false(any(attr(out, "corrected")))
  # 7157 has four digits, so only 42248 looks like a date serial
  expect_equal(attr(out, "skipped"), "42248")
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("a skipped list reports its tokens in a message", {
  x <- c("42248", "43352", "7157", "672")
  expect_message(
    out <- correct_genenames(x, mapping = hs_now()),
    "2 left alone: mostly numbers"
  )
  expect_equal(attr(out, "skipped"), c("42248", "43352"))
  expect_message(correct_genenames(x, mapping = hs_now()), "e.g. 42248", fixed = TRUE)
  expect_silent(correct_genenames(x, mapping = hs_now(), quiet = TRUE))
})

test_that("a list with no damage has an empty skipped attribute", {
  out <- correct_genenames(c("7157", "672", "1956"), mapping = hs_now(), quiet = TRUE)
  expect_identical(attr(out, "skipped"), character())
  out <- correct_genenames(c("TP53", "MYC"), quiet = TRUE)
  expect_identical(attr(out, "skipped"), character())
})

test_that("a numeric list needs no annotation", {
  local_mocked_bindings(.locate_annotation = function(...) stop("must not fetch"))
  x <- c("42248", "7157", "672")
  out <- correct_genenames(x, quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_equal(attr(out, "skipped"), "42248")
})

test_that("a lone serial among symbols is still repaired", {
  # 42248 is 2015-09-01, and 43352 is 2018-09-09
  out <- run("42248")
  expect_equal(as.vector(out)[1], "SEPTIN1")
  expect_true(attr(out, "corrected")[1])
  expect_length(attr(out, "skipped"), 0L)
  expect_equal(fixed(c("TP53", "43352")), c("TP53", "SEPTIN9"))
})

test_that("one token alone is an explicit request", {
  expect_equal(as.vector(correct_genenames("42248", mapping = hs_now(), quiet = TRUE)), "SEPTIN1")
  # NA and empty tokens do not make it a list
  x <- c("42248", NA, "")
  expect_equal(as.vector(correct_genenames(x, mapping = hs_now(), quiet = TRUE)), c("SEPTIN1", NA, ""))
})

test_that("the rule needs more than half numeric tokens", {
  numeric_list <- genevintage:::.numeric_list
  expect_false(numeric_list("42248"))
  expect_false(numeric_list(c("42248", "TP53")))
  expect_false(numeric_list(c("42248", "7157", "TP53", "MYC")))
  expect_true(numeric_list(c("42248", "7157", "TP53")))
  expect_true(numeric_list(c("42248", "7157", "672", "TP53", "MYC")))
  expect_false(numeric_list(character()))
  # NA and empty tokens are not counted
  expect_false(numeric_list(c("42248", NA, "", "TP53")))
  expect_true(numeric_list(c("42248", "7157", NA, "", "TP53")))
})

test_that("exactly half numeric is still a gene list", {
  expect_equal(
    as.vector(correct_genenames(c("42248", "TP53"), mapping = hs_now(), quiet = TRUE)),
    c("SEPTIN1", "TP53")
  )
})

test_that("array coordinates stay unchanged", {
  x <- c("3/1/10", "3/1/11", "3/2/10", "1.2.3", "10-2-3", "TP53")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_false(any(attr(out, "corrected")))
  expect_true(all(c("3/1/10", "3/1/11", "3/2/10") %in% attr(out, "skipped")))
  expect_length(attr(out, "ambiguous"), 0L)
})

test_that("instrument timestamps stay unchanged", {
  x <- c("2018/09/01 10:32:31", "2018/09/01 10:33:02", "2018/09/02 08:00:00", "2018/09/02 08:01:15")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_false(any(attr(out, "corrected")))
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("a clear month token is repaired inside a numeric list", {
  x <- c("9-Sep", "42248", "7157", "672")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), c("SEPTIN9", "42248", "7157", "672"))
  expect_equal(attr(out, "corrected"), c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(attr(out, "skipped"), "42248")
})

test_that("scientific notation is not a numeric token", {
  x <- c("2.310009E+19", "42248", "7157")
  out <- correct_genenames(x, mapping = nom("2310009E13"), quiet = TRUE)
  expect_equal(as.vector(out), c("2310009E13", "42248", "7157"))
  expect_equal(attr(out, "skipped"), "42248")
})

test_that("placeholders are not numbers", {
  x <- c(".", "-", "--", "..", "TP53", "43352")
  expect_false(genevintage:::.numeric_list(x))
  expect_equal(
    as.vector(correct_genenames(x, mapping = hs_now(), quiet = TRUE)),
    c(".", "-", "--", "..", "TP53", "SEPTIN9")
  )
})

test_that("skipped keeps the raw token", {
  x <- c("'42248'", "7157", "672")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_equal(attr(out, "skipped"), "'42248'")
  # the time part counts as part of the raw token too
  y <- c("2019-09-09 00:00:00", "7157", "672")
  out <- correct_genenames(y, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), y)
  expect_equal(attr(out, "skipped"), y[1])
})

test_that("skipped keeps duplicates, order and alignment", {
  x <- c("43352", NA, "42248", "42248", "", "7157")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_length(attr(out, "corrected"), length(x))
  expect_equal(attr(out, "skipped"), c("43352", "42248", "42248"))
})

test_that("a matrix of Entrez row names keeps them", {
  ids <- c("42248", "7157", "672", "1956")
  out <- correct_genenames(mk(ids), mapping = hs_now(), quiet = TRUE)
  expect_equal(rownames(out), ids)
})

test_that("a plain number among symbols keeps a serial from becoming a gene", {
  x <- c("42248", "7157", "TP53", "MYC", "GAPDH")
  expect_false(genevintage:::.numeric_list(x))
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_equal(attr(out, "skipped"), "42248")
  expect_message(correct_genenames(x, mapping = hs_now()), "the list has plain numbers")
  # a month token is still damage
  y <- c("9-Sep", x)
  expect_equal(as.vector(correct_genenames(y, mapping = hs_now(), quiet = TRUE))[1], "SEPTIN9")
  # a digit-only symbol in the annotation is no evidence, and neither is a second serial
  m <- nom("SEPTIN1", "5", "TP53")
  expect_equal(fixed(c("42248", "5"), m), c("SEPTIN1", "5"))
  expect_equal(fixed(c("42248", "43352")), c("SEPTIN1", "SEPTIN9"))
})
