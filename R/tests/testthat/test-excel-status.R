# The status attributes, and what must stay untouched.

status <- c("corrected", "ambiguous", "unresolved", "skipped")

test_that("every status attribute is always present", {
  for (x in list(character(), "TP53", "9-Sep", c("42248", "7157"), c(NA, ""))) {
    out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
    expect_true(all(status %in% names(attributes(out))))
    expect_type(attr(out, "corrected"), "logical")
    expect_type(attr(out, "ambiguous"), "character")
    expect_type(attr(out, "unresolved"), "character")
    expect_type(attr(out, "skipped"), "character")
    expect_length(attr(out, "corrected"), length(x))
    expect_length(out, length(x))
  }
})

test_that("empty input gives empty attributes", {
  out <- correct_genenames(character(), quiet = TRUE)
  expect_identical(attr(out, "corrected"), logical())
  expect_identical(attr(out, "ambiguous"), character())
  expect_identical(attr(out, "unresolved"), character())
  expect_identical(attr(out, "skipped"), character())
})

test_that("a mixed list keeps order, length and alignment", {
  x <- c(
    NA, "", "9-Sep", "TP53", "1-Mar", "30-Sep", "9-Sep", "SEPTIN2", "MYC",
    "3-Mar", "ACTB", "EGFR", "BRCA1", "KRAS", NA
  )
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  want <- x
  want[c(3, 7)] <- "SEPTIN9"
  want[10] <- "MARCHF3"
  expect_equal(as.vector(out), want)
  expect_equal(attr(out, "corrected"), seq_along(x) %in% c(3, 7, 10))
  expect_equal(attr(out, "ambiguous"), "1-Mar")
  expect_equal(attr(out, "unresolved"), "30-Sep")
  expect_length(attr(out, "skipped"), 0L)
})

test_that("each damaged token is in exactly one status", {
  x <- c("9-Sep", "1-Mar", "30-Sep", "3-Mar", "15-Sep", "42248", plain)
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  damaged <- sum(!is.na(genevintage:::.excel_kind(x)))
  total <- sum(attr(out, "corrected")) + length(attr(out, "ambiguous")) +
    length(attr(out, "unresolved")) + length(attr(out, "skipped"))
  expect_equal(total, damaged)
})

test_that("a skipped token is not also in another status", {
  x <- c("42248", "1-Mar", "7157", "672")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(attr(out, "skipped"), "42248")
  expect_equal(attr(out, "ambiguous"), "1-Mar")
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("NA and empty strings stay untouched", {
  x <- c(NA, "", "   ", NA_character_, "9-Sep")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), c(NA, "", "   ", NA, "SEPTIN9"))
  expect_equal(attr(out, "corrected"), c(FALSE, FALSE, FALSE, FALSE, TRUE))
})

test_that("only NA and empty input needs no annotation", {
  local_mocked_bindings(.locate_annotation = function(...) stop("must not fetch"))
  out <- correct_genenames(c(NA, "", NA), quiet = TRUE)
  expect_equal(as.vector(out), c(NA, "", NA))
  expect_length(attr(out, "corrected"), 3L)
})

test_that("each duplicate is repaired and reported", {
  expect_equal(fixed(c("9-Sep", "9-Sep")), c("SEPTIN9", "SEPTIN9"))
  x <- c("1-Mar", "1-Mar", "30-Sep", "30-Sep")
  out <- run(x)
  expect_equal(as.vector(out)[1:4], x)
  expect_equal(attr(out, "ambiguous"), x[1:2])
  expect_equal(attr(out, "unresolved"), x[3:4])
})

test_that("case does not change the reading", {
  for (tok in c("9-sep", "9-SEP", "9-sEp", "SEP-09", "sep-09", "SEPT-09", "9-SEPTEMBER")) {
    expect_equal(fixed(tok), "SEPTIN9", label = tok)
  }
  for (tok in c("3-MAR", "3-mar", "MAR-03", "MARCH-03")) {
    expect_equal(fixed(tok), "MARCHF3", label = tok)
  }
})

test_that("case does not change a matching mouse symbol", {
  out <- run("9-SEP", mm_now(), species = "mouse")
  expect_equal(as.vector(out)[1], "Septin9")
})

test_that("real symbols that resemble dates stay untouched", {
  x <- c("SEPT9", "MARCH1", "Sept9", "MAR1", "Oct4", "DEC1", "NKX2-1", "HLA-A", "MT-ND1", "H2-K1")
  expect_true(all(is.na(genevintage:::.excel_kind(x))))
  local_mocked_bindings(.locate_annotation = function(...) stop("must not fetch"))
  out <- correct_genenames(x, quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_false(any(attr(out, "corrected")))
})

test_that("a worm symbol shaped like a date stays untouched", {
  x <- c("sep-1", "ago-1", "mars-1", "oct-1", "mai-2", "dec-1")
  out <- correct_genenames(x, mapping = nom(x), species = "worm", quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_false(any(attr(out, "corrected")))
  expect_length(attr(out, "unresolved"), 0L)
  expect_length(attr(out, "ambiguous"), 0L)
})

test_that("a symbol token matches its annotation in any case", {
  out <- run("SEP-1", nom("sep-1"))
  expect_equal(as.vector(out)[1], "SEP-1")
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("a valid gene in the list removes that candidate", {
  expect_equal(fixed(c("1-Mar", "MARCHF1")), c("MTARC1", "MARCHF1"))
  expect_equal(fixed(c("1-Mar", "MTARC1")), c("MARCHF1", "MTARC1"))
  expect_equal(fixed(c("2-Mar", "MARCHF2")), c("MTARC2", "MARCHF2"))
  expect_equal(fixed(c("MARCHF2", "Mar-02")), c("MARCHF2", "MTARC2"))
})

test_that("a valid gene counts in any case and under its old spelling", {
  expect_equal(fixed(c("1-Mar", "marchf1")), c("MTARC1", "marchf1"))
  expect_equal(fixed(c("1-Mar", "MARCH1")), c("MTARC1", "MARCH1"))
  expect_equal(fixed(c("1-Mar", "MARC1")), c("MARCHF1", "MARC1"))
})

test_that("a valid gene outside the candidates changes nothing", {
  out <- run(c("1-Mar", "SEPTIN2", "MARCHF2"))
  expect_equal(as.vector(out)[1], "1-Mar")
  expect_equal(attr(out, "ambiguous"), "1-Mar")
})

test_that("both candidates present leaves the token ambiguous", {
  x <- c("1-Mar", "MARCHF1", "MTARC1")
  out <- run(x)
  expect_equal(as.vector(out)[1:3], x)
  expect_equal(attr(out, "ambiguous"), "1-Mar")
})

test_that("mouse Septin2 in the list points 2-Sep to the old Sep2", {
  # Septin2 is present, so Apoa1 remains
  m <- mm_now()
  expect_equal(fixed(c("2-Sep", "Septin2"), m, species = "mouse"), c("Apoa1", "Septin2"))
})

test_that("reversing the input reverses the output", {
  x <- c("TP53", "9-Sep", NA, "3-Mar", "", "1-Dec", "MYC", "GAPDH", "ACTB")
  a <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  b <- correct_genenames(rev(x), mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(b), rev(as.vector(a)))
  expect_equal(attr(b, "corrected"), rev(attr(a, "corrected")))
})

test_that("a factor gets the same result as its labels", {
  x <- c("9-Sep", "TP53", "1-Mar")
  a <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  b <- correct_genenames(factor(x), mapping = hs_now(), quiet = TRUE)
  expect_equal(b, a)
})

test_that("unique keeps the attributes aligned", {
  x <- c("9-Sep", "SEPTIN9", "1-Mar", "1-Mar", "TP53")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE, unique = TRUE)
  expect_length(out, length(x))
  expect_equal(anyDuplicated(as.vector(out)), 0L)
  expect_equal(as.vector(out)[1], "SEPTIN9_9-Sep")
  expect_equal(attr(out, "corrected"), c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(attr(out, "ambiguous"), c("1-Mar", "1-Mar"))
})

test_that("messages count only the tokens they describe", {
  x <- c("9-Sep", "30-Sep", "42248", "TP53", "MYC", "ACTB")
  msg <- capture_messages(correct_genenames(x, mapping = hs_now()))
  expect_match(msg, "repaired 2 of 3 damaged names", all = FALSE)
  expect_match(msg, "1 left alone: no symbol in this annotation produces 30-Sep", all = FALSE)
})

test_that("quiet silences every report", {
  x <- c("9-Sep", "30-Sep", "1-Mar", "TP53")
  expect_silent(correct_genenames(x, mapping = hs_now(), quiet = TRUE))
  expect_silent(correct_genenames(c("42248", "7157"), mapping = hs_now(), quiet = TRUE))
})
