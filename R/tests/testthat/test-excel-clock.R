# Only midnight is a date export. Any other time marks a real timestamp.

midnights <- c(
  "0:00", "00:00", "0:00:00", "00:00:00", "12:00 AM", "12:00:00 AM",
  "12:00AM", "12:00 am", "12:00:00 am"
)

test_that("every spelling of midnight is stripped", {
  for (t in midnights) {
    for (d in c("9-Sep-2019", "9/9/2019", "2019-09-09", "2019/9/9")) {
      tok <- paste(d, t)
      expect_equal(
        as.vector(correct_genenames(tok, mapping = hs_now(), quiet = TRUE)),
        "SEPTIN9",
        label = tok
      )
    }
  }
})

test_that("the ISO T form of midnight is stripped", {
  for (tok in c("2019-09-09T00:00:00", "2019-09-09T00:00")) {
    expect_equal(
      as.vector(correct_genenames(tok, mapping = hs_now(), quiet = TRUE)),
      "SEPTIN9",
      label = tok
    )
  }
})

test_that("the cleaner removes midnight and nothing else", {
  clean <- genevintage:::.excel_clean
  expect_equal(clean("2019-09-09 00:00:00"), "2019-09-09")
  expect_equal(clean("2019-09-09T00:00:00"), "2019-09-09")
  expect_equal(clean("9/9/2019 12:00 AM"), "9/9/2019")
  expect_equal(clean("2018/09/01 10:32:31"), "2018/09/01 10:32:31")
  expect_equal(clean("9/9/2019 12:00 PM"), "9/9/2019 12:00 PM")
  expect_equal(clean("TP53"), "TP53")
  expect_equal(clean(NA_character_), NA_character_)
})

test_that("an instrument timestamp is not repaired and is not midnight", {
  tok <- "2018/09/01 10:32:31"
  expect_true(is.na(genevintage:::.excel_kind(tok)))
  out <- correct_genenames(c(tok, "TP53"), mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), c(tok, "TP53"))
  expect_false(any(attr(out, "corrected")))
  # a timestamp is not damage, so no status lists it
  expect_length(attr(out, "unresolved"), 0L)
  expect_length(attr(out, "ambiguous"), 0L)
  expect_length(attr(out, "skipped"), 0L)
})

test_that("times near midnight are not midnight", {
  toks <- c(
    "2019-09-09 00:00:01", "2019-09-09 00:01", "2019-09-09 00:30:00",
    "2019-09-09 01:00", "2019-09-09 10:32:31", "2019-09-09 23:59:59",
    "9/9/2019 12:30 AM", "9/9/2019 1:00 AM", "2019-09-09T10:32:31"
  )
  for (tok in toks) {
    out <- correct_genenames(tok, mapping = hs_now(), quiet = TRUE)
    expect_equal(as.vector(out), tok, label = tok)
    expect_false(attr(out, "corrected"), label = tok)
    expect_length(attr(out, "unresolved"), 0L)
  }
})

test_that("noon is not midnight", {
  for (tok in c("9/9/2019 12:00 PM", "9/9/2019 12:00:00 PM", "2019-09-09 12:00:00", "2019-09-09 12:00")) {
    out <- correct_genenames(tok, mapping = hs_now(), quiet = TRUE)
    expect_equal(as.vector(out), tok, label = tok)
  }
})

test_that("a real timestamp needs no annotation", {
  local_mocked_bindings(.locate_annotation = function(...) stop("must not fetch"))
  x <- c("2018/09/01 10:32:31", "2018/09/02 08:00:00", "TP53")
  expect_equal(as.vector(correct_genenames(x, quiet = TRUE)), x)
})

test_that("a repaired midnight token keeps the raw token in its status", {
  tok <- "1-Mar 00:00"
  out <- correct_genenames(tok, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), tok)
  expect_equal(attr(out, "ambiguous"), tok)
})
