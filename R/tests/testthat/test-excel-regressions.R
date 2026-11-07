# False positives, locale months and the four-digit year rule.

test_that("labels that only resemble dates are not damage", {
  x <- c("ago.1", "Mars.3", "mars.1", "ago.2", "Sep 4, 2023", "September 4, 2023", "7183.9.15")
  expect_true(all(is.na(genevintage:::.excel_kind(x))))
  # the annotation has the gene that a false reading would give
  out <- run(x, nom("AGO1", "AGO2", "MARS1", "MARS3", "MARCHF3", "SEPTIN4", "SELENOF"))
  expect_equal(as.vector(out)[seq_along(x)], x)
  expect_false(any(attr(out, "corrected")))
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("a four-digit number must be a 19xx or 20xx year", {
  expect_true(is.na(genevintage:::.excel_kind("7183.9.15")))
  expect_true(is.na(genevintage:::.excel_kind("3609.1.15")))
  expect_true(is.na(genevintage:::.excel_kind("2119.9.15")))
})

test_that("a dotted date with a real year is read as a date", {
  # NEWS: a four-digit year must start with 19 or 20, and the rule allows dots
  expect_equal(genevintage:::.excel_kind("2008.09.15"), "date_day")
  expect_equal(fixed("2008.09.15"), "SELENOF")
})

test_that("placeholders are not damage", {
  x <- c(".", "..", "...", "-", "--", "---", "/", "NA", "n/a", "N/A", "?", " ", "TP53")
  expect_true(all(is.na(genevintage:::.excel_kind(x))))
  local_mocked_bindings(.locate_annotation = function(...) stop("must not fetch"))
  out <- correct_genenames(x, quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_false(any(attr(out, "corrected")))
})

test_that("a letter that only folds to a month does not stop the call", {
  long_s <- "ſep-1"
  expect_no_error(out <- correct_genenames(long_s, mapping = hs_now(), quiet = TRUE))
  expect_equal(as.vector(out), long_s)
  expect_false(attr(out, "corrected"))
})

test_that("a folded letter leaves the rest of the list intact", {
  x <- c("ſep-1", "ſept-09", "oKt-09", "9-Sep", "TP53")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), c(x[1:3], "SEPTIN9", "TP53"))
  expect_equal(attr(out, "corrected"), c(FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("mars-02 has the March and the MARS candidates", {
  cand <- genevintage:::.excel_candidates("mars-02", "homo_sapiens")[[1]]
  expect_true(all(c("MARCH2", "MARCHF2", "MTARC2") %in% cand))
  expect_true("MARS2" %in% toupper(cand))
})

test_that("mars-02 is settled only by what is present", {
  m <- nom("MARCHF2", "MARS2")
  out <- run("mars-02", m)
  expect_equal(as.vector(out)[1], "mars-02")
  expect_equal(attr(out, "ambiguous"), "mars-02")
  expect_equal(fixed(c("mars-02", "MARS2"), m), c("MARCHF2", "MARS2"))
  expect_equal(fixed(c("mars-02", "MARCHF2"), m), c("MARS2", "MARCHF2"))
  # only one of them exists in this annotation
  expect_equal(fixed("mars-02", nom("MARS2")), "MARS2")
  expect_equal(fixed("mars-02", nom("MARCHF2")), "MARCHF2")
})

test_that("ago-01 can be AGO1 or a date in August", {
  expect_equal(fixed("ago-01", nom("AGO1")), "AGO1")
  for (tok in c("AGO-01", "Ago-01", "ago 01", "ago.01", "ago/01")) {
    expect_equal(fixed(tok, nom("AGO1")), "AGO1", label = tok)
  }
  out <- run("ago-01", nom("AGO1", "AUG1"))
  expect_equal(attr(out, "ambiguous"), "ago-01")
})

test_that("locale months with a two-digit day map to the right month", {
  cases <- c(
    "Dic-01" = "DELEC1", "dez-01" = "DELEC1", "mrt-03" = "MARCHF3",
    "maart-03" = "MARCHF3", "agosto-09" = "AUG9", "okt-04" = "OCT4",
    "oktober-04" = "OCT4", "mei-07" = "MAY7", "mai-07" = "MAY7"
  )
  m <- nom("DELEC1", "MARCHF3", "AUG9", "OCT4", "MAY7")
  for (tok in names(cases)) {
    expect_equal(fixed(tok, m), cases[[tok]], label = tok)
  }
})

test_that("a locale month that is a symbol prefix gives that symbol", {
  expect_equal(fixed("mei-01", nom("MEI1")), "MEI1")
  expect_equal(fixed("MARS-02", nom("MARS2")), "MARS2")
})

test_that("a locale month with no gene is reported", {
  out <- run("mei-01", nom("TP53"))
  expect_equal(as.vector(out)[1], "mei-01")
  expect_equal(attr(out, "unresolved"), "mei-01")
})

test_that("a locale month first with one digit is a label", {
  x <- c("Dic-1", "AGO-1", "MARS-1", "ago-1", "mei.1", "okt 1")
  expect_true(all(is.na(genevintage:::.excel_kind(x))))
  out <- run(x, nom("DELEC1", "AGO1", "MARS1", "MARCHF1", "MEI1"))
  expect_equal(as.vector(out)[seq_along(x)], x)
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("an English month first with one digit is a date", {
  cases <- c(
    "Sep-1" = "SEPTIN1", "sep.9" = "SEPTIN9", "SEP 9" = "SEPTIN9", "Sept-9" = "SEPTIN9",
    "Mar-3" = "MARCHF3", "MAR-3" = "MARCHF3", "March-3" = "MARCHF3", "Dec-1" = "DELEC1"
  )
  for (tok in names(cases)) {
    expect_equal(fixed(tok), cases[[tok]], label = tok)
  }
})

test_that("a day-first locale month needs no padding", {
  cases <- c("1-Dic" = "DELEC1", "1-dez" = "DELEC1", "1-ago" = "AGO1", "3-mars" = "MARCHF3", "03-maart" = "MARCHF3")
  for (tok in names(cases)) {
    expect_equal(fixed(tok), cases[[tok]], label = tok)
  }
  expect_equal(attr(run("2-mars"), "ambiguous"), "2-mars")
})

test_that("an ISO date with slashes and no padding is read", {
  expect_equal(genevintage:::.excel_kind("2019/3/1"), "date_day")
  expect_equal(fixed("2019/3/1", nom("MARCHF1")), "MARCHF1")
  expect_equal(fixed("2019/3/1", nom("MTARC1")), "MTARC1")
  expect_equal(fixed(c("2019/3/3", "2019/9/9", "2019/09/09", "2019.9.9"))[1:4], c("MARCHF3", rep("SEPTIN9", 3)))
})

test_that("2019/3/1 has both March genes, so it is ambiguous", {
  out <- run("2019/3/1")
  expect_equal(as.vector(out)[1], "2019/3/1")
  expect_equal(attr(out, "ambiguous"), "2019/3/1")
})

test_that("the year rule accepts 19xx and 20xx only", {
  expect_equal(fixed("1999/9/9"), "SEPTIN9")
  expect_equal(fixed("2099/9/9"), "SEPTIN9")
  for (tok in c("1899/9/9", "2119/9/9", "7183/9/15", "0219/9/9")) {
    expect_true(is.na(genevintage:::.excel_kind(tok)), label = tok)
    expect_equal(fixed(tok), tok, label = tok)
  }
})

test_that("a month or day that cannot exist gets no candidate", {
  x <- c("2019/13/1", "2019/9/32", "2019/0/9", "2019/9/0")
  out <- run(x)
  expect_equal(as.vector(out)[seq_along(x)], x)
  expect_false(any(attr(out, "corrected")))
  expect_equal(attr(out, "unresolved"), x)
})

test_that("each local month word has its languages", {
  expect_setequal(names(genevintage:::.LOCALE_MONTHS), names(genevintage:::.LOCALE_LANGUAGES))
})
