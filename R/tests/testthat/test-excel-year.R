# The gene number can sit in the year: SEPT3 stored as 2003-09-01.

test_that("two years in one month put the gene number in the year", {
  out <- run(c("2003-09-01", "2005-09-01"))
  expect_equal(as.vector(out)[1:2], c("SEPTIN3", "SEPTIN5"))
  expect_equal(attr(out, "corrected")[1:2], c(TRUE, TRUE))
  expect_length(attr(out, "ambiguous"), 0L)
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("a single date gives no sign of the year style", {
  expect_equal(fixed("2003-09-01"), "SEPTIN1")
  # the same date twice is still one year
  expect_equal(fixed(c("2003-09-01", "2003-09-01")), c("SEPTIN1", "SEPTIN1"))
  # a date of another shape is not a second year
  expect_equal(fixed(c("2003-09-01", "9-Sep")), c("SEPTIN1", "SEPTIN9"))
})

test_that("a one-year list keeps the day reading", {
  days <- sprintf("2020-09-%02d", 1:14)
  want <- c(paste0("SEPTIN", 1:12), "SEPTIN7P2", "SEPTIN14")
  expect_equal(fixed(days), want)
  expect_equal(fixed(days)[1:2], c("SEPTIN1", "SEPTIN2"))
  # two months in one year are not two years
  expect_equal(
    fixed(c("2020-09-01", "2020-03-03", "2020-03-05")),
    c("SEPTIN1", "MARCHF3", "MARCHF5")
  )
})

test_that("the March family reads the same way", {
  x <- sprintf("20%02d-03-01", 3:11)
  expect_equal(fixed(x), paste0("MARCHF", 3:11))
  # the signal can come from another month
  expect_equal(
    fixed(c("2003-03-01", "2003-09-01", "2005-09-01")),
    c("MARCHF3", "SEPTIN3", "SEPTIN5")
  )
})

test_that("MARCH1 and MARCH2 stay ambiguous in the year style", {
  x <- c("2001-03-01", "2002-03-01", "2003-03-01")
  out <- run(x)
  expect_equal(as.vector(out)[1:3], c("2001-03-01", "2002-03-01", "MARCHF3"))
  expect_equal(attr(out, "ambiguous"), x[1:2])
  expect_equal(attr(out, "corrected")[1:3], c(FALSE, FALSE, TRUE))
})

test_that("an undamaged partner settles the year-style March tokens", {
  expect_equal(
    fixed(c("2001-03-01", "2005-03-01", "MTARC1")),
    c("MARCHF1", "MARCHF5", "MTARC1")
  )
  expect_equal(
    fixed(c("2001-03-01", "2005-03-01", "MARCHF1")),
    c("MTARC1", "MARCHF5", "MARCHF1")
  )
  # the partner can carry its old spelling
  # MARS2 is also a candidate here, so this needs a list without it
  expect_equal(
    fixed(c("2002-03-01", "2005-03-01", "MARC2"), nom("MARCHF2", "MARCHF5", "MTARC2")),
    c("MARCHF2", "MARCHF5", "MARC2")
  )
})

test_that("the high years reach the renamed genes", {
  x <- c("2012-09-01", "2013-09-01", "2014-09-01", "2015-09-01")
  expect_equal(fixed(x), c("SEPTIN12", "SEPTIN7P2", "SEPTIN14", "SELENOF"))
  # MGI has SEP15 too
  mouse <- fixed(c("2014-09-01", "2015-09-01"), mm_now(), species = "mouse")
  expect_equal(mouse, c("Septin14", "Selenof"))
})

test_that("day 1 of every month is rewritten in a year-style list", {
  x <- c(sprintf("2003-%02d-01", 1:12), "2005-09-01")
  expect_equal(
    genevintage:::.month_year(x),
    c(paste0("3-", month.abb), "5-Sep")
  )
  # DEC1 is now DELEC1, and FEB3 is an old name of SCN1A
  expect_equal(
    fixed(c("2001-12-01", "2003-02-01", "2003-09-01", "2005-09-01")),
    c("DELEC1", "SCN1A", "SEPTIN3", "SEPTIN5")
  )
})

test_that("a rewritten month with no gene is left alone and reported", {
  x <- c("2001-08-01", "2003-09-01", "2005-09-01")
  out <- run(x, nom("SEPTIN3", "SEPTIN5", "TP53"))
  expect_equal(as.vector(out)[1:3], c("2001-08-01", "SEPTIN3", "SEPTIN5"))
  expect_equal(attr(out, "unresolved"), "2001-08-01")
  # ago is August in Spanish and Italian, so AGO1 is a candidate
  expect_equal(fixed(x)[1], "AGO1")
})

test_that("a year with no gene number is not read as day 1", {
  # no gene has number 0 or 45, and 2020 would be SEPT20
  x <- c("2000-09-01", "2045-09-01", "2020-09-01", "2003-09-01")
  out <- run(x)
  expect_equal(as.vector(out)[1:4], c(x[1:3], "SEPTIN3"))
  expect_equal(attr(out, "unresolved"), x[1:3])
})

test_that("the time part does not hide the year style", {
  expect_equal(
    fixed(c("2003-09-01 00:00:00", "2005-09-01 00:00:00")),
    c("SEPTIN3", "SEPTIN5")
  )
  # some rows have the time part and some do not
  expect_equal(
    fixed(c("2003-09-01 00:00:00", "2005-09-01", "2007-09-01 0:00")),
    c("SEPTIN3", "SEPTIN5", "SEPTIN7")
  )
  expect_equal(fixed("2003-09-01 00:00:00"), "SEPTIN1")
})

test_that("a time other than midnight is not year-style evidence", {
  x <- c("2003-09-01 10:00:00", "2005-09-01 10:00:00")
  out <- run(x)
  expect_equal(as.vector(out)[1:2], x)
  expect_false(any(attr(out, "corrected")))
})

test_that("other separators and an unpadded month are read", {
  forms <- list(
    c("2003/09/01", "2005/09/01"),
    c("2003.09.01", "2005.09.01"),
    c("2003-9-01", "2005-9-01"),
    c("2003/9/1", "2005/9/1"),
    c("2003.9.1", "2005.9.1"),
    c("2003-09-1", "2005-09-1"),
    # the rows can differ
    c("2003/09/01", "2005-9-1")
  )
  for (x in forms) {
    expect_equal(fixed(x), c("SEPTIN3", "SEPTIN5"), label = x[1])
  }
})

test_that("a two-digit month is not taken for a padded one", {
  expect_equal(
    genevintage:::.month_year(c("2003-10-01", "2005-10-01", "2003-12-01")),
    c("3-Oct", "5-Oct", "3-Dec")
  )
})

test_that("years outside 20xx are not rewritten", {
  x <- c("1903-09-01", "1905-09-01")
  expect_identical(genevintage:::.month_year(x), x)
  # the day reading still applies, and the day is 1
  expect_equal(fixed(x), c("SEPTIN1", "SEPTIN1"))
  # this year is not a date at all
  y <- c("2103-09-01", "2105-09-01")
  expect_identical(genevintage:::.month_year(y), y)
  out <- run(y)
  expect_equal(as.vector(out)[1:2], y)
  expect_false(any(attr(out, "corrected")))
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("rows that are not day 1 are not rewritten", {
  x <- c("2003-09-01", "2005-09-01", "2007-09-02", "2020-09-14")
  expect_equal(
    genevintage:::.month_year(x),
    c("3-Sep", "5-Sep", "2007-09-02", "2020-09-14")
  )
  expect_equal(fixed(x), c("SEPTIN3", "SEPTIN5", "SEPTIN2", "SEPTIN14"))
  # two years on day 2 are not a signal
  y <- c("2003-09-02", "2005-09-02")
  expect_identical(genevintage:::.month_year(y), y)
  expect_equal(fixed(y), c("SEPTIN2", "SEPTIN2"))
})

test_that("two years in different months are not a signal", {
  x <- c("2003-09-01", "2005-03-03")
  expect_identical(genevintage:::.month_year(x), x)
  y <- c("2003-09-01", "2005-03-01")
  expect_identical(genevintage:::.month_year(y), y)
  expect_equal(fixed(y)[1], "SEPTIN1")
})

test_that("the year reading ignores row order", {
  x <- c("2003-09-01", "TP53", "2005-09-01", "2007-03-01")
  want <- c("SEPTIN3", "TP53", "SEPTIN5", "MARCHF7")
  expect_equal(fixed(x), want)
  expect_equal(fixed(rev(x)), rev(want))
})

test_that("the rewrite keeps other tokens and empty input", {
  expect_identical(genevintage:::.month_year(character()), character())
  x <- c(NA, "", "TP53", "9-Sep", "2003-09-01", "2005-09-01")
  expect_identical(
    genevintage:::.month_year(x),
    c(NA, "", "TP53", "9-Sep", "3-Sep", "5-Sep")
  )
})

test_that("a quoted date counts as a year-style date", {
  expect_equal(fixed(c("'2003-09-01'", "\"2005-09-01\"")), c("SEPTIN3", "SEPTIN5"))
})

test_that("the original token is what the attributes report", {
  x <- c("2001-03-01 00:00:00", "2005-03-01 00:00:00")
  out <- run(x)
  expect_equal(attr(out, "ambiguous"), x[1])
  expect_equal(as.vector(out)[1], x[1])
})

test_that("a short list of dates is a numeric list and stays unchanged", {
  # two of the three tokens are numeric, so the numeric rule comes first
  x <- c("2003-09-01", "2005-09-01", "TP53")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), x)
  expect_equal(attr(out, "skipped"), x[1:2])
})

test_that("row names of a matrix get the year reading", {
  ids <- c("2003-09-01", "2005-09-01", "2001-03-01", plain)
  out <- correct_genenames(mk(ids), mapping = hs_now(), quiet = TRUE)
  expect_equal(rownames(out)[1:3], c("SEPTIN3", "SEPTIN5", "2001-03-01"))
  expect_equal(dim(out), c(length(ids), 2L))
  expect_equal(unname(out[, "s1"]), seq_along(ids))
})

test_that("row names of a data frame get the year reading", {
  ids <- c("2003-09-01", "2005-09-01", plain)
  df <- data.frame(s1 = seq_along(ids), row.names = ids)
  out <- correct_genenames(df, mapping = hs_now(), quiet = TRUE)
  expect_equal(rownames(out)[1:2], c("SEPTIN3", "SEPTIN5"))
  expect_equal(out$s1, df$s1)
})

test_that("year-style repairs that collide get unique row names", {
  ids <- c("2003-09-01", "2005-09-01", "3-Sep", plain)
  out <- correct_genenames(mk(ids), mapping = hs_now(), quiet = TRUE)
  expect_equal(anyDuplicated(rownames(out)), 0L)
  expect_true(all(startsWith(rownames(out)[c(1, 3)], "SEPTIN3")))
  expect_equal(rownames(out)[2], "SEPTIN5")
})

test_that("a month with other days keeps the day reading in a merged list", {
  m <- nom("SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN12", "SEPTIN14")
  x <- c("2012-09-01", "2012-09-02", "2012-09-03", "2014-09-01", "2014-09-02", "2014-09-03")
  # Two files from two years give two copies of each gene.
  out <- run(x, m)
  expect_equal(attr(out, "corrected")[1:6], rep(TRUE, 6))
  expect_equal(as.vector(out)[1:6], rep(c("SEPTIN1", "SEPTIN2", "SEPTIN3"), 2))
})

test_that("day-style months from years that have no such gene stay on the day", {
  m <- nom("SEPTIN1", "SEPTIN2")
  x <- c("2016-09-01", "2016-09-02", "2017-09-01", "2017-09-02")
  expect_equal(fixed(x, m)[1:4], rep(c("SEPTIN1", "SEPTIN2"), 2))
})

test_that("a year and month with other days keep the day reading", {
  m <- nom("MARCHF3", "MARCHF5", "SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN5")
  x <- c("2003-03-01", "2005-03-01", "2003-09-01", "2005-09-01", "2005-09-02")
  # 2005-09 has day 2 (day-style); 2003-09 has only day 1 (year-style)
  expect_equal(fixed(x, m)[1:5], c("MARCHF3", "MARCHF5", "SEPTIN3", "SEPTIN1", "SEPTIN2"))
})

test_that("a merge of a year-style file and a day-style file reads each correctly", {
  m <- nom("SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN5", "SEPTIN14")
  x <- c("2003-09-01", "2005-09-01", "2020-09-01", "2020-09-02", "2020-09-14")
  expect_equal(fixed(x, m), c("SEPTIN3", "SEPTIN5", "SEPTIN1", "SEPTIN2", "SEPTIN14"))
})

test_that("a year-style list with only day 1 still gets the year reading", {
  m <- nom("SEPTIN3", "SEPTIN5", "SEPTIN6", "MARCHF7")
  x <- c("2003-09-01", "2005-09-01", "2006-09-01", "2007-03-01")
  expect_equal(fixed(x, m)[1:4], c("SEPTIN3", "SEPTIN5", "SEPTIN6", "MARCHF7"))
})

test_that("the year reading offers other-language month names", {
  m <- nom("MARCHF2", "MARCHF3", "MARS2", "SEPTIN3", "SEPTIN5")
  # MARS2 is a real gene behind 2002-03-01 in GSE202940 and GSE213252.
  x <- c("2002-03-01", "2003-09-01", "2005-09-01")
  out <- run(x, m)
  expect_true("2002-03-01" %in% attr(out, "ambiguous"))
  expect_equal(fixed(c(x, "MARCHF2"), m)[1], "MARS2")
  expect_equal(fixed(c(x, "MARS2"), m)[1], "MARCHF2")
})

test_that("the year reading finds a gene that only the local month name gives", {
  m <- nom("MEI1", "MEI4", "SEPTIN3", "SEPTIN5")
  x <- c("2001-05-01", "2004-05-01", "2003-09-01", "2005-09-01")
  expect_equal(fixed(x, m)[1:4], c("MEI1", "MEI4", "SEPTIN3", "SEPTIN5"))
  m <- nom("AGO2", "SEPTIN3", "SEPTIN5")
  expect_equal(fixed(c("2002-08-01", "2003-09-01", "2005-09-01"), m)[1], "AGO2")
})

test_that("the day reading does not add other-language month names", {
  m <- nom("MARCHF2", "MARS2")
  # one year puts the number in the day, and no token says mars
  expect_equal(fixed(c("2020-03-02"), m)[1], "MARCHF2")
})

test_that("a year reading that names no gene falls back to the day", {
  x <- c("2016-09-01", "2017-09-01")
  out <- run(x)
  expect_equal(as.vector(out)[1:2], c("SEPTIN1", "SEPTIN1"))
  expect_length(attr(out, "unresolved"), 0L)
  # one gene from the year reading is enough to keep it
  y <- c(x, "2003-09-01")
  expect_equal(fixed(y), c(x, "SEPTIN3"))
})
