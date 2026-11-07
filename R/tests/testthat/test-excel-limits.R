# Known limits. Each test states the wanted result and is skipped.

test_that("a timezone-shifted midnight is read as its date", {
  skip("known limit: timezone-shifted times are not read as midnight")
  # one hour west of the file's zone
  expect_equal(fixed("2019-09-08 23:00:00"), "SEPTIN9")
  expect_equal(fixed("2019-09-09 01:00:00"), "SEPTIN9")
})

test_that("a UTC suffix on midnight is read as its date", {
  skip("known limit: timezone-shifted times are not read as midnight")
  expect_equal(fixed("2019-09-09T00:00:00Z"), "SEPTIN9")
  expect_equal(fixed("2019-09-09 00:00:00+00:00"), "SEPTIN9")
})

test_that("unpadded locale tokens are read", {
  skip("known limit: a locale month first needs a two-digit day")
  expect_equal(fixed("Dic-1"), "DELEC1")
  expect_equal(fixed("AGO-1"), "AGO1")
  expect_equal(fixed("MARS-1", nom("MARS1")), "MARS1")
})

test_that("other-language months are read", {
  skip("known limit: other-language months outside the locale table")
  # German, French and Italian
  expect_equal(fixed("märz-03"), "MARCHF3")
  expect_equal(fixed("déc-01"), "DELEC1")
  expect_equal(fixed("set-09"), "SEPTIN9")
})

test_that("a Dutch May in the year style reaches MEI1 and MEI4", {
  # real rows: the IDs give Mei1 and Mei4
  x <- c("2001-05-01 00:00:00", "2004-05-01 00:00:00", "2003-09-01", "2005-09-01")
  expect_equal(fixed(x, mm_now(), species = "mouse")[1:2], c("Mei1", "Mei4"))
})

test_that("the month-year rule reads day-first forms", {
  skip("known limit: the month-year rule covers ISO forms only")
  expect_equal(fixed(c("1-Sep-2003", "1-Sep-2005")), c("SEPTIN3", "SEPTIN5"))
  expect_equal(fixed(c("01-Sep-03", "01-Sep-05")), c("SEPTIN3", "SEPTIN5"))
})

test_that("the month-year rule reads slash and dot forms", {
  skip("known limit: the month-year rule covers ISO forms only")
  expect_equal(fixed(c("01/09/2003", "01/09/2005")), c("SEPTIN3", "SEPTIN5"))
  expect_equal(fixed(c("9/1/2003", "9/1/2005")), c("SEPTIN3", "SEPTIN5"))
  expect_equal(fixed(c("01.09.2003", "01.09.2005")), c("SEPTIN3", "SEPTIN5"))
})
