# The category and locale attributes describe each damaged token.

cat_of <- function(x, m = hs_now(), ...) attr(run(x, m, ...), "category")[seq_along(x)]
loc_of <- function(x, m = hs_now(), ...) attr(run(x, m, ...), "locale")[seq_along(x)]

test_that("attributes have one value per input token", {
  x <- c("2-Sep", "TP53", "Mar-03")
  out <- run(x)
  expect_length(attr(out, "category"), length(out))
  expect_length(attr(out, "locale"), length(out))
  expect_type(attr(out, "category"), "character")
  expect_type(attr(out, "locale"), "character")
})

test_that("an undamaged token has no category and no locale", {
  expect_equal(cat_of(c("TP53", "SEPTIN2")), c(NA_character_, NA_character_))
  expect_equal(loc_of(c("TP53", "SEPTIN2")), c(NA_character_, NA_character_))
})

test_that("each form of date has its category", {
  expect_equal(cat_of("2-Sep"), "day_month")
  expect_equal(cat_of("2/Sep"), "day_month")
  expect_equal(cat_of("Sep-02"), "month_day")
  expect_equal(cat_of("Sep 2"), "month_day")
  expect_equal(cat_of("2-Sep-2019"), "day_month_year")
  expect_equal(cat_of("9/2/2019"), "numeric_date")
  expect_equal(cat_of("2019/9/2"), "date_day")
  expect_equal(cat_of("2020-09-02 00:00:00"), "date_day")
})

test_that("serials and scientific notation have their category", {
  m <- nom("SEPTIN1", "2310009E13", "2310009E13Rik")
  expect_equal(cat_of("42248", m), "serial")
  expect_equal(cat_of("2.31E+13", m), "scientific")
})

test_that("the year reading has its own category", {
  m <- nom("SEPTIN3", "SEPTIN5", "SEPTIN1", "SEPTIN2")
  expect_equal(cat_of(c("2003-09-01", "2005-09-01"), m), c("date_year", "date_year"))
  expect_equal(cat_of("2003-09-01", m), "date_day")
})

test_that("a merged list labels each token by its own reading", {
  m <- nom("SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN5")
  x <- c("2003-09-01", "2005-09-01", "2020-09-01", "2020-09-02")
  expect_equal(cat_of(x, m), c("date_year", "date_year", "date_day", "date_day"))
})

test_that("a damaged token left alone keeps its category", {
  out <- run(c("1-Mar", "2001-08-01"), nom("MARCHF1", "MTARC1"))
  expect_equal(attr(out, "category")[1:2], c("day_month", "date_day"))
  expect_true("1-Mar" %in% attr(out, "ambiguous"))
})

test_that("a skipped token keeps its category", {
  x <- c("42248", "7157", "672", "2-Sep")
  out <- correct_genenames(x, mapping = hs_now(), quiet = TRUE)
  expect_equal(attr(out, "category")[1], "serial")
})

test_that("a real symbol that looks like a date has no category", {
  worm <- nom("ago-1", "sep-1")
  out <- correct_genenames(c("ago-1", "sep-1"), mapping = worm, quiet = TRUE)
  expect_equal(attr(out, "category"), c(NA_character_, NA_character_))
  expect_equal(attr(out, "locale"), c(NA_character_, NA_character_))
})

test_that("a local month word in the token gives its languages", {
  expect_equal(loc_of("mars-02"), "fr,sv,nb")
  expect_equal(loc_of("mei-01", nom("MEI1")), "nl")
  expect_equal(loc_of("ago-01"), "es,it,pt")
  expect_equal(loc_of("1-mrt"), "nl")
  expect_equal(loc_of("dic-01"), "es,it")
})

test_that("an English month gives no locale", {
  expect_equal(loc_of("2-Sep"), NA_character_)
  expect_equal(loc_of("Mar-03"), NA_character_)
  expect_equal(loc_of("2020-09-02"), NA_character_)
})

test_that("the year reading gets a locale only from the repaired gene", {
  x <- c("2001-05-01", "2004-05-01", "2003-09-01", "2005-09-01")
  expect_equal(loc_of(x, mm_now(), species = "mouse"), c("nl", "nl", NA, NA))
  m <- nom("MARCHF2", "MARS2", "SEPTIN3", "SEPTIN5")
  expect_equal(loc_of(c("2002-03-01", "2003-09-01", "2005-09-01", "MARCHF2"), m)[1], "fr,sv,nb")
  # an ambiguous token has no repaired gene, so it proves no language
  expect_equal(loc_of(c("2002-03-01", "2003-09-01", "2005-09-01"), m)[1], NA_character_)
})

test_that("the attributes do not change the repair", {
  x <- c("2-Sep", "Mar-03", "2003-09-01", "2005-09-01", "mars-02")
  out <- run(x)
  expect_equal(as.vector(out)[1:2], c("SEPTIN2", "MARCHF3"))
  expect_length(attr(out, "corrected"), length(out))
})

test_that("matrices and data frames keep none of the attributes", {
  m <- matrix(1:6, 3, dimnames = list(c("2-Sep", "TP53", "MYC"), NULL))
  out <- correct_genenames(m, mapping = nom("SEPTIN2", "TP53", "MYC"), quiet = TRUE)
  expect_null(attr(out, "category"))
  expect_null(attr(out, "locale"))
})

test_that("unique keeps the attributes aligned", {
  out <- run(c("2-Sep", "Sep-02"), nom("SEPTIN2", "TP53"), unique = TRUE)
  expect_equal(attr(out, "category")[1:2], c("day_month", "month_day"))
})
