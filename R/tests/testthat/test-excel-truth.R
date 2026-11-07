# Repairs checked against GEO rows whose Ensembl ID gives the true gene.
# Each list element is one kind of file, because a list decides the year reading.
# The mouse truth is upper case, as the ID lookup gave it.

truth_human <- list(
  day_first = c(
    "1-Dec" = "DELEC1", "1-Sep" = "SEPTIN1", "2-Sep" = "SEPTIN2", "3-Mar" = "MARCHF3", "3-Sep" = "SEPTIN3",
    "4-Mar" = "MARCHF4", "4-Sep" = "SEPTIN4", "5-Mar" = "MARCHF5", "5-Sep" = "SEPTIN5", "6-Mar" = "MARCHF6",
    "6-Sep" = "SEPTIN6", "7-Mar" = "MARCHF7", "7-Sep" = "SEPTIN7", "8-Mar" = "MARCHF8", "8-Sep" = "SEPTIN8",
    "9-Mar" = "MARCHF9", "9-Sep" = "SEPTIN9", "10-Mar" = "MARCHF10", "10-Sep" = "SEPTIN10", "11-Mar" = "MARCHF11",
    "11-Sep" = "SEPTIN11", "12-Sep" = "SEPTIN12", "14-Sep" = "SEPTIN14", "15-Sep" = "SELENOF"
  ),
  month_first = c(
    "Dec-01" = "DELEC1", "Mar-03" = "MARCHF3", "Mar-04" = "MARCHF4", "Mar-05" = "MARCHF5", "Mar-06" = "MARCHF6",
    "Mar-07" = "MARCHF7", "Mar-08" = "MARCHF8", "Mar-09" = "MARCHF9", "Mar-10" = "MARCHF10", "Mar-11" = "MARCHF11",
    "Sep-01" = "SEPTIN1", "Sep-02" = "SEPTIN2", "Sep-03" = "SEPTIN3", "Sep-04" = "SEPTIN4", "Sep-05" = "SEPTIN5",
    "Sep-06" = "SEPTIN6", "Sep-07" = "SEPTIN7", "Sep-08" = "SEPTIN8", "Sep-09" = "SEPTIN9", "Sep-10" = "SEPTIN10",
    "Sep-11" = "SEPTIN11", "Sep-12" = "SEPTIN12", "Sep-14" = "SEPTIN14", "Sep-15" = "SELENOF"
  ),
  y2017 = c(
    "2017-03-03 00:00:00" = "MARCHF3", "2017-03-05 00:00:00" = "MARCHF5", "2017-03-06 00:00:00" = "MARCHF6",
    "2017-03-07 00:00:00" = "MARCHF7", "2017-03-08 00:00:00" = "MARCHF8", "2017-03-09 00:00:00" = "MARCHF9",
    "2017-09-01 00:00:00" = "SEPTIN1", "2017-09-02 00:00:00" = "SEPTIN2", "2017-09-03 00:00:00" = "SEPTIN3",
    "2017-09-04 00:00:00" = "SEPTIN4", "2017-09-05 00:00:00" = "SEPTIN5", "2017-09-06 00:00:00" = "SEPTIN6",
    "2017-09-07 00:00:00" = "SEPTIN7", "2017-09-08 00:00:00" = "SEPTIN8", "2017-09-09 00:00:00" = "SEPTIN9",
    "2017-09-10 00:00:00" = "SEPTIN10", "2017-09-11 00:00:00" = "SEPTIN11", "2017-09-15 00:00:00" = "SELENOF"
  ),
  y2023 = c(
    "2023-03-03 00:00:00" = "MARCHF3", "2023-03-04 00:00:00" = "MARCHF4", "2023-03-05 00:00:00" = "MARCHF5",
    "2023-03-06 00:00:00" = "MARCHF6", "2023-03-07 00:00:00" = "MARCHF7", "2023-03-08 00:00:00" = "MARCHF8",
    "2023-03-09 00:00:00" = "MARCHF9", "2023-03-10 00:00:00" = "MARCHF10", "2023-03-11 00:00:00" = "MARCHF11",
    "2023-09-01 00:00:00" = "SEPTIN1", "2023-09-02 00:00:00" = "SEPTIN2", "2023-09-03 00:00:00" = "SEPTIN3",
    "2023-09-04 00:00:00" = "SEPTIN4", "2023-09-05 00:00:00" = "SEPTIN5", "2023-09-06 00:00:00" = "SEPTIN6",
    "2023-09-07 00:00:00" = "SEPTIN7", "2023-09-08 00:00:00" = "SEPTIN8", "2023-09-09 00:00:00" = "SEPTIN9",
    "2023-09-10 00:00:00" = "SEPTIN10", "2023-09-11 00:00:00" = "SEPTIN11", "2023-09-12 00:00:00" = "SEPTIN12",
    "2023-09-14 00:00:00" = "SEPTIN14", "2023-12-01 00:00:00" = "DELEC1"
  ),
  y2025 = c(
    "2025-03-03 00:00:00" = "MARCHF3", "2025-03-04 00:00:00" = "MARCHF4", "2025-03-05 00:00:00" = "MARCHF5",
    "2025-03-06 00:00:00" = "MARCHF6", "2025-03-07 00:00:00" = "MARCHF7", "2025-03-08 00:00:00" = "MARCHF8",
    "2025-03-09 00:00:00" = "MARCHF9", "2025-03-10 00:00:00" = "MARCHF10", "2025-03-11 00:00:00" = "MARCHF11",
    "2025-09-01 00:00:00" = "SEPTIN1", "2025-09-02 00:00:00" = "SEPTIN2", "2025-09-03 00:00:00" = "SEPTIN3",
    "2025-09-04 00:00:00" = "SEPTIN4", "2025-09-05 00:00:00" = "SEPTIN5", "2025-09-06 00:00:00" = "SEPTIN6",
    "2025-09-07 00:00:00" = "SEPTIN7", "2025-09-08 00:00:00" = "SEPTIN8", "2025-09-09 00:00:00" = "SEPTIN9",
    "2025-09-10 00:00:00" = "SEPTIN10", "2025-09-11 00:00:00" = "SEPTIN11", "2025-09-12 00:00:00" = "SEPTIN12",
    "2025-09-14 00:00:00" = "SEPTIN14", "2025-12-01 00:00:00" = "DELEC1"
  )
)

truth_mouse <- list(
  day_first = c(
    "1-Sep" = "SEPTIN1", "2-Sep" = "SEPTIN2", "3-Mar" = "MARCHF3", "3-Sep" = "SEPTIN3", "4-Mar" = "MARCHF4",
    "4-Sep" = "SEPTIN4", "5-Mar" = "MARCHF5", "5-Sep" = "SEPTIN5", "6-Mar" = "MARCHF6", "6-Sep" = "SEPTIN6",
    "7-Mar" = "MARCHF7", "7-Sep" = "SEPTIN7", "8-Mar" = "MARCHF8", "8-Sep" = "SEPTIN8", "9-Mar" = "MARCHF9",
    "9-Sep" = "SEPTIN9", "10-Mar" = "MARCHF10", "10-Sep" = "SEPTIN10", "11-Mar" = "MARCHF11", "11-Sep" = "SEPTIN11",
    "12-Sep" = "SEPTIN12", "14-Sep" = "SEPTIN14", "15-Sep" = "SELENOF"
  ),
  month_first = c(
    "Mar-03" = "MARCHF3", "Mar-04" = "MARCHF4", "Mar-05" = "MARCHF5", "Mar-06" = "MARCHF6", "Mar-07" = "MARCHF7",
    "Mar-08" = "MARCHF8", "Mar-09" = "MARCHF9", "Mar-10" = "MARCHF10", "Mar-11" = "MARCHF11", "Sep-01" = "SEPTIN1",
    "Sep-02" = "SEPTIN2", "Sep-03" = "SEPTIN3", "Sep-04" = "SEPTIN4", "Sep-05" = "SEPTIN5", "Sep-06" = "SEPTIN6",
    "Sep-07" = "SEPTIN7", "Sep-08" = "SEPTIN8", "Sep-09" = "SEPTIN9", "Sep-10" = "SEPTIN10", "Sep-11" = "SEPTIN11",
    "Sep-12" = "SEPTIN12", "Sep-14" = "SEPTIN14", "Sep-15" = "SELENOF"
  ),
  y2018 = c(
    "2018-03-03 00:00:00" = "MARCHF3", "2018-03-04 00:00:00" = "MARCHF4", "2018-03-05 00:00:00" = "MARCHF5",
    "2018-03-06 00:00:00" = "MARCHF6", "2018-03-07 00:00:00" = "MARCHF7", "2018-03-08 00:00:00" = "MARCHF8",
    "2018-03-09 00:00:00" = "MARCHF9", "2018-03-10 00:00:00" = "MARCHF10", "2018-03-11 00:00:00" = "MARCHF11",
    "2018-09-01 00:00:00" = "SEPTIN1", "2018-09-02 00:00:00" = "SEPTIN2", "2018-09-03 00:00:00" = "SEPTIN3",
    "2018-09-04 00:00:00" = "SEPTIN4", "2018-09-05 00:00:00" = "SEPTIN5", "2018-09-06 00:00:00" = "SEPTIN6",
    "2018-09-07 00:00:00" = "SEPTIN7", "2018-09-08 00:00:00" = "SEPTIN8", "2018-09-09 00:00:00" = "SEPTIN9",
    "2018-09-10 00:00:00" = "SEPTIN10", "2018-09-11 00:00:00" = "SEPTIN11", "2018-09-12 00:00:00" = "SEPTIN12",
    "2018-09-14 00:00:00" = "SEPTIN14"
  ),
  y2019 = c(
    "2019-03-03 00:00:00" = "MARCHF3", "2019-03-04 00:00:00" = "MARCHF4", "2019-03-05 00:00:00" = "MARCHF5",
    "2019-03-06 00:00:00" = "MARCHF6", "2019-03-07 00:00:00" = "MARCHF7", "2019-03-08 00:00:00" = "MARCHF8",
    "2019-03-09 00:00:00" = "MARCHF9", "2019-03-10 00:00:00" = "MARCHF10", "2019-03-11 00:00:00" = "MARCHF11",
    "2019-09-01 00:00:00" = "SEPTIN1", "2019-09-02 00:00:00" = "SEPTIN2", "2019-09-03 00:00:00" = "SEPTIN3",
    "2019-09-04 00:00:00" = "SEPTIN4", "2019-09-05 00:00:00" = "SEPTIN5", "2019-09-06 00:00:00" = "SEPTIN6",
    "2019-09-07 00:00:00" = "SEPTIN7", "2019-09-08 00:00:00" = "SEPTIN8", "2019-09-09 00:00:00" = "SEPTIN9",
    "2019-09-10 00:00:00" = "SEPTIN10", "2019-09-11 00:00:00" = "SEPTIN11", "2019-09-12 00:00:00" = "SEPTIN12",
    "2019-09-14 00:00:00" = "SEPTIN14", "2019-09-15 00:00:00" = "SELENOF"
  ),
  y2020 = c(
    "2020-03-03 00:00:00" = "MARCHF3", "2020-03-04 00:00:00" = "MARCHF4", "2020-03-05 00:00:00" = "MARCHF5",
    "2020-03-06 00:00:00" = "MARCHF6", "2020-03-07 00:00:00" = "MARCHF7", "2020-03-08 00:00:00" = "MARCHF8",
    "2020-03-09 00:00:00" = "MARCHF9", "2020-03-10 00:00:00" = "MARCHF10", "2020-03-11 00:00:00" = "MARCHF11",
    "2020-09-01 00:00:00" = "SEPTIN1", "2020-09-02 00:00:00" = "SEPTIN2", "2020-09-03 00:00:00" = "SEPTIN3",
    "2020-09-04 00:00:00" = "SEPTIN4", "2020-09-05 00:00:00" = "SEPTIN5", "2020-09-06 00:00:00" = "SEPTIN6",
    "2020-09-07 00:00:00" = "SEPTIN7", "2020-09-08 00:00:00" = "SEPTIN8", "2020-09-09 00:00:00" = "SEPTIN9",
    "2020-09-10 00:00:00" = "SEPTIN10", "2020-09-11 00:00:00" = "SEPTIN11", "2020-09-12 00:00:00" = "SEPTIN12",
    "2020-09-14 00:00:00" = "SEPTIN14", "2020-09-15 00:00:00" = "SELENOF"
  ),
  y2022 = c(
    "2022-03-03 00:00:00" = "MARCHF3", "2022-03-04 00:00:00" = "MARCHF4", "2022-03-05 00:00:00" = "MARCHF5",
    "2022-03-06 00:00:00" = "MARCHF6", "2022-03-07 00:00:00" = "MARCHF7", "2022-03-08 00:00:00" = "MARCHF8",
    "2022-03-09 00:00:00" = "MARCHF9", "2022-03-10 00:00:00" = "MARCHF10", "2022-03-11 00:00:00" = "MARCHF11",
    "2022-09-01 00:00:00" = "SEPTIN1", "2022-09-02 00:00:00" = "SEPTIN2", "2022-09-03 00:00:00" = "SEPTIN3",
    "2022-09-04 00:00:00" = "SEPTIN4", "2022-09-05 00:00:00" = "SEPTIN5", "2022-09-06 00:00:00" = "SEPTIN6",
    "2022-09-07 00:00:00" = "SEPTIN7", "2022-09-08 00:00:00" = "SEPTIN8", "2022-09-09 00:00:00" = "SEPTIN9",
    "2022-09-10 00:00:00" = "SEPTIN10", "2022-09-11 00:00:00" = "SEPTIN11", "2022-09-12 00:00:00" = "SEPTIN12",
    "2022-09-15 00:00:00" = "SELENOF"
  ),
  y2024 = c(
    "2024-03-03 00:00:00" = "MARCHF3", "2024-03-04 00:00:00" = "MARCHF4", "2024-03-05 00:00:00" = "MARCHF5",
    "2024-03-06 00:00:00" = "MARCHF6", "2024-03-07 00:00:00" = "MARCHF7", "2024-03-08 00:00:00" = "MARCHF8",
    "2024-03-09 00:00:00" = "MARCHF9", "2024-03-10 00:00:00" = "MARCHF10", "2024-03-11 00:00:00" = "MARCHF11",
    "2024-09-01 00:00:00" = "SEPTIN1", "2024-09-02 00:00:00" = "SEPTIN2", "2024-09-03 00:00:00" = "SEPTIN3",
    "2024-09-04 00:00:00" = "SEPTIN4", "2024-09-05 00:00:00" = "SEPTIN5", "2024-09-06 00:00:00" = "SEPTIN6",
    "2024-09-07 00:00:00" = "SEPTIN7", "2024-09-08 00:00:00" = "SEPTIN8", "2024-09-09 00:00:00" = "SEPTIN9",
    "2024-09-10 00:00:00" = "SEPTIN10", "2024-09-11 00:00:00" = "SEPTIN11", "2024-09-12 00:00:00" = "SEPTIN12",
    "2024-09-14 00:00:00" = "SEPTIN14", "2024-09-15 00:00:00" = "SELENOF"
  ),
  year_style = c(
    "2001-09-01 00:00:00" = "SEPTIN1", "2002-09-01 00:00:00" = "SEPTIN2", "2003-03-01 00:00:00" = "MARCHF3",
    "2003-09-01 00:00:00" = "SEPTIN3", "2004-03-01 00:00:00" = "MARCHF4", "2004-09-01 00:00:00" = "SEPTIN4",
    "2005-03-01 00:00:00" = "MARCHF5", "2005-09-01 00:00:00" = "SEPTIN5", "2006-03-01 00:00:00" = "MARCHF6",
    "2006-09-01 00:00:00" = "SEPTIN6", "2007-03-01 00:00:00" = "MARCHF7", "2007-09-01 00:00:00" = "SEPTIN7",
    "2008-03-01 00:00:00" = "MARCHF8", "2008-09-01 00:00:00" = "SEPTIN8", "2009-03-01 00:00:00" = "MARCHF9",
    "2009-09-01 00:00:00" = "SEPTIN9", "2010-03-01 00:00:00" = "MARCHF10", "2010-09-01 00:00:00" = "SEPTIN10",
    "2011-03-01 00:00:00" = "MARCHF11", "2011-09-01 00:00:00" = "SEPTIN11", "2012-09-01 00:00:00" = "SEPTIN12",
    "2014-09-01 00:00:00" = "SEPTIN14", "2015-09-01 00:00:00" = "SELENOF"
  )
)

not_recognised <- c(
  "ATRIP.1", "ATRIP.2", "CASR.1", "CASR.2", "CATX-2", "ERVE-1",
  "ERVFRD-1", "ERVFRD-3", "ERVH-1", "ERVK-28", "ERVV-1", "ERVV-2",
  "ERVW-1", "GHRLOS.1", "GHRLOS.2", "HULC.1", "HULC.2", "HYMAI.1",
  "HYMAI.2", "IGHVIII-44", "IGHVIII-82", "IGLVI-20", "IGLVI-38", "IGLVI-42",
  "IGLVI-56", "IGLVI-63", "IGLVI-68", "IGLVI-70", "IGLVIV-53", "IGLVIV-59",
  "IGLVIV-64", "IGLVIV-65", "IGLVV-58", "IGLVV-66", "ITCH.1", "ITCH.2",
  "RMRP.1", "RMRP.2", "Tdh.1", "Vault.1", "Vault.2", "Vault.3",
  "Vault.4", "pRNA.1", "pRNA.2", "pRNA.3", "pRNA.4", "pRNA.5",
  "pRNA.6", "pRNA.7", "pRNA.8", "pRNA.9"
)

# Rows with two possible genes, with the genes the IDs gave.
ambiguous_human <- list(
  "1-Mar" = c("MARCHF1", "MTARC1"),
  "2-Mar" = c("MARCHF2", "MTARC2"),
  "Mar-01" = c("MARCHF1", "MTARC1"),
  "Mar-02" = c("MARCHF2", "MTARC2"),
  "2023-03-01 00:00:00" = c("MARCHF1", "MTARC1"),
  "2023-03-02 00:00:00" = c("MARCHF2", "MTARC2"),
  "2025-03-01 00:00:00" = c("MARCHF1", "MTARC1"),
  "2025-03-02 00:00:00" = c("MARCHF2", "MTARC2")
)
ambiguous_mouse <- list(
  "1-Mar" = c("MARCHF1", "MTARC1"),
  "2-Mar" = c("MARCHF2", "MTARC2"),
  "Mar-01" = c("MARCHF1", "MTARC1"),
  "Mar-02" = c("MARCHF2", "MTARC2"),
  "2020-03-01 00:00:00" = c("MARCHF1", "MTARC1"),
  "2020-03-02 00:00:00" = c("MARCHF2", "MTARC2")
)

# real mouse lists have Apoa1 undamaged
mouse_file <- function(tokens, mapping = mm_now()) {
  as.vector(correct_genenames(c(tokens, "Apoa1", plain),
    mapping = mapping, species = "mouse", quiet = TRUE
  ))[seq_along(tokens)]
}

test_that("human rows repair to the gene their Ensembl ID names", {
  for (file in names(truth_human)) {
    want <- truth_human[[file]]
    out <- run(names(want))
    expect_equal(as.vector(out)[seq_along(want)], unname(want), label = file)
    expect_true(all(attr(out, "corrected")[seq_along(want)]), label = file)
    expect_length(attr(out, "ambiguous"), 0L)
    expect_length(attr(out, "unresolved"), 0L)
  }
})

test_that("mouse rows repair to the gene their Ensembl ID names", {
  for (file in names(truth_mouse)) {
    want <- truth_mouse[[file]]
    expect_equal(toupper(mouse_file(names(want))), unname(want), label = file)
  }
})

test_that("mouse repairs use the mouse spelling", {
  expect_equal(
    mouse_file(c("1-Sep", "3-Mar", "15-Sep")),
    c("Septin1", "Marchf3", "Selenof")
  )
})

test_that("each token gives its gene outside its own list too", {
  # year-style rows need their list, so they stay out of this loop
  for (want in truth_human[c("day_first", "month_first")]) {
    for (tok in names(want)) {
      expect_equal(
        as.vector(correct_genenames(tok, mapping = hs_now(), quiet = TRUE)),
        want[[tok]],
        label = tok
      )
    }
  }
})

test_that("a token with two real genes stays ambiguous in human", {
  for (tok in names(ambiguous_human)) {
    out <- run(tok)
    expect_equal(as.vector(out)[1], tok, label = tok)
    expect_false(attr(out, "corrected")[1], label = tok)
    expect_equal(attr(out, "ambiguous"), tok, label = tok)
  }
})

test_that("a token with two real genes stays ambiguous in mouse", {
  for (tok in names(ambiguous_mouse)) {
    out <- run(tok, mm_now(), species = "mouse")
    expect_equal(as.vector(out)[1], tok, label = tok)
    expect_equal(attr(out, "ambiguous"), tok, label = tok)
  }
})

test_that("both genes behind an ambiguous token are candidates", {
  cand <- function(tok, sp) {
    toupper(genevintage:::.excel_candidates(genevintage:::.excel_clean(tok), sp)[[1]])
  }
  for (tok in names(ambiguous_human)) {
    expect_true(all(ambiguous_human[[tok]] %in% cand(tok, "homo_sapiens")), label = tok)
  }
  for (tok in names(ambiguous_mouse)) {
    expect_true(all(ambiguous_mouse[[tok]] %in% cand(tok, "mus_musculus")), label = tok)
  }
})

test_that("year-style March rows match the table", {
  # these rows are ambiguous in the table, with MARS2 behind one of them
  x <- c("2001-03-01 00:00:00", "2002-03-01 00:00:00", "2003-03-01 00:00:00", "2005-09-01 00:00:00")
  out <- run(x, mm_now(), species = "mouse")
  expect_equal(as.vector(out)[1:4], c(x[1:2], "Marchf3", "Septin5"))
  expect_equal(attr(out, "ambiguous"), x[1:2])
})

test_that("rows the table corrected with the partner present still resolve", {
  expect_equal(fixed(c("1-Mar", "MTARC1")), c("MARCHF1", "MTARC1"))
  expect_equal(fixed(c("2-Mar", "MTARC2")), c("MARCHF2", "MTARC2"))
  expect_equal(mouse_file(c("2-Mar", "Mtarc2")), c("Marchf2", "Mtarc2"))
})

test_that("mouse 2-Sep needs Apoa1 in the list to be settled", {
  # MGI lists Sep2 as an old name of Apoa1, so 2-Sep has two candidates
  out <- run("2-Sep", mm_now(), species = "mouse")
  expect_equal(as.vector(out)[1], "2-Sep")
  expect_equal(attr(out, "ambiguous"), "2-Sep")
  expect_equal(mouse_file("2-Sep"), "Septin2")
  # the year form needs a second year in its list
  x <- c("2002-09-01 00:00:00", "2003-09-01 00:00:00")
  out <- run(x, mm_now(), species = "mouse")
  expect_equal(as.vector(out)[1:2], c(x[1], "Septin3"))
  expect_equal(attr(out, "ambiguous"), x[1])
  expect_equal(mouse_file(x), c("Septin2", "Septin3"))
})

test_that("the human annotation has no Apoa1 rule", {
  expect_equal(fixed("2-Sep", nom("SEPTIN2", "APOA1")), "SEPTIN2")
})

test_that("rows the table did not read as damage stay unchanged", {
  out <- correct_genenames(not_recognised, mapping = hs_now(), quiet = TRUE)
  expect_equal(as.vector(out), not_recognised)
  expect_false(any(attr(out, "corrected")))
  expect_length(attr(out, "unresolved"), 0L)
  expect_length(attr(out, "ambiguous"), 0L)
  expect_length(attr(out, "skipped"), 0L)
})

test_that("Dutch May rows match their Ensembl IDs", {
  # the IDs give Mei1 and Mei4
  x <- c("2001-05-01 00:00:00", "2004-05-01 00:00:00", "2003-09-01 00:00:00", "2005-09-01 00:00:00")
  out <- run(x, mm_now(), species = "mouse")
  expect_equal(as.vector(out)[1:4], c("Mei1", "Mei4", "Septin3", "Septin5"))
  expect_length(attr(out, "unresolved"), 0L)
})

test_that("the pairs hold against the real human annotation", {
  skip_net()
  m <- fetch_mapping("human", release = 116)
  for (file in names(truth_human)) {
    want <- truth_human[[file]]
    expect_equal(fixed(names(want), m, species = "human"), unname(want), label = file)
  }
  for (tok in c("1-Mar", "2-Mar")) {
    expect_equal(attr(run(tok, m, species = "human"), "ambiguous"), tok, label = tok)
  }
})

test_that("the pairs hold against the real mouse annotation", {
  skip_net()
  m <- fetch_mapping("mouse", release = 116)
  for (file in names(truth_mouse)) {
    want <- truth_mouse[[file]]
    expect_equal(toupper(mouse_file(names(want), m)), unname(want), label = file)
  }
  for (tok in c("1-Mar", "2-Mar")) {
    expect_equal(attr(run(tok, m, species = "mouse"), "ambiguous"), tok, label = tok)
  }
})
