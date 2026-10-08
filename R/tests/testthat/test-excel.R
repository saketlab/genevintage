# corrections depend on vintage, so the annotation decides them

# an annotation of the old nomenclature, and of the 2019 renaming
old_nom <- data.frame(
  id = paste0("g", 1:5),
  name = c("SEPT9", "MARCH1", "DEC1", "TP53", "SELENOF"),
  chr = "1", biotype = "protein_coding", stringsAsFactors = FALSE
)
new_nom <- data.frame(
  id = paste0("g", 1:6),
  name = c("SEPTIN9", "MARCHF1", "DELEC1", "TP53", "SELENOF", "SEPTIN1"),
  chr = "1", biotype = "protein_coding", stringsAsFactors = FALSE
)

test_that("the same token corrects differently against different vintages", {
  # renamed in 2019, so no fixed table fits both old and current matrices
  tok <- c("9-Sep", "1-Mar", "1-Dec")
  expect_equal(
    as.vector(correct_genenames(tok, mapping = old_nom, quiet = TRUE)),
    c("SEPT9", "MARCH1", "DEC1")
  )
  expect_equal(
    as.vector(correct_genenames(tok, mapping = new_nom, quiet = TRUE)),
    c("SEPTIN9", "MARCHF1", "DELEC1")
  )
})

test_that("undamaged names are returned unchanged", {
  out <- correct_genenames(c("TP53", "SEPTIN9", "9-Sep"), mapping = new_nom, quiet = TRUE)
  expect_equal(as.vector(out), c("TP53", "SEPTIN9", "SEPTIN9"))
  expect_equal(attr(out, "corrected"), c(FALSE, FALSE, TRUE))
  expect_length(out, 3L)
})

test_that("date-like spelling variants are recognised", {
  for (tok in c(
    "9-Sep", "Sep-09", "09-Sep", "sept-09", "9/Sep", "09.Sep",
    "Sep 9", "9-Sep-2019", "'9-Sep'", "2019-09-09",
    "9/9/2019", "09.09.2019", "9-Sep-2019 10:30"
  )) {
    expect_equal(as.vector(correct_genenames(tok, mapping = new_nom, quiet = TRUE)),
      "SEPTIN9",
      label = tok
    )
  }
})

test_that("a five-digit serial is read as the date Excel stored", {
  # 42248 is 2015-09-01, so the day is 1 and the symbol is SEPTIN1
  expect_equal(
    as.vector(correct_genenames("42248", mapping = new_nom, quiet = TRUE)),
    "SEPTIN1"
  )
})

test_that("scientific notation recovers the exponent the symbol had", {
  # Excel shows E+19; seven mantissa digits shift the point six, so E13
  riken <- data.frame(
    id = "g1", name = "2310009E13Rik", chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_equal(
    as.vector(correct_genenames("2.310009E+19", mapping = riken, quiet = TRUE)),
    "2310009E13Rik"
  )
  bare <- data.frame(
    id = "g1", name = "2310009E13", chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_equal(
    as.vector(correct_genenames("2.310009E+19", mapping = bare, quiet = TRUE)),
    "2310009E13"
  )
})

test_that("three significant figures are not enough to recover a symbol", {
  # default display drops the digits; do not guess
  bare <- data.frame(
    id = "g1", name = "2310009E13", chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  out <- correct_genenames("2.31E+19", mapping = bare, quiet = TRUE)
  expect_equal(as.vector(out), "2.31E+19")
  expect_false(attr(out, "corrected"))
})

test_that("the annotation resolves which way round a slash date was written", {
  # 11/9 is 11 September or 9 November depending on locale; only one is a gene
  expect_equal(as.vector(correct_genenames("11/9/2020", mapping = data.frame(
    id = "g1", name = "SEPTIN11", chr = "1", biotype = "protein_coding",
    stringsAsFactors = FALSE
  ), quiet = TRUE)), "SEPTIN11")
})

test_that("a token matching two real symbols is left alone and reported", {
  both <- data.frame(
    id = c("g1", "g2"), name = c("SEPT9", "SEPTIN9"), chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_warning(
    out <- suppressMessages(correct_genenames("9-Sep", mapping = both)),
    "match more than one symbol"
  )
  expect_equal(as.vector(out), "9-Sep")
  expect_equal(attr(out, "ambiguous"), "9-Sep")
})

test_that("a damaged token with no surviving symbol is left alone", {
  # no SEPTIN30 in this annotation, so nothing produces 30-Sep
  out <- correct_genenames("30-Sep", mapping = new_nom, quiet = TRUE)
  expect_equal(as.vector(out), "30-Sep")
  expect_false(attr(out, "corrected"))
})

test_that("an annotation's own symbols are never mistaken for damage", {
  # a false positive silently renames a gene
  skip_net()
  for (sp in c("human", "mouse", "zebrafish", "yeast")) {
    m <- fetch_mapping(sp, release = 116)
    expect_equal(sum(!is.na(genevintage:::.excel_kind(stats::na.omit(m$name)))), 0L,
      label = sp
    )
  }
})

test_that("nothing damaged means nothing fetched", {
  local_mocked_bindings(.locate_annotation = function(...) stop("must not fetch"))
  out <- correct_genenames(c("TP53", "BRCA1"), quiet = TRUE)
  expect_equal(as.vector(out), c("TP53", "BRCA1"))
  expect_false(any(attr(out, "corrected")))
})

test_that("the repair is reported unless silenced", {
  expect_message(correct_genenames("9-Sep", mapping = new_nom), "repaired 1 of 1")
  expect_silent(correct_genenames("9-Sep", mapping = new_nom, quiet = TRUE))
})

test_that("a quoted or padded token is matched and extracted from the same string", {
  # regression: regexec matched the cleaned token, regmatches the raw one
  riken <- data.frame(
    id = "g1", name = "2310009E13", chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  for (tok in c("'2.310009E+19'", " 2.310009E+19 ", "\"2.310009E+19\"")) {
    expect_equal(as.vector(correct_genenames(tok, mapping = riken, quiet = TRUE)),
      "2310009E13",
      label = tok
    )
  }
})

test_that("a day outside a month is not treated as damage", {
  # SEPT0 exists in this annotation, so a lax day rule would "repair" 0-Sep to it
  m <- data.frame(
    id = c("g1", "g2", "g3"), name = c("SEPTIN9", "SEPT0", "MARCH32"),
    chr = "1", biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_equal(as.vector(correct_genenames("0-Sep", mapping = m, quiet = TRUE)), "0-Sep")
  expect_equal(as.vector(correct_genenames("32-Mar", mapping = m, quiet = TRUE)), "32-Mar")
  expect_equal(as.vector(correct_genenames("0/9/2020", mapping = m, quiet = TRUE)), "0/9/2020")
  expect_equal(as.vector(correct_genenames("9-Sep", mapping = m, quiet = TRUE)), "SEPTIN9")
})

test_that("serials inside Excel's 1900 leap bug are refused", {
  # serials up to 60 hit Excel's phantom 29 Feb 1900
  m <- data.frame(
    id = c("g1", "g2"), name = c("FEB27", "FEB28"), chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  for (tok in c("00059", "00060")) {
    expect_equal(as.vector(correct_genenames(tok, mapping = m, quiet = TRUE)), tok)
  }
})

test_that("an exponent too large for an integer does not abort the vector", {
  riken <- data.frame(
    id = "g1", name = "2310009E13", chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  out <- correct_genenames(c("2.310009E+9999999999", "2.310009E+19"),
    mapping = riken, quiet = TRUE
  )
  expect_equal(as.vector(out), c("2.310009E+9999999999", "2310009E13"))
})

test_that("damage left alone is reported, not silently passed through", {
  expect_message(out <- correct_genenames("30-Sep", mapping = new_nom), "left alone")
  expect_equal(attr(out, "unresolved"), "30-Sep")
  expect_length(attr(out, "ambiguous"), 0L)
})

test_that("one mixed vector keeps its order and its attributes", {
  m <- data.frame(
    id = paste0("g", 1:4),
    name = c("SEPTIN9", "MARCHF1", "SEPTIN1", "TP53"), chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  x <- c("TP53", "9-Sep", NA, "1-Mar", "", "42248", "15-Sep", "not-a-date")
  out <- suppressMessages(correct_genenames(x, mapping = m))
  expect_length(out, length(x))
  expect_equal(
    as.vector(out),
    c("TP53", "SEPTIN9", NA, "MARCHF1", "", "SEPTIN1", "15-Sep", "not-a-date")
  )
  expect_equal(
    attr(out, "corrected"),
    c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)
  )
  expect_equal(attr(out, "unresolved"), "15-Sep")
})

test_that("a slash date resolves whichever way round the annotation allows", {
  only_nov <- data.frame(
    id = "g1", name = "NOV9", chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  # 11/9 is 11 September or 9 November; here only the second is a gene
  expect_equal(
    as.vector(correct_genenames("11/9/2020", mapping = only_nov, quiet = TRUE)),
    "NOV9"
  )

  both <- data.frame(
    id = c("g1", "g2"), name = c("SEPTIN11", "NOV9"), chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_warning(
    o <- suppressMessages(correct_genenames("11/9/2020", mapping = both)),
    "more than one"
  )
  expect_equal(as.vector(o), "11/9/2020")

  neither <- data.frame(
    id = "g1", name = "TP53", chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_equal(
    as.vector(correct_genenames("11/9/2020", mapping = neither, quiet = TRUE)),
    "11/9/2020"
  )

  # 13 cannot be a month, so only one reading survives
  # HGNC renamed SEPT13 to SEPTIN7P2
  expect_equal(as.vector(correct_genenames("13/9/2020", mapping = nom("SEPTIN7P2"), quiet = TRUE)), "SEPTIN7P2")
})

test_that("input with nothing to repair never reaches the annotation", {
  local_mocked_bindings(.locate_annotation = function(...) stop("must not fetch"))
  for (x in list(character(), NA_character_, c("", ""), c("TP53", NA))) {
    out <- correct_genenames(x, quiet = TRUE)
    expect_length(out, length(x))
    expect_false(any(attr(out, "corrected")))
  }
})

test_that("the long month name is recognised", {
  expect_equal(
    as.vector(correct_genenames(c("September-9", "9-September"),
      mapping = new_nom, quiet = TRUE
    )),
    c("SEPTIN9", "SEPTIN9")
  )
})

test_that("quiet silences the ambiguity warning too", {
  expect_silent(out <- correct_genenames("9-Sep", mapping = nom("SEPT9", "SEPTIN9"), quiet = TRUE))
  expect_equal(attr(out, "ambiguous"), "9-Sep")
})

test_that("a duplicated annotation row is not a second candidate", {
  out <- correct_genenames("9-Sep", mapping = nom("SEPTIN9", "SEPTIN9"), quiet = TRUE)
  expect_equal(as.vector(out), "SEPTIN9")
  expect_length(attr(out, "ambiguous"), 0L)
})

test_that("a damaged token without mapping or species is refused", {
  expect_error(correct_genenames("9-Sep"), "species")
  expect_error(correct_genenames("9-Sep", release = 116), "species")
  # undamaged input needs no annotation
  expect_silent(correct_genenames("TP53", quiet = TRUE))
})

test_that("a species without a release corrects against the newest one", {
  skip_net()
  out <- correct_genenames(c("9-Sep", "1-Dec"), species = "human", quiet = TRUE)
  expect_equal(as.vector(out), c("SEPTIN9", "DELEC1"))
})

test_that("every HGNC rename of a date-like symbol is recovered from its date", {
  koh <- data.frame(
    old = c(
      "DEC1", "MARC1", "MARC2", paste0("MARCH", 1:11), paste0("SEPT", 1:14), "SEP15"
    ),
    tok = c(
      "1-Dec", "1-Mar", "2-Mar", paste0(1:11, "-Mar"), paste0(1:14, "-Sep"), "15-Sep"
    ),
    want = c(
      "DELEC1", "MTARC1", "MTARC2", paste0("MARCHF", 1:11),
      paste0("SEPTIN", 1:12), "SEPTIN7P2", "SEPTIN14", "SELENOF"
    )
  )
  for (i in seq_len(nrow(koh))) {
    expect_equal(
      as.vector(correct_genenames(koh$tok[i], mapping = nom(koh$want[i]), quiet = TRUE)),
      koh$want[i],
      label = koh$old[i]
    )
  }
})

test_that("renames come from the species' own authority", {
  # ZFIN: sept15 -> septin15 and sep15 -> selenof both became 15-Sep; humans have no septin 15
  zf <- function(...) correct_genenames(..., species = "zebrafish", quiet = TRUE)
  expect_equal(attr(suppressWarnings(zf("15-Sep", mapping = nom("septin15", "selenof"))), "ambiguous"), "15-Sep")
  expect_equal(as.vector(zf("15-Sep", mapping = nom("septin15"))), "septin15")
  # ZFIN renamed sept9 to septin9a, not septin9
  expect_equal(as.vector(zf("9-Sep", mapping = nom("septin9a"))), "septin9a")
  # MGI's Sep2 -> Apoa1 is mouse nomenclature and must not reach human
  expect_equal(
    as.vector(correct_genenames("2-Sep", mapping = nom("SEPTIN2", "APOA1"), species = "human", quiet = TRUE)),
    "SEPTIN2"
  )
})

test_that("HGNC previous symbols Ensembl never shipped still count", {
  # HGNC lists FEB3 (febrile convulsions 3) as a previous symbol of SCN1A
  expect_equal(as.vector(correct_genenames("3-Feb", mapping = nom("SCN1A"), quiet = TRUE)), "SCN1A")
})

test_that("MARCH1 and MARC1 collide on 1-Mar, and an undamaged row resolves it", {
  m <- nom("MARCHF1", "MTARC1")
  expect_warning(
    out <- correct_genenames("Mar-01", mapping = m, quiet = FALSE),
    "more than one"
  )
  expect_equal(as.vector(out), "Mar-01")
  expect_equal(
    as.vector(correct_genenames(c("Mar-01", "MTARC1"), mapping = m, quiet = TRUE)),
    c("MARCHF1", "MTARC1")
  )
  # the undamaged row may carry the old spelling
  expect_equal(
    as.vector(correct_genenames(c("Mar-01", "MARC1"), mapping = m, quiet = TRUE)),
    c("MARCHF1", "MARC1")
  )
  # both damaged: two identical tokens, nothing tells them apart
  out <- suppressWarnings(correct_genenames(c("1-Mar", "1-Mar"), mapping = m, quiet = TRUE))
  expect_equal(as.vector(out), c("1-Mar", "1-Mar"))
})
