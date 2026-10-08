test_that("the corpus is present and diverse", {
  skip_if_not(dir.exists(fx_dir))
  expect_gt(nrow(manifest), 100)
  expect_gte(length(unique(manifest$species)), 12)
  expect_gte(length(unique(manifest$scheme)), 6)
  expect_true(any(manifest$frac_versioned > 0.9)) # versioned IDs are covered
  expect_true(any(manifest$frac_versioned < 0.01)) # so are unversioned
})

test_that("guess_scheme reproduces the manifest on every fixture", {
  skip_if_not(dir.exists(fx_dir))
  got <- vapply(manifest$file, function(f) guess_scheme(read_fx(f))$scheme, character(1))
  expect_equal(unname(got), manifest$scheme)
})

test_that("parse_ens splits versioned and unversioned IDs alike", {
  p <- parse_ens(c(
    "ENSG00000141510", "ENSG00000141510.16", "ENSMUSG00000051951.6",
    " ENSDARG00000102141 ", "ENSG00000182378.14_PAR_Y",
    "TP53", "AT1G01010", ""
  ))
  expect_equal(p$prefix, c("ENS", "ENS", "ENSMUS", "ENSDAR", "ENS", NA, NA, NA))
  expect_equal(p$feature, c("G", "G", "G", "G", "G", NA, NA, NA))
  expect_equal(p$num, c(141510, 141510, 51951, 102141, 182378, NA, NA, NA))
  expect_equal(p$id_version, c(NA, 16L, 6L, NA, 14L, NA, NA, NA))
  expect_equal(p$tag, c(NA, NA, NA, NA, "_PAR_Y", NA, NA, NA))
})

test_that("version suffixes parse on real data, not just handmade input", {
  skip_if_not(dir.exists(fx_dir))
  f <- manifest$file[manifest$frac_versioned > 0.9 & manifest$scheme == "ensembl"]
  skip_if(length(f) == 0)
  for (i in f) {
    expect_gt(mean(grepl("\\.\\d+$", read_fx(i))), 0.9, label = i)
  }
})

test_that("species is recoverable from the Ensembl prefix", {
  skip_if_not(dir.exists(fx_dir))
  ens <- manifest[manifest$scheme == "ensembl" & manifest$purity > 0.95, ]
  # only the species stem is wanted, so match it directly rather than
  # parsing every field of two million identifiers
  pfx <- vapply(ens$file, function(f) {
    m <- regmatches(
      read_fx(f),
      regexpr("^ENS[A-Z]{0,4}?(?=[EGTPR][0-9]{11})", read_fx(f), perl = TRUE)
    )
    names(sort(table(m), decreasing = TRUE))[1]
  }, character(1))
  # one prefix per species -- except mmulatta, which holds a known mislabelled study
  tab <- split(unname(pfx), ens$species)
  expect_equal(unique(tab$homo_sapiens), "ENS")
  expect_equal(unique(tab$mus_musculus), "ENSMUS")
  expect_equal(unique(tab$danio_rerio), "ENSDAR")
  expect_true("ENSPPAG" %in% paste0(tab$macaca_mulatta, "G")) # GSE127774 is bonobo
})

test_that(".id_number takes only a clean numeric tail", {
  x <- c(
    "ENSG00000141510.16", "WBGene00000003", "FBgn0267431",
    "si.dkey.219e24.7", # a real zebrafish clone symbol
    "AT1G01010", "YDR387C", "1e5", "TP53"
  )
  # anything but PREFIX+digits is NA -- "219e24" must never read as 2.19e26,
  # which would make every release infeasible in detect_release()
  # TP53 is letters-then-digits like any Ensembl id, so it yields 53 -- this is
  # a lexical helper; deciding whether that is a gene id is guess_scheme()'s job
  expect_equal(.id_number(x), c(141510, 3, 267431, NA, NA, NA, NA, 53))
})

test_that("parse_ens survives NA and empty input", {
  p <- parse_ens(c(NA_character_, "ENSG00000141510", ""))
  expect_equal(p$num, c(NA, 141510, NA))
  expect_equal(nrow(parse_ens(character())), 0L)
})

test_that("SGD covers mitochondrial and plasmid systematic names", {
  expect_equal(guess_scheme(c("YDR387C", "Q0010", "R0010W"))$scheme, "sgd")
})
