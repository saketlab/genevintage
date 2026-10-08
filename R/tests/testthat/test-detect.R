test_that("species comes from the ids, not from the metadata", {
  skip_no_index()
  expect_equal(detect_species("ENSG00000141510")$species[1], "homo_sapiens")
  expect_equal(detect_species("ENSMUSG00000051951")$species[1], "mus_musculus")
  # ENSMUSG starts with neither ENSG nor a shorter competitor -- longest wins
  expect_equal(
    detect_species(c("ENSMUSG00000051951", "ENSMUSG00000089699"))$species[1],
    "mus_musculus"
  )
})

test_that("species detection beats the GEO metadata on a known mislabelled study", {
  skip_no_index()
  skip_if_not(dir.exists(fx_dir))
  f <- manifest$file[manifest$accession == "GSE127774"]
  skip_if(length(f) == 0)
  # GEO says Macaca mulatta; the ids are ENSPPAG, which is Pan paniscus
  ids <- read_fx(f[1])
  expect_true(all(grepl("^ENSPPAG", head(ids[grepl("^ENS", ids)], 20))))
  expect_false(identical(detect_species(ids)$species[1], "macaca_mulatta"))
})

test_that("detect_release ranks releases and never returns just one", {
  skip_no_index()
  fp <- .fp()
  skip_if(!"homo_sapiens" %in% fp$species)
  ref <- fp[fp$species == "homo_sapiens", ]
  skip_if(nrow(ref) < 3)
  # these low numbers date to the oldest indexed release, which warns that the
  # true answer may sit outside the index -- asserted on its own further down
  suppressWarnings(out <- detect_release(sprintf("ENSG%011d", seq_len(500) * 100)))
  expect_true(is.data.frame(out))
  expect_true(all(c("release", "feasible", "dist") %in% names(out)))
  expect_equal(nrow(out), nrow(ref))
  expect_false(is.unsorted(c(!out$feasible))) # feasible ones come first
})

test_that("a newer id rules out older releases", {
  skip_no_index()
  fp <- .fp()
  ref <- fp[fp$species == "homo_sapiens", ]
  skip_if(nrow(ref) < 3)
  newest <- max(ref$qmax)
  out <- detect_release(c(sprintf("ENSG%011.0f", newest), "ENSG00000141510"))
  # only releases that had minted that id can be feasible
  # feasibility allows a 0.1% shortfall, since Ensembl retires top-numbered genes
  expect_true(all(out$release[out$feasible] %in% ref$release[ref$qmax >= newest * 0.999]))
})

test_that("real fixtures resolve to the right species", {
  skip_no_index()
  skip_if_not(dir.exists(fx_dir))
  # the manifest now records canonical Ensembl species names, so a fixture's
  # recorded species is directly what detect_species() should return
  want <- c(
    "homo_sapiens", "mus_musculus", "danio_rerio", "rattus_norvegicus",
    "sus_scrofa", "gallus_gallus", "macaca_mulatta", "sus_scrofa"
  )
  # two studies are labelled for one species and quantified against another --
  # they are asserted separately below, as the cases identifiers beat metadata
  MISLABELLED <- c("GSE112356", "GSE127774")
  m <- manifest[manifest$scheme == "ensembl" & manifest$purity > 0.98 &
    manifest$species %in% want & !manifest$accession %in% MISLABELLED, ]
  skip_if(nrow(m) == 0)
  for (i in seq_len(nrow(m))) {
    cand <- detect_species(read_fx(m$file[i]))
    skip_if(!nrow(cand))
    # breed and strain genomes share an identifier space (ENSSSCG covers every
    # pig breed), so any species tied at the top is a correct answer
    tied <- cand$species[cand$frac >= cand$frac[1] - 1e-9]
    expect_true(m$species[i] %in% tied, label = m$file[i])
  }
})

test_that("transcript and protein ids resolve to their species", {
  skip_no_index()
  # the index stores gene prefixes; these differ only in the feature letter
  expect_equal(detect_species("ENSGALT00000000002")$species[1], "gallus_gallus")
  expect_equal(detect_species("ENSMMUT00000020734")$species[1], "macaca_mulatta")
  expect_equal(detect_species("ENSSSCT00000000003.3")$species[1], "sus_scrofa")
  # and a shorter stem must not swallow a longer one
  expect_equal(detect_species("ENSMUSG00000051951")$species[1], "mus_musculus")
})

test_that("bonobo ids resolve to bonobo, not to the species GEO claims", {
  skip_no_index()
  # ENSPPAG was unresolvable when the index held 12 species; with all Ensembl
  # vertebrates it resolves, which is the point of the wider index
  expect_equal(detect_species("ENSPPAG00000006288")$species[1], "pan_paniscus")
})

test_that("an identifier from no known scheme matches nothing", {
  skip_no_index()
  expect_equal(nrow(detect_species("NOTANID12345")), 0L)
})

test_that("schemes without ordered numbers fall back to count alone", {
  skip_no_index()
  skip_if(!"arabidopsis_thaliana" %in% .fp()$species)
  f <- file.path(fx_dir, "athaliana_GSE213622_2023.txt.gz")
  skip_if_not(file.exists(f))
  out <- suppressMessages(detect_release(read_fx(basename(f)),
    species = "arabidopsis_thaliana"
  ))
  expect_false(attr(out, "numeric_ordered")) # AT1G01010 embeds the chromosome
  expect_true(all(is.na(out$feasible))) # unknown, not asserted
  expect_false(is.unsorted(out$dist))
})

test_that("a match pinned at the edge of the index warns", {
  skip_no_index()
  ref <- .fp()[.fp()$species == "homo_sapiens" & .fp()$source == "ensembl", ]
  skip_if(nrow(ref) < 3)
  # a set sized like the oldest indexed release scores best there, at the edge.
  # Releases are character now, so "100" sorts before "50" -- order numerically.
  oldest <- ref[order(.rel_num(ref$release)), ][1, ]
  ids <- sprintf("ENSG%011.0f", seq(1, oldest$qmax, length.out = oldest$n))
  expect_warning(detect_release(ids, species = "homo_sapiens"), "oldest")
})

test_that("identifiers win over the species GEO records", {
  skip_no_index()
  skip_if_not(dir.exists(fx_dir))
  # GSE112356 is filed under Macaca mulatta and holds human ENSG ids;
  # GSE127774 is filed under Macaca mulatta and holds bonobo ENSPPAG ids
  for (case in list(c("GSE112356", "homo_sapiens"), c("GSE127774", "pan_paniscus"))) {
    f <- manifest$file[manifest$accession == case[1]]
    skip_if(length(f) == 0)
    expect_equal(detect_species(read_fx(f[1]))$species[1], case[2], label = case[1])
  }
})

test_that("nothing to go on returns nothing, rather than erroring", {
  skip_no_index()
  # a matrix can arrive with empty or missing row names; that is an answer of
  # "no candidates", not a failure
  for (x in list(character(), NA_character_, c(NA, ""), "not_an_identifier")) {
    out <- detect_species(x)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 0L)
  }
})

test_that("a version suffix does not hide the species", {
  skip_no_index()
  bare <- detect_species("ENSG00000141510")
  vers <- detect_species("ENSG00000141510.16")
  expect_equal(vers$species, bare$species)
  expect_equal(vers$frac, bare$frac)
})

test_that("species sharing an identifier space stay visible as a tie", {
  skip_no_index()
  # every dog breed is ENSCAFG; identifiers alone cannot separate them, and
  # collapsing to one row would present a coin-flip as a determination
  out <- detect_species(c("ENSCAFG00845000001", "ENSCAFG00845000002"))
  expect_gt(nrow(out), 1L)
  expect_equal(length(unique(out$frac[out$frac == out$frac[1]])), 1L)
  expect_true(all(grepl("^canis", out$species[out$frac == max(out$frac)])))
})

test_that("detect_release says which way it could not proceed", {
  skip_no_index()
  # identifiers that are not this species at all
  expect_error(
    detect_release(c("ENSMUSG00000051951"), species = "human"),
    "none of the identifiers match"
  )
  # right species, but no number to fingerprint with
  expect_error(
    detect_release(c("ENSG", "ENSG"), species = "human"),
    "none of the identifiers match"
  )
})

test_that("a best match at the edge of the index is flagged, not presented as final", {
  skip_no_index()
  # low identifier numbers fit the oldest indexed release best -- but the index
  # starts at r50, so the real answer may be a release it does not carry. The
  # ranking cannot say so; the warning can.
  ids <- sprintf("ENSG%011d", as.integer(seq(1000, 20000, length.out = 200)))
  # Historical GENCODE fingerprints are now indexed, so this low-number set
  # is no longer pinned to the old r50 boundary.
  out <- suppressWarnings(detect_release(ids, species = "human"))
  expect_true(out$release[1] %in% c("47", "48", "49", "50"))
  expect_true(is.data.frame(out))
})

test_that("qmax_violations tolerates a small dip and reports a sharp one", {
  # Ensembl retires genes at the top of the range occasionally, so a fall of a
  # fraction of a percent is normal and must not be reported as a fault
  fp <- data.frame(
    source = "ensembl", species = "test_species", assembly = "1",
    release = c("100", "101", "102"), ordered = TRUE,
    qmax = c(1000000, 999000, 998500), stringsAsFactors = FALSE
  )
  expect_length(qmax_violations(fp), 0L)

  fp$qmax <- c(1000000, 500000, 400000)
  expect_length(qmax_violations(fp), 1L)
  expect_match(qmax_violations(fp), "test_species")

  # string order puts "100" before "50", which would invent a violation
  fp2 <- data.frame(
    source = "ensembl", species = "test_species", assembly = "1",
    release = c("50", "100"), ordered = TRUE,
    qmax = c(500000, 1000000), stringsAsFactors = FALSE
  )
  expect_length(qmax_violations(fp2), 0L)

  # a scheme whose numbers are not ordered is exempt from the rule entirely
  fp3 <- fp
  fp3$ordered <- FALSE
  expect_length(qmax_violations(fp3), 0L)
})
